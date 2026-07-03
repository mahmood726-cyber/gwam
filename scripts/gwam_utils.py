#!/usr/bin/env python3
"""Shared GWAM utility logic for filtering and classification."""

from __future__ import annotations

import math
import re
import sys
import time
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

if sys.version_info < (3, 10):
    raise SystemExit("GWAM requires Python >= 3.10 (uses X | Y union type syntax).")

# Small curated alias map for common intervention names.
INTERVENTION_ALIASES: dict[str, list[str]] = {
    "escitalopram": ["lexapro", "cipralex"],
    "sertraline": ["zoloft", "lustral"],
    "duloxetine": ["cymbalta"],
    "quetiapine": ["seroquel"],
    "pregabalin": ["lyrica"],
    "oseltamivir": ["tamiflu"],
}


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def safe_float(value: Any) -> float:
    """Parse a value to float, returning NaN on failure or non-finite input."""
    if value is None:
        return float("nan")
    try:
        out = float(value)
        return out if math.isfinite(out) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def normal_cdf(x: float) -> float:
    """Standard normal CDF (scalar, via scipy.stats.norm.cdf)."""
    from scipy.stats import norm
    return float(norm.cdf(x))


_CSV_FORMULA_CHARS = frozenset("=+@\t\r")


def sanitize_csv_cell(value: str) -> str:
    """Prepend apostrophe to cells that could trigger spreadsheet formula injection.

    Note: minus sign (-) is intentionally NOT escaped to preserve negative numeric
    values in medical data (e.g. "-0.5 mmHg").
    """
    if value and value[0] in _CSV_FORMULA_CHARS:
        return "'" + value
    return value


def sanitize_path_component(text: str) -> str:
    """Remove path traversal and OS-invalid characters from a path component."""
    # Strip path separators and parent-directory references
    text = text.replace("/", "_").replace("\\", "_").replace("..", "_")
    # Remove null bytes and Windows-invalid filename characters
    text = re.sub(r'[\0<>:"|?*]', "", text)
    return text


def normal_quantile(p: float) -> float:
    """Inverse normal CDF (exact, via scipy.stats.norm.ppf)."""
    from scipy.stats import norm
    return float(norm.ppf(p))


def normalize_text(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def phrase_in_text(phrase: str, text: str) -> bool:
    phrase_n = normalize_text(phrase)
    text_n = normalize_text(text)
    if not phrase_n or not text_n:
        return False
    return re.search(rf"\b{re.escape(phrase_n)}\b", text_n) is not None


def parse_timeframe_values_to_weeks(timeframe: str) -> list[float]:
    """Parse all numeric timeframe mentions into weeks (best-effort)."""
    if not timeframe:
        return []
    text = timeframe.lower()
    matches: list[tuple[float, str]] = []

    for value_raw, unit in re.findall(
        r"(\d+(?:\.\d+)?)\s*(day|days|week|weeks|month|months|year|years)\b",
        text,
    ):
        matches.append((float(value_raw), unit))
    for unit, value_raw in re.findall(
        r"\b(day|days|week|weeks|month|months|year|years)\s*(\d+(?:\.\d+)?)",
        text,
    ):
        matches.append((float(value_raw), unit))
    if not matches:
        return []
    values_weeks: list[float] = []
    for value, unit in matches:
        if unit.startswith("day"):
            values_weeks.append(value / 7.0)
        elif unit.startswith("week"):
            values_weeks.append(value)
        elif unit.startswith("month"):
            values_weeks.append(value * 4.348)  # 365.25/12/7
        elif unit.startswith("year"):
            values_weeks.append(value * 52.179)  # 365.25/7
    return values_weeks


def parse_timeframe_to_weeks(timeframe: str) -> float | None:
    """Backward-compatible helper returning the furthest parsed horizon."""
    values_weeks = parse_timeframe_values_to_weeks(timeframe)
    if not values_weeks:
        return None
    return max(values_weeks)


def _canonicalize_phase(phase: str) -> str:
    return re.sub(r"\s+", "", str(phase).upper().strip())


def trial_passes_estimand_filters(
    protocol: dict,
    *,
    allowed_phases: Iterable[str] | None,
    required_outcome_keywords: Iterable[str] | None,
    outcome_timeframe_min_weeks: float | None,
    outcome_timeframe_max_weeks: float | None,
    missing_outcome_policy: str = "include_as_unknown",
) -> bool:
    if missing_outcome_policy not in {"exclude", "include_as_unknown"}:
        raise ValueError(f"Unsupported missing_outcome_policy: {missing_outcome_policy}")
    phases = protocol.get("designModule", {}).get("phases", []) or []
    allowed_phase_set = {
        _canonicalize_phase(p)
        for p in (allowed_phases or [])
        if str(p).strip()
    }
    if allowed_phase_set:
        protocol_phases = {_canonicalize_phase(p) for p in phases if str(p).strip()}
        if not protocol_phases.intersection(allowed_phase_set):
            return False

    need_outcome_keyword = any(str(k).strip() for k in (required_outcome_keywords or []))
    need_timeframe = outcome_timeframe_min_weeks is not None or outcome_timeframe_max_weeks is not None
    if not (need_outcome_keyword or need_timeframe):
        return True

    outcomes = protocol.get("outcomesModule", {}).get("primaryOutcomes", []) or []
    if not outcomes:
        return missing_outcome_policy == "include_as_unknown"

    keywords = [str(k).strip() for k in (required_outcome_keywords or []) if str(k).strip()]
    keyword_matched_outcomes: list[dict] = []
    for outcome in outcomes:
        text = " ".join(
            [
                str(outcome.get("measure", "")),
                str(outcome.get("description", "")),
            ]
        )
        if not keywords:
            keyword_matched_outcomes.append(outcome)
        elif any(phrase_in_text(keyword, text) for keyword in keywords):
            keyword_matched_outcomes.append(outcome)

    if need_outcome_keyword and not keyword_matched_outcomes:
        return False

    if not need_timeframe:
        return True

    timeframe_candidates = keyword_matched_outcomes if keyword_matched_outcomes else outcomes
    timeframe_values: list[float] = []
    for outcome in timeframe_candidates:
        timeframe = str(outcome.get("timeFrame", ""))
        timeframe_values.extend(parse_timeframe_values_to_weeks(timeframe))
    if not timeframe_values:
        return missing_outcome_policy == "include_as_unknown"

    min_w = outcome_timeframe_min_weeks
    max_w = outcome_timeframe_max_weeks
    for weeks in timeframe_values:
        if min_w is not None and weeks < min_w:
            continue
        if max_w is not None and weeks > max_w:
            continue
        return True
    return False


def intervention_aliases_for(intervention: str, extra_aliases: Iterable[str] | None = None) -> list[str]:
    key = normalize_text(intervention)
    aliases = [intervention]
    aliases.extend(INTERVENTION_ALIASES.get(key, []))
    if extra_aliases:
        aliases.extend([a for a in extra_aliases if a and a.strip()])
    # Preserve order but deduplicate normalized forms.
    seen: set[str] = set()
    out: list[str] = []
    for alias in aliases:
        norm = normalize_text(alias)
        if norm and norm not in seen:
            seen.add(norm)
            out.append(alias)
    return out


def intervention_matches(
    intervention_names: Iterable[str],
    *,
    intervention: str,
    extra_aliases: Iterable[str] | None = None,
) -> bool:
    names = [n for n in intervention_names if n and str(n).strip()]
    if not names:
        return False
    aliases = intervention_aliases_for(intervention, extra_aliases=extra_aliases)
    for name in names:
        if any(phrase_in_text(alias, str(name)) for alias in aliases):
            return True
    return False


def has_placebo_marker(protocol: dict) -> bool:
    arms = protocol.get("armsInterventionsModule", {}).get("armGroups", []) or []
    interventions = protocol.get("armsInterventionsModule", {}).get("interventions", []) or []
    texts: list[str] = []
    for arm in arms:
        arm_type = str(arm.get("type", "")).upper()
        if "PLACEBO" in arm_type:
            return True
        texts.append(str(arm.get("label", "")))
        texts.append(str(arm.get("description", "")))
        texts.append(str(arm.get("interventionNames", "")))
    for intervention in interventions:
        int_type = str(intervention.get("type", "")).upper()
        if "PLACEBO" in int_type:
            return True
        texts.append(str(intervention.get("name", "")))
        texts.append(str(intervention.get("description", "")))
        texts.append(str(intervention.get("type", "")))
    return any(phrase_in_text("placebo", text) for text in texts)


def has_active_comparator_marker(protocol: dict) -> bool:
    """Detect if a trial has an active comparator arm.

    NOTE: The fallback heuristic (non-placebo DRUG/BIOLOGICAL interventions) may
    over-include some trials. Combined with the ``has_placebo`` exclusion in
    ``trial_passes_design_filters``, this works correctly for the intended use case
    (finding head-to-head trials that lack a placebo arm). The combined filter
    ``has_active AND NOT has_placebo`` is the operative gate.
    """
    arms = protocol.get("armsInterventionsModule", {}).get("armGroups", []) or []
    interventions = protocol.get("armsInterventionsModule", {}).get("interventions", []) or []

    non_placebo_arms = 0
    for arm in arms:
        arm_type = str(arm.get("type", "")).upper()
        if "ACTIVE" in arm_type and "COMPARATOR" in arm_type:
            return True
        if arm_type == "COMPARATOR":
            return True
        label = str(arm.get("label", ""))
        desc = str(arm.get("description", ""))
        if phrase_in_text("active comparator", f"{label} {desc}"):
            return True
        if not phrase_in_text("placebo", f"{label} {desc}"):
            non_placebo_arms += 1

    for intervention in interventions:
        int_type = str(intervention.get("type", "")).upper()
        if "PLACEBO" in int_type:
            continue
        if int_type in {"DRUG", "BIOLOGICAL", "DEVICE", "PROCEDURE"}:
            return True

    # Fallback: two-arm non-placebo interventional design is likely an active-comparator setup.
    return non_placebo_arms >= 2


def extract_pmids(
    refs: Iterable[dict],
    *,
    publication_linkage: str,
) -> tuple[list[str], list[str], bool]:
    pmids_any: list[str] = []
    pmids_results: list[str] = []
    seen_any: set[str] = set()
    seen_results: set[str] = set()

    for ref in refs:
        pmid = str(ref.get("pmid", "")).strip()
        if not pmid:
            continue
        if pmid not in seen_any:
            seen_any.add(pmid)
            pmids_any.append(pmid)

        ref_type = str(ref.get("type", "")).upper().strip()
        if publication_linkage == "results_pmid_strict":
            linked_tokens = ("RESULT", "PRIMARY")
        elif publication_linkage == "results_pmid":
            # Broader linkage mode including DERIVED references. Note: DERIVED
            # references on CT.gov may point to secondary publications (subgroup
            # analyses, cost-effectiveness) rather than primary outcome reports.
            # Use "results_pmid_strict" for primary analysis to avoid over-counting.
            linked_tokens = ("RESULT", "PRIMARY", "DERIVED")
        elif publication_linkage == "any_pmid":
            linked_tokens = tuple()
        else:
            raise ValueError(f"Unsupported publication_linkage: {publication_linkage}")

        is_results_linked = any(token in ref_type for token in linked_tokens)
        if is_results_linked and pmid not in seen_results:
            seen_results.add(pmid)
            pmids_results.append(pmid)

    if publication_linkage in {"results_pmid_strict", "results_pmid"}:
        has_linked_publication = bool(pmids_results)
    elif publication_linkage == "any_pmid":
        has_linked_publication = bool(pmids_any)
    else:
        raise ValueError(f"Unsupported publication_linkage: {publication_linkage}")

    return pmids_any, pmids_results, has_linked_publication


def trial_passes_design_filters(
    protocol: dict,
    *,
    require_randomized: bool,
    require_treatment_purpose: bool,
    comparator_scope: str,
    allowed_phases: Iterable[str] | None = None,
    required_outcome_keywords: Iterable[str] | None = None,
    outcome_timeframe_min_weeks: float | None = None,
    outcome_timeframe_max_weeks: float | None = None,
    missing_outcome_policy: str = "include_as_unknown",
) -> bool:
    design = protocol.get("designModule", {}) or {}
    if str(design.get("studyType", "")).upper() != "INTERVENTIONAL":
        return False

    design_info = design.get("designInfo", {}) or {}
    allocation = str(design_info.get("allocation", "")).upper()
    primary_purpose = str(design_info.get("primaryPurpose", "")).upper()
    if require_randomized and allocation != "RANDOMIZED":
        return False
    if require_treatment_purpose and primary_purpose != "TREATMENT":
        return False

    arm_groups = protocol.get("armsInterventionsModule", {}).get("armGroups", []) or []
    has_placebo = has_placebo_marker(protocol)
    has_active = has_active_comparator_marker(protocol)
    if comparator_scope == "active":
        if len(arm_groups) < 2:
            return False
        if not has_active:
            return False
        if has_placebo:
            return False
    elif comparator_scope == "placebo":
        if len(arm_groups) < 2:
            return False
        if not has_placebo:
            return False
    elif comparator_scope != "any":
        raise ValueError(f"Unknown comparator_scope: {comparator_scope}")

    if not trial_passes_estimand_filters(
        protocol,
        allowed_phases=allowed_phases,
        required_outcome_keywords=required_outcome_keywords,
        outcome_timeframe_min_weeks=outcome_timeframe_min_weeks,
        outcome_timeframe_max_weeks=outcome_timeframe_max_weeks,
        missing_outcome_policy=missing_outcome_policy,
    ):
        return False
    return True


def classify_publication_status(
    *,
    has_pmid: bool,
    has_results: bool,
    has_any_pmid: bool = False,
    ghost_require_no_any_pmid: bool = False,
    ghost_definition: str,
) -> tuple[bool, bool]:
    """Return (is_unpublished_with_results, is_ghost).

    ghost_definition controls classification:
      - "no_pmid_no_results" (recommended for primary analysis): only trials with
        neither PMID nor posted results are ghosts. Results-only trials contribute
        observed data and are kept separate.
      - "no_pmid" (sensitivity analysis only): all trials without PMID are ghosts,
        including those with posted results. This overstates ghost prevalence and
        should only be used for worst-case sensitivity bounds.
    """
    is_unpublished_with_results = (not has_pmid) and has_results
    if ghost_definition == "no_pmid_no_results":
        is_ghost = (not has_pmid) and (not has_results)
    elif ghost_definition == "no_pmid":
        is_ghost = not has_pmid
    else:
        raise ValueError(f"Unsupported ghost_definition: {ghost_definition}")
    # Prevent classifying protocols with any PMID metadata as fully latent ghosts
    # unless explicitly requested. This avoids over-calling missingness when strict
    # linkage rules exclude non-results references.
    if is_ghost and ghost_require_no_any_pmid and has_any_pmid:
        is_ghost = False
    return is_unpublished_with_results, is_ghost


def build_environment_metadata() -> dict[str, str]:
    """Build environment metadata dict for JSON output traceability."""
    import datetime
    import platform
    import sys

    meta: dict[str, str] = {
        "python_version": sys.version.split()[0],
        "numpy_version": __import__("numpy").__version__,
        "scipy_version": __import__("scipy").__version__,
        "platform": platform.platform(),
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    try:
        # Try package-relative import first (when used as package)
        from . import __version__ as _gwam_ver
        meta["gwam_version"] = _gwam_ver
    except (ImportError, SystemError):
        try:
            # Fallback for standalone script execution
            import importlib
            _init = importlib.import_module("__init__")
            meta["gwam_version"] = _init.__version__
        except (ImportError, AttributeError):
            meta["gwam_version"] = "unknown"
    try:
        meta["pandas_version"] = __import__("pandas").__version__
    except ImportError:
        pass
    try:
        meta["requests_version"] = __import__("requests").__version__
    except ImportError:
        pass
    try:
        meta["pyreadr_version"] = __import__("pyreadr").__version__
    except ImportError:
        pass
    return meta


def request_json_with_retry(
    session: Any,
    *,
    url: str,
    params: dict[str, str | int],
    timeout: int,
    attempts: int = 5,
) -> dict:
    """GET a JSON endpoint with exponential-backoff retry."""
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            last_error = exc
            if attempt == attempts:
                break
            sleep_s = min(12.0, 1.5 * (2 ** (attempt - 1)))
            time.sleep(sleep_s)
    if last_error is None:
        raise RuntimeError("request_json_with_retry: no attempts made")
    raise last_error


# ---------------------------------------------------------------------------
# Reusable GWAM correction API
# ---------------------------------------------------------------------------
#
# The deterministic core of GWAM is a weighted correction: studies are
# partitioned into three integrity strata (linked-PMID, results-only, and
# fully-latent "ghost" protocols) using registry enrollment as the weight, and
# the corrected pooled mean is the enrollment-weighted average of each
# stratum's assumed mean. Historically this math lived only inside
# ``model_gwam.main()`` (CSV in / JSON out). The functions below expose it as
# an importable, side-effect-free API so callers can run the correction
# directly on in-memory rows. The arithmetic is intentionally identical to
# ``model_gwam.py`` so refactoring introduces no numeric drift.


@dataclass(frozen=True)
class GwamResult:
    """Result of a deterministic GWAM correction.

    Attributes:
        lambda_pmid_only: Enrollment-weighted fraction of studies with a linked
            results PMID (strictest integrity ratio).
        lambda_non_ghost: Enrollment-weighted fraction of non-ghost studies
            (linked-PMID plus results-only).
        mu_published: The published pooled effect passed in.
        mu_gwam_null_point: Deterministic GWAM-corrected pooled mean under the
            supplied stratum means (default: ghost/results-only mean = 0).
        weight_pmid_only: Total enrollment weight of the linked-PMID stratum.
        weight_results_only: Total enrollment weight of the results-only
            (no PMID) stratum.
        weight_ghost: Total enrollment weight of the ghost stratum.
        weight_total: Sum of all stratum weights.
        n_pmid_only: Count of studies in the linked-PMID stratum.
        n_results_only: Count of studies in the results-only stratum.
        n_ghost: Count of studies in the ghost stratum.
    """

    lambda_pmid_only: float
    lambda_non_ghost: float
    mu_published: float
    mu_gwam_null_point: float
    weight_pmid_only: float
    weight_results_only: float
    weight_ghost: float
    weight_total: float
    n_pmid_only: int
    n_results_only: int
    n_ghost: int

    def to_dict(self) -> dict[str, float | int]:
        """Return a plain-dict view suitable for JSON serialisation."""
        return {
            "lambda_pmid_only": self.lambda_pmid_only,
            "lambda_non_ghost": self.lambda_non_ghost,
            "mu_published": self.mu_published,
            "mu_gwam_null_point": self.mu_gwam_null_point,
            "weight_pmid_only": self.weight_pmid_only,
            "weight_results_only": self.weight_results_only,
            "weight_ghost": self.weight_ghost,
            "weight_total": self.weight_total,
            "n_pmid_only": self.n_pmid_only,
            "n_results_only": self.n_results_only,
            "n_ghost": self.n_ghost,
        }


def partition_registry_weights(
    rows: Sequence[Mapping[str, Any]],
    *,
    weight_column: str = "weight_n",
) -> tuple[list[float], list[float], list[float]]:
    """Partition registry rows into (pmid_only, results_only, ghost) weight lists.

    Mirrors the stratification logic in ``model_gwam.main()``:
    ghost protocols form their own stratum; non-ghost rows are split by
    ``has_pmid``; non-ghost rows lacking a PMID (with or without posted
    results) fall into the results-only stratum.

    Args:
        rows: Iterable of mapping-like registry rows (e.g. ``csv.DictReader``).
        weight_column: Column holding a positive, finite numeric weight.

    Returns:
        Three lists of floats: ``(pmid_weights, results_only_weights,
        ghost_weights)``.

    Raises:
        ValueError: If ``rows`` is empty, the weight column is missing, or any
            weight is non-numeric, non-finite, or <= 0.
    """
    if not rows:
        raise ValueError("Registry rows are empty.")
    if weight_column not in rows[0]:
        raise ValueError(f"Weight column '{weight_column}' not found in rows.")

    pmid_weights: list[float] = []
    results_only_weights: list[float] = []
    ghost_weights: list[float] = []
    for row_idx, row in enumerate(rows, start=2):  # row 2 = first data row after header
        raw_weight = row.get(weight_column)
        try:
            weight = float(raw_weight)  # type: ignore[arg-type]
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Non-numeric weight '{raw_weight}' in row {row_idx}, "
                f"column '{weight_column}'."
            ) from exc
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError(
                f"Invalid weight {weight} in row {row_idx}, column '{weight_column}'."
            )
        is_ghost = parse_bool(row.get("is_ghost_protocol", ""))
        has_pmid = parse_bool(row.get("has_pmid", ""))
        if is_ghost:
            ghost_weights.append(weight)
        elif has_pmid:
            pmid_weights.append(weight)
        else:
            # Non-ghost rows without a PMID (results-only or missing flags)
            # are treated as observed-unknown.
            results_only_weights.append(weight)
    return pmid_weights, results_only_weights, ghost_weights


def gwam_correct(
    rows: Sequence[Mapping[str, Any]],
    *,
    published_mu: float,
    weight_column: str = "weight_n",
    ghost_mu: float = 0.0,
    results_only_mode: str = "as_unknown",
    results_only_mu: float = 0.0,
) -> GwamResult:
    """Compute the deterministic GWAM correction from registry rows.

    This is the importable equivalent of the point-estimate layer in
    ``model_gwam.py``. It computes the two integrity ratios and the
    enrollment-weighted corrected pooled mean::

        mu_gwam = (w_pmid * mu_published
                   + w_results_only * mu_results_only
                   + w_ghost * ghost_mu) / w_total

    where ``mu_results_only == mu_published`` under ``results_only_mode
    == "as_observed"`` and ``results_only_mu`` otherwise.

    Args:
        rows: Registry rows (mapping-like) with ``is_ghost_protocol``,
            ``has_pmid`` and a positive weight column.
        published_mu: Published pooled effect (e.g. Hedges g). Must be finite.
        weight_column: Column holding the study weight (default ``weight_n``).
        ghost_mu: Assumed mean effect for ghost studies (default 0.0).
        results_only_mode: ``"as_unknown"`` (use ``results_only_mu``) or
            ``"as_observed"`` (treat results-only studies at ``published_mu``).
        results_only_mu: Assumed mean for results-only studies under
            ``as_unknown`` (default 0.0).

    Returns:
        A :class:`GwamResult`.

    Raises:
        ValueError: On empty input, missing/invalid weights, non-finite
            ``published_mu``, unknown ``results_only_mode``, or a non-positive
            total weight.
    """
    if not math.isfinite(float(published_mu)):
        raise ValueError(f"published_mu must be finite, got {published_mu!r}.")
    if results_only_mode not in {"as_unknown", "as_observed"}:
        raise ValueError(
            f"results_only_mode must be 'as_unknown' or 'as_observed', "
            f"got {results_only_mode!r}."
        )

    pmid_weights, results_only_weights, ghost_weights = partition_registry_weights(
        rows, weight_column=weight_column
    )

    w_pmid = float(sum(pmid_weights))
    w_results_only = float(sum(results_only_weights))
    w_ghost = float(sum(ghost_weights))
    w_total = w_pmid + w_results_only + w_ghost
    if not math.isfinite(w_total) or w_total <= 0:
        raise ValueError(f"Total weight is invalid ({w_total}); check input data.")

    lambda_pmid_only = w_pmid / w_total
    lambda_non_ghost = (w_pmid + w_results_only) / w_total

    mu_published = float(published_mu)
    mu_results_only = (
        mu_published if results_only_mode == "as_observed" else float(results_only_mu)
    )
    mu_gwam_null = (
        w_pmid * mu_published
        + w_results_only * mu_results_only
        + w_ghost * float(ghost_mu)
    ) / w_total

    return GwamResult(
        lambda_pmid_only=lambda_pmid_only,
        lambda_non_ghost=lambda_non_ghost,
        mu_published=mu_published,
        mu_gwam_null_point=mu_gwam_null,
        weight_pmid_only=w_pmid,
        weight_results_only=w_results_only,
        weight_ghost=w_ghost,
        weight_total=w_total,
        n_pmid_only=len(pmid_weights),
        n_results_only=len(results_only_weights),
        n_ghost=len(ghost_weights),
    )
