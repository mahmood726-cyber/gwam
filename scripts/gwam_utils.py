#!/usr/bin/env python3
"""Shared GWAM utility logic for filtering and classification."""

from __future__ import annotations

import math
import re
import sys
import time
from typing import Any, Iterable

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
