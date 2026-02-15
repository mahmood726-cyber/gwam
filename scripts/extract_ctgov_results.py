#!/usr/bin/env python3
"""Extract effect sizes from ClinicalTrials.gov posted results for results-only studies.

For trials with has_results=True but no PMID, fetches the results section from
CT.gov API v2 and attempts to extract:
  - Binary outcomes: 2x2 table -> log(OR) + SE
  - Continuous outcomes: means/SDs/Ns -> SMD + SE
  - Direct effect estimates from analyses section (OR/RR + CI)

Output: enriched registry CSV with results_effect, results_se, results_effect_type columns.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import time
from pathlib import Path
from typing import Any

import requests

from gwam_utils import normal_quantile, parse_bool, phrase_in_text, safe_float, sanitize_csv_cell

API_BASE = "https://clinicaltrials.gov/api/v2/studies"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-csv", type=Path, required=True, help="Input registry CSV.")
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Output enriched CSV. Defaults to <input>_with_results.csv.",
    )
    parser.add_argument(
        "--outcome-keywords",
        default="",
        help="Comma-separated keywords to match primary outcomes (default: match any).",
    )
    parser.add_argument(
        "--sleep-sec",
        type=float,
        default=0.3,
        help="Sleep between API requests (seconds).",
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=30.0,
        help="HTTP timeout per request.",
    )
    parser.add_argument(
        "--max-trials",
        type=int,
        default=None,
        help="Optional cap on number of trials to process.",
    )
    parser.add_argument(
        "--require-completed",
        action="store_true",
        default=True,
        help="Only extract results from COMPLETED trials (default: True, this flag is "
             "implicit). Use --no-require-completed to include TERMINATED/SUSPENDED.",
    )
    parser.add_argument(
        "--no-require-completed",
        action="store_false",
        dest="require_completed",
    )
    return parser.parse_args()


_NCT_FORMAT = re.compile(r"^NCT\d{8}$", re.IGNORECASE)

# CT.gov arm types used to distinguish treatment from control
_CONTROL_TYPES = {"PLACEBO_COMPARATOR", "ACTIVE_COMPARATOR", "SHAM_COMPARATOR", "NO_INTERVENTION"}
_TREATMENT_TYPES = {"EXPERIMENTAL"}


def _order_groups_by_arm_type(
    group_ids: list[str],
    outcome_groups: list[dict],
    arm_groups: list[dict],
) -> tuple[str, str]:
    """Return (treatment_gid, control_gid) using armGroups[].type when available.

    Falls back to first-two ordering when arm types cannot resolve the pair.
    """
    if len(group_ids) < 2:
        raise ValueError("Need at least 2 group IDs")

    # Build a map from outcome-group title to armGroup type
    # Outcome groups have 'id' and 'title'; armGroups have 'label' and 'type'
    arm_type_by_label: dict[str, str] = {}
    for ag in arm_groups:
        label = str(ag.get("label", "")).strip().upper()
        atype = str(ag.get("type", "")).strip().upper()
        if label and atype:
            arm_type_by_label[label] = atype

    gid_type: dict[str, str] = {}
    for og in outcome_groups:
        gid = og.get("id", "")
        title = str(og.get("title", "")).strip().upper()
        if gid in group_ids and title in arm_type_by_label:
            gid_type[gid] = arm_type_by_label[title]

    treatment_ids = [g for g in group_ids if gid_type.get(g) in _TREATMENT_TYPES]
    control_ids = [g for g in group_ids if gid_type.get(g) in _CONTROL_TYPES]

    if treatment_ids and control_ids:
        return treatment_ids[0], control_ids[0]

    # Fallback: first two in original order
    return group_ids[0], group_ids[1]


def fetch_results_section(
    nct_id: str, timeout_sec: float, *, attempts: int = 3,
    session: requests.Session | None = None,
) -> dict | None:
    """Fetch the results section for a trial from CT.gov API v2 (with retry)."""
    if not _NCT_FORMAT.match(nct_id):
        return None
    url = f"{API_BASE}/{nct_id}"
    params = {
        "fields": (
            "resultsSection.outcomeMeasuresModule,"
            "protocolSection.armsInterventionsModule,"
            "protocolSection.statusModule.overallStatus"
        ),
    }
    getter = session or requests
    for attempt in range(1, attempts + 1):
        try:
            resp = getter.get(url, params=params, timeout=timeout_sec)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, ValueError) as exc:
            # ValueError covers json.JSONDecodeError (subclass) for non-JSON responses
            if attempt == attempts:
                print(f"  [WARN] {nct_id}: request failed after {attempts} attempts ({type(exc).__name__}: {exc})")
                return None
            sleep_s = min(8.0, 1.5 * (2 ** (attempt - 1)))
            time.sleep(sleep_s)


def _match_outcome(title: str, keywords: list[str]) -> bool:
    """Check if outcome title matches any keyword (word-boundary matching)."""
    if not keywords:
        return True
    return any(phrase_in_text(kw, title) for kw in keywords)


def _build_denom_map(denoms: list[dict]) -> dict[str, int]:
    """Parse denominator counts from CT.gov outcome denoms list."""
    denom_map: dict[str, int] = {}
    for denom in denoms:
        for count in denom.get("counts", []):
            gid = count.get("groupId", "")
            val = safe_float(count.get("value"))
            if math.isfinite(val) and val > 0:
                denom_map[gid] = int(val)
    return denom_map


def extract_binary_from_measurements(
    outcome: dict,
    arm_groups: list[dict],
) -> tuple[float, float] | None:
    """Extract log(OR) + SE from 2x2 table in outcome measurements.

    Returns (log_or, se_log_or) or None.
    """
    classes = outcome.get("classes", [])
    if not classes:
        return None

    # Try to identify treatment and control groups
    groups = outcome.get("groups", [])
    if len(groups) < 2:
        return None

    denom_map = _build_denom_map(outcome.get("denoms", []))

    # Look for event counts across groups
    for cls in classes:
        categories = cls.get("categories", [])
        for cat in categories:
            measurements = cat.get("measurements", [])
            if len(measurements) < 2:
                continue

            events: dict[str, float] = {}
            for m in measurements:
                gid = m.get("groupId", "")
                val = safe_float(m.get("value"))
                if math.isfinite(val):
                    events[gid] = val

            if len(events) < 2:
                continue

            # Identify treatment vs control using arm type metadata
            gids = list(events.keys())
            if len(gids) < 2:
                continue
            g1, g2 = _order_groups_by_arm_type(gids[:2], groups, arm_groups)
            n1 = denom_map.get(g1, 0)
            n2 = denom_map.get(g2, 0)
            if n1 <= 0 or n2 <= 0:
                continue

            e1 = events[g1]
            e2 = events[g2]
            if not (0 <= e1 <= n1 and 0 <= e2 <= n2):
                continue

            # Double-zero: both arms zero events or both arms all events — OR undefined
            if (e1 == 0 and e2 == 0) or (e1 == n1 and e2 == n2):
                continue

            # Haldane-Anscombe correction: apply only when at least one cell is zero
            has_zero = (e1 == 0 or e2 == 0 or e1 == n1 or e2 == n2)
            cc = 0.5 if has_zero else 0.0
            a = e1 + cc
            b = (n1 - e1) + cc
            c = e2 + cc
            d = (n2 - e2) + cc
            if a <= 0 or b <= 0 or c <= 0 or d <= 0:
                continue
            log_or = math.log((a * d) / (b * c))
            se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
            return log_or, se

    return None


def extract_continuous_from_measurements(
    outcome: dict,
    arm_groups: list[dict] | None = None,
) -> tuple[float, float] | None:
    """Extract SMD + SE from means/SDs in outcome measurements.

    Returns (smd, se_smd) or None.

    CT.gov API v2 structure: dispersionType is at the **outcome** level,
    while the numeric spread value is at the measurement level.
    """
    classes = outcome.get("classes", [])
    if not classes:
        return None

    groups = outcome.get("groups", [])
    if len(groups) < 2:
        return None

    # Dispersion type is at outcome level (not measurement level)
    # CT.gov API v2 values: "Standard Deviation", "Standard Error of the Mean", etc.
    dispersion_type = str(outcome.get("dispersionType", "")).upper()
    is_sd = "STANDARD DEVIATION" in dispersion_type
    is_se = "STANDARD ERROR" in dispersion_type
    if not (is_sd or is_se):
        return None

    denom_map = _build_denom_map(outcome.get("denoms", []))

    for cls in classes:
        categories = cls.get("categories", [])
        for cat in categories:
            measurements = cat.get("measurements", [])
            if len(measurements) < 2:
                continue

            means: dict[str, float] = {}
            spreads: dict[str, float] = {}
            for m in measurements:
                gid = m.get("groupId", "")
                val = safe_float(m.get("value"))
                spread = safe_float(m.get("spread"))

                if math.isfinite(val):
                    means[gid] = val
                if math.isfinite(spread) and spread > 0:
                    spreads[gid] = spread

            if len(means) < 2 or len(spreads) < 2:
                continue

            gids = list(means.keys())
            if len(gids) < 2:
                continue
            g1, g2 = _order_groups_by_arm_type(gids[:2], groups, arm_groups or [])
            if g1 not in spreads or g2 not in spreads:
                continue

            n1 = denom_map.get(g1, 0)
            n2 = denom_map.get(g2, 0)
            if n1 <= 1 or n2 <= 1:
                continue

            mean1, mean2 = means[g1], means[g2]
            spread1, spread2 = spreads[g1], spreads[g2]

            # Convert SE to SD if dispersion type is Standard Error
            if is_se:
                sd1 = spread1 * math.sqrt(n1)
                sd2 = spread2 * math.sqrt(n2)
            else:
                sd1 = spread1
                sd2 = spread2

            if sd1 <= 0 or sd2 <= 0:
                continue

            # Pooled SD (Hedges' g denominator)
            df = n1 + n2 - 2
            s_pooled = math.sqrt(((n1 - 1) * sd1 ** 2 + (n2 - 1) * sd2 ** 2) / df)
            if s_pooled <= 0:
                continue

            smd = (mean1 - mean2) / s_pooled
            # Hedges' correction
            j = 1.0 - 3.0 / (4.0 * df - 1.0) if df > 1 else 1.0
            smd_corrected = smd * j
            # SE of Hedges' g using Hedges & Olkin (1985) convention: df = n1+n2-2.
            # Borenstein et al. (2009) uses n1+n2 in the denominator — both are valid
            # large-sample approximations; the H&O form is slightly more precise.
            se = math.sqrt((n1 + n2) / (n1 * n2) + smd_corrected ** 2 / (2.0 * df))
            return smd_corrected, se

    return None


def _ci_z_divisor(analysis: dict) -> float:
    """Compute the z-divisor for CI-to-SE conversion based on ciPctValue and ciNumSides.

    Reads ciPctValue (e.g., 95, 90, 99) and ciNumSides (ONE_SIDED, TWO_SIDED)
    from the analysis dict. Defaults to 95% two-sided (z_divisor = 2*1.96).

    NOTE on one-sided CIs: When ciNumSides=ONE_SIDED, the interpretation of
    (ciLowerLimit, ciUpperLimit) is ambiguous. CT.gov rarely reports one-sided
    CIs. If encountered, this function computes z from the full one-sided alpha
    and returns 2*z, which assumes the width spans both bounds. This may
    overestimate SE for truly one-sided intervals. In practice, the CI bracket
    validation in callers (ci_lower <= param <= ci_upper) catches most invalid
    one-sided entries.
    """
    ci_pct = safe_float(analysis.get("ciPctValue"))
    if not math.isfinite(ci_pct) or ci_pct <= 0 or ci_pct >= 100:
        ci_pct = 95.0

    num_sides = str(analysis.get("ciNumSides", "")).upper()
    if "ONE" in num_sides:
        # One-sided CI: alpha = (100 - ci_pct)/100
        alpha = (100.0 - ci_pct) / 100.0
    else:
        # Two-sided CI (default): alpha/2
        alpha = (100.0 - ci_pct) / 200.0

    z = normal_quantile(1.0 - alpha)
    return 2.0 * z


def extract_from_analyses(outcome: dict) -> tuple[float, float, str] | None:
    """Extract effect estimate from analyses section.

    Returns (effect, se, effect_type) or None.
    """
    analyses = outcome.get("analyses", [])
    for analysis in analyses:
        param_type = str(analysis.get("paramType", "")).strip()
        param_value = safe_float(analysis.get("paramValue"))
        ci_lower = safe_float(analysis.get("ciLowerLimit"))
        ci_upper = safe_float(analysis.get("ciUpperLimit"))

        if not math.isfinite(param_value):
            continue

        # Validate CI brackets: lower < upper
        if math.isfinite(ci_lower) and math.isfinite(ci_upper) and ci_lower >= ci_upper:
            continue

        z_divisor = _ci_z_divisor(analysis)
        if z_divisor <= 0:
            continue

        param_upper = param_type.upper()
        # Ratio measures: convert to log scale
        if any(kw in param_upper for kw in ["ODDS RATIO", "RISK RATIO", "HAZARD RATIO"]):
            if param_value <= 0:
                continue
            log_effect = math.log(param_value)
            if math.isfinite(ci_lower) and math.isfinite(ci_upper) and ci_lower > 0 and ci_upper > 0:
                # Validate bracket: ci_lower <= param <= ci_upper
                if not (ci_lower <= param_value <= ci_upper):
                    continue
                se = (math.log(ci_upper) - math.log(ci_lower)) / z_divisor
            else:
                continue
            effect_type = "log_" + param_upper.replace(" ", "_").lower()
            return log_effect, se, effect_type

        # Difference measures: use directly
        if any(kw in param_upper for kw in ["MEAN DIFFERENCE", "DIFFERENCE", "RISK DIFFERENCE"]):
            if math.isfinite(ci_lower) and math.isfinite(ci_upper):
                # Validate bracket: ci_lower <= param <= ci_upper
                if not (ci_lower <= param_value <= ci_upper):
                    continue
                se = (ci_upper - ci_lower) / z_divisor
            else:
                continue
            effect_type = param_upper.replace(" ", "_").lower()
            return param_value, se, effect_type

    return None


def extract_outcome_effect(
    payload: dict,
    outcome_keywords: list[str],
) -> tuple[float, float, str] | None:
    """Extract the best available effect estimate from CT.gov results.

    Priority:
    1. Analyses section (pre-computed OR/RR + CI)
    2. Binary 2x2 table from measurements
    3. Continuous SMD from measurements

    Returns (effect, se, effect_type) or None.
    """
    results_section = payload.get("resultsSection", {})
    outcome_module = results_section.get("outcomeMeasuresModule", {})
    outcomes = outcome_module.get("outcomeMeasures", [])
    if not outcomes:
        return None

    arms_module = payload.get("protocolSection", {}).get("armsInterventionsModule", {})
    arm_groups = arms_module.get("armGroups", []) if arms_module else []

    # Filter to primary outcomes matching keywords
    primary_outcomes = [
        o for o in outcomes
        if str(o.get("type", "")).upper() == "PRIMARY"
        and _match_outcome(str(o.get("title", "")), outcome_keywords)
    ]

    # Fall back to all matching outcomes if no primary
    if not primary_outcomes:
        primary_outcomes = [
            o for o in outcomes
            if _match_outcome(str(o.get("title", "")), outcome_keywords)
        ]

    if not primary_outcomes:
        primary_outcomes = outcomes  # last resort: try all

    for outcome in primary_outcomes:
        # Try analyses section first (most reliable)
        result = extract_from_analyses(outcome)
        if result is not None:
            return result

        # Try binary extraction
        binary = extract_binary_from_measurements(outcome, arm_groups)
        if binary is not None:
            return binary[0], binary[1], "log_or_2x2"

        # Try continuous extraction
        continuous = extract_continuous_from_measurements(outcome, arm_groups)
        if continuous is not None:
            return continuous[0], continuous[1], "smd_hedges_g"

    return None


def main() -> int:
    args = parse_args()
    if args.output_csv is None:
        stem = args.registry_csv.stem
        args.output_csv = args.registry_csv.parent / f"{stem}_with_results.csv"
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    outcome_keywords = [kw.strip() for kw in args.outcome_keywords.split(",") if kw.strip()]

    with args.registry_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if not rows:
        raise ValueError("Registry CSV is empty.")

    # Add new columns
    new_cols = ["results_effect", "results_se", "results_effect_type"]
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    n_attempted = 0
    n_extracted = 0
    n_failed = 0

    session = requests.Session()
    session.headers["User-Agent"] = "GWAM-pipeline/1.0 (clinicaltrials.gov registry analysis)"

    for row in rows:
        row.setdefault("results_effect", "")
        row.setdefault("results_se", "")
        row.setdefault("results_effect_type", "")

        has_results = parse_bool(row.get("has_results", ""))
        has_pmid = parse_bool(row.get("has_pmid", ""))
        is_ghost = parse_bool(row.get("is_ghost_protocol", ""))

        # Only extract for results-only studies (has_results but no PMID)
        if not (has_results and not has_pmid and not is_ghost):
            continue

        if args.max_trials is not None and n_attempted >= args.max_trials:
            break

        nct_id = str(row.get("nct_id", "")).strip()
        if not nct_id:
            continue

        n_attempted += 1

        payload = fetch_results_section(nct_id, args.timeout_sec, session=session)
        if payload is None:
            n_failed += 1
            if args.sleep_sec > 0:
                time.sleep(args.sleep_sec)
            continue

        # Filter by trial completion status
        if args.require_completed:
            status = (
                payload.get("protocolSection", {})
                .get("statusModule", {})
                .get("overallStatus", "")
            )
            if str(status).upper() != "COMPLETED":
                n_failed += 1
                if args.sleep_sec > 0:
                    time.sleep(args.sleep_sec)
                continue

        result = extract_outcome_effect(payload, outcome_keywords)
        if result is not None:
            effect, se, effect_type = result
            row["results_effect"] = str(effect)
            row["results_se"] = str(se)
            row["results_effect_type"] = effect_type
            n_extracted += 1
            print(f"  {nct_id}: {effect_type} = {effect:.4f} (SE={se:.4f})")
        else:
            n_failed += 1

        if args.sleep_sec > 0:
            time.sleep(args.sleep_sec)

    # Sanitize string cells against spreadsheet formula injection
    for row in rows:
        for key, val in row.items():
            if isinstance(val, str):
                row[key] = sanitize_csv_cell(val)

    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output_csv}")
    print(
        f"Results extraction: attempted={n_attempted}, "
        f"extracted={n_extracted}, failed={n_failed}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
