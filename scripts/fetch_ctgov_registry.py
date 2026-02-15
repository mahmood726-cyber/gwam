#!/usr/bin/env python3
"""Download completed ClinicalTrials.gov records for one intervention/condition pair.

Output is a flat trial-level CSV with ghost-protocol classification fields.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import time
from pathlib import Path

import requests

from gwam_utils import (
    classify_publication_status,
    extract_pmids,
    intervention_matches,
    parse_bool,
    sanitize_csv_cell,
    sanitize_path_component,
    trial_passes_design_filters,
)

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"


def parse_int(value: object) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intervention", required=True, help="Intervention name.")
    parser.add_argument("--condition", required=True, help="Condition name.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path. Defaults to data/raw/<intervention>__<condition>.csv",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=100,
        help="API page size (max 100).",
    )
    parser.add_argument(
        "--ghost-definition",
        choices=["no_pmid", "no_pmid_no_results"],
        default="no_pmid_no_results",
        help="How to define ghost protocols for modeling.",
    )
    parser.add_argument(
        "--ghost-require-no-any-pmid",
        action="store_true",
        default=True,
        help=(
            "Require no PMID at all when classifying ghost protocols. "
            "This avoids labeling trials with non-results PMID references as fully latent ghosts."
        ),
    )
    parser.add_argument(
        "--allow-any-pmid-ghost",
        action="store_false",
        dest="ghost_require_no_any_pmid",
        help="Allow ghost classification even if a non-linked PMID exists.",
    )
    parser.add_argument(
        "--intervention-alias",
        action="append",
        default=[],
        help="Additional alias for intervention matching; can be passed multiple times.",
    )
    parser.add_argument(
        "--comparator-scope",
        choices=["any", "active", "placebo"],
        default="active",
        help="Comparator filtering heuristic.",
    )
    parser.add_argument(
        "--publication-linkage",
        choices=["results_pmid_strict", "results_pmid", "any_pmid"],
        default="results_pmid_strict",
        help="Publication linkage rule used for has_pmid and ghost classification.",
    )
    parser.add_argument(
        "--require-randomized",
        action="store_true",
        default=True,
        help="Require randomized allocation.",
    )
    parser.add_argument(
        "--no-require-randomized",
        action="store_false",
        dest="require_randomized",
        help="Allow non-randomized studies.",
    )
    parser.add_argument(
        "--require-treatment-purpose",
        action="store_true",
        default=True,
        help="Require primary purpose == TREATMENT.",
    )
    parser.add_argument(
        "--no-require-treatment-purpose",
        action="store_false",
        dest="require_treatment_purpose",
        help="Allow non-treatment primary purpose.",
    )
    parser.add_argument(
        "--allowed-phase",
        action="append",
        default=[],
        help="Allowed trial phase (repeatable), e.g. PHASE3.",
    )
    parser.add_argument(
        "--required-outcome-keyword",
        action="append",
        default=[],
        help="Required primary-outcome keyword for estimand alignment (repeatable).",
    )
    parser.add_argument(
        "--outcome-timeframe-min-weeks",
        type=float,
        default=None,
        help="Minimum allowed primary-outcome timeframe in weeks.",
    )
    parser.add_argument(
        "--outcome-timeframe-max-weeks",
        type=float,
        default=None,
        help="Maximum allowed primary-outcome timeframe in weeks.",
    )
    parser.add_argument(
        "--missing-outcome-policy",
        choices=["exclude", "include_as_unknown"],
        default="include_as_unknown",
        help="How to handle trials with missing/non-parseable primary outcome metadata.",
    )
    parser.add_argument(
        "--extract-results",
        action="store_true",
        help=(
            "After building the registry, extract effect sizes from CT.gov posted "
            "results for results-only trials (has_results=True, no PMID). "
            "Adds results_effect, results_se, results_effect_type columns."
        ),
    )
    parser.add_argument(
        "--results-outcome-keywords",
        default="",
        help="Comma-separated keywords for outcome matching during results extraction.",
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=100,
        help="Maximum number of API pages to fetch (safety limit against runaway queries).",
    )
    parser.add_argument(
        "--include-terminated",
        action="store_true",
        help=(
            "Include TERMINATED and WITHDRAWN trials (classified as ghost protocols). "
            "NOTE: TERMINATED trials may have been stopped early for efficacy or harm, "
            "so the default ghost prior mu=0 may not be appropriate. Consider sensitivity "
            "analysis with non-zero ghost_mu when this flag is used."
        ),
    )
    return parser.parse_args()


def request_json_with_retry(
    session: requests.Session,
    *,
    url: str,
    params: dict[str, str | int],
    timeout: int,
    attempts: int = 5,
) -> dict:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as exc:
            last_error = exc
            if attempt == attempts:
                break
            sleep_s = min(12.0, 1.5 * (2 ** (attempt - 1)))
            time.sleep(sleep_s)
    if last_error is None:
        raise RuntimeError("request_json_with_retry: no attempts made")
    raise last_error


def main() -> int:
    args = parse_args()
    if args.max_pages < 1:
        raise ValueError("--max-pages must be >= 1.")
    if args.outcome_timeframe_min_weeks is not None and args.outcome_timeframe_max_weeks is not None:
        if args.outcome_timeframe_min_weeks > args.outcome_timeframe_max_weeks:
            raise ValueError("--outcome-timeframe-min-weeks cannot exceed --outcome-timeframe-max-weeks.")
    safe_intr = sanitize_path_component(args.intervention.strip().lower().replace(" ", "_"))
    safe_cond = sanitize_path_component(args.condition.strip().lower().replace(" ", "_"))
    safe_name = f"{safe_intr}__{safe_cond}.csv"
    output_path = args.output or Path("data/raw") / safe_name
    output_path.parent.mkdir(parents=True, exist_ok=True)

    statuses = ["COMPLETED"]
    if args.include_terminated:
        statuses.extend(["TERMINATED", "WITHDRAWN"])

    params: dict[str, str | int] = {
        "format": "json",
        "pageSize": args.page_size,
        "filter.overallStatus": ",".join(statuses),
        "query.intr": args.intervention,
        "query.cond": args.condition,
    }

    rows: list[dict[str, str | int | float]] = []
    seen_nct_ids: set[str] = set()
    next_page_token: str | None = None
    pages_fetched = 0
    with requests.Session() as session:
        session.headers["User-Agent"] = "GWAM-pipeline/1.0 (clinicaltrials.gov registry analysis)"
        while pages_fetched < args.max_pages:
            if next_page_token:
                params["pageToken"] = next_page_token
            else:
                params.pop("pageToken", None)

            payload = request_json_with_retry(session, url=BASE_URL, params=params, timeout=90)
            pages_fetched += 1
            studies = payload.get("studies", [])
            if not studies:
                break

            for study in studies:
                protocol = study.get("protocolSection", {})
                ident = protocol.get("identificationModule", {})
                status = protocol.get("statusModule", {})
                sponsor = protocol.get("sponsorCollaboratorsModule", {})
                design = protocol.get("designModule", {})
                refs = protocol.get("referencesModule", {}).get("references", [])

                if not trial_passes_design_filters(
                    protocol,
                    require_randomized=args.require_randomized,
                    require_treatment_purpose=args.require_treatment_purpose,
                    comparator_scope=args.comparator_scope,
                    allowed_phases=args.allowed_phase,
                    required_outcome_keywords=args.required_outcome_keyword,
                    outcome_timeframe_min_weeks=args.outcome_timeframe_min_weeks,
                    outcome_timeframe_max_weeks=args.outcome_timeframe_max_weeks,
                    missing_outcome_policy=args.missing_outcome_policy,
                ):
                    continue

                # Deduplicate NCT IDs across API pages
                nct_id = ident.get("nctId", "")
                if nct_id in seen_nct_ids:
                    continue
                seen_nct_ids.add(nct_id)

                interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                intervention_names = [
                    str(item.get("name", "")).strip().lower() for item in interventions if item.get("name")
                ]
                if not intervention_matches(
                    intervention_names,
                    intervention=args.intervention,
                    extra_aliases=args.intervention_alias,
                ):
                    continue

                pmids_any, pmids_results, has_pmid = extract_pmids(
                    refs,
                    publication_linkage=args.publication_linkage,
                )
                has_any_pmid = bool(pmids_any)
                enrollment = design.get("enrollmentInfo", {})
                enrollment_count = parse_int(enrollment.get("count"))
                phases = design.get("phases", [])
                has_results = bool(study.get("hasResults"))
                is_unpublished_with_results, is_ghost = classify_publication_status(
                    has_pmid=has_pmid,
                    has_results=has_results,
                    has_any_pmid=has_any_pmid,
                    ghost_require_no_any_pmid=args.ghost_require_no_any_pmid,
                    ghost_definition=args.ghost_definition,
                )

                rows.append(
                    {
                        "nct_id": ident.get("nctId", ""),
                        "brief_title": ident.get("briefTitle", ""),
                        "overall_status": status.get("overallStatus", ""),
                        "completion_date": status.get("completionDateStruct", {}).get("date", ""),
                        "study_first_post_date": status.get("studyFirstPostDateStruct", {}).get("date", ""),
                        "phase": "|".join(phases) if phases else "",
                        "has_results": has_results,
                        "enrollment_count": enrollment_count if enrollment_count is not None else "",
                        "enrollment_type": enrollment.get("type", ""),
                        "lead_sponsor_name": sponsor.get("leadSponsor", {}).get("name", ""),
                        "lead_sponsor_class": sponsor.get("leadSponsor", {}).get("class", ""),
                        "pmid_count_any": len(pmids_any),
                        "pmids_any": ";".join(pmids_any),
                        "pmid_count_results": len(pmids_results),
                        "pmids_results": ";".join(pmids_results),
                        "publication_linkage_rule": args.publication_linkage,
                        "has_any_pmid": has_any_pmid,
                        "has_pmid": has_pmid,
                        "is_unlinked_any_pmid": (not has_pmid) and has_any_pmid,
                        "is_unpublished_with_results": is_unpublished_with_results,
                        "is_ghost_protocol": is_ghost,
                    }
                )

            next_page_token = payload.get("nextPageToken")
            if not next_page_token:
                break

    # Impute missing enrollment with median for weighting workflows.
    non_missing_n = sorted(
        [int(r["enrollment_count"]) for r in rows if str(r["enrollment_count"]).strip() != ""]
    )
    if non_missing_n:
        mid = len(non_missing_n) // 2
        if len(non_missing_n) % 2 == 0:
            median_n = (non_missing_n[mid - 1] + non_missing_n[mid]) // 2
        else:
            median_n = non_missing_n[mid]
    else:
        median_n = 100

    query_date = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for row in rows:
        n_str = str(row["enrollment_count"]).strip()
        weight_n = int(n_str) if n_str else median_n
        row["weight_n"] = weight_n
        row["weight_n_imputation"] = "none" if n_str else "median"
        row["query_date_utc"] = query_date

    fieldnames = [
        "nct_id",
        "brief_title",
        "overall_status",
        "completion_date",
        "study_first_post_date",
        "phase",
        "has_results",
        "enrollment_count",
        "enrollment_type",
        "lead_sponsor_name",
        "lead_sponsor_class",
        "pmid_count_any",
        "pmids_any",
        "pmid_count_results",
        "pmids_results",
        "publication_linkage_rule",
        "has_any_pmid",
        "has_pmid",
        "is_unlinked_any_pmid",
        "is_unpublished_with_results",
        "is_ghost_protocol",
        "weight_n",
        "weight_n_imputation",
        "query_date_utc",
    ]
    # Sanitize string cells against spreadsheet formula injection (P1-5)
    for row in rows:
        for key, val in row.items():
            if isinstance(val, str):
                row[key] = sanitize_csv_cell(val)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    with_pmid = sum(1 for row in rows if row["has_pmid"])
    unpublished_with_results = sum(1 for row in rows if row["is_unpublished_with_results"])
    ghost = sum(1 for row in rows if row["is_ghost_protocol"])
    print(f"Wrote {output_path}")
    print(
        f"completed={total}, with_pmid={with_pmid}, ghost={ghost}, "
        f"unpublished_with_results={unpublished_with_results}, "
        f"ghost_ratio={ghost / total if total else 0:.3f}, "
        f"median_weight_n={median_n}, ghost_definition={args.ghost_definition}, "
        f"ghost_require_no_any_pmid={args.ghost_require_no_any_pmid}, "
        f"publication_linkage={args.publication_linkage}, "
        f"comparator_scope={args.comparator_scope}, "
        f"require_randomized={args.require_randomized}, "
        f"require_treatment_purpose={args.require_treatment_purpose}, "
        f"allowed_phase={args.allowed_phase or 'ANY'}, "
        f"required_outcome_keyword={args.required_outcome_keyword or 'ANY'}, "
        f"outcome_timeframe_min_weeks={args.outcome_timeframe_min_weeks}, "
        f"outcome_timeframe_max_weeks={args.outcome_timeframe_max_weeks}, "
        f"missing_outcome_policy={args.missing_outcome_policy}"
    )

    # Phase 2: optionally extract effect sizes from CT.gov posted results
    if args.extract_results:
        from extract_ctgov_results import (
            extract_outcome_effect,
            fetch_results_section,
        )

        outcome_keywords = [
            kw.strip()
            for kw in args.results_outcome_keywords.split(",")
            if kw.strip()
        ]

        # Re-read the CSV we just wrote so we can add new columns
        with output_path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            enriched_fieldnames = list(reader.fieldnames or [])
            enriched_rows = list(reader)

        new_cols = ["results_effect", "results_se", "results_effect_type"]
        for col in new_cols:
            if col not in enriched_fieldnames:
                enriched_fieldnames.append(col)

        n_results_attempted = 0
        n_results_extracted = 0
        for erow in enriched_rows:
            erow.setdefault("results_effect", "")
            erow.setdefault("results_se", "")
            erow.setdefault("results_effect_type", "")

            is_results_only = (
                parse_bool(erow.get("has_results", ""))
                and not parse_bool(erow.get("has_pmid", ""))
                and not parse_bool(erow.get("is_ghost_protocol", ""))
            )
            if not is_results_only:
                continue

            nct_id = str(erow.get("nct_id", "")).strip()
            if not nct_id:
                continue

            n_results_attempted += 1
            result_payload = fetch_results_section(nct_id, timeout_sec=30.0)
            if result_payload is None:
                continue

            effect_result = extract_outcome_effect(result_payload, outcome_keywords)
            if effect_result is not None:
                effect, se, effect_type = effect_result
                erow["results_effect"] = str(effect)
                erow["results_se"] = str(se)
                erow["results_effect_type"] = effect_type
                n_results_extracted += 1
                print(f"  Results: {nct_id}: {effect_type} = {effect:.4f} (SE={se:.4f})")

            time.sleep(0.3)  # rate limit

        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=enriched_fieldnames)
            writer.writeheader()
            writer.writerows(enriched_rows)

        print(
            f"Results extraction: attempted={n_results_attempted}, "
            f"extracted={n_results_extracted}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
