#!/usr/bin/env python3
"""Scan candidate intervention/condition pairs on ClinicalTrials.gov.

The script counts completed studies under analysis filters and reports:
- linked-PMID studies
- unpublished-with-results studies
- strict ghost protocols
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from gwam_utils import (
    classify_publication_status,
    extract_pmids,
    intervention_matches,
    sanitize_csv_cell,
    trial_passes_design_filters,
)

BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

DEFAULT_CANDIDATES: list[tuple[str, str]] = [
    ("sertraline", "depression"),
    ("escitalopram", "depression"),
    ("duloxetine", "depression"),
    ("quetiapine", "depression"),
    ("pregabalin", "neuropathic pain"),
    ("oseltamivir", "influenza"),
]


def parse_candidate_token(token: str) -> tuple[str, str]:
    if "::" not in token:
        raise ValueError(f"Invalid --candidate value '{token}'. Use format: intervention::condition")
    intervention, condition = token.split("::", 1)
    intervention = intervention.strip()
    condition = condition.strip()
    if not intervention or not condition:
        raise ValueError(f"Invalid --candidate value '{token}'. Both intervention and condition are required.")
    return intervention, condition


@dataclass
class CandidateStats:
    intervention: str
    condition: str
    total_completed: int
    with_pmid: int
    unpublished_with_results: int
    ghost_protocols: int
    ghost_ratio: float
    with_results: int


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


def fetch_candidate_stats(
    session: requests.Session,
    intervention: str,
    condition: str,
    *,
    page_size: int = 100,
    timeout: int = 60,
    intervention_aliases: list[str] | None = None,
    comparator_scope: str = "active",
    require_randomized: bool = True,
    require_treatment_purpose: bool = True,
    ghost_definition: str = "no_pmid_no_results",
    ghost_require_no_any_pmid: bool = True,
    publication_linkage: str = "results_pmid",
    allowed_phases: list[str] | None = None,
    required_outcome_keywords: list[str] | None = None,
    outcome_timeframe_min_weeks: float | None = None,
    outcome_timeframe_max_weeks: float | None = None,
    missing_outcome_policy: str = "include_as_unknown",
    max_pages: int = 100,
) -> CandidateStats:
    params: dict[str, str | int] = {
        "format": "json",
        "pageSize": page_size,
        "filter.overallStatus": "COMPLETED",
        "query.intr": intervention,
        "query.cond": condition,
    }

    total_completed = 0
    with_pmid = 0
    unpublished_with_results = 0
    ghost_protocols = 0
    with_results = 0
    next_page_token: str | None = None
    pages_fetched = 0

    while True:
        if next_page_token:
            params["pageToken"] = next_page_token
        else:
            params.pop("pageToken", None)

        pages_fetched += 1
        if pages_fetched > max_pages:
            break
        payload = request_json_with_retry(session, url=BASE_URL, params=params, timeout=timeout)
        studies = payload.get("studies", [])
        if not studies:
            break

        for study in studies:
            protocol = study.get("protocolSection", {})
            if not trial_passes_design_filters(
                protocol,
                require_randomized=require_randomized,
                require_treatment_purpose=require_treatment_purpose,
                comparator_scope=comparator_scope,
                allowed_phases=allowed_phases,
                required_outcome_keywords=required_outcome_keywords,
                outcome_timeframe_min_weeks=outcome_timeframe_min_weeks,
                outcome_timeframe_max_weeks=outcome_timeframe_max_weeks,
                missing_outcome_policy=missing_outcome_policy,
            ):
                continue

            interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
            intervention_names = [
                str(item.get("name", "")).strip().lower() for item in interventions if item.get("name")
            ]
            if not intervention_matches(
                intervention_names,
                intervention=intervention,
                extra_aliases=intervention_aliases,
            ):
                continue

            total_completed += 1
            refs = protocol.get("referencesModule", {}).get("references", [])
            pmids_any, _, has_pmid = extract_pmids(refs, publication_linkage=publication_linkage)
            has_any_pmid = bool(pmids_any)
            has_results = bool(study.get("hasResults"))
            if has_pmid:
                with_pmid += 1
            if has_results:
                with_results += 1
            is_unpublished_with_results, is_ghost = classify_publication_status(
                has_pmid=has_pmid,
                has_results=has_results,
                has_any_pmid=has_any_pmid,
                ghost_require_no_any_pmid=ghost_require_no_any_pmid,
                ghost_definition=ghost_definition,
            )
            if is_unpublished_with_results:
                unpublished_with_results += 1
            if is_ghost:
                ghost_protocols += 1

        next_page_token = payload.get("nextPageToken")
        if not next_page_token:
            break

    ghost_ratio = (ghost_protocols / total_completed) if total_completed else 0.0
    return CandidateStats(
        intervention=intervention,
        condition=condition,
        total_completed=total_completed,
        with_pmid=with_pmid,
        unpublished_with_results=unpublished_with_results,
        ghost_protocols=ghost_protocols,
        ghost_ratio=ghost_ratio,
        with_results=with_results,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/registry_candidate_summary.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--ghost-definition",
        choices=["no_pmid", "no_pmid_no_results"],
        default="no_pmid_no_results",
        help="How to define ghost protocols in summary.",
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
        "--max-pages",
        type=int,
        default=100,
        help="Maximum number of API pages to fetch per candidate (safety limit).",
    )
    parser.add_argument(
        "--candidate",
        action="append",
        default=[],
        help="Optional candidate pair in format intervention::condition. Repeatable.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.outcome_timeframe_min_weeks is not None and args.outcome_timeframe_max_weeks is not None:
        if args.outcome_timeframe_min_weeks > args.outcome_timeframe_max_weeks:
            raise ValueError("--outcome-timeframe-min-weeks cannot exceed --outcome-timeframe-max-weeks.")

    candidates = DEFAULT_CANDIDATES
    if args.candidate:
        candidates = [parse_candidate_token(token) for token in args.candidate]

    run_date = dt.date.today().isoformat()
    rows: list[CandidateStats] = []
    with requests.Session() as session:
        for intervention, condition in candidates:
            rows.append(
                fetch_candidate_stats(
                    session,
                    intervention,
                    condition,
                    ghost_definition=args.ghost_definition,
                    ghost_require_no_any_pmid=args.ghost_require_no_any_pmid,
                    comparator_scope=args.comparator_scope,
                    require_randomized=args.require_randomized,
                    require_treatment_purpose=args.require_treatment_purpose,
                    publication_linkage=args.publication_linkage,
                    allowed_phases=args.allowed_phase,
                    required_outcome_keywords=args.required_outcome_keyword,
                    outcome_timeframe_min_weeks=args.outcome_timeframe_min_weeks,
                    outcome_timeframe_max_weeks=args.outcome_timeframe_max_weeks,
                    missing_outcome_policy=args.missing_outcome_policy,
                    max_pages=args.max_pages,
                )
            )

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "run_date",
                "intervention",
                "condition",
                "total_completed",
                "with_pmid",
                "unpublished_with_results",
                "ghost_protocols",
                "ghost_ratio",
                "with_results",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    run_date,
                    sanitize_csv_cell(row.intervention),
                    sanitize_csv_cell(row.condition),
                    row.total_completed,
                    row.with_pmid,
                    row.unpublished_with_results,
                    row.ghost_protocols,
                    f"{row.ghost_ratio:.6f}",
                    row.with_results,
                ]
            )

    print(f"Wrote {args.output}")
    print(
        f"filters: ghost_definition={args.ghost_definition}, "
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
    for row in rows:
        print(
            f"{row.intervention}/{row.condition}: "
            f"completed={row.total_completed}, with_pmid={row.with_pmid}, "
            f"unpublished_with_results={row.unpublished_with_results}, "
            f"ghost={row.ghost_protocols} ({row.ghost_ratio:.3f}), "
            f"has_results={row.with_results}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
