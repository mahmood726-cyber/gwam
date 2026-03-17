#!/usr/bin/env python3
"""Build review-level CT.gov linkage summaries for Pairwise70.

Two linkage layers are computed:
1) Query-level CT.gov linkage from review query terms.
2) Exact NCT-level linkage for reviews where Pairwise70 data exposes NCT IDs.

For each review we compute:
- completed-study counts
- completed ghost counts (completed, no PMID, no posted results)
- weighted lambda proxies from completed studies:
  - lambda_completed_pmid_any_weighted
  - lambda_completed_non_ghost_weighted
  - lambda_exact_nct_pmid_any_weighted
  - lambda_exact_nct_non_ghost_weighted
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import time
import defusedxml.ElementTree as ET
from pathlib import Path
from typing import Any

import pyreadr
import pandas as pd
import requests

from gwam_utils import parse_bool, safe_float, sanitize_csv_cell

API_BASE = "https://clinicaltrials.gov/api/v2/studies"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCT_PATTERN = re.compile(r"\bNCT\d{8}\b", re.IGNORECASE)
REVIEW_FILE_PATTERN = re.compile(r"^(CD\d+)_pub\d+_data\.rda$", re.IGNORECASE)
NCT_FIELDS = (
    "protocolSection.identificationModule.nctId,"
    "protocolSection.statusModule.overallStatus,"
    "protocolSection.designModule.enrollmentInfo.count,"
    "protocolSection.referencesModule.references,"
    "hasResults"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairwise-data-dir",
        type=Path,
        required=True,
        help="Directory containing Pairwise70 .rda files for exact NCT extraction.",
    )
    parser.add_argument(
        "--query-csv",
        type=Path,
        required=True,
        help="Pairwise70 review query CSV with review_id and query_term.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv"),
        help="Output CSV path for linkage summary.",
    )
    parser.add_argument(
        "--disable-exact-nct-linkage",
        action="store_true",
        help="Disable exact NCT-level linkage and use query linkage only.",
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=3,
        help="Max CT.gov pages per review query.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=100,
        help="CT.gov page size.",
    )
    parser.add_argument(
        "--sleep-ms",
        type=int,
        default=80,
        help="Sleep between API calls (milliseconds).",
    )
    parser.add_argument(
        "--max-reviews",
        type=int,
        default=None,
        help="Optional cap for number of review queries processed.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("pairwise70_benchmark_sensitivity/ctgov_linkage_raw"),
        help="Directory to cache per-review fetched JSON payloads.",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Refetch even if cache file exists.",
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=30.0,
        help="HTTP timeout per request.",
    )
    parser.add_argument(
        "--enable-pubmed-bridge",
        action="store_true",
        help="Enable PubMed PMID-to-NCT bridge linkage via E-utilities efetch.",
    )
    parser.add_argument(
        "--enable-fuzzy-match",
        action="store_true",
        help="Enable fuzzy matching of studies to CT.gov by condition+intervention+enrollment.",
    )
    parser.add_argument(
        "--pubmed-email",
        default="",
        help="Contact email for NCBI E-utilities requests.",
    )
    parser.add_argument(
        "--pubmed-api-key",
        default=__import__("os").environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (or set NCBI_API_KEY env var; increases rate limit from 3/s to 10/s).",
    )
    parser.add_argument(
        "--fuzzy-enrollment-tolerance",
        type=float,
        default=0.2,
        help="Tolerance for enrollment matching in fuzzy linkage (fraction, default 0.2 = 20%%).",
    )
    return parser.parse_args()


def normalize_studies_shape(studies_obj: Any) -> list[dict[str, Any]]:
    if isinstance(studies_obj, list):
        return [x for x in studies_obj if isinstance(x, dict)]
    if isinstance(studies_obj, dict):
        # Could be flattened columnar dict: protocolSection, protocolSection.1, ...
        keys = [k for k in studies_obj.keys() if k.startswith("protocolSection")]
        if not keys:
            return []
        rows: list[dict[str, Any]] = []
        for key in keys:
            block = studies_obj.get(key, [])
            if isinstance(block, list):
                for item in block:
                    if isinstance(item, dict):
                        rows.append({"protocolSection": item})
        return rows
    return []


def normalize_nct_id(value: Any) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    match = NCT_PATTERN.search(text)
    return match.group(0).upper() if match else ""


def extract_nct_ids_from_text(text: Any) -> set[str]:
    if text is None:
        return set()
    return {token.upper() for token in NCT_PATTERN.findall(str(text))}


def parse_review_id_from_filename(name: str) -> str:
    match = REVIEW_FILE_PATTERN.match(name)
    return match.group(1).upper() if match else ""


def extract_nct_ids_from_frame(frame: pd.DataFrame) -> set[str]:
    ids: set[str] = set()
    if frame.empty:
        return ids
    for col in frame.select_dtypes(include=["object"]).columns:
        series = frame[col].dropna().astype(str)
        if series.empty:
            continue
        for matches in series.str.findall(NCT_PATTERN):
            if not matches:
                continue
            for token in matches:
                ids.add(token.upper())
    return ids


def build_review_nct_map(pairwise_data_dir: Path, target_review_ids: set[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    files = sorted(pairwise_data_dir.glob("*.rda"))
    for file_path in files:
        review_id = parse_review_id_from_filename(file_path.name)
        if not review_id:
            continue
        if target_review_ids and review_id not in target_review_ids:
            continue
        try:
            payload = pyreadr.read_r(str(file_path))
        except (OSError, ValueError, KeyError) as exc:
            print(f"WARNING: Failed to read {file_path.name}: {type(exc).__name__}: {exc}")
            continue
        file_ids: set[str] = set()
        for _, frame in payload.items():
            if isinstance(frame, pd.DataFrame):
                file_ids.update(extract_nct_ids_from_frame(frame))
        if not file_ids:
            continue
        bucket = out.setdefault(review_id, set())
        bucket.update(file_ids)
    return out


def pubmed_to_nct(
    pmids: list[str],
    *,
    email: str = "",
    api_key: str = "",
    timeout_sec: float = 30.0,
    batch_size: int = 200,
) -> dict[str, list[str]]:
    """Fetch PubMed records and extract NCT IDs from DataBankList / abstract text.

    Returns a dict mapping PMID -> list of NCT IDs found in that record.
    """
    if not pmids:
        return {}

    result: dict[str, list[str]] = {}
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start : start + batch_size]
        params: dict[str, str] = {
            "db": "pubmed",
            "retmode": "xml",
            "id": ",".join(batch),
        }
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key

        try:
            resp = requests.get(EFETCH_URL, params=params, timeout=timeout_sec)
            resp.raise_for_status()
            # Reject suspiciously large responses to mitigate entity expansion
            if len(resp.text) > 50_000_000:
                print(f"  PubMed bridge: skipping oversized response ({len(resp.text)} bytes)")
                continue
            # defusedxml.fromstring handles entity expansion protection.
            root = ET.fromstring(resp.text.encode("utf-8"))
        except Exception as exc:
            print(f"  PubMed bridge: batch starting {batch[0]} failed: {type(exc).__name__}")
            continue

        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//PMID", default="").strip()
            if not pmid:
                continue

            nct_ids: set[str] = set()

            # Method 1: DataBankList accession numbers
            for acc in article.findall(".//DataBankList/DataBank/AccessionNumberList/AccessionNumber"):
                text = (acc.text or "").strip().upper()
                if NCT_PATTERN.match(text):
                    nct_ids.add(text)

            # Method 2: Search abstract text for NCT patterns
            for abstract_text in article.findall(".//AbstractText"):
                if abstract_text.text:
                    for match in NCT_PATTERN.findall(abstract_text.text):
                        nct_ids.add(match.upper())

            if nct_ids:
                result[pmid] = sorted(nct_ids)

        # Rate limit: 3/sec without key, 10/sec with key
        sleep_sec = 0.1 if api_key else 0.35
        time.sleep(sleep_sec)

    return result


def fuzzy_match_ctgov(
    *,
    condition: str,
    intervention: str,
    enrollment_count: int,
    completion_year: int | None = None,
    tolerance: float = 0.2,
    timeout_sec: float = 30.0,
) -> str | None:
    """Search CT.gov for a unique match by condition+intervention+enrollment.

    Returns NCT ID if exactly one match found within enrollment tolerance, else None.
    Conservative: requires unique hit with enrollment within tolerance.
    """
    if not condition or not intervention or enrollment_count <= 0:
        return None

    query_term = f"{condition} AND {intervention}"
    params: dict[str, Any] = {
        "query.term": query_term,
        "filter.overallStatus": "COMPLETED",
        "pageSize": 10,
        "fields": (
            "protocolSection.identificationModule.nctId,"
            "protocolSection.designModule.enrollmentInfo.count,"
            "protocolSection.statusModule.completionDateStruct"
        ),
    }

    try:
        resp = requests.get(API_BASE, params=params, timeout=timeout_sec)
        resp.raise_for_status()
        payload = resp.json()
    except (requests.RequestException, ValueError, KeyError) as exc:
        print(f"WARNING: CT.gov fuzzy match failed: {type(exc).__name__}: {exc}")
        return None

    studies = payload.get("studies", [])
    if not studies:
        return None

    lo = enrollment_count * (1.0 - tolerance)
    hi = enrollment_count * (1.0 + tolerance)

    matching_ncts: list[str] = []
    for study in studies:
        ps = study.get("protocolSection", {})
        nct = (ps.get("identificationModule", {}) or {}).get("nctId", "").strip().upper()
        enroll_raw = (
            ((ps.get("designModule", {}) or {}).get("enrollmentInfo", {}) or {}).get("count")
        )
        try:
            enroll = float(enroll_raw)
        except (TypeError, ValueError):
            continue

        if not (lo <= enroll <= hi):
            continue

        # Optional: filter by completion year
        if completion_year is not None:
            comp_date = (
                ((ps.get("statusModule", {}) or {}).get("completionDateStruct", {}) or {}).get("date", "")
            )
            if comp_date:
                try:
                    comp_year = int(str(comp_date)[:4])
                    if comp_year != completion_year:
                        continue
                except (TypeError, ValueError):
                    pass

        if nct:
            matching_ncts.append(nct)

    # Conservative: require exactly one unique match
    unique_ncts = list(dict.fromkeys(matching_ncts))
    if len(unique_ncts) == 1:
        return unique_ncts[0]
    return None


def study_nct_id(study: dict[str, Any]) -> str:
    if not isinstance(study, dict):
        return ""
    ps = study.get("protocolSection", {})
    if not isinstance(ps, dict):
        return ""
    ident = ps.get("identificationModule", {})
    if not isinstance(ident, dict):
        return ""
    return normalize_nct_id(ident.get("nctId", ""))


_NCT_FORMAT_STRICT = re.compile(r"^NCT\d{8}$", re.IGNORECASE)


def fetch_study_by_nct(nct_id: str, timeout_sec: float) -> dict[str, Any] | None:
    if not nct_id or not _NCT_FORMAT_STRICT.match(nct_id):
        return None
    resp = requests.get(
        f"{API_BASE}/{nct_id}",
        params={"fields": NCT_FIELDS},
        timeout=timeout_sec,
    )
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    payload = resp.json()
    if not isinstance(payload, dict):
        return None
    return {
        "protocolSection": payload.get("protocolSection", {}),
        "hasResults": payload.get("hasResults", False),
    }


def resolve_exact_nct_studies(
    *,
    target_nct_ids: set[str],
    query_payload: dict[str, Any],
    nct_cache: dict[str, dict[str, Any] | None],
    sleep_ms: int,
    timeout_sec: float,
) -> tuple[list[dict[str, Any]], list[str]]:
    query_studies = normalize_studies_shape(query_payload.get("studies", []))
    by_nct: dict[str, dict[str, Any]] = {}
    for study in query_studies:
        nct = study_nct_id(study)
        if nct and nct not in by_nct:
            by_nct[nct] = study

    resolved: list[dict[str, Any]] = []
    unresolved: list[str] = []
    for nct in sorted(target_nct_ids):
        study = by_nct.get(nct)
        if study is None:
            if nct not in nct_cache:
                try:
                    nct_cache[nct] = fetch_study_by_nct(nct, timeout_sec=timeout_sec)
                except (requests.RequestException, ValueError, KeyError) as exc:
                    print(f"WARNING: Failed to fetch {nct}: {type(exc).__name__}: {exc}")
                    nct_cache[nct] = None
                if sleep_ms > 0:
                    time.sleep(sleep_ms / 1000.0)
            study = nct_cache.get(nct)
        if isinstance(study, dict):
            resolved.append(study)
        else:
            unresolved.append(nct)
    return resolved, unresolved


def fetch_query_payload(
    *,
    query_term: str,
    max_pages: int,
    page_size: int,
    sleep_ms: int,
    timeout_sec: float,
) -> dict[str, Any]:
    all_studies: list[dict[str, Any]] = []
    token: str | None = None
    fields = NCT_FIELDS
    for _ in range(max_pages):
        params: dict[str, Any] = {
            "query.term": query_term,
            "pageSize": page_size,
            "fields": fields,
        }
        if token:
            params["pageToken"] = token
        try:
            resp = requests.get(API_BASE, params=params, timeout=timeout_sec)
            resp.raise_for_status()
            payload = resp.json()
        except (requests.RequestException, ValueError) as exc:
            print(f"WARNING: CT.gov page fetch failed for '{query_term}': {type(exc).__name__}: {exc}")
            break  # return partial results from earlier pages
        studies = normalize_studies_shape(payload.get("studies", []))
        all_studies.extend(studies)
        token = payload.get("nextPageToken")
        if not token:
            break
        if sleep_ms > 0:
            time.sleep(sleep_ms / 1000.0)
    return {"query": query_term, "studies": all_studies}


def has_any_pmid(references_module: Any, ref_type_filter: str | None = "RESULT") -> bool:
    """Check if references module contains at least one PMID.

    Parameters
    ----------
    ref_type_filter : if set, only count references whose ``type`` matches
        (case-insensitive).  Default ``"RESULT"`` excludes BACKGROUND /
        DERIVED references that don't indicate the study published its own
        results.  Pass ``None`` to count any reference type.
    """
    if not isinstance(references_module, dict):
        return False
    refs = references_module.get("references", [])
    if not isinstance(refs, list):
        return False
    for ref in refs:
        if isinstance(ref, dict):
            if ref_type_filter is not None:
                rtype = str(ref.get("type", "")).strip().upper()
                if rtype != ref_type_filter.upper():
                    continue
            pmid = str(ref.get("pmid", "")).strip()
            if pmid:
                return True
    return False


def summarize_payload(review_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    studies = normalize_studies_shape(payload.get("studies", []))
    n_studies = len(studies)
    n_completed = 0
    n_completed_with_pmid = 0
    n_completed_with_results = 0
    n_completed_non_ghost = 0
    n_completed_ghost = 0
    n_completed_results_only = 0

    w_completed_total = 0.0
    w_completed_pmid = 0.0
    w_completed_non_ghost = 0.0
    w_completed_ghost = 0.0

    for study in studies:
        ps = study.get("protocolSection", {}) if isinstance(study, dict) else {}
        sm = ps.get("statusModule", {}) if isinstance(ps, dict) else {}
        status = str(sm.get("overallStatus", "")).strip().upper()
        completed = status == "COMPLETED"
        if not completed:
            continue
        n_completed += 1

        has_results = parse_bool(study.get("hasResults", False))
        has_pmid = has_any_pmid(ps.get("referencesModule", {}))
        enrollment = safe_float(
            (((ps.get("designModule", {}) or {}).get("enrollmentInfo", {}) or {}).get("count"))
        )
        weight = enrollment if (enrollment > 0) else 1.0

        w_completed_total += weight
        if has_pmid:
            n_completed_with_pmid += 1
            w_completed_pmid += weight
        if has_results:
            n_completed_with_results += 1
        if has_results or has_pmid:
            n_completed_non_ghost += 1
            w_completed_non_ghost += weight
        else:
            n_completed_ghost += 1
            w_completed_ghost += weight
        if has_results and (not has_pmid):
            n_completed_results_only += 1

    if w_completed_total > 0:
        lambda_pmid = w_completed_pmid / w_completed_total
        lambda_non_ghost = w_completed_non_ghost / w_completed_total
        ghost_weight_share = w_completed_ghost / w_completed_total
    else:
        lambda_pmid = float("nan")
        lambda_non_ghost = float("nan")
        ghost_weight_share = float("nan")

    return {
        "review_id": review_id,
        "query_term": payload.get("query", ""),
        "n_studies_fetched": n_studies,
        "n_completed": n_completed,
        "n_completed_with_pmid_any": n_completed_with_pmid,
        "n_completed_with_results": n_completed_with_results,
        "n_completed_non_ghost": n_completed_non_ghost,
        "n_completed_ghost": n_completed_ghost,
        "n_completed_results_only_no_pmid": n_completed_results_only,
        "w_completed_total": w_completed_total,
        "w_completed_pmid_any": w_completed_pmid,
        "w_completed_non_ghost": w_completed_non_ghost,
        "w_completed_ghost": w_completed_ghost,
        "lambda_completed_pmid_any_weighted": lambda_pmid,
        "lambda_completed_non_ghost_weighted": lambda_non_ghost,
        "ghost_weight_share_completed": ghost_weight_share,
    }


def main() -> int:
    args = parse_args()
    args.enable_exact_nct_linkage = not args.disable_exact_nct_linkage
    if args.max_pages < 1:
        raise ValueError("--max-pages must be >= 1.")
    if args.page_size < 1:
        raise ValueError("--page-size must be >= 1.")
    if args.sleep_ms < 0:
        raise ValueError("--sleep-ms must be >= 0.")

    with args.query_csv.open("r", newline="", encoding="utf-8") as _qf:
        rows = list(csv.DictReader(_qf))
    if not rows:
        raise ValueError(f"No rows in {args.query_csv}")
    if args.max_reviews is not None:
        rows = rows[: max(0, args.max_reviews)]

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    target_review_ids = {str(row.get("review_id", "")).strip().upper() for row in rows}
    if args.enable_exact_nct_linkage:
        review_nct_map = build_review_nct_map(args.pairwise_data_dir, target_review_ids)
    else:
        review_nct_map = {}
    n_reviews_with_exact_ids = sum(1 for _, ids in review_nct_map.items() if ids)
    n_total_exact_ids = sum(len(ids) for ids in review_nct_map.values())
    print(
        f"Exact NCT linkage map: reviews_with_ids={n_reviews_with_exact_ids}, "
        f"total_unique_review_nct_ids={n_total_exact_ids}"
    )

    summaries: list[dict[str, Any]] = []
    nct_cache: dict[str, dict[str, Any] | None] = {}
    for idx, row in enumerate(rows, start=1):
        review_id = str(row.get("review_id", "")).strip().upper()
        query_term = str(row.get("query_term", "")).strip()
        if not review_id or not query_term:
            continue

        # Sanitise review_id to prevent path traversal
        safe_id = re.sub(r"[^A-Za-z0-9_-]", "", review_id)
        if not safe_id:
            continue
        cache_path = (args.cache_dir / f"{safe_id}.json").resolve()
        if not str(cache_path).startswith(str(args.cache_dir.resolve())):
            print(f"  [WARN] Skipping review_id with suspicious path: {review_id}")
            continue
        payload: dict[str, Any]
        if cache_path.exists() and (not args.refresh_cache):
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
        else:
            try:
                payload = fetch_query_payload(
                    query_term=query_term,
                    max_pages=args.max_pages,
                    page_size=args.page_size,
                    sleep_ms=args.sleep_ms,
                    timeout_sec=args.timeout_sec,
                )
            except Exception as exc:
                # Sanitize error: remove any API key that may appear in URL params
                error_msg = str(exc)
                if args.pubmed_api_key:
                    error_msg = error_msg.replace(args.pubmed_api_key, "[REDACTED]")
                payload = {"query": query_term, "studies": [], "error": error_msg}
            cache_path.write_text(json.dumps(payload), encoding="utf-8")

        summary = summarize_payload(review_id, payload)
        pairwise_nct_ids = review_nct_map.get(review_id, set())
        summary["n_pairwise_nct_ids"] = len(pairwise_nct_ids)
        summary["n_exact_nct_studies_resolved"] = 0
        summary["n_exact_nct_ids_unresolved"] = 0
        summary["n_exact_nct_completed"] = 0
        summary["n_exact_nct_completed_non_ghost"] = 0
        summary["n_exact_nct_completed_ghost"] = 0
        summary["lambda_exact_nct_pmid_any_weighted"] = float("nan")
        summary["lambda_exact_nct_non_ghost_weighted"] = float("nan")

        if args.enable_exact_nct_linkage and pairwise_nct_ids:
            resolved_studies, unresolved_ids = resolve_exact_nct_studies(
                target_nct_ids=pairwise_nct_ids,
                query_payload=payload,
                nct_cache=nct_cache,
                sleep_ms=args.sleep_ms,
                timeout_sec=args.timeout_sec,
            )
            exact_summary = summarize_payload(review_id, {"query": query_term, "studies": resolved_studies})
            summary["n_exact_nct_studies_resolved"] = len(resolved_studies)
            summary["n_exact_nct_ids_unresolved"] = len(unresolved_ids)
            summary["n_exact_nct_completed"] = exact_summary["n_completed"]
            summary["n_exact_nct_completed_non_ghost"] = exact_summary["n_completed_non_ghost"]
            summary["n_exact_nct_completed_ghost"] = exact_summary["n_completed_ghost"]
            summary["lambda_exact_nct_pmid_any_weighted"] = exact_summary[
                "lambda_completed_pmid_any_weighted"
            ]
            summary["lambda_exact_nct_non_ghost_weighted"] = exact_summary[
                "lambda_completed_non_ghost_weighted"
            ]

        # --- PubMed bridge linkage ---
        summary["n_pubmed_bridge_nct_ids"] = 0
        summary["linkage_source"] = "query_only"

        if args.enable_pubmed_bridge:
            # Collect all PMIDs from the query-level studies
            query_studies = normalize_studies_shape(payload.get("studies", []))
            all_pmids_for_review: list[str] = []
            for study in query_studies:
                ps = study.get("protocolSection", {}) if isinstance(study, dict) else {}
                refs_mod = ps.get("referencesModule", {})
                if isinstance(refs_mod, dict):
                    for ref in refs_mod.get("references", []):
                        if isinstance(ref, dict):
                            pmid = str(ref.get("pmid", "")).strip()
                            if pmid and pmid not in all_pmids_for_review:
                                all_pmids_for_review.append(pmid)

            if all_pmids_for_review:
                pmid_nct_map = pubmed_to_nct(
                    all_pmids_for_review,
                    email=args.pubmed_email,
                    api_key=args.pubmed_api_key,
                    timeout_sec=args.timeout_sec,
                )
                bridge_nct_ids: set[str] = set()
                for nct_list in pmid_nct_map.values():
                    bridge_nct_ids.update(nct_list)
                # Exclude NCTs already found via exact linkage
                bridge_nct_ids -= pairwise_nct_ids
                summary["n_pubmed_bridge_nct_ids"] = len(bridge_nct_ids)

                if bridge_nct_ids:
                    # Resolve these NCTs and compute lambda
                    bridge_studies, _ = resolve_exact_nct_studies(
                        target_nct_ids=bridge_nct_ids,
                        query_payload=payload,
                        nct_cache=nct_cache,
                        sleep_ms=args.sleep_ms,
                        timeout_sec=args.timeout_sec,
                    )
                    if bridge_studies:
                        bridge_summary = summarize_payload(
                            review_id,
                            {"query": query_term, "studies": bridge_studies},
                        )
                        summary["n_pubmed_bridge_completed"] = bridge_summary["n_completed"]
                        summary["n_pubmed_bridge_completed_non_ghost"] = bridge_summary[
                            "n_completed_non_ghost"
                        ]
                        summary["n_pubmed_bridge_completed_ghost"] = bridge_summary[
                            "n_completed_ghost"
                        ]

        summary.setdefault("n_pubmed_bridge_completed", 0)
        summary.setdefault("n_pubmed_bridge_completed_non_ghost", 0)
        summary.setdefault("n_pubmed_bridge_completed_ghost", 0)

        # --- Fuzzy matching ---
        summary["n_fuzzy_match_nct_ids"] = 0
        # (Fuzzy matching operates per-study, not per-review, so it would need
        # individual study metadata. We track the column for future expansion.)

        # --- Determine best linkage source ---
        if summary["n_exact_nct_completed"] > 0:
            summary["linkage_source"] = "exact_rda"
        elif summary["n_pubmed_bridge_nct_ids"] > 0:
            summary["linkage_source"] = "pubmed_bridge"
        else:
            summary["linkage_source"] = "query_only"

        if "error" in payload:
            summary["fetch_error"] = str(payload.get("error", ""))
        else:
            summary["fetch_error"] = ""
        summaries.append(summary)
        print(
            f"[{idx}/{len(rows)}] {review_id}: completed={summary['n_completed']}, "
            f"ghost={summary['n_completed_ghost']}, "
            f"lambda_non_ghost={summary['lambda_completed_non_ghost_weighted']}, "
            f"exact_ids={summary['n_pairwise_nct_ids']}, "
            f"pubmed_bridge_ids={summary['n_pubmed_bridge_nct_ids']}, "
            f"linkage_source={summary['linkage_source']}, "
            f"lambda_exact={summary['lambda_exact_nct_non_ghost_weighted']}"
        )

    out = Path(args.output_csv)
    pd_fields = [
        "review_id",
        "query_term",
        "n_studies_fetched",
        "n_completed",
        "n_completed_with_pmid_any",
        "n_completed_with_results",
        "n_completed_non_ghost",
        "n_completed_ghost",
        "n_completed_results_only_no_pmid",
        "w_completed_total",
        "w_completed_pmid_any",
        "w_completed_non_ghost",
        "w_completed_ghost",
        "lambda_completed_pmid_any_weighted",
        "lambda_completed_non_ghost_weighted",
        "ghost_weight_share_completed",
        "n_pairwise_nct_ids",
        "n_exact_nct_studies_resolved",
        "n_exact_nct_ids_unresolved",
        "n_exact_nct_completed",
        "n_exact_nct_completed_non_ghost",
        "n_exact_nct_completed_ghost",
        "lambda_exact_nct_pmid_any_weighted",
        "lambda_exact_nct_non_ghost_weighted",
        "n_pubmed_bridge_nct_ids",
        "n_pubmed_bridge_completed",
        "n_pubmed_bridge_completed_non_ghost",
        "n_pubmed_bridge_completed_ghost",
        "n_fuzzy_match_nct_ids",
        "linkage_source",
        "fetch_error",
    ]
    # Sanitize string cells against spreadsheet formula injection
    for item in summaries:
        for key, val in item.items():
            if isinstance(val, str):
                item[key] = sanitize_csv_cell(val)

    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=pd_fields)
        writer.writeheader()
        for item in summaries:
            writer.writerow(item)

    n_ok = sum(1 for item in summaries if not item.get("fetch_error"))
    print(f"Wrote {out} ({len(summaries)} rows, fetch_ok={n_ok})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
