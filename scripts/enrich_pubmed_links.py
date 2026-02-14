#!/usr/bin/env python3
"""Enrich trial-level registry PMIDs with PubMed metadata.

Reads a CSV from fetch_ctgov_registry.py and writes a trial-publication map.
"""

from __future__ import annotations

import argparse
import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

from gwam_utils import sanitize_csv_cell

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-csv", type=Path, required=True, help="Input registry CSV.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path. Defaults to data/raw/pubmed_links_<stem>.csv",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Contact email for NCBI API request etiquette.",
    )
    parser.add_argument("--api-key", default="", help="Optional NCBI API key.")
    parser.add_argument(
        "--pmid-column",
        default="pmids_results",
        help="CSV column containing semicolon-separated PMIDs (default: pmids_results).",
    )
    return parser.parse_args()


def chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[idx : idx + size] for idx in range(0, len(values), size)]


def get_text_with_retry(
    session: requests.Session,
    *,
    url: str,
    params: dict[str, str],
    timeout: int,
    attempts: int = 5,
) -> str:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.text
        except requests.RequestException as exc:
            last_error = exc
            if attempt == attempts:
                break
            sleep_s = min(12.0, 1.5 * (2 ** (attempt - 1)))
            time.sleep(sleep_s)
    assert last_error is not None
    raise last_error


def fetch_pubmed_metadata(
    session: requests.Session, pmids: list[str], *, email: str = "", api_key: str = ""
) -> dict[str, dict[str, str]]:
    if not pmids:
        return {}

    params: dict[str, str] = {
        "db": "pubmed",
        "retmode": "xml",
        "id": ",".join(pmids),
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    xml_text = get_text_with_retry(session, url=EFETCH_URL, params=params, timeout=60)
    if len(xml_text) > 50_000_000:  # 50 MB safety cap
        raise ValueError(f"PubMed XML response too large ({len(xml_text)} bytes); refusing to parse.")
    parser = ET.XMLParser(encoding="utf-8")
    root = ET.fromstring(xml_text.encode("utf-8"), parser=parser)

    out: dict[str, dict[str, str]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID", default="").strip()
        title = article.findtext(".//ArticleTitle", default="").strip()
        journal = article.findtext(".//Journal/Title", default="").strip()
        year = article.findtext(".//JournalIssue/PubDate/Year", default="").strip()
        if not year:
            medline_date = article.findtext(".//JournalIssue/PubDate/MedlineDate", default="").strip()
            year = medline_date[:4] if medline_date else ""
        out[pmid] = {
            "pubmed_title": title,
            "pubmed_journal": journal,
            "pubmed_year": year,
        }
    return out


def main() -> int:
    args = parse_args()
    default_output = Path("data/raw") / f"pubmed_links_{args.registry_csv.stem}.csv"
    output = args.output or default_output
    output.parent.mkdir(parents=True, exist_ok=True)

    with args.registry_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    trial_pub_rows: list[dict[str, str]] = []
    unique_pmids: set[str] = set()
    pmid_column = args.pmid_column
    if rows and pmid_column not in rows[0]:
        if "pmids" in rows[0]:
            pmid_column = "pmids"
        elif "pmids_any" in rows[0]:
            pmid_column = "pmids_any"
        else:
            raise ValueError(f"PMID column '{args.pmid_column}' not found in registry CSV.")
    for row in rows:
        pmids = [x.strip() for x in row.get(pmid_column, "").split(";") if x.strip()]
        for pmid in pmids:
            unique_pmids.add(pmid)
            trial_pub_rows.append(
                {
                    "nct_id": row.get("nct_id", ""),
                    "pmid": pmid,
                }
            )

    meta: dict[str, dict[str, str]] = {}
    with requests.Session() as session:
        for pmid_batch in chunked(sorted(unique_pmids), 200):
            try:
                meta.update(
                    fetch_pubmed_metadata(
                        session,
                        pmid_batch,
                        email=args.email,
                        api_key=args.api_key,
                    )
                )
            except Exception as exc:
                # Sanitize exception message to avoid leaking API key
                msg = str(exc)
                if args.api_key and args.api_key in msg:
                    msg = msg.replace(args.api_key, "[REDACTED]")
                print(f"WARNING: PubMed batch fetch failed ({len(pmid_batch)} PMIDs): {msg}")
                continue

    output_rows: list[dict[str, str]] = []
    for item in trial_pub_rows:
        details = meta.get(item["pmid"], {})
        output_rows.append(
            {
                "nct_id": item["nct_id"],
                "pmid": item["pmid"],
                "pubmed_year": sanitize_csv_cell(details.get("pubmed_year", "")),
                "pubmed_journal": sanitize_csv_cell(details.get("pubmed_journal", "")),
                "pubmed_title": sanitize_csv_cell(details.get("pubmed_title", "")),
            }
        )

    fieldnames = ["nct_id", "pmid", "pubmed_year", "pubmed_journal", "pubmed_title"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {output}")
    print(
        f"pmid_column={pmid_column}, unique_pmids={len(unique_pmids)}, "
        f"trial_publication_links={len(output_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
