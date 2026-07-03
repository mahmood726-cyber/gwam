#!/usr/bin/env python3
"""Self-contained, reproducible GWAM benchmark on bundled registry data.

Unlike the Pairwise70 benchmark (which requires an external dataset path),
this benchmark runs entirely on the registry CSVs shipped in ``data/raw/``.
It exercises the reusable ``gwam_correct()`` API on every bundled
intervention/condition pair and reports the structural GWAM quantities:

    - lambda_pmid_only   : enrollment-weighted fraction with a linked results PMID
    - lambda_non_ghost   : enrollment-weighted non-ghost fraction
    - attenuation        : mu_gwam / mu_published under the null-ghost assumption
                           (equals lambda_pmid_only when ghost/results-only
                           means are held at 0)
    - study counts per integrity stratum

These quantities depend only on the registry weight structure, NOT on the
magnitude of the published effect, so the cross-dataset comparison is
reproducible without per-dataset literature estimates. For datasets whose
published pooled effect is documented in this repo, the corrected pooled mean
is also reported; otherwise a normalised ``published_mu = 1.0`` illustrative
correction is shown (clearly labelled -- it is a structural attenuation, not a
scientific effect claim).

Usage::

    python scripts/run_bundled_benchmark.py
    python scripts/run_bundled_benchmark.py --output-csv reports/bundled_benchmark.csv

Output is deterministic (no Monte-Carlo): the same inputs always produce the
same table.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from gwam_utils import build_environment_metadata, gwam_correct, sanitize_csv_cell

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"

# Published pooled effects that are documented and verified within this repo.
# Keyed by the registry CSV stem. Values NOT listed here are benchmarked on a
# normalised published_mu = 1.0 (structural attenuation only). Do not add a
# value here without a documented source in the repository.
DOCUMENTED_PUBLISHED_MU: dict[str, float] = {
    # log(1.29); see reports/final_run_report_escitalopram_depression.md
    "escitalopram__depression": 0.25464221837358075,
}


def discover_registry_csvs() -> list[Path]:
    """Return bundled registry CSVs that carry the required GWAM columns."""
    csvs: list[Path] = []
    for path in sorted(RAW_DIR.glob("*.csv")):
        if "pubmed" in path.name.lower():
            continue
        try:
            with path.open("r", newline="", encoding="utf-8") as handle:
                header = next(csv.reader(handle), [])
        except (OSError, StopIteration):
            continue
        if "is_ghost_protocol" in header and "weight_n" in header:
            csvs.append(path)
    return csvs


def _load_rows(path: Path) -> list[dict]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def benchmark_dataset(path: Path) -> dict:
    """Run the deterministic GWAM correction on one registry CSV."""
    rows = _load_rows(path)
    stem = path.stem
    documented_mu = DOCUMENTED_PUBLISHED_MU.get(stem)
    published_mu = documented_mu if documented_mu is not None else 1.0
    result = gwam_correct(rows, published_mu=published_mu)

    # Under the default null-ghost assumption the correction is a pure
    # multiplicative attenuation of the published effect.
    attenuation = (
        result.mu_gwam_null_point / published_mu if published_mu != 0 else float("nan")
    )
    return {
        "dataset": stem,
        "n_trials": len(rows),
        "n_pmid_only": result.n_pmid_only,
        "n_results_only": result.n_results_only,
        "n_ghost": result.n_ghost,
        "lambda_pmid_only": result.lambda_pmid_only,
        "lambda_non_ghost": result.lambda_non_ghost,
        "attenuation": attenuation,
        "published_mu_source": "documented" if documented_mu is not None else "normalised(1.0)",
        "published_mu": published_mu,
        "mu_gwam_null_point": result.mu_gwam_null_point,
    }


def format_table(records: list[dict]) -> str:
    """Render an aligned, monospace text table."""
    headers = [
        ("dataset", 34, "s"),
        ("k", 4, "d"),
        ("pmid", 5, "d"),
        ("res", 4, "d"),
        ("ghost", 6, "d"),
        ("lam_pmid", 9, "f"),
        ("lam_ngh", 8, "f"),
        ("attenu", 7, "f"),
        ("mu_src", 15, "s"),
    ]
    keymap = {
        "dataset": "dataset",
        "k": "n_trials",
        "pmid": "n_pmid_only",
        "res": "n_results_only",
        "ghost": "n_ghost",
        "lam_pmid": "lambda_pmid_only",
        "lam_ngh": "lambda_non_ghost",
        "attenu": "attenuation",
        "mu_src": "published_mu_source",
    }
    lines = []
    head = " ".join(f"{name:<{w}}" for name, w, _ in headers)
    lines.append(head)
    lines.append("-" * len(head))
    for rec in records:
        cells = []
        for name, w, fmt in headers:
            val = rec[keymap[name]]
            if fmt == "s":
                cells.append(f"{str(val):<{w}.{w}}")
            elif fmt == "d":
                cells.append(f"{int(val):<{w}d}")
            else:  # "f" -> fixed 4-decimal float, left-padded to width
                cells.append(f"{float(val):<{w}.4f}")
        lines.append(" ".join(cells))
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Optional path to write the per-dataset benchmark table as CSV.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to write full results (with environment metadata) as JSON.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    csvs = discover_registry_csvs()
    if not csvs:
        print(f"ERROR: no bundled registry CSVs found under {RAW_DIR}", file=sys.stderr)
        return 1

    records = [benchmark_dataset(path) for path in csvs]

    print(f"GWAM bundled benchmark -- {len(records)} dataset(s) under {RAW_DIR}")
    print()
    print(format_table(records))
    print()
    print(
        "Notes: lam_pmid/lam_ngh are enrollment-weighted integrity ratios; "
        "attenu = mu_gwam/mu_published (structural, published_mu-independent). "
        "mu_src=normalised(1.0) rows report structural attenuation only."
    )

    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = list(records[0].keys())
        with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for rec in records:
                writer.writerow({k: sanitize_csv_cell(str(v)) for k, v in rec.items()})
        print(f"Wrote {args.output_csv}")

    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {"metadata": build_environment_metadata(), "results": records}
        with args.output_json.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
        print(f"Wrote {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
