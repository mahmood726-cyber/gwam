#!/usr/bin/env python3
"""Compare multiple Pairwise70 benchmark summary.json files side-by-side."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scenario",
        action="append",
        required=True,
        help="Scenario in the format name=path/to/summary.json (repeatable).",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("pairwise70_benchmark_sensitivity/side_by_side_summary.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--output-md",
        type=Path,
        default=Path("reports/pairwise70_lambda_sensitivity_comparison.md"),
        help="Output markdown report path.",
    )
    return parser.parse_args()


def parse_named_paths(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --scenario entry '{item}'. Expected name=path.")
        name, raw_path = item.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"Invalid --scenario entry '{item}': missing name.")
        path = Path(raw_path.strip())
        if not path.exists():
            raise FileNotFoundError(f"Scenario summary not found: {path}")
        out[name] = path
    if len(out) < 2:
        raise ValueError("Provide at least two --scenario entries.")
    return out


def get_nested(data: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def load_row(name: str, path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    assumptions = get_nested(payload, "assumptions", default={}) or {}
    headline = get_nested(payload, "headline_counts", default={}) or {}
    dist = get_nested(payload, "distribution", default={}) or {}
    cross_010_count = headline.get("re_to_gwam_proxy_crossed_below_0p10")
    cross_020_count = headline.get("re_to_gwam_proxy_crossed_below_0p20")
    cross_010_pct_all = headline.get("re_to_gwam_proxy_crossed_below_0p10_pct_all")
    cross_020_pct_all = headline.get("re_to_gwam_proxy_crossed_below_0p20_pct_all")
    row = {
        "scenario": name,
        "summary_path": str(path),
        "n_files_scanned": get_nested(payload, "n_files_scanned"),
        "n_files_with_valid_binary_analyses": get_nested(payload, "n_files_with_valid_binary_analyses"),
        "n_analyses_estimable": get_nested(payload, "n_analyses_estimable"),
        "n_pet_peese_estimable": get_nested(payload, "n_pet_peese_estimable"),
        "n_selection_ipw_estimable": get_nested(payload, "n_selection_ipw_estimable"),
        "default_lambda_proxy": assumptions.get("default_lambda_proxy"),
        "lambda_proxy_lower": assumptions.get("lambda_proxy_lower"),
        "lambda_proxy_upper": assumptions.get("lambda_proxy_upper"),
        "median_lambda_proxy": dist.get("median_lambda_proxy"),
        "median_abs_shift_proxy": dist.get("median_abs_shift_proxy"),
        "median_re": dist.get("median_re"),
        "median_gwam_proxy": dist.get("median_gwam_proxy"),
        "lambda_proxy_at_lower_bound_count": dist.get("lambda_proxy_at_lower_bound_count"),
        "lambda_proxy_at_upper_bound_count": dist.get("lambda_proxy_at_upper_bound_count"),
        "lambda_proxy_at_lower_bound_pct": dist.get("lambda_proxy_at_lower_bound_pct"),
        "lambda_proxy_at_upper_bound_pct": dist.get("lambda_proxy_at_upper_bound_pct"),
        "re_abs_effect_ge_0p10": headline.get("re_abs_effect_ge_0p10"),
        "gwam_proxy_abs_effect_ge_0p10": headline.get("gwam_proxy_abs_effect_ge_0p10"),
        "re_to_gwam_proxy_crossed_below_0p10": cross_010_count,
        "re_to_gwam_proxy_crossed_below_0p10_pct_all": cross_010_pct_all,
        "re_abs_effect_ge_0p20": headline.get("re_abs_effect_ge_0p20"),
        "gwam_proxy_abs_effect_ge_0p20": headline.get("gwam_proxy_abs_effect_ge_0p20"),
        "re_to_gwam_proxy_crossed_below_0p20": cross_020_count,
        "re_to_gwam_proxy_crossed_below_0p20_pct_all": cross_020_pct_all,
    }
    n = row.get("n_analyses_estimable")
    try:
        n_float = float(n)
    except (TypeError, ValueError):
        n_float = float("nan")
    if cross_010_pct_all is not None and cross_020_pct_all is not None:
        row["crossed_below_0p10_pct"] = float(cross_010_pct_all)
        row["crossed_below_0p20_pct"] = float(cross_020_pct_all)
    elif n_float > 0:
        row["crossed_below_0p10_pct"] = float(row["re_to_gwam_proxy_crossed_below_0p10"]) / n_float
        row["crossed_below_0p20_pct"] = float(row["re_to_gwam_proxy_crossed_below_0p20"]) / n_float
    else:
        row["crossed_below_0p10_pct"] = float("nan")
        row["crossed_below_0p20_pct"] = float("nan")
    return row


def fmt_float(value: Any, digits: int = 3) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_pct(value: Any, digits: int = 1) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(x):
        return "NA"
    return f"{100.0 * x:.{digits}f}%"


def write_markdown(df: pd.DataFrame, output_md: Path) -> None:
    output_md.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Pairwise70 Lambda Sensitivity Comparison")
    lines.append("")
    lines.append("| Scenario | Lambda bounds | Default lambda | Median lambda | Lower-bound saturation | Crossed <0.10 | Crossed <0.20 | Median |d| (RE->GWAM) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for _, row in df.iterrows():
        bounds = f"[{fmt_float(row['lambda_proxy_lower'],2)}, {fmt_float(row['lambda_proxy_upper'],2)}]"
        lines.append(
            "| "
            + str(row["scenario"])
            + " | "
            + bounds
            + " | "
            + fmt_float(row["default_lambda_proxy"], 2)
            + " | "
            + fmt_float(row["median_lambda_proxy"], 3)
            + " | "
            + fmt_pct(row["lambda_proxy_at_lower_bound_pct"], 1)
            + " | "
            + fmt_pct(row["crossed_below_0p10_pct"], 1)
            + " | "
            + fmt_pct(row["crossed_below_0p20_pct"], 1)
            + " | "
            + fmt_float(row["median_abs_shift_proxy"], 3)
            + " |"
        )
    lines.append("")
    lines.append("Notes:")
    lines.append("- `Crossed <0.10` = analyses with `|mu_RE|>=0.10` and `|mu_GWAM_proxy|<0.10`.")
    lines.append("- `Crossed <0.20` = analyses with `|mu_RE|>=0.20` and `|mu_GWAM_proxy|<0.20`.")
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    named = parse_named_paths(args.scenario)
    rows = [load_row(name, path) for name, path in named.items()]
    df = pd.DataFrame(rows).sort_values("scenario").reset_index(drop=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    write_markdown(df, args.output_md)

    print(f"Wrote {args.output_csv}")
    print(f"Wrote {args.output_md}")
    print(df[["scenario", "median_lambda_proxy", "crossed_below_0p10_pct", "crossed_below_0p20_pct"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
