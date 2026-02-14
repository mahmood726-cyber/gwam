#!/usr/bin/env python3
"""Run a real-data benchmark on Pairwise70 using GWAM-style sensitivity outputs.

This script:
1) Loads Pairwise70 .rda datasets.
2) Builds binary-outcome analysis units from study-level event counts.
3) Runs RE, selection-IPW RE, and PET-PEESE.
4) Applies GWAM shrinkage as a sensitivity grid and a conservative CT.gov proxy lambda.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyreadr

from scipy.special import ndtr as normal_cdf

from gwam_utils import build_environment_metadata, normal_quantile, safe_float
from grey_meta_v8 import GRMA
from model_gwam_bayesian import PriorSpec, bayesian_gwam_posterior
from simulate_gwam_vs_re import pet_peese_estimate, random_effects_dl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairwise-data-dir",
        type=Path,
        required=True,
        help="Directory containing Pairwise70 .rda files.",
    )
    parser.add_argument(
        "--ctgov-covariates-csv",
        type=Path,
        required=True,
        help="Review-level CT.gov covariates used to build lambda proxies.",
    )
    parser.add_argument(
        "--ctgov-linkage-csv",
        type=Path,
        default=Path("pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv"),
        help=(
            "Review-level CT.gov linkage summary built by "
            "scripts/build_pairwise70_ctgov_linkage_summary.py."
        ),
    )
    parser.add_argument(
        "--lambda-source-mode",
        choices=[
            "auto",
            "transport_proxy",
            "linkage_pmid",
            "linkage_non_ghost",
            "exact_nct_pmid",
            "exact_nct_non_ghost",
            "pubmed_bridge_non_ghost",
        ],
        default="auto",
        help=(
            "Source for lambda proxies. auto prefers exact_nct_non_ghost, then "
            "pubmed_bridge_non_ghost, then linkage_non_ghost, "
            "then transport proxy, then default lambda fallback."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("pairwise70_benchmark"),
        help="Output directory for benchmark artifacts.",
    )
    parser.add_argument("--min-k", type=int, default=2, help="Minimum studies required for RE estimation.")
    parser.add_argument(
        "--min-k-pet",
        type=int,
        default=3,
        help="Minimum studies required for PET-PEESE estimation.",
    )
    parser.add_argument(
        "--lambda-grid",
        default="0.2,0.4,0.6,0.8",
        help="Comma-separated GWAM lambda sensitivity values.",
    )
    parser.add_argument(
        "--default-lambda-proxy",
        type=float,
        default=0.7,
        help="Fallback lambda used when CT.gov covariates are unavailable.",
    )
    parser.add_argument(
        "--lambda-proxy-lower",
        type=float,
        default=0.2,
        help="Lower bound for CT.gov-derived lambda proxy.",
    )
    parser.add_argument(
        "--lambda-proxy-upper",
        type=float,
        default=0.95,
        help="Upper bound for CT.gov-derived lambda proxy.",
    )
    parser.add_argument(
        "--p-sig",
        type=float,
        default=0.95,
        help="Assumed publication probability for significant positive studies in IPW model.",
    )
    parser.add_argument(
        "--lambda-uncertainty-model",
        choices=["deterministic", "beta"],
        default="deterministic",
        help=(
            "Lambda treatment: deterministic uses point lambda; beta uses a Beta posterior "
            "from completed-study counts (or pseudo-Beta fallback) and propagates uncertainty."
        ),
    )
    parser.add_argument(
        "--lambda-posterior-draws",
        type=int,
        default=400,
        help="Monte Carlo draws per analysis when --lambda-uncertainty-model=beta.",
    )
    parser.add_argument(
        "--lambda-beta-prior-alpha",
        type=float,
        default=1.0,
        help="Alpha prior for count-based Beta lambda posterior.",
    )
    parser.add_argument(
        "--lambda-beta-prior-beta",
        type=float,
        default=1.0,
        help="Beta prior for count-based Beta lambda posterior.",
    )
    parser.add_argument(
        "--lambda-pseudo-concentration",
        type=float,
        default=20.0,
        help=(
            "Effective sample size for pseudo-Beta lambda posterior when completed-study counts "
            "are unavailable (e.g., transport/default sources)."
        ),
    )
    parser.add_argument(
        "--disable-lambda-clipping",
        action="store_true",
        help="Disable hard lower/upper clipping of deterministic lambda proxies.",
    )
    parser.add_argument(
        "--review-bootstrap-iters",
        type=int,
        default=400,
        help="Cluster bootstrap iterations for review-level uncertainty summaries.",
    )
    parser.add_argument(
        "--review-bootstrap-seed",
        type=int,
        default=42,
        help="Seed for review-level cluster bootstrap.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Optional cap for number of .rda files processed (debug mode).",
    )
    return parser.parse_args()


def parse_lambda_grid(text: str) -> list[float]:
    values: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        val = float(token)
        if not (0 < val <= 1):
            raise ValueError("All lambda-grid values must be in (0, 1].")
        values.append(val)
    if not values:
        raise ValueError("--lambda-grid produced no values.")
    return values


# normal_cdf is imported from scipy.special.ndtr (truly vectorized, no Python loop)


def extract_review_id(file_name: str) -> str:
    stem = Path(file_name).stem
    match = re.match(r"^(CD\d+)(?:_pub\d+)?_data$", stem)
    if match:
        return match.group(1)
    return stem


def make_analysis_key(df: pd.DataFrame) -> pd.Series:
    pieces = []
    for col in ["Analysis.group", "Analysis.number", "Analysis.name", "Subgroup"]:
        if col in df.columns:
            pieces.append(df[col].fillna("").astype(str).str.strip())
        else:
            pieces.append(pd.Series([""] * len(df), index=df.index))
    return pieces[0] + "|" + pieces[1] + "|" + pieces[2] + "|" + pieces[3]


def load_ctgov_covariates(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["review_id", "trial_count", "enrollment_mean"])
    frame = pd.read_csv(path)
    keep = ["review_id", "trial_count", "enrollment_mean"]
    for col in keep:
        if col not in frame.columns:
            frame[col] = np.nan
    frame = frame[keep].copy()
    frame["trial_count"] = pd.to_numeric(frame["trial_count"], errors="coerce")
    frame["enrollment_mean"] = pd.to_numeric(frame["enrollment_mean"], errors="coerce")
    return frame


def load_ctgov_linkage(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=[
                "review_id",
                "n_completed",
                "n_completed_with_pmid_any",
                "lambda_completed_pmid_any_weighted",
                "lambda_completed_non_ghost_weighted",
                "lambda_exact_nct_pmid_any_weighted",
                "lambda_exact_nct_non_ghost_weighted",
                "n_completed_ghost",
                "n_completed_non_ghost",
                "n_pairwise_nct_ids",
                "n_exact_nct_studies_resolved",
                "n_pubmed_bridge_nct_ids",
                "n_pubmed_bridge_completed",
                "n_pubmed_bridge_completed_non_ghost",
                "n_pubmed_bridge_completed_ghost",
                "n_fuzzy_match_nct_ids",
                "linkage_source",
                "n_exact_nct_ids_unresolved",
                "n_exact_nct_completed",
                "n_exact_nct_completed_non_ghost",
                "n_exact_nct_completed_ghost",
            ]
        )
    frame = pd.read_csv(path)
    keep = [
        "review_id",
        "n_completed",
        "n_completed_with_pmid_any",
        "lambda_completed_pmid_any_weighted",
        "lambda_completed_non_ghost_weighted",
        "lambda_exact_nct_pmid_any_weighted",
        "lambda_exact_nct_non_ghost_weighted",
        "n_completed_ghost",
        "n_completed_non_ghost",
        "n_pairwise_nct_ids",
        "n_exact_nct_studies_resolved",
        "n_exact_nct_ids_unresolved",
        "n_exact_nct_completed",
        "n_exact_nct_completed_non_ghost",
        "n_exact_nct_completed_ghost",
        "n_pubmed_bridge_nct_ids",
        "n_pubmed_bridge_completed",
        "n_pubmed_bridge_completed_non_ghost",
        "n_pubmed_bridge_completed_ghost",
        "n_fuzzy_match_nct_ids",
        "linkage_source",
    ]
    for col in keep:
        if col not in frame.columns:
            frame[col] = np.nan if col != "linkage_source" else ""
    frame = frame[keep].copy()
    # Coerce numeric columns (skip review_id and linkage_source)
    numeric_cols = [c for c in keep[1:] if c != "linkage_source"]
    for col in numeric_cols:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def calibrate_p_nonsig_from_lambda(lambda_target: float, q_sig: float, p_sig: float) -> float:
    q = float(min(max(q_sig, 0.0), 1.0))
    lam = float(min(max(lambda_target, 0.01), p_sig))
    if q >= 0.999:
        return float(min(max(lam, 0.01), p_sig))
    denom = max(1e-9, 1.0 - q)
    p_nonsig = (lam - (q * p_sig)) / denom
    return float(min(max(p_nonsig, 0.01), p_sig))


def apply_lambda_bounds(
    value: float,
    *,
    lower: float,
    upper: float,
    clipping_enabled: bool,
) -> float:
    if not math.isfinite(value):
        return float("nan")
    if clipping_enabled:
        return float(min(max(value, lower), upper))
    return float(value)


def _bounded_unit(value: float, eps: float = 1e-6) -> float:
    return float(min(max(value, eps), 1.0 - eps))


def build_lambda_posterior_params(
    *,
    lambda_source: str,
    lambda_center: float,
    review_link: dict[str, Any],
    prior_alpha: float,
    prior_beta: float,
    pseudo_concentration: float,
) -> dict[str, float | str]:
    success = float("nan")
    failure = float("nan")
    kind = "pseudo_beta"

    if lambda_source == "pairwise_exact_nct_non_ghost":
        success = safe_float(review_link.get("n_exact_nct_completed_non_ghost", float("nan")))
        failure = safe_float(review_link.get("n_exact_nct_completed_ghost", float("nan")))
        kind = "count_beta"
    elif lambda_source == "pubmed_bridge_non_ghost":
        success = safe_float(review_link.get("n_pubmed_bridge_completed_non_ghost", float("nan")))
        failure = safe_float(review_link.get("n_pubmed_bridge_completed_ghost", float("nan")))
        kind = "count_beta"
    elif lambda_source == "ctgov_linkage_non_ghost":
        success = safe_float(review_link.get("n_completed_non_ghost", float("nan")))
        failure = safe_float(review_link.get("n_completed_ghost", float("nan")))
        kind = "count_beta"
    elif lambda_source == "ctgov_linkage_pmid":
        completed = safe_float(review_link.get("n_completed", float("nan")))
        with_pmid = safe_float(review_link.get("n_completed_with_pmid_any", float("nan")))
        if math.isfinite(completed) and math.isfinite(with_pmid):
            success = max(0.0, with_pmid)
            failure = max(0.0, completed - with_pmid)
            kind = "count_beta"

    if (
        kind == "count_beta"
        and math.isfinite(success)
        and math.isfinite(failure)
        and success >= 0
        and failure >= 0
        and (success + failure) > 0
    ):
        alpha = float(prior_alpha + success)
        beta = float(prior_beta + failure)
    else:
        kind = "pseudo_beta"
        m = _bounded_unit(lambda_center if math.isfinite(lambda_center) else 0.7)
        eff_n = max(float(pseudo_concentration), 2.0)
        alpha = float(max(prior_alpha, m * eff_n))
        beta = float(max(prior_beta, (1.0 - m) * eff_n))
        success = float("nan")
        failure = float("nan")

    total = alpha + beta
    mean = float(alpha / total)
    var = float((alpha * beta) / ((total**2) * (total + 1.0)))
    sd = float(math.sqrt(max(var, 0.0)))
    return {
        "kind": kind,
        "alpha": alpha,
        "beta": beta,
        "mean": mean,
        "sd": sd,
        "effective_n": float(total),
        "source_success": float(success) if math.isfinite(success) else float("nan"),
        "source_failure": float(failure) if math.isfinite(failure) else float("nan"),
    }


def build_review_cluster_summary(
    df: pd.DataFrame,
    *,
    bootstrap_iters: int,
    bootstrap_seed: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    if df.empty:
        empty = {
            "n_reviews": 0,
            "n_analyses": 0,
            "mean_analyses_per_review": float("nan"),
            "median_analyses_per_review": float("nan"),
            "crossed_below_0p10_pct_all_weighted": float("nan"),
            "crossed_below_0p20_pct_all_weighted": float("nan"),
            "crossed_below_0p10_pct_all_review_mean": float("nan"),
            "crossed_below_0p20_pct_all_review_mean": float("nan"),
            "crossed_below_0p10_pct_all_weighted_boot_ci_lo": float("nan"),
            "crossed_below_0p10_pct_all_weighted_boot_ci_hi": float("nan"),
            "crossed_below_0p20_pct_all_weighted_boot_ci_lo": float("nan"),
            "crossed_below_0p20_pct_all_weighted_boot_ci_hi": float("nan"),
            "bootstrap_iters": int(bootstrap_iters),
        }
        return empty, pd.DataFrame()

    grouped = (
        df.groupby("review_id", dropna=False)
        .agg(
            n_analyses=("review_id", "size"),
            crossed_010=("cross_prob_010", "sum"),
            crossed_020=("cross_prob_020", "sum"),
        )
        .reset_index()
    )
    grouped["crossed_010_rate_within_review"] = grouped["crossed_010"] / grouped["n_analyses"]
    grouped["crossed_020_rate_within_review"] = grouped["crossed_020"] / grouped["n_analyses"]

    total_n = float(grouped["n_analyses"].sum())
    weighted_010 = float(grouped["crossed_010"].sum() / total_n)
    weighted_020 = float(grouped["crossed_020"].sum() / total_n)
    review_mean_010 = float(grouped["crossed_010_rate_within_review"].mean())
    review_mean_020 = float(grouped["crossed_020_rate_within_review"].mean())

    ci_010_lo = float("nan")
    ci_010_hi = float("nan")
    ci_020_lo = float("nan")
    ci_020_hi = float("nan")

    if bootstrap_iters > 0 and len(grouped) > 1:
        rng = np.random.default_rng(int(bootstrap_seed))
        arr_n = grouped["n_analyses"].to_numpy(dtype=float)
        arr_010 = grouped["crossed_010"].to_numpy(dtype=float)
        arr_020 = grouped["crossed_020"].to_numpy(dtype=float)
        k = len(grouped)
        boot_010: list[float] = []
        boot_020: list[float] = []
        for _ in range(int(bootstrap_iters)):
            idx = rng.integers(0, k, size=k)
            denom = float(np.sum(arr_n[idx]))
            if denom <= 0:
                continue
            boot_010.append(float(np.sum(arr_010[idx]) / denom))
            boot_020.append(float(np.sum(arr_020[idx]) / denom))
        if boot_010 and boot_020:
            ci_010_lo = float(np.quantile(boot_010, 0.025))
            ci_010_hi = float(np.quantile(boot_010, 0.975))
            ci_020_lo = float(np.quantile(boot_020, 0.025))
            ci_020_hi = float(np.quantile(boot_020, 0.975))

    summary = {
        "n_reviews": int(len(grouped)),
        "n_analyses": int(total_n),
        "mean_analyses_per_review": float(np.mean(grouped["n_analyses"])),
        "median_analyses_per_review": float(np.median(grouped["n_analyses"])),
        "crossed_below_0p10_pct_all_weighted": weighted_010,
        "crossed_below_0p20_pct_all_weighted": weighted_020,
        "crossed_below_0p10_pct_all_review_mean": review_mean_010,
        "crossed_below_0p20_pct_all_review_mean": review_mean_020,
        "crossed_below_0p10_pct_all_weighted_boot_ci_lo": ci_010_lo,
        "crossed_below_0p10_pct_all_weighted_boot_ci_hi": ci_010_hi,
        "crossed_below_0p20_pct_all_weighted_boot_ci_lo": ci_020_lo,
        "crossed_below_0p20_pct_all_weighted_boot_ci_hi": ci_020_hi,
        "bootstrap_iters": int(bootstrap_iters),
    }
    return summary, grouped


def parse_pairwise_object(file_path: Path) -> pd.DataFrame:
    payload = pyreadr.read_r(str(file_path))
    if not payload:
        return pd.DataFrame()
    frame = next(iter(payload.values()))
    if not isinstance(frame, pd.DataFrame):
        return pd.DataFrame()
    return frame.copy()


def lambda_column_suffix(value: float) -> str:
    return str(value).replace(".", "p")


def build_attenuation_metrics(
    df: pd.DataFrame,
    *,
    mu_gwam_col: str = "mu_gwam_proxy",
    cross_prob_010_col: str = "cross_prob_010",
    cross_prob_020_col: str = "cross_prob_020",
) -> dict[str, float]:
    n = int(len(df))
    if n == 0:
        return {
            "re_abs_effect_ge_0p10": 0.0,
            "gwam_proxy_abs_effect_ge_0p10": 0.0,
            "re_to_gwam_proxy_crossed_below_0p10": 0.0,
            "re_to_gwam_proxy_crossed_below_0p10_pct_all": float("nan"),
            "re_to_gwam_proxy_crossed_below_0p10_pct_within_re_ge": float("nan"),
            "re_abs_effect_ge_0p20": 0.0,
            "gwam_proxy_abs_effect_ge_0p20": 0.0,
            "re_to_gwam_proxy_crossed_below_0p20": 0.0,
            "re_to_gwam_proxy_crossed_below_0p20_pct_all": float("nan"),
            "re_to_gwam_proxy_crossed_below_0p20_pct_within_re_ge": float("nan"),
        }

    abs_re = np.abs(df["mu_re"].to_numpy(dtype=float))
    abs_gwam = np.abs(df[mu_gwam_col].to_numpy(dtype=float))

    re_ge_010 = int(np.sum(abs_re >= 0.10))
    gwam_ge_010 = int(np.sum(abs_gwam >= 0.10))
    if cross_prob_010_col in df.columns:
        crossed_010 = float(np.sum(pd.to_numeric(df[cross_prob_010_col], errors="coerce").fillna(0.0)))
    else:
        crossed_010 = float(np.sum((abs_re >= 0.10) & (abs_gwam < 0.10)))

    re_ge_020 = int(np.sum(abs_re >= 0.20))
    gwam_ge_020 = int(np.sum(abs_gwam >= 0.20))
    if cross_prob_020_col in df.columns:
        crossed_020 = float(np.sum(pd.to_numeric(df[cross_prob_020_col], errors="coerce").fillna(0.0)))
    else:
        crossed_020 = float(np.sum((abs_re >= 0.20) & (abs_gwam < 0.20)))

    return {
        "re_abs_effect_ge_0p10": float(re_ge_010),
        "gwam_proxy_abs_effect_ge_0p10": float(gwam_ge_010),
        "re_to_gwam_proxy_crossed_below_0p10": crossed_010,
        "re_to_gwam_proxy_crossed_below_0p10_pct_all": float(crossed_010 / n),
        "re_to_gwam_proxy_crossed_below_0p10_pct_within_re_ge": (
            float(crossed_010 / re_ge_010) if re_ge_010 > 0 else float("nan")
        ),
        "re_abs_effect_ge_0p20": float(re_ge_020),
        "gwam_proxy_abs_effect_ge_0p20": float(gwam_ge_020),
        "re_to_gwam_proxy_crossed_below_0p20": crossed_020,
        "re_to_gwam_proxy_crossed_below_0p20_pct_all": float(crossed_020 / n),
        "re_to_gwam_proxy_crossed_below_0p20_pct_within_re_ge": (
            float(crossed_020 / re_ge_020) if re_ge_020 > 0 else float("nan")
        ),
    }


def build_distribution_metrics(
    df: pd.DataFrame,
    *,
    lambda_proxy_lower: float,
    lambda_proxy_upper: float,
    clipping_enabled: bool,
) -> dict[str, float]:
    n = int(len(df))
    if n == 0:
        return {
            "median_k": float("nan"),
            "median_lambda_proxy": float("nan"),
            "median_abs_shift_proxy": float("nan"),
            "median_re": float("nan"),
            "median_gwam_proxy": float("nan"),
            "lambda_proxy_at_lower_bound_count": 0.0,
            "lambda_proxy_at_upper_bound_count": 0.0,
            "lambda_proxy_at_lower_bound_pct": float("nan"),
            "lambda_proxy_at_upper_bound_pct": float("nan"),
        }
    if clipping_enabled:
        lambda_at_lower = int(np.sum(np.isclose(df["lambda_proxy"], lambda_proxy_lower)))
        lambda_at_upper = int(np.sum(np.isclose(df["lambda_proxy"], lambda_proxy_upper)))
        lambda_lower_pct = float(lambda_at_lower / n)
        lambda_upper_pct = float(lambda_at_upper / n)
    else:
        lambda_at_lower = 0
        lambda_at_upper = 0
        lambda_lower_pct = float("nan")
        lambda_upper_pct = float("nan")
    return {
        "median_k": float(np.nanmedian(df["k"])),
        "median_lambda_proxy": float(np.nanmedian(df["lambda_proxy"])),
        "median_abs_shift_proxy": float(np.nanmedian(np.abs(df["mu_re"] - df["mu_gwam_proxy"]))),
        "median_re": float(np.nanmedian(df["mu_re"])),
        "median_gwam_proxy": float(np.nanmedian(df["mu_gwam_proxy"])),
        "lambda_proxy_at_lower_bound_count": float(lambda_at_lower),
        "lambda_proxy_at_upper_bound_count": float(lambda_at_upper),
        "lambda_proxy_at_lower_bound_pct": lambda_lower_pct,
        "lambda_proxy_at_upper_bound_pct": lambda_upper_pct,
    }


def build_subset_summary(
    df: pd.DataFrame,
    *,
    lambda_proxy_lower: float,
    lambda_proxy_upper: float,
    clipping_enabled: bool,
    mu_gwam_col: str = "mu_gwam_proxy",
    cross_prob_010_col: str = "cross_prob_010",
    cross_prob_020_col: str = "cross_prob_020",
) -> dict[str, Any]:
    n = int(len(df))
    return {
        "n_analyses": n,
        "distribution": build_distribution_metrics(
            df,
            lambda_proxy_lower=lambda_proxy_lower,
            lambda_proxy_upper=lambda_proxy_upper,
            clipping_enabled=clipping_enabled,
        ),
        "attenuation": build_attenuation_metrics(
            df,
            mu_gwam_col=mu_gwam_col,
            cross_prob_010_col=cross_prob_010_col,
            cross_prob_020_col=cross_prob_020_col,
        ),
    }


def flatten_subset_row(name: str, subset_summary: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {"subset": name, "n_analyses": subset_summary["n_analyses"]}
    for k, v in subset_summary["distribution"].items():
        out[k] = v
    for k, v in subset_summary["attenuation"].items():
        out[k] = v
    return out


def _build_grma_summary(results: pd.DataFrame) -> dict[str, Any]:
    """Build GRMA vs RE comparison metrics for summary.json."""
    grma_valid = results.loc[np.isfinite(results["mu_grma"])].copy()
    n_grma = len(grma_valid)
    if n_grma == 0:
        return {"n_grma_estimable": 0}
    abs_shift = np.abs(grma_valid["mu_grma"].values - grma_valid["mu_re"].values)
    rel_shift = abs_shift / np.maximum(np.abs(grma_valid["mu_re"].values), 1e-9)
    n_grma_pos_sig = int(np.sum(grma_valid["grma_positive_sig"]))
    n_re_pos_sig_in_grma = int(np.sum(grma_valid["re_positive_sig"]))
    n_re_to_grma_sig_drop = int(
        np.sum(grma_valid["re_positive_sig"] & (~grma_valid["grma_positive_sig"]))
    )
    n_grma_to_re_sig_gain = int(
        np.sum((~grma_valid["re_positive_sig"]) & grma_valid["grma_positive_sig"])
    )
    return {
        "n_grma_estimable": n_grma,
        "n_grma_positive_significant": n_grma_pos_sig,
        "n_re_positive_significant_where_grma_estimable": n_re_pos_sig_in_grma,
        "re_to_grma_significance_drops": n_re_to_grma_sig_drop,
        "grma_to_re_significance_gains": n_grma_to_re_sig_gain,
        "median_abs_shift_grma_vs_re": float(np.median(abs_shift)),
        "mean_abs_shift_grma_vs_re": float(np.mean(abs_shift)),
        "median_rel_shift_grma_vs_re": float(np.median(rel_shift)),
        "median_w_max_grma": float(np.median(grma_valid["w_max_grma"].values)),
        "median_n_eff_grma": float(np.median(grma_valid["n_eff_grma"].values)),
        "correlation_grma_re": float(np.corrcoef(grma_valid["mu_grma"].values, grma_valid["mu_re"].values)[0, 1]),
        "by_k_band": {
            "k_3": _grma_k_band_summary(grma_valid, 3, 3),
            "k_4_to_9": _grma_k_band_summary(grma_valid, 4, 9),
            "k_ge_10": _grma_k_band_summary(grma_valid, 10, 999999),
        },
    }


def _grma_k_band_summary(df: pd.DataFrame, k_lo: int, k_hi: int) -> dict[str, Any]:
    sub = df.loc[(df["k"] >= k_lo) & (df["k"] <= k_hi)]
    if len(sub) == 0:
        return {"n": 0}
    abs_shift = np.abs(sub["mu_grma"].values - sub["mu_re"].values)
    return {
        "n": int(len(sub)),
        "median_abs_shift": float(np.median(abs_shift)),
        "mean_abs_shift": float(np.mean(abs_shift)),
        "correlation_grma_re": float(np.corrcoef(sub["mu_grma"].values, sub["mu_re"].values)[0, 1]),
        "median_w_max": float(np.median(sub["w_max_grma"].values)),
        "median_n_eff": float(np.median(sub["n_eff_grma"].values)),
    }


def main() -> int:
    args = parse_args()
    if not (0 < args.default_lambda_proxy <= 1):
        raise ValueError("--default-lambda-proxy must be in (0, 1].")
    if not (0 < args.lambda_proxy_lower <= 1):
        raise ValueError("--lambda-proxy-lower must be in (0, 1].")
    if not (0 < args.lambda_proxy_upper <= 1):
        raise ValueError("--lambda-proxy-upper must be in (0, 1].")
    if args.lambda_proxy_lower > args.lambda_proxy_upper:
        raise ValueError("--lambda-proxy-lower cannot exceed --lambda-proxy-upper.")
    if not (0 < args.p_sig <= 1):
        raise ValueError("--p-sig must be in (0, 1].")
    if args.lambda_posterior_draws < 50:
        raise ValueError("--lambda-posterior-draws must be >= 50.")
    if args.lambda_beta_prior_alpha <= 0 or args.lambda_beta_prior_beta <= 0:
        raise ValueError("--lambda-beta-prior-alpha and --lambda-beta-prior-beta must be > 0.")
    if args.lambda_pseudo_concentration <= 0:
        raise ValueError("--lambda-pseudo-concentration must be > 0.")
    if args.review_bootstrap_iters < 0:
        raise ValueError("--review-bootstrap-iters must be >= 0.")
    if args.min_k < 2:
        raise ValueError("--min-k must be >= 2.")
    if args.min_k_pet < 3:
        raise ValueError("--min-k-pet must be >= 3.")

    lambda_grid = parse_lambda_grid(args.lambda_grid)
    clipping_enabled = not args.disable_lambda_clipping
    rng_lambda = np.random.default_rng(int(args.review_bootstrap_seed) + 1001)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cov = load_ctgov_covariates(args.ctgov_covariates_csv)
    cov_map = cov.set_index("review_id").to_dict(orient="index")
    linkage = load_ctgov_linkage(args.ctgov_linkage_csv)
    linkage_map = linkage.set_index("review_id").to_dict(orient="index")

    files = sorted(args.pairwise_data_dir.glob("*.rda"))
    if args.max_files is not None:
        files = files[: max(0, args.max_files)]
    if not files:
        raise FileNotFoundError(f"No .rda files found in {args.pairwise_data_dir}")

    all_rows: list[dict[str, Any]] = []
    files_with_valid_binary = 0

    for file_path in files:
        review_id = extract_review_id(file_path.name)
        df = parse_pairwise_object(file_path)
        if df.empty:
            continue

        required = ["Experimental.cases", "Experimental.N", "Control.cases", "Control.N"]
        if not all(col in df.columns for col in required):
            continue

        for col in required:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df["analysis_key"] = make_analysis_key(df)

        valid_file_binary = False
        grouped = df.groupby("analysis_key", sort=False)
        for analysis_key, sub in grouped:
            clean = sub.dropna(subset=required).copy()
            if clean.empty:
                continue

            e_t = clean["Experimental.cases"].to_numpy(dtype=float)
            n_t = clean["Experimental.N"].to_numpy(dtype=float)
            e_c = clean["Control.cases"].to_numpy(dtype=float)
            n_c = clean["Control.N"].to_numpy(dtype=float)

            valid_mask = (
                (e_t >= 0)
                & (e_c >= 0)
                & (n_t > 0)
                & (n_c > 0)
                & (e_t <= n_t)
                & (e_c <= n_c)
            )
            if not np.any(valid_mask):
                continue
            e_t = e_t[valid_mask]
            n_t = n_t[valid_mask]
            e_c = e_c[valid_mask]
            n_c = n_c[valid_mask]
            clean = clean.loc[clean.index[valid_mask], :]

            k = int(e_t.size)
            if k < args.min_k:
                continue

            valid_file_binary = True

            # Exclude double-zero cells (both arms zero events or all events) — uninformative
            double_zero = ((e_t == 0) & (e_c == 0)) | ((e_t == n_t) & (e_c == n_c))
            if np.any(double_zero):
                keep = ~double_zero
                e_t = e_t[keep]
                n_t = n_t[keep]
                e_c = e_c[keep]
                n_c = n_c[keep]
                clean = clean.loc[clean.index[keep], :]
                k = int(e_t.size)
                if k < args.min_k:
                    continue

            # Haldane-Anscombe continuity correction: add 0.5 ONLY to zero-cell studies
            has_zero = (e_t == 0) | ((n_t - e_t) == 0) | (e_c == 0) | ((n_c - e_c) == 0)
            cc = np.where(has_zero, 0.5, 0.0)
            a = e_t + cc
            b = (n_t - e_t) + cc
            c = e_c + cc
            d = (n_c - e_c) + cc
            yi = np.log((a * d) / (b * c))
            vi = (1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d)

            re_fit = random_effects_dl(yi, vi)
            if re_fit is None:
                continue
            mu_re, se_re, tau2_re = re_fit
            if not (math.isfinite(mu_re) and math.isfinite(se_re) and se_re > 0):
                continue

            se_i = np.sqrt(vi)
            z = yi / se_i
            p_two = 2.0 * (1.0 - normal_cdf(np.abs(z)))
            positive_sig = (p_two < 0.05) & (yi > 0)
            q_sig = float(np.mean(positive_sig))

            review_cov = cov_map.get(review_id, {})
            review_link = linkage_map.get(review_id, {})
            trial_count = safe_float(review_cov.get("trial_count", float("nan")))
            enrollment_mean = safe_float(review_cov.get("enrollment_mean", float("nan")))
            published_weight = float(np.sum(n_t + n_c))
            registry_weight_proxy = (
                trial_count * enrollment_mean
                if math.isfinite(trial_count) and math.isfinite(enrollment_mean) and trial_count > 0 and enrollment_mean > 0
                else float("nan")
            )
            lambda_transport_raw = (
                (published_weight / registry_weight_proxy)
                if (math.isfinite(registry_weight_proxy) and registry_weight_proxy > 0)
                else float("nan")
            )
            lambda_linkage_pmid_raw = safe_float(review_link.get("lambda_completed_pmid_any_weighted", float("nan")))
            lambda_linkage_non_ghost_raw = safe_float(
                review_link.get("lambda_completed_non_ghost_weighted", float("nan"))
            )
            lambda_exact_nct_pmid_raw = safe_float(
                review_link.get("lambda_exact_nct_pmid_any_weighted", float("nan"))
            )
            lambda_exact_nct_non_ghost_raw = safe_float(
                review_link.get("lambda_exact_nct_non_ghost_weighted", float("nan"))
            )
            n_pairwise_nct_ids = safe_float(review_link.get("n_pairwise_nct_ids", float("nan")))
            n_exact_nct_resolved = safe_float(
                review_link.get("n_exact_nct_studies_resolved", float("nan"))
            )
            n_pubmed_bridge_nct_ids = safe_float(
                review_link.get("n_pubmed_bridge_nct_ids", float("nan"))
            )

            # Compute PubMed bridge lambda: combine exact + bridge resolved studies
            n_bridge_non_ghost = safe_float(
                review_link.get("n_pubmed_bridge_completed_non_ghost", float("nan"))
            )
            n_bridge_ghost = safe_float(
                review_link.get("n_pubmed_bridge_completed_ghost", float("nan"))
            )
            if (
                math.isfinite(n_bridge_non_ghost)
                and math.isfinite(n_bridge_ghost)
                and (n_bridge_non_ghost + n_bridge_ghost) > 0
            ):
                lambda_pubmed_bridge_non_ghost_raw = n_bridge_non_ghost / (
                    n_bridge_non_ghost + n_bridge_ghost
                )
            else:
                lambda_pubmed_bridge_non_ghost_raw = float("nan")

            lambda_selected_raw = float("nan")
            lambda_source = "default"
            if args.lambda_source_mode == "transport_proxy":
                if math.isfinite(lambda_transport_raw):
                    lambda_selected_raw = lambda_transport_raw
                    lambda_source = "ctgov_transport_proxy"
            elif args.lambda_source_mode == "linkage_pmid":
                if math.isfinite(lambda_linkage_pmid_raw):
                    lambda_selected_raw = lambda_linkage_pmid_raw
                    lambda_source = "ctgov_linkage_pmid"
            elif args.lambda_source_mode == "linkage_non_ghost":
                if math.isfinite(lambda_linkage_non_ghost_raw):
                    lambda_selected_raw = lambda_linkage_non_ghost_raw
                    lambda_source = "ctgov_linkage_non_ghost"
            elif args.lambda_source_mode == "exact_nct_pmid":
                if math.isfinite(lambda_exact_nct_pmid_raw):
                    lambda_selected_raw = lambda_exact_nct_pmid_raw
                    lambda_source = "pairwise_exact_nct_pmid"
            elif args.lambda_source_mode == "exact_nct_non_ghost":
                if math.isfinite(lambda_exact_nct_non_ghost_raw):
                    lambda_selected_raw = lambda_exact_nct_non_ghost_raw
                    lambda_source = "pairwise_exact_nct_non_ghost"
            elif args.lambda_source_mode == "pubmed_bridge_non_ghost":
                if math.isfinite(lambda_pubmed_bridge_non_ghost_raw):
                    lambda_selected_raw = lambda_pubmed_bridge_non_ghost_raw
                    lambda_source = "pubmed_bridge_non_ghost"
            else:  # auto — 4-tier hierarchy
                if math.isfinite(lambda_exact_nct_non_ghost_raw):
                    lambda_selected_raw = lambda_exact_nct_non_ghost_raw
                    lambda_source = "pairwise_exact_nct_non_ghost"
                elif math.isfinite(lambda_pubmed_bridge_non_ghost_raw):
                    lambda_selected_raw = lambda_pubmed_bridge_non_ghost_raw
                    lambda_source = "pubmed_bridge_non_ghost"
                elif math.isfinite(lambda_linkage_non_ghost_raw):
                    lambda_selected_raw = lambda_linkage_non_ghost_raw
                    lambda_source = "ctgov_linkage_non_ghost"
                elif math.isfinite(lambda_transport_raw):
                    lambda_selected_raw = lambda_transport_raw
                    lambda_source = "ctgov_transport_proxy"

            if math.isfinite(lambda_selected_raw):
                lambda_raw = float(lambda_selected_raw)
            else:
                lambda_raw = float("nan")
                lambda_source = "default"
            lambda_unclipped = lambda_raw if math.isfinite(lambda_raw) else float(args.default_lambda_proxy)
            if not math.isfinite(lambda_unclipped) or (lambda_unclipped <= 0):
                lambda_unclipped = float(args.default_lambda_proxy)

            posterior = {
                "kind": "deterministic",
                "alpha": float("nan"),
                "beta": float("nan"),
                "mean": float("nan"),
                "sd": float("nan"),
                "effective_n": float("nan"),
                "source_success": float("nan"),
                "source_failure": float("nan"),
            }

            if args.lambda_uncertainty_model == "beta":
                posterior = build_lambda_posterior_params(
                    lambda_source=lambda_source,
                    lambda_center=lambda_unclipped,
                    review_link=review_link,
                    prior_alpha=float(args.lambda_beta_prior_alpha),
                    prior_beta=float(args.lambda_beta_prior_beta),
                    pseudo_concentration=float(args.lambda_pseudo_concentration),
                )
                lambda_draws_unclipped = rng_lambda.beta(
                    float(posterior["alpha"]),
                    float(posterior["beta"]),
                    size=int(args.lambda_posterior_draws),
                )
                lambda_draws = (
                    np.clip(lambda_draws_unclipped, args.lambda_proxy_lower, args.lambda_proxy_upper)
                    if clipping_enabled
                    else lambda_draws_unclipped.copy()
                )
                lambda_proxy = float(np.mean(lambda_draws))
                lambda_proxy_unclipped = float(np.mean(lambda_draws_unclipped))
                lambda_draw_mean = float(np.mean(lambda_draws))
                lambda_draw_sd = float(np.std(lambda_draws, ddof=1))
                lambda_draw_mean_unclipped = float(np.mean(lambda_draws_unclipped))
                lambda_draw_sd_unclipped = float(np.std(lambda_draws_unclipped, ddof=1))
            else:
                lambda_proxy = apply_lambda_bounds(
                    lambda_unclipped,
                    lower=args.lambda_proxy_lower,
                    upper=args.lambda_proxy_upper,
                    clipping_enabled=clipping_enabled,
                )
                lambda_proxy_unclipped = float(lambda_unclipped)
                lambda_draw_mean = float("nan")
                lambda_draw_sd = float("nan")
                lambda_draw_mean_unclipped = float("nan")
                lambda_draw_sd_unclipped = float("nan")

            if not math.isfinite(lambda_proxy) or lambda_proxy <= 0:
                lambda_proxy = float(args.default_lambda_proxy)
                lambda_source = "default"
            if not math.isfinite(lambda_proxy_unclipped) or lambda_proxy_unclipped <= 0:
                lambda_proxy_unclipped = float(args.default_lambda_proxy)

            p_nonsig = calibrate_p_nonsig_from_lambda(lambda_proxy, q_sig=q_sig, p_sig=args.p_sig)
            pub_prob = np.where(positive_sig, args.p_sig, p_nonsig)
            ipw_fit = random_effects_dl(yi, vi, weight_multiplier=(1.0 / np.maximum(pub_prob, 1e-9)))
            if ipw_fit is None:
                mu_ipw = float("nan")
                se_ipw = float("nan")
                tau2_ipw = float("nan")
            else:
                mu_ipw, se_ipw, tau2_ipw = ipw_fit

            if k >= args.min_k_pet:
                pet_fit = pet_peese_estimate(yi, vi)
            else:
                pet_fit = None
            if pet_fit is None:
                mu_pet = float("nan")
                se_pet = float("nan")
            else:
                mu_pet, se_pet = pet_fit

            # --- GRMA (Grey Relational Meta-Analysis) ---
            min_k_grma = 3
            if k >= min_k_grma:
                try:
                    grma_obj = GRMA(effect_guard=True)
                    grma_fit = grma_obj.fit(yi, vi)
                    mu_grma = float(grma_fit["estimate"])
                    # Jackknife SE: leave-one-out
                    loo_ests = np.empty(k)
                    for jj in range(k):
                        idx_loo = np.concatenate([np.arange(jj), np.arange(jj + 1, k)])
                        loo_fit = grma_obj._core(yi[idx_loo], vi[idx_loo])
                        loo_ests[jj] = loo_fit["estimate"]
                    loo_mean = np.mean(loo_ests)
                    se_grma = float(np.sqrt((k - 1) / k * np.sum((loo_ests - loo_mean) ** 2)))
                    w_max_grma = float(grma_fit["w_max"])
                    n_eff_grma = float(grma_fit["n_eff"])
                except Exception:
                    mu_grma = float("nan")
                    se_grma = float("nan")
                    w_max_grma = float("nan")
                    n_eff_grma = float("nan")
            else:
                mu_grma = float("nan")
                se_grma = float("nan")
                w_max_grma = float("nan")
                n_eff_grma = float("nan")

            # --- Bayesian GWAM ---
            # Use the published weight from this analysis's studies and estimated ghost weight
            # from the lambda proxy. lambda = w_pub / w_total => w_total = w_pub / lambda
            # w_ghost = w_total - w_pub = w_pub * (1/lambda - 1)
            bayesian_prior = PriorSpec()
            if lambda_proxy > 0 and lambda_proxy < 1:
                w_pub_bayesian = published_weight
                w_total_bayesian = w_pub_bayesian / lambda_proxy
                w_ghost_bayesian = w_total_bayesian - w_pub_bayesian
                # Split ghost weight across hypothetical ghost studies (use 1 ghost "pool")
                w_ghost_arr = np.array([w_ghost_bayesian]) if w_ghost_bayesian > 0 else np.array([], dtype=float)
                bayesian_result = bayesian_gwam_posterior(
                    w_pub=w_pub_bayesian,
                    mu_pub=mu_re,
                    se_pub=se_re,
                    w_results_only_individual=np.array([], dtype=float),
                    w_ghost_individual=w_ghost_arr,
                    w_total=w_total_bayesian,
                    prior=bayesian_prior,
                    ghost_sigma=0.1,
                    results_only_sigma=0.1,
                )
                mu_bayesian_gwam = bayesian_result.posterior_mean
                se_bayesian_gwam = bayesian_result.posterior_sd
                ci_bayesian_gwam_lo = bayesian_result.cri_lo
                ci_bayesian_gwam_hi = bayesian_result.cri_hi
            else:
                mu_bayesian_gwam = float(mu_re)
                se_bayesian_gwam = float(se_re)
                z_crit = normal_quantile(0.975)
                ci_bayesian_gwam_lo = mu_bayesian_gwam - z_crit * se_bayesian_gwam
                ci_bayesian_gwam_hi = mu_bayesian_gwam + z_crit * se_bayesian_gwam

            mu_gwam_proxy = float(lambda_proxy * mu_re)
            se_gwam_proxy = float(lambda_proxy * se_re)
            mu_gwam_proxy_unclipped = float(lambda_proxy_unclipped * mu_re)
            se_gwam_proxy_unclipped = float(lambda_proxy_unclipped * se_re)
            abs_mu_re = abs(float(mu_re))
            if abs_mu_re >= 0.10:
                if args.lambda_uncertainty_model == "beta":
                    cross_prob_010 = float(np.mean(np.abs(lambda_draws * mu_re) < 0.10))
                    cross_prob_010_unclipped = float(np.mean(np.abs(lambda_draws_unclipped * mu_re) < 0.10))
                else:
                    cross_prob_010 = float(abs(mu_gwam_proxy) < 0.10)
                    cross_prob_010_unclipped = float(abs(mu_gwam_proxy_unclipped) < 0.10)
            else:
                cross_prob_010 = 0.0
                cross_prob_010_unclipped = 0.0

            if abs_mu_re >= 0.20:
                if args.lambda_uncertainty_model == "beta":
                    cross_prob_020 = float(np.mean(np.abs(lambda_draws * mu_re) < 0.20))
                    cross_prob_020_unclipped = float(np.mean(np.abs(lambda_draws_unclipped * mu_re) < 0.20))
                else:
                    cross_prob_020 = float(abs(mu_gwam_proxy) < 0.20)
                    cross_prob_020_unclipped = float(abs(mu_gwam_proxy_unclipped) < 0.20)
            else:
                cross_prob_020 = 0.0
                cross_prob_020_unclipped = 0.0

            analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
            subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
            analysis_group = (
                str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
            )
            analysis_number = (
                str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
            )
            review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""

            row: dict[str, Any] = {
                "dataset": file_path.stem,
                "review_id": review_id,
                "review_doi": review_doi,
                "analysis_key": analysis_key,
                "analysis_group": analysis_group,
                "analysis_number": analysis_number,
                "analysis_name": analysis_name,
                "subgroup": subgroup,
                "k": k,
                "published_weight_n_total": published_weight,
                "mu_re": float(mu_re),
                "se_re": float(se_re),
                "tau2_re": float(tau2_re),
                "mu_ipw_re": float(mu_ipw) if math.isfinite(mu_ipw) else float("nan"),
                "se_ipw_re": float(se_ipw) if math.isfinite(se_ipw) else float("nan"),
                "tau2_ipw_re": float(tau2_ipw) if math.isfinite(tau2_ipw) else float("nan"),
                "mu_pet_peese": float(mu_pet) if math.isfinite(mu_pet) else float("nan"),
                "se_pet_peese": float(se_pet) if math.isfinite(se_pet) else float("nan"),
                "mu_grma": float(mu_grma) if math.isfinite(mu_grma) else float("nan"),
                "se_grma": float(se_grma) if math.isfinite(se_grma) else float("nan"),
                "w_max_grma": float(w_max_grma) if math.isfinite(w_max_grma) else float("nan"),
                "n_eff_grma": float(n_eff_grma) if math.isfinite(n_eff_grma) else float("nan"),
                "q_sig_positive": q_sig,
                "p_sig_assumed": float(args.p_sig),
                "p_nonsig_calibrated": float(p_nonsig),
                "lambda_proxy_source": lambda_source,
                "lambda_source_mode_requested": args.lambda_source_mode,
                "lambda_uncertainty_model": args.lambda_uncertainty_model,
                "lambda_clipping_enabled": bool(clipping_enabled),
                "lambda_proxy_raw": float(lambda_raw) if math.isfinite(lambda_raw) else float("nan"),
                "lambda_proxy": float(lambda_proxy),
                "lambda_proxy_unclipped": float(lambda_proxy_unclipped),
                "lambda_posterior_kind": str(posterior["kind"]),
                "lambda_posterior_alpha": float(posterior["alpha"]),
                "lambda_posterior_beta": float(posterior["beta"]),
                "lambda_posterior_mean": float(posterior["mean"]),
                "lambda_posterior_sd": float(posterior["sd"]),
                "lambda_posterior_effective_n": float(posterior["effective_n"]),
                "lambda_posterior_source_success": float(posterior["source_success"]),
                "lambda_posterior_source_failure": float(posterior["source_failure"]),
                "lambda_draw_mean": float(lambda_draw_mean),
                "lambda_draw_sd": float(lambda_draw_sd),
                "lambda_draw_mean_unclipped": float(lambda_draw_mean_unclipped),
                "lambda_draw_sd_unclipped": float(lambda_draw_sd_unclipped),
                "lambda_transport_raw": (
                    float(lambda_transport_raw) if math.isfinite(lambda_transport_raw) else float("nan")
                ),
                "lambda_linkage_pmid_raw": (
                    float(lambda_linkage_pmid_raw) if math.isfinite(lambda_linkage_pmid_raw) else float("nan")
                ),
                "lambda_linkage_non_ghost_raw": (
                    float(lambda_linkage_non_ghost_raw) if math.isfinite(lambda_linkage_non_ghost_raw) else float("nan")
                ),
                "lambda_exact_nct_pmid_raw": (
                    float(lambda_exact_nct_pmid_raw) if math.isfinite(lambda_exact_nct_pmid_raw) else float("nan")
                ),
                "lambda_exact_nct_non_ghost_raw": (
                    float(lambda_exact_nct_non_ghost_raw)
                    if math.isfinite(lambda_exact_nct_non_ghost_raw)
                    else float("nan")
                ),
                "n_pairwise_nct_ids": (
                    float(n_pairwise_nct_ids) if math.isfinite(n_pairwise_nct_ids) else float("nan")
                ),
                "n_exact_nct_studies_resolved": (
                    float(n_exact_nct_resolved) if math.isfinite(n_exact_nct_resolved) else float("nan")
                ),
                "ctgov_trial_count": float(trial_count) if math.isfinite(trial_count) else float("nan"),
                "ctgov_enrollment_mean": float(enrollment_mean) if math.isfinite(enrollment_mean) else float("nan"),
                "ctgov_registry_weight_proxy": (
                    float(registry_weight_proxy) if math.isfinite(registry_weight_proxy) else float("nan")
                ),
                "mu_gwam_proxy": mu_gwam_proxy,
                "se_gwam_proxy": se_gwam_proxy,
                "mu_gwam_proxy_unclipped": mu_gwam_proxy_unclipped,
                "se_gwam_proxy_unclipped": se_gwam_proxy_unclipped,
                "mu_bayesian_gwam": float(mu_bayesian_gwam),
                "se_bayesian_gwam": float(se_bayesian_gwam),
                "ci_bayesian_gwam_lo": float(ci_bayesian_gwam_lo),
                "ci_bayesian_gwam_hi": float(ci_bayesian_gwam_hi),
                "cross_prob_010": cross_prob_010,
                "cross_prob_020": cross_prob_020,
                "cross_prob_010_unclipped": cross_prob_010_unclipped,
                "cross_prob_020_unclipped": cross_prob_020_unclipped,
                "re_ge_010_indicator": bool(abs_mu_re >= 0.10),
                "re_ge_020_indicator": bool(abs_mu_re >= 0.20),
            }

            for lam in lambda_grid:
                suffix = lambda_column_suffix(lam)
                row[f"mu_gwam_lambda_{suffix}"] = float(lam * mu_re)
                row[f"se_gwam_lambda_{suffix}"] = float(lam * se_re)

            z_ci = normal_quantile(0.975)
            row["re_ci_lo"] = float(mu_re - (z_ci * se_re))
            row["re_ci_hi"] = float(mu_re + (z_ci * se_re))
            row["gwam_proxy_ci_lo"] = float(row["mu_gwam_proxy"] - (z_ci * row["se_gwam_proxy"]))
            row["gwam_proxy_ci_hi"] = float(row["mu_gwam_proxy"] + (z_ci * row["se_gwam_proxy"]))
            row["re_positive_sig"] = bool(row["re_ci_lo"] > 0)
            row["gwam_proxy_positive_sig"] = bool(row["gwam_proxy_ci_lo"] > 0)
            row["ipw_positive_sig"] = bool(
                math.isfinite(row["mu_ipw_re"])
                and math.isfinite(row["se_ipw_re"])
                and ((row["mu_ipw_re"] - z_ci * row["se_ipw_re"]) > 0)
            )
            row["pet_positive_sig"] = bool(
                math.isfinite(row["mu_pet_peese"])
                and math.isfinite(row["se_pet_peese"])
                and ((row["mu_pet_peese"] - z_ci * row["se_pet_peese"]) > 0)
            )
            row["grma_positive_sig"] = bool(
                math.isfinite(row["mu_grma"])
                and math.isfinite(row["se_grma"])
                and row["se_grma"] > 0
                and ((row["mu_grma"] - z_ci * row["se_grma"]) > 0)
            )

            all_rows.append(row)

        if valid_file_binary:
            files_with_valid_binary += 1

    if not all_rows:
        raise RuntimeError("No valid binary analyses were estimable.")

    results = pd.DataFrame(all_rows)
    results = results.sort_values(["review_id", "dataset", "analysis_key"]).reset_index(drop=True)
    results_path = output_dir / "analysis_results.csv"
    results.to_csv(results_path, index=False)

    top_shifted = results.copy()
    top_shifted["abs_gwam_shift_proxy"] = (top_shifted["mu_re"] - top_shifted["mu_gwam_proxy"]).abs()
    top_shifted = top_shifted.sort_values("abs_gwam_shift_proxy", ascending=False).head(200)
    top_shifted_path = output_dir / "top_shifted_analyses.csv"
    top_shifted.to_csv(top_shifted_path, index=False)

    n_total = int(len(results))
    n_pet = int(np.sum(np.isfinite(results["mu_pet_peese"])))
    n_ipw = int(np.sum(np.isfinite(results["mu_ipw_re"])))
    n_re_pos = int(np.sum(results["re_positive_sig"]))
    n_gwam_pos = int(np.sum(results["gwam_proxy_positive_sig"]))
    n_drop = int(np.sum(results["re_positive_sig"] & (~results["gwam_proxy_positive_sig"])))

    attenuation_metrics = build_attenuation_metrics(
        results,
        mu_gwam_col="mu_gwam_proxy",
        cross_prob_010_col="cross_prob_010",
        cross_prob_020_col="cross_prob_020",
    )
    attenuation_metrics_unclipped = build_attenuation_metrics(
        results,
        mu_gwam_col="mu_gwam_proxy_unclipped",
        cross_prob_010_col="cross_prob_010_unclipped",
        cross_prob_020_col="cross_prob_020_unclipped",
    )
    distribution_metrics = build_distribution_metrics(
        results,
        lambda_proxy_lower=args.lambda_proxy_lower,
        lambda_proxy_upper=args.lambda_proxy_upper,
        clipping_enabled=clipping_enabled,
    )

    by_lambda_source: dict[str, Any] = {}
    for source in sorted(results["lambda_proxy_source"].dropna().unique()):
        source_df = results.loc[results["lambda_proxy_source"] == source].copy()
        by_lambda_source[str(source)] = build_subset_summary(
            source_df,
            lambda_proxy_lower=args.lambda_proxy_lower,
            lambda_proxy_upper=args.lambda_proxy_upper,
            clipping_enabled=clipping_enabled,
        )

    by_k_band = {
        "k_eq_2": build_subset_summary(
            results.loc[results["k"] == 2].copy(),
            lambda_proxy_lower=args.lambda_proxy_lower,
            lambda_proxy_upper=args.lambda_proxy_upper,
            clipping_enabled=clipping_enabled,
        ),
        "k_eq_3": build_subset_summary(
            results.loc[results["k"] == 3].copy(),
            lambda_proxy_lower=args.lambda_proxy_lower,
            lambda_proxy_upper=args.lambda_proxy_upper,
            clipping_enabled=clipping_enabled,
        ),
        "k_ge_4": build_subset_summary(
            results.loc[results["k"] >= 4].copy(),
            lambda_proxy_lower=args.lambda_proxy_lower,
            lambda_proxy_upper=args.lambda_proxy_upper,
            clipping_enabled=clipping_enabled,
        ),
    }

    robust_threshold = 4.0
    robust_df = results.loc[np.abs(results["mu_re"]) <= robust_threshold].copy()
    n_excluded = int(len(results) - len(robust_df))
    full_subset = build_subset_summary(
        results,
        lambda_proxy_lower=args.lambda_proxy_lower,
        lambda_proxy_upper=args.lambda_proxy_upper,
        clipping_enabled=clipping_enabled,
    )
    robust_subset = build_subset_summary(
        robust_df,
        lambda_proxy_lower=args.lambda_proxy_lower,
        lambda_proxy_upper=args.lambda_proxy_upper,
        clipping_enabled=clipping_enabled,
    )
    review_cluster_summary, review_cluster_table = build_review_cluster_summary(
        results,
        bootstrap_iters=int(args.review_bootstrap_iters),
        bootstrap_seed=int(args.review_bootstrap_seed),
    )
    robust_delta = {
        "crossed_below_0p10_pct_all_filtered_minus_full": (
            robust_subset["attenuation"]["re_to_gwam_proxy_crossed_below_0p10_pct_all"]
            - full_subset["attenuation"]["re_to_gwam_proxy_crossed_below_0p10_pct_all"]
        ),
        "crossed_below_0p20_pct_all_filtered_minus_full": (
            robust_subset["attenuation"]["re_to_gwam_proxy_crossed_below_0p20_pct_all"]
            - full_subset["attenuation"]["re_to_gwam_proxy_crossed_below_0p20_pct_all"]
        ),
        "median_abs_shift_filtered_minus_full": (
            robust_subset["distribution"]["median_abs_shift_proxy"]
            - full_subset["distribution"]["median_abs_shift_proxy"]
        ),
    }

    gwam_grid_sig: dict[str, int] = {}
    for lam in lambda_grid:
        suffix = lambda_column_suffix(lam)
        mu_col = f"mu_gwam_lambda_{suffix}"
        se_col = f"se_gwam_lambda_{suffix}"
        ci_lo = results[mu_col] - (normal_quantile(0.975) * results[se_col])
        gwam_grid_sig[str(lam)] = int(np.sum(ci_lo > 0))

    summary = {
        "metadata": build_environment_metadata(),
        "pairwise_data_dir": str(args.pairwise_data_dir),
        "ctgov_covariates_csv": str(args.ctgov_covariates_csv),
        "ctgov_linkage_csv": str(args.ctgov_linkage_csv),
        "n_files_scanned": int(len(files)),
        "n_files_with_valid_binary_analyses": int(files_with_valid_binary),
        "n_analyses_estimable": n_total,
        "n_pet_peese_estimable": n_pet,
        "n_grma_estimable": int(np.sum(np.isfinite(results["mu_grma"]))),
        "n_selection_ipw_estimable": n_ipw,
        "lambda_grid": lambda_grid,
        "assumptions": {
            "min_k": args.min_k,
            "min_k_pet": args.min_k_pet,
            "default_lambda_proxy": args.default_lambda_proxy,
            "lambda_proxy_lower": args.lambda_proxy_lower,
            "lambda_proxy_upper": args.lambda_proxy_upper,
            "lambda_source_mode": args.lambda_source_mode,
            "p_sig": args.p_sig,
            "lambda_uncertainty_model": args.lambda_uncertainty_model,
            "lambda_posterior_draws": int(args.lambda_posterior_draws),
            "lambda_beta_prior_alpha": float(args.lambda_beta_prior_alpha),
            "lambda_beta_prior_beta": float(args.lambda_beta_prior_beta),
            "lambda_pseudo_concentration": float(args.lambda_pseudo_concentration),
            "lambda_clipping_enabled": bool(clipping_enabled),
            "review_bootstrap_iters": int(args.review_bootstrap_iters),
            "ipw_p_nonsig_calibration": "Solved from (q_sig*p_sig + (1-q_sig)*p_nonsig = lambda_proxy), then clipped.",
            "gwam_mode": "Sensitivity grid + proxy lambda shrinkage on RE estimate.",
            "gwam_significance_note": (
                "GWAM proxy scales both mean and SE by lambda; significance-drop counts are diagnostic and expected to be invariant."
            ),
            "proxy_linkage_note": (
                "Auto lambda now prefers exact Pairwise70 NCT linkage (when available), then query-derived CT.gov linkage/proxy."
            ),
        },
        "headline_counts": {
            **attenuation_metrics,
        },
        "headline_counts_unclipped_lambda_sensitivity": {
            **attenuation_metrics_unclipped,
        },
        "diagnostic_significance_invariance": {
            "re_positive_significant": n_re_pos,
            "gwam_proxy_positive_significant": n_gwam_pos,
            "re_to_gwam_proxy_significance_drops": n_drop,
            "re_to_gwam_proxy_drop_rate_within_re_positive": (
                float(n_drop / n_re_pos) if n_re_pos > 0 else float("nan")
            ),
            "all_equal_re_vs_gwam_positive_significant": bool(
                np.all(results["re_positive_sig"] == results["gwam_proxy_positive_sig"])
            ),
        },
        "clustered_by_review": review_cluster_summary,
        "distribution": distribution_metrics,
        "bound_free_lambda_diagnostics": {
            "n_lambda_proxy_gt_1": int(np.sum(results["lambda_proxy"] > 1.0)),
            "pct_lambda_proxy_gt_1": float(np.mean(results["lambda_proxy"] > 1.0)),
            "n_lambda_proxy_unclipped_gt_1": int(np.sum(results["lambda_proxy_unclipped"] > 1.0)),
            "pct_lambda_proxy_unclipped_gt_1": float(np.mean(results["lambda_proxy_unclipped"] > 1.0)),
            "n_lambda_proxy_unclipped_le_0": int(np.sum(results["lambda_proxy_unclipped"] <= 0.0)),
            "pct_lambda_proxy_unclipped_le_0": float(np.mean(results["lambda_proxy_unclipped"] <= 0.0)),
        },
        "stratified": {
            "by_lambda_proxy_source": by_lambda_source,
            "by_k_band": by_k_band,
        },
        "robust_sensitivity_excluding_abs_mu_re_gt_4": {
            "threshold_abs_mu_re": robust_threshold,
            "n_excluded": n_excluded,
            "excluded_pct": float(n_excluded / n_total),
            "full": full_subset,
            "filtered": robust_subset,
            "delta_filtered_minus_full": robust_delta,
        },
        "grma_comparison": _build_grma_summary(results),
        "gwam_grid_positive_significant_counts": gwam_grid_sig,
        "outputs": {
            "analysis_results_csv": str(results_path),
            "top_shifted_analyses_csv": str(top_shifted_path),
        },
    }

    summary_path = output_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    strat_lambda_rows = [
        flatten_subset_row(name, subset) for name, subset in summary["stratified"]["by_lambda_proxy_source"].items()
    ]
    strat_k_rows = [
        flatten_subset_row(name, subset) for name, subset in summary["stratified"]["by_k_band"].items()
    ]
    strat_lambda_path = output_dir / "stratified_by_lambda_source.csv"
    strat_k_path = output_dir / "stratified_by_k_band.csv"
    pd.DataFrame(strat_lambda_rows).to_csv(strat_lambda_path, index=False)
    pd.DataFrame(strat_k_rows).to_csv(strat_k_path, index=False)
    review_cluster_path = output_dir / "clustered_by_review.csv"
    review_cluster_table.to_csv(review_cluster_path, index=False)

    robust_path = output_dir / "robust_sensitivity_excluding_abs_mu_re_gt_4.json"
    robust_path.write_text(
        json.dumps(summary["robust_sensitivity_excluding_abs_mu_re_gt_4"], indent=2),
        encoding="utf-8",
    )
    summary["outputs"]["stratified_by_lambda_source_csv"] = str(strat_lambda_path)
    summary["outputs"]["stratified_by_k_band_csv"] = str(strat_k_path)
    summary["outputs"]["robust_sensitivity_json"] = str(robust_path)
    summary["outputs"]["clustered_by_review_csv"] = str(review_cluster_path)

    # Rewrite summary after output paths are appended.
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Wrote {results_path}")
    print(f"Wrote {top_shifted_path}")
    print(f"Wrote {strat_lambda_path}")
    print(f"Wrote {strat_k_path}")
    print(f"Wrote {review_cluster_path}")
    print(f"Wrote {robust_path}")
    print(f"Wrote {summary_path}")
    print(
        f"n_analyses={n_total}, median_lambda_proxy={summary['distribution']['median_lambda_proxy']:.3f}, "
        f"crossed_0p10={summary['headline_counts']['re_to_gwam_proxy_crossed_below_0p10_pct_all']:.3f}, "
        f"crossed_0p20={summary['headline_counts']['re_to_gwam_proxy_crossed_below_0p20_pct_all']:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
