#!/usr/bin/env python3
"""Run GWAM corrections from a registry CSV.

This script supports two layers:
1) Deterministic correction: mu_gwam = lambda * mu_published when ghost mean = 0.
2) Sensitivity simulation: ghost study effects sampled from Normal(mu_ghost, sd_ghost).

Input registry CSV is expected from scripts/fetch_ctgov_registry.py.

NOTE: Weights are based on enrollment counts (weight_n column), NOT inverse-variance.
This is by design: GWAM uses registry enrollment as the "universe" denominator for
the integrity ratio lambda, since variance estimates are unavailable for ghost protocols
and many results-only studies. This design choice is discussed in the manuscript.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from gwam_utils import (
    build_environment_metadata,
    gwam_correct,
    parse_bool,
    partition_registry_weights,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-csv", type=Path, required=True, help="Trial-level registry CSV.")
    parser.add_argument(
        "--published-mu",
        type=float,
        required=True,
        help="Published pooled effect size (e.g., Hedges g from literature meta-analysis).",
    )
    parser.add_argument(
        "--weight-column",
        default="weight_n",
        help="Column to use as study weights (default: weight_n).",
    )
    parser.add_argument("--sim-n", type=int, default=5000, help="Monte Carlo simulation count (max 10_000_000).")
    parser.add_argument("--ghost-mu", type=float, default=0.0, help="Mean effect for ghost studies.")
    parser.add_argument("--ghost-sd", type=float, default=0.10, help="SD for ghost effect prior/sensitivity.")
    parser.add_argument(
        "--results-only-mode",
        choices=["as_unknown", "as_observed"],
        default="as_unknown",
        help="How to treat studies with no PMID but has posted registry results.",
    )
    parser.add_argument(
        "--results-only-mu",
        type=float,
        default=0.0,
        help="Mean effect used when --results-only-mode=as_unknown.",
    )
    parser.add_argument(
        "--results-only-sd",
        type=float,
        default=0.10,
        help="SD used when --results-only-mode=as_unknown.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("data/analysis/gwam_summary.json"),
        help="Output JSON summary path.",
    )
    return parser.parse_args()


def quantile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values, q, method="linear"))


def main() -> int:
    args = parse_args()
    if args.sim_n < 2:
        raise ValueError(f"--sim-n={args.sim_n} must be at least 2 for meaningful Monte Carlo.")
    if args.sim_n > 10_000_000:
        raise ValueError(f"--sim-n={args.sim_n} exceeds maximum of 10,000,000.")
    args.output_json.parent.mkdir(parents=True, exist_ok=True)

    with args.registry_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError("Registry CSV is empty.")
    if args.weight_column not in rows[0]:
        raise ValueError(f"Weight column '{args.weight_column}' not found in CSV.")

    # Diagnostic-only counters that the reusable API does not track. The
    # weight partitioning and lambda/mu arithmetic below is delegated to
    # gwam_utils.gwam_correct() so the CLI and the importable API stay in sync.
    n_ghost_no_pmid_results = 0
    for row in rows:
        if parse_bool(row.get("is_ghost_protocol", "")):
            if (not parse_bool(row.get("has_pmid", ""))) and parse_bool(
                row.get("has_results", "")
            ):
                n_ghost_no_pmid_results += 1

    result = gwam_correct(
        rows,
        published_mu=float(args.published_mu),
        weight_column=args.weight_column,
        ghost_mu=args.ghost_mu,
        results_only_mode=args.results_only_mode,
        results_only_mu=args.results_only_mu,
    )

    # Per-study weight lists (original registry order) for the Monte-Carlo
    # simulation path below. partition_registry_weights() preserves the same
    # iteration order the summation loop used, so simulation draws are
    # unchanged.
    _pmid_weights, results_only_weights, ghost_weights = partition_registry_weights(
        rows, weight_column=args.weight_column
    )
    w_pmid = result.weight_pmid_only
    w_results_only = result.weight_results_only
    w_ghost = result.weight_ghost
    w_total = result.weight_total
    n_with_pmid = result.n_pmid_only
    n_results_only = result.n_results_only

    lambda_pmid_only = result.lambda_pmid_only
    lambda_non_ghost = result.lambda_non_ghost
    if lambda_pmid_only == 0:
        print("WARNING: lambda_pmid_only=0 -- all studies are ghost/results-only. GWAM correction is undefined.")
    if lambda_non_ghost == 0:
        print("WARNING: lambda_non_ghost=0 -- all studies are ghosts. GWAM correction is undefined.")
    mu_published = result.mu_published
    mu_gwam_null = result.mu_gwam_null_point

    rng = np.random.default_rng(args.seed)
    mu_sim = np.full(args.sim_n, (w_pmid * mu_published) / w_total, dtype=float)

    if args.results_only_mode == "as_observed":
        mu_sim += (w_results_only * mu_published) / w_total
    elif results_only_weights:
        rw = np.asarray(results_only_weights, dtype=float)
        if args.results_only_sd > 0:
            ro_effects = rng.normal(
                loc=args.results_only_mu,
                scale=args.results_only_sd,
                size=(args.sim_n, rw.size),
            )
            mu_sim += (ro_effects @ rw) / w_total
        else:
            mu_sim += (w_results_only * args.results_only_mu) / w_total

    if ghost_weights:
        gw = np.asarray(ghost_weights, dtype=float)
        if args.ghost_sd > 0:
            ghost_effects = rng.normal(loc=args.ghost_mu, scale=args.ghost_sd, size=(args.sim_n, gw.size))
            mu_sim += (ghost_effects @ gw) / w_total
        else:
            mu_sim += (w_ghost * args.ghost_mu) / w_total

    summary = {
        "metadata": build_environment_metadata(),
        "registry_csv": str(args.registry_csv),
        "n_trials_total": len(rows),
        "n_trials_with_pmid": n_with_pmid,
        "n_trials_results_only_no_pmid": n_results_only,
        "n_trials_ghost_with_results_no_pmid": n_ghost_no_pmid_results,
        "n_trials_observed_non_ghost_proxy": n_with_pmid + n_results_only,
        "n_trials_pmid_only": result.n_pmid_only,
        "n_trials_results_only": result.n_results_only,
        "n_trials_ghost_proxy": result.n_ghost,
        "weights": {
            "pmid_only": w_pmid,
            "results_only_no_pmid": w_results_only,
            "observed_non_ghost": w_pmid + w_results_only,
            "ghost": w_ghost,
            "total": w_total,
            "integrity_ratio_lambda_pmid_only": lambda_pmid_only,
            "integrity_ratio_lambda_non_ghost": lambda_non_ghost,
        },
        "inputs": {
            "published_mu": mu_published,
            "sim_n": args.sim_n,
            "results_only_mode": args.results_only_mode,
            "results_only_mu": args.results_only_mu,
            "results_only_sd": args.results_only_sd,
            "ghost_mu": args.ghost_mu,
            "ghost_sd": args.ghost_sd,
            "seed": args.seed,
        },
        "estimates": {
            "mu_random_effects_published_only": mu_published,
            "mu_gwam_null_point": mu_gwam_null,
            "mu_gwam_sim_mean": float(np.mean(mu_sim)),
            "mu_gwam_sim_sd": float(np.std(mu_sim, ddof=1)),
            "mu_gwam_sim_q025": quantile(mu_sim, 0.025),
            "mu_gwam_sim_q975": quantile(mu_sim, 0.975),
        },
    }

    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Wrote {args.output_json}")
    print(
        f"lambda_pmid_only={lambda_pmid_only:.3f}, "
        f"lambda_non_ghost={lambda_non_ghost:.3f}, "
        f"published_mu={mu_published:.3f}, "
        f"gwam_null_mu={mu_gwam_null:.3f}, "
        f"sim_mean={summary['estimates']['mu_gwam_sim_mean']:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
