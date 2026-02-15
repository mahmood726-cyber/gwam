#!/usr/bin/env python3
"""Bayesian hierarchical GWAM correction with analytical Normal conjugacy posterior.

Replaces the fixed ghost=0 point mass with a proper hierarchical model:

Level 1: theta_gi | mu_G, sigma_G ~ N(mu_G, sigma_G^2)   [individual ghost effects]
Level 2: mu_G ~ N(m_G, s_mG^2)                            [hyperprior on ghost mean]
         sigma_G from sensitivity grid                      [within-ghost heterogeneity]

Similarly for results-only studies without observed effects.

Because ghost effects are unobserved, the posterior for mu_GWAM is analytically
tractable (Normal conjugacy):

  mu_post = (w_pub * mu_pub + sum(w_Rj * theta_Rj) + w_G * m_G) / w_total
  sigma_post^2 = (w_pub^2 * se_pub^2
                  + sum(w_Rj_obs^2 * se_Rj^2)          [observed RO: measurement SE]
                  + sum(w_Rj_unobs^2) * sigma_R^2       [unobserved RO: Level 1]
                  + w_R_unobs^2 * s_mR^2                [unobserved RO: Level 2]
                  + sum(w_Gg^2) * sigma_G^2 + w_G^2 * s_mG^2
                 ) / w_total^2

When --use-observed-results is set, results-only studies with finite effects
contribute their actual SE (measurement uncertainty) rather than sigma_R.
This gives principled CrI without empirical multipliers.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from gwam_utils import build_environment_metadata, normal_cdf, normal_quantile, parse_bool


@dataclass
class PriorSpec:
    """Specification for hierarchical priors on ghost and results-only effects."""

    # Ghost hyperprior: mu_G ~ N(m_G, s_mG^2)
    # NOTE: Default mu=0 assumes ghost protocols have null effects on average.
    # This is the conservative assumption for the GWAM correction (unpublished
    # trials are assumed to show no benefit). Sensitivity to this assumption
    # is explored via the sigma grid. The manuscript must justify this choice.
    ghost_mu_mean: float = 0.0
    ghost_mu_sd: float = 0.15

    # Within-ghost heterogeneity sensitivity grid
    ghost_sigma_grid: list[float] = field(
        default_factory=lambda: [0.05, 0.1, 0.2, 0.3]
    )

    # Results-only hyperprior: mu_R ~ N(m_R, s_mR^2)
    results_only_mu_mean: float = 0.0
    results_only_mu_sd: float = 0.15

    # Within-results-only heterogeneity sensitivity grid
    results_only_sigma_grid: list[float] = field(
        default_factory=lambda: [0.05, 0.1, 0.2, 0.3]
    )


@dataclass
class BayesianResult:
    """Result of a single Bayesian GWAM posterior computation."""

    posterior_mean: float
    posterior_sd: float
    cri_lo: float  # credible interval lower bound
    cri_hi: float  # credible interval upper bound
    pr_positive: float  # Pr(mu_GWAM > 0)
    ghost_sigma: float
    results_only_sigma: float


def bayesian_gwam_posterior(
    *,
    w_pub: float,
    mu_pub: float,
    se_pub: float,
    w_results_only_individual: np.ndarray,
    w_ghost_individual: np.ndarray,
    w_total: float,
    prior: PriorSpec,
    ghost_sigma: float,
    results_only_sigma: float,
    results_only_observed: np.ndarray | None = None,
    results_only_observed_se: np.ndarray | None = None,
    alpha: float = 0.05,
) -> BayesianResult:
    """Compute the analytical Bayesian GWAM posterior for one prior configuration.

    Parameters
    ----------
    w_pub : total weight of published (PMID-linked) studies
    mu_pub : published pooled effect estimate
    se_pub : SE of published pooled effect
    w_results_only_individual : per-study weights for results-only studies
    w_ghost_individual : per-study weights for ghost studies
    w_total : total weight across all strata
    prior : PriorSpec with hyperprior parameters
    ghost_sigma : within-ghost heterogeneity (Level 1 SD)
    results_only_sigma : within-results-only heterogeneity (Level 1 SD)
    results_only_observed : optional array of observed effects for results-only studies
        (same length as w_results_only_individual); NaN entries fall back to prior
    results_only_observed_se : optional array of SEs for observed results-only studies
        (same length as w_results_only_individual); NaN entries treated as exact
    alpha : significance level for credible interval (default 0.05 -> 95% CrI)

    Returns
    -------
    BayesianResult with posterior mean, SD, credible interval, and Pr(positive).
    """
    if w_total <= 0:
        raise ValueError("w_total must be > 0")

    if results_only_observed is not None and w_results_only_individual.size > 0:
        if results_only_observed.size != w_results_only_individual.size:
            raise ValueError(
                f"results_only_observed length ({results_only_observed.size}) "
                f"must match w_results_only_individual length ({w_results_only_individual.size})"
            )

    w_R = float(np.sum(w_results_only_individual))
    w_G = float(np.sum(w_ghost_individual))

    # --- Posterior mean ---
    pub_mean_component = w_pub * mu_pub

    # Results-only: use observed effect when available, prior mean otherwise
    if results_only_observed is not None and w_results_only_individual.size > 0:
        ro_effects = np.where(
            np.isfinite(results_only_observed),
            results_only_observed,
            prior.results_only_mu_mean,
        )
        ro_mean_component = float(np.sum(w_results_only_individual * ro_effects))
    else:
        ro_mean_component = w_R * prior.results_only_mu_mean

    # Ghost: always prior (unobserved)
    ghost_mean_component = w_G * prior.ghost_mu_mean

    mu_post = (pub_mean_component + ro_mean_component + ghost_mean_component) / w_total

    # --- Posterior variance ---
    # Published: w_pub^2 * se_pub^2
    pub_var = (w_pub * se_pub) ** 2

    # Results-only uncertainty
    if results_only_observed is not None and w_results_only_individual.size > 0:
        observed_mask = np.isfinite(results_only_observed)
        # Observed studies: contribute their own SE² (measurement uncertainty)
        # If SE is unknown (NaN or not provided), treat as exact (zero variance)
        ro_obs_var = 0.0
        if observed_mask.any():
            w_obs = w_results_only_individual[observed_mask]
            if results_only_observed_se is not None:
                se_obs = results_only_observed_se[observed_mask]
                # Only add variance for finite SEs
                finite_se = np.where(np.isfinite(se_obs), se_obs, 0.0)
                ro_obs_var = float(np.sum((w_obs * finite_se) ** 2))
        # Unobserved: Level 1 (within-study) + Level 2 (hyperprior)
        w_unobs = w_results_only_individual[~observed_mask]
        ro_level1_var = float(np.sum(w_unobs ** 2)) * (results_only_sigma ** 2) + ro_obs_var
        w_R_unobs = float(np.sum(w_unobs))
        ro_level2_var = (w_R_unobs ** 2) * (prior.results_only_mu_sd ** 2)
    else:
        # All results-only are unobserved
        ro_level1_var = float(np.sum(w_results_only_individual ** 2)) * (results_only_sigma ** 2)
        ro_level2_var = (w_R ** 2) * (prior.results_only_mu_sd ** 2)

    ro_var = ro_level1_var + ro_level2_var

    # Ghost uncertainty (always unobserved)
    ghost_level1_var = float(np.sum(w_ghost_individual ** 2)) * (ghost_sigma ** 2)
    ghost_level2_var = (w_G ** 2) * (prior.ghost_mu_sd ** 2)
    ghost_var = ghost_level1_var + ghost_level2_var

    sigma_post_sq = (pub_var + ro_var + ghost_var) / (w_total ** 2)
    sigma_post = math.sqrt(max(sigma_post_sq, 0.0))

    # --- Credible interval ---
    z = normal_quantile(1.0 - alpha / 2.0)
    cri_lo = mu_post - z * sigma_post
    cri_hi = mu_post + z * sigma_post

    # --- Pr(mu_GWAM > 0) ---
    if sigma_post > 0:
        pr_positive = normal_cdf(mu_post / sigma_post)
    else:
        pr_positive = 1.0 if mu_post > 0 else (0.0 if mu_post < 0 else 0.5)

    return BayesianResult(
        posterior_mean=mu_post,
        posterior_sd=sigma_post,
        cri_lo=cri_lo,
        cri_hi=cri_hi,
        pr_positive=pr_positive,
        ghost_sigma=ghost_sigma,
        results_only_sigma=results_only_sigma,
    )


def run_sensitivity_grid(
    *,
    w_pub: float,
    mu_pub: float,
    se_pub: float,
    w_results_only_individual: np.ndarray,
    w_ghost_individual: np.ndarray,
    w_total: float,
    prior: PriorSpec,
    results_only_observed: np.ndarray | None = None,
    results_only_observed_se: np.ndarray | None = None,
    alpha: float = 0.05,
) -> list[BayesianResult]:
    """Run Bayesian GWAM across the full sensitivity grid of sigma values."""
    if not prior.ghost_sigma_grid or not prior.results_only_sigma_grid:
        raise ValueError("Sensitivity grid requires at least one value in each sigma grid.")
    results: list[BayesianResult] = []
    for ghost_sigma in prior.ghost_sigma_grid:
        for ro_sigma in prior.results_only_sigma_grid:
            result = bayesian_gwam_posterior(
                w_pub=w_pub,
                mu_pub=mu_pub,
                se_pub=se_pub,
                w_results_only_individual=w_results_only_individual,
                w_ghost_individual=w_ghost_individual,
                w_total=w_total,
                prior=prior,
                ghost_sigma=ghost_sigma,
                results_only_sigma=ro_sigma,
                results_only_observed=results_only_observed,
                results_only_observed_se=results_only_observed_se,
                alpha=alpha,
            )
            results.append(result)
    return results


def load_registry_weights(
    path: Path,
    weight_column: str,
    *,
    collect_observed: bool = False,
) -> tuple[float, np.ndarray, np.ndarray, float, dict[str, int],
           np.ndarray | None, np.ndarray | None]:
    """Load registry CSV and return weight components.

    Parameters
    ----------
    collect_observed : if True, also extract ``results_effect`` and
        ``results_se`` columns for results-only studies (avoids a second CSV
        read).

    Returns
    -------
    (w_pub, w_results_only_individual, w_ghost_individual, w_total, counts,
     results_only_observed, results_only_observed_se)

    ``results_only_observed`` and ``results_only_observed_se`` are ``None``
    when ``collect_observed`` is False or no effects are present.
    """
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise ValueError("Registry CSV is empty.")
    if weight_column not in rows[0]:
        raise ValueError(f"Weight column '{weight_column}' not found in CSV.")

    pmid_weights: list[float] = []
    results_only_weights: list[float] = []
    ghost_weights: list[float] = []
    ro_effects: list[float] = []
    ro_ses: list[float] = []

    for row_idx, row in enumerate(rows, start=2):  # row 2 = first data row (after header)
        raw_weight = row[weight_column]
        try:
            weight = float(raw_weight)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Non-numeric weight '{raw_weight}' in row {row_idx}, column '{weight_column}'."
            ) from exc
        if not math.isfinite(weight) or weight <= 0:
            raise ValueError(
                f"Invalid weight {weight} in row {row_idx}, column '{weight_column}'."
            )
        is_ghost = parse_bool(row.get("is_ghost_protocol", ""))
        has_pmid = parse_bool(row.get("has_pmid", ""))

        if is_ghost:
            ghost_weights.append(weight)
        elif has_pmid:
            pmid_weights.append(weight)
        else:
            results_only_weights.append(weight)
            if collect_observed:
                effect_str = row.get("results_effect", "")
                se_str = row.get("results_se", "")
                try:
                    ro_effects.append(float(effect_str) if effect_str else float("nan"))
                except (TypeError, ValueError):
                    ro_effects.append(float("nan"))
                try:
                    ro_ses.append(float(se_str) if se_str else float("nan"))
                except (TypeError, ValueError):
                    ro_ses.append(float("nan"))

    w_pub = float(np.sum(pmid_weights))
    w_ro = np.asarray(results_only_weights, dtype=float)
    w_g = np.asarray(ghost_weights, dtype=float)
    w_total = w_pub + float(np.sum(w_ro)) + float(np.sum(w_g))

    counts = {
        "n_total": len(rows),
        "n_with_pmid": len(pmid_weights),
        "n_results_only": len(results_only_weights),
        "n_ghost": len(ghost_weights),
    }

    ro_obs: np.ndarray | None = None
    ro_obs_se: np.ndarray | None = None
    if collect_observed and ro_effects:
        ro_obs = np.asarray(ro_effects, dtype=float)
        ro_obs_se = np.asarray(ro_ses, dtype=float)
        # If all NaN, treat as None
        if not np.any(np.isfinite(ro_obs)):
            ro_obs = None
            ro_obs_se = None

    return w_pub, w_ro, w_g, w_total, counts, ro_obs, ro_obs_se


def parse_float_list(text: str) -> list[float]:
    """Parse comma-separated floats (rejects NaN/Inf)."""
    out: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if token:
            val = float(token)
            if not math.isfinite(val):
                raise ValueError(f"Non-finite value in float list: {token!r}")
            out.append(val)
    return out


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
        "--published-se",
        type=float,
        required=True,
        help="SE of the published pooled effect size.",
    )
    parser.add_argument(
        "--weight-column",
        default="weight_n",
        help="Column to use as study weights (default: weight_n).",
    )
    parser.add_argument(
        "--prior-ghost-mu-mean",
        type=float,
        default=0.0,
        help="Hyperprior mean for ghost effect mean (default: 0.0).",
    )
    parser.add_argument(
        "--prior-ghost-mu-sd",
        type=float,
        default=0.15,
        help="Hyperprior SD for ghost effect mean (default: 0.15).",
    )
    parser.add_argument(
        "--prior-ghost-sigma-grid",
        default="0.05,0.1,0.2,0.3",
        help="Comma-separated within-ghost heterogeneity values for sensitivity grid.",
    )
    parser.add_argument(
        "--prior-ro-mu-mean",
        type=float,
        default=0.0,
        help="Hyperprior mean for results-only effect mean (default: 0.0).",
    )
    parser.add_argument(
        "--prior-ro-mu-sd",
        type=float,
        default=0.15,
        help="Hyperprior SD for results-only effect mean (default: 0.15).",
    )
    parser.add_argument(
        "--prior-ro-sigma-grid",
        default="0.05,0.1,0.2,0.3",
        help="Comma-separated within-results-only heterogeneity values for sensitivity grid.",
    )
    parser.add_argument(
        "--use-observed-results",
        action="store_true",
        help="Use observed effects from results_effect column for results-only studies.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance level for credible interval (default: 0.05 -> 95%% CrI).",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("data/analysis/gwam_bayesian_summary.json"),
        help="Output JSON summary path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_json.parent.mkdir(parents=True, exist_ok=True)

    prior = PriorSpec(
        ghost_mu_mean=args.prior_ghost_mu_mean,
        ghost_mu_sd=args.prior_ghost_mu_sd,
        ghost_sigma_grid=parse_float_list(args.prior_ghost_sigma_grid),
        results_only_mu_mean=args.prior_ro_mu_mean,
        results_only_mu_sd=args.prior_ro_mu_sd,
        results_only_sigma_grid=parse_float_list(args.prior_ro_sigma_grid),
    )

    if not math.isfinite(args.published_se) or args.published_se <= 0:
        raise ValueError("--published-se must be a finite positive number.")
    if not (0 < args.alpha < 1):
        raise ValueError("--alpha must be in (0, 1).")
    if not prior.ghost_sigma_grid:
        raise ValueError("--prior-ghost-sigma-grid must have at least one value.")
    if any(v < 0 for v in prior.ghost_sigma_grid):
        raise ValueError("--prior-ghost-sigma-grid values must be >= 0.")
    if not prior.results_only_sigma_grid:
        raise ValueError("--prior-ro-sigma-grid must have at least one value.")
    if any(v < 0 for v in prior.results_only_sigma_grid):
        raise ValueError("--prior-ro-sigma-grid values must be >= 0.")

    w_pub, w_ro, w_g, w_total, counts, ro_obs, ro_obs_se = load_registry_weights(
        args.registry_csv, args.weight_column,
        collect_observed=args.use_observed_results,
    )

    if not math.isfinite(w_total) or w_total <= 0:
        raise ValueError(f"Total weight is invalid ({w_total}); check input data.")

    # Run sensitivity grid
    grid_results = run_sensitivity_grid(
        w_pub=w_pub,
        mu_pub=args.published_mu,
        se_pub=args.published_se,
        w_results_only_individual=w_ro,
        w_ghost_individual=w_g,
        w_total=w_total,
        prior=prior,
        results_only_observed=ro_obs,
        results_only_observed_se=ro_obs_se,
        alpha=args.alpha,
    )

    # Find reference result (median SD across grid)
    sds = [r.posterior_sd for r in grid_results]
    median_sd = float(np.median(sds))
    reference_idx = int(np.argmin([abs(r.posterior_sd - median_sd) for r in grid_results]))
    reference = grid_results[reference_idx]

    # Build sensitivity table
    sensitivity_table: list[dict[str, float]] = []
    for r in grid_results:
        sensitivity_table.append({
            "ghost_sigma": r.ghost_sigma,
            "results_only_sigma": r.results_only_sigma,
            "posterior_mean": r.posterior_mean,
            "posterior_sd": r.posterior_sd,
            "cri_lo": r.cri_lo,
            "cri_hi": r.cri_hi,
            "cri_width": r.cri_hi - r.cri_lo,
            "pr_positive": r.pr_positive,
        })

    lambda_pmid_only = w_pub / w_total
    lambda_non_ghost = (w_pub + float(np.sum(w_ro))) / w_total

    n_observed_results = 0
    if ro_obs is not None:
        n_observed_results = int(np.sum(np.isfinite(ro_obs)))

    summary = {
        "metadata": build_environment_metadata(),
        "registry_csv": str(args.registry_csv),
        "n_trials_total": counts["n_total"],
        "n_trials_with_pmid": counts["n_with_pmid"],
        "n_trials_results_only": counts["n_results_only"],
        "n_trials_ghost": counts["n_ghost"],
        "n_results_only_with_observed_effect": n_observed_results,
        "weights": {
            "pmid_only": w_pub,
            "results_only": float(np.sum(w_ro)),
            "ghost": float(np.sum(w_g)),
            "total": w_total,
            "integrity_ratio_lambda_pmid_only": lambda_pmid_only,
            "integrity_ratio_lambda_non_ghost": lambda_non_ghost,
        },
        "inputs": {
            "published_mu": args.published_mu,
            "published_se": args.published_se,
            "alpha": args.alpha,
            "use_observed_results": args.use_observed_results,
        },
        "prior_spec": {
            "ghost_mu_mean": prior.ghost_mu_mean,
            "ghost_mu_sd": prior.ghost_mu_sd,
            "ghost_sigma_grid": prior.ghost_sigma_grid,
            "results_only_mu_mean": prior.results_only_mu_mean,
            "results_only_mu_sd": prior.results_only_mu_sd,
            "results_only_sigma_grid": prior.results_only_sigma_grid,
        },
        "bayesian_gwam": {
            "reference_posterior_mean": reference.posterior_mean,
            "reference_posterior_sd": reference.posterior_sd,
            "reference_cri_lo": reference.cri_lo,
            "reference_cri_hi": reference.cri_hi,
            "reference_cri_width": reference.cri_hi - reference.cri_lo,
            "reference_pr_positive": reference.pr_positive,
            "reference_ghost_sigma": reference.ghost_sigma,
            "reference_results_only_sigma": reference.results_only_sigma,
            "grid_posterior_mean_min": float(min(r.posterior_mean for r in grid_results)),
            "grid_posterior_mean_max": float(max(r.posterior_mean for r in grid_results)),
            "grid_posterior_sd_min": float(min(r.posterior_sd for r in grid_results)),
            "grid_posterior_sd_max": float(max(r.posterior_sd for r in grid_results)),
            "grid_cri_width_min": float(min(r.cri_hi - r.cri_lo for r in grid_results)),
            "grid_cri_width_max": float(max(r.cri_hi - r.cri_lo for r in grid_results)),
        },
        "sensitivity_table": sensitivity_table,
    }

    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Wrote {args.output_json}")
    print(
        f"lambda_pmid_only={lambda_pmid_only:.3f}, "
        f"lambda_non_ghost={lambda_non_ghost:.3f}, "
        f"published_mu={args.published_mu:.3f}, "
        f"bayesian_gwam_mean={reference.posterior_mean:.4f}, "
        f"bayesian_gwam_sd={reference.posterior_sd:.4f}, "
        f"cri=[{reference.cri_lo:.4f}, {reference.cri_hi:.4f}], "
        f"pr_positive={reference.pr_positive:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
