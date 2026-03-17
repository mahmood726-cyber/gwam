#!/usr/bin/env python3
"""Compare GWAM vs published-only random effects under publication bias.

The simulation uses real registry-derived study weights from one intervention
domain. Publication bias is modeled by selective publication probabilities.

Note: The IPW comparator uses *oracle* (true) publication probabilities, not
estimated ones. This provides an idealized upper bound for what IPW could
achieve if the selection function were perfectly known.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from scipy.special import ndtr as _array_normal_cdf

from gwam_utils import build_environment_metadata, normal_quantile, parse_bool
from model_gwam_bayesian import PriorSpec, bayesian_gwam_posterior


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry-csv", type=Path, required=True, help="Input registry CSV.")
    parser.add_argument(
        "--weight-column",
        default="weight_n",
        help="Study-weight source column in registry CSV.",
    )
    parser.add_argument(
        "--mu-true-list",
        default="0.0,0.1,0.2",
        help="Comma-separated true pooled effects to simulate.",
    )
    parser.add_argument("--tau2", type=float, default=0.1, help="Between-study variance.")
    parser.add_argument("--n-meta", type=int, default=5000, help="Replicates per scenario.")
    parser.add_argument(
        "--p-sig",
        type=float,
        default=0.95,
        help="Publication probability for statistically significant positive results.",
    )
    parser.add_argument(
        "--p-nonsig",
        type=float,
        default=0.20,
        help="Publication probability for non-significant/negative results.",
    )
    parser.add_argument(
        "--calibrate-nonsig",
        action="store_true",
        help="Calibrate p-nonsig to match observed registry lambda (weight-level).",
    )
    parser.add_argument(
        "--calibration-runs",
        type=int,
        default=1200,
        help="Replicates used per lambda-calibration candidate evaluation.",
    )
    parser.add_argument(
        "--calibration-grid-size",
        type=int,
        default=25,
        help="Candidate count per p-nonsig calibration grid pass.",
    )
    parser.add_argument(
        "--calibration-refine-rounds",
        type=int,
        default=3,
        help="Local refinement rounds around the best p-nonsig candidate.",
    )
    parser.add_argument(
        "--calibration-max-abs-error",
        type=float,
        default=0.03,
        help="Advisory max absolute lambda-calibration error reported in output diagnostics.",
    )
    parser.add_argument(
        "--calibration-enforce-tolerance",
        action="store_true",
        help="Fail run when lambda calibration absolute error exceeds --calibration-max-abs-error.",
    )
    parser.add_argument(
        "--winsorize-pct",
        type=float,
        default=99.0,
        help="Winsorize weight distribution at this percentile (100 disables).",
    )
    parser.add_argument(
        "--control-event-rate",
        type=float,
        default=0.35,
        help="Assumed control-group event rate for binary log(OR) simulation.",
    )
    parser.add_argument(
        "--control-event-sd-logit",
        type=float,
        default=0.35,
        help="Between-study SD of baseline control risk on logit scale (0 disables).",
    )
    parser.add_argument(
        "--allocation-min-frac-treat",
        type=float,
        default=0.45,
        help="Minimum treatment-arm fraction used to sample per-study allocation ratio.",
    )
    parser.add_argument(
        "--allocation-max-frac-treat",
        type=float,
        default=0.55,
        help="Maximum treatment-arm fraction used to sample per-study allocation ratio.",
    )
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
    parser.add_argument(
        "--ghost-sigma",
        type=float,
        default=0.10,
        help="Assumed SD of ghost effects for GWAM SE computation (0 disables ghost variance).",
    )
    parser.add_argument(
        "--lambda-target",
        choices=["pmid_only", "non_ghost"],
        default="non_ghost",
        help="Lambda definition used for calibration summary keys.",
    )
    parser.add_argument(
        "--calibration-objective",
        choices=["auto", "pmid_only", "non_ghost"],
        default="auto",
        help=(
            "Objective lambda used to calibrate p-nonsig. "
            "auto uses pmid_only when results-only pool exists, otherwise follows --lambda-target."
        ),
    )
    parser.add_argument(
        "--ci-calibration",
        choices=["none", "empirical"],
        default="empirical",
        help="How to calibrate CI width under selection bias in simulation summaries.",
    )
    parser.add_argument(
        "--ci-target-coverage",
        type=float,
        default=0.95,
        help="Target two-sided coverage used for empirical CI calibration.",
    )
    parser.add_argument(
        "--ci-calibration-runs",
        type=int,
        default=400,
        help="Independent replicates used to estimate empirical CI multipliers.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--calibration-mu",
        type=float,
        default=None,
        help="Mu value used to calibrate p_nonsig. Defaults to first mu_true_list value.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("data/analysis/simulation_summary.json"),
        help="Output summary JSON path.",
    )
    return parser.parse_args()


def logistic(x: np.ndarray) -> np.ndarray:
    from scipy.special import expit  # numerically stable sigmoid
    return expit(x)


def random_effects_dl(
    y: np.ndarray,
    v: np.ndarray,
    weight_multiplier: np.ndarray | None = None,
) -> tuple[float, float, float] | None:
    k = y.size
    if k < 2:
        return None

    if weight_multiplier is None:
        mult = np.ones_like(v, dtype=float)
    else:
        mult = np.asarray(weight_multiplier, dtype=float)
        if mult.shape != v.shape:
            return None
        if np.any(mult <= 0):
            return None

    if np.any(v <= 0):
        return None

    w = mult / v
    w_sum = np.sum(w)
    if w_sum <= 0:
        return None

    mu_fixed = float(np.sum(w * y) / w_sum)
    q = float(np.sum(w * (y - mu_fixed) ** 2))
    c = float(w_sum - (np.sum(w**2) / w_sum))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0

    w_star = mult / (v + tau2)
    w_star_sum = np.sum(w_star)
    if w_star_sum <= 0:
        return None

    mu_re = float(np.sum(w_star * y) / w_star_sum)
    se_re = float(math.sqrt(1.0 / w_star_sum))
    return mu_re, se_re, tau2


def weighted_linear_intercept(
    *,
    y: np.ndarray,
    x: np.ndarray,
    w: np.ndarray,
) -> tuple[float, float] | None:
    if y.size < 3:
        return None
    if not (y.size == x.size == w.size):
        return None
    if np.any(w <= 0):
        return None

    xmat = np.column_stack([np.ones_like(x), x])
    wx = xmat * w[:, None]
    xtwx = xmat.T @ wx
    xtwy = xmat.T @ (w * y)
    try:
        beta = np.linalg.solve(xtwx, xtwy)
        xtwx_inv = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        return None

    resid = y - (xmat @ beta)
    dof = max(1, y.size - 2)
    sigma2 = float(np.sum(w * (resid**2)) / dof)
    cov = xtwx_inv * sigma2
    se0 = float(math.sqrt(max(cov[0, 0], 0.0)))
    return float(beta[0]), se0


def pet_peese_estimate(y: np.ndarray, v: np.ndarray) -> tuple[float, float] | None:
    if y.size < 3:
        return None
    if np.any(v <= 0):
        return None
    w = 1.0 / v
    se = np.sqrt(v)

    pet = weighted_linear_intercept(y=y, x=se, w=w)
    if pet is None:
        return None
    pet_mu, pet_se = pet
    if pet_se <= 0:
        return None

    # Standard PET-PEESE decision rule (Stanley & Doucouliagos 2014):
    # If the PET intercept is significant at the two-sided alpha=0.10 level
    # (i.e., |z| > 1.645), switch to the PEESE model; otherwise keep PET.
    use_peese = abs(pet_mu / pet_se) > 1.645
    if use_peese:
        peese = weighted_linear_intercept(y=y, x=v, w=w)
        if peese is None:
            return pet_mu, pet_se
        return peese
    return pet_mu, pet_se


def simulate_one_meta(
    *,
    weights_selection: np.ndarray,
    results_only_weights: np.ndarray,
    total_weight: float,
    mu_true: float,
    tau2: float,
    p_sig: float,
    p_nonsig: float,
    control_event_rate: float,
    control_event_sd_logit: float,
    allocation_min_frac_treat: float,
    allocation_max_frac_treat: float,
    ghost_sigma: float = 0.10,
    results_only_mode: str,
    results_only_mu: float,
    results_only_sd: float,
    rng: np.random.Generator,
    max_attempts: int = 50,
) -> dict[str, float] | None:
    n_total = np.maximum(np.round(weights_selection).astype(int), 4)
    frac_treat = rng.uniform(allocation_min_frac_treat, allocation_max_frac_treat, size=n_total.size)
    n_t = np.maximum(np.round(n_total * frac_treat).astype(int), 1)
    n_t = np.minimum(n_t, n_total - 1)
    n_t = np.maximum(n_t, 1)
    n_c = np.maximum(n_total - n_t, 1)
    logit_p0 = math.log(control_event_rate / (1.0 - control_event_rate))
    k = weights_selection.size
    w_results_only = float(np.sum(results_only_weights))
    all_weight_sum = float(total_weight)

    for _ in range(max_attempts):
        theta = rng.normal(loc=mu_true, scale=math.sqrt(tau2), size=k)
        if control_event_sd_logit > 0:
            logit_p0_i = rng.normal(loc=logit_p0, scale=control_event_sd_logit, size=k)
            p0_i = logistic(logit_p0_i)
        else:
            logit_p0_i = np.full(k, logit_p0, dtype=float)
            p0_i = np.full(k, control_event_rate, dtype=float)
        p_t = logistic(logit_p0_i + theta)

        events_t = rng.binomial(n_t, p_t)
        events_c = rng.binomial(n_c, p0_i)
        non_events_t = n_t - events_t
        non_events_c = n_c - events_c

        # Note: Double-zero studies (both arms zero events, or both arms all events)
        # are NOT excluded here, unlike the benchmark script.  Under the simulation's
        # parameter range (control_event_rate >= 0.10, per-arm n >= 2), such studies
        # are statistically negligible, and the continuity correction makes them finite.
        #
        # Haldane-Anscombe continuity correction: apply only to studies with
        # at least one zero cell to avoid unnecessary bias toward the null.
        has_zero = (events_t == 0) | (non_events_t == 0) | (events_c == 0) | (non_events_c == 0)
        cc = np.where(has_zero, 0.5, 0.0)
        a = events_t + cc
        b = non_events_t + cc
        c = events_c + cc
        d = non_events_c + cc
        observed = np.log((a * d) / (b * c))
        variances = (1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d)
        ses = np.sqrt(variances)
        z_vals = observed / ses
        p_vals = 2.0 * (1.0 - _array_normal_cdf(np.abs(z_vals)))
        positive_sig = (p_vals < 0.05) & (observed > 0)
        pub_probs = np.where(positive_sig, p_sig, p_nonsig)
        is_published = rng.random(k) < pub_probs
        if np.sum(is_published) >= 2:
            break
    else:
        return None

    y_pub = observed[is_published]
    v_pub = variances[is_published]
    w_pub = weights_selection[is_published]
    re = random_effects_dl(y_pub, v_pub)
    if re is None:
        return None

    mu_re, se_re, _ = re
    # Oracle IPW-RE: uses the *true* (known) publication probabilities for
    # inverse-probability weighting.  In practice these are unknown; this serves
    # as an idealized upper-bound comparator for GWAM.  The SE from
    # random_effects_dl uses the standard 1/sqrt(sum(w*)) formula, which does
    # not account for the IPW multiplier and may underestimate true variance.
    # A sandwich estimator would be more appropriate but is omitted here for
    # comparability with the standard DL framework.
    pub_positive_sig = positive_sig[is_published]
    pub_prob = np.where(pub_positive_sig, p_sig, p_nonsig)
    pub_prob = np.maximum(pub_prob, 1e-6)
    ipw_re = random_effects_dl(y_pub, v_pub, weight_multiplier=(1.0 / pub_prob))
    if ipw_re is None:
        mu_ipw = float("nan")
        se_ipw = float("nan")
    else:
        mu_ipw, se_ipw, _ = ipw_re
    pet_peese = pet_peese_estimate(y_pub, v_pub)
    if pet_peese is None:
        mu_pet = float("nan")
        se_pet = float("nan")
    else:
        mu_pet, se_pet = pet_peese

    w_pub_sum = float(np.sum(w_pub))
    lambda_pmid_only = w_pub_sum / all_weight_sum if all_weight_sum > 0 else 0.0
    lambda_non_ghost = (w_pub_sum + w_results_only) / all_weight_sum if all_weight_sum > 0 else 0.0

    if results_only_mode == "as_observed":
        gwam_scale = lambda_non_ghost
        mu_gwam = gwam_scale * mu_re
        # Ghost SE: per-study squared weights (consistent with Bayesian code)
        ghost_w_ind = weights_selection[~is_published]
        ghost_se_sq = float(np.sum((ghost_w_ind / all_weight_sum) ** 2) * (ghost_sigma ** 2))
        se_gwam = math.sqrt((gwam_scale * se_re) ** 2 + ghost_se_sq)
    else:
        mu_results_component = 0.0
        se_results_component_sq = 0.0
        if results_only_weights.size:
            if results_only_sd > 0:
                ro_effects = rng.normal(results_only_mu, results_only_sd, size=results_only_weights.size)
                mu_results_component = float(np.sum(ro_effects * results_only_weights) / all_weight_sum)
                se_results_component_sq = float(
                    np.sum(((results_only_weights / all_weight_sum) ** 2) * (results_only_sd**2))
                )
            else:
                mu_results_component = float(np.sum(results_only_weights) * results_only_mu / all_weight_sum)

        # Ghost SE: use per-study squared weights (consistent with Bayesian code)
        ghost_weights_individual = weights_selection[~is_published]
        se_ghost_component_sq = float(
            np.sum((ghost_weights_individual / all_weight_sum) ** 2) * (ghost_sigma ** 2)
        )

        mu_gwam = (lambda_pmid_only * mu_re) + mu_results_component
        se_gwam = math.sqrt((lambda_pmid_only * se_re) ** 2 + se_results_component_sq + se_ghost_component_sq)

    # --- Bayesian GWAM (5th comparator) ---
    # Compute analytical posterior with default priors
    bayesian_prior = PriorSpec()
    bayesian_result = bayesian_gwam_posterior(
        w_pub=w_pub_sum,
        mu_pub=mu_re,
        se_pub=se_re,
        w_results_only_individual=results_only_weights,
        w_ghost_individual=weights_selection[~is_published],
        w_total=all_weight_sum,
        prior=bayesian_prior,
        ghost_sigma=ghost_sigma,
        results_only_sigma=ghost_sigma,
    )
    mu_bayesian_gwam = bayesian_result.posterior_mean
    se_bayesian_gwam = bayesian_result.posterior_sd

    z_crit = normal_quantile(0.975)  # two-sided 95% CI
    ci_re_lo = mu_re - z_crit * se_re
    ci_re_hi = mu_re + z_crit * se_re
    ci_gwam_lo = mu_gwam - z_crit * se_gwam
    ci_gwam_hi = mu_gwam + z_crit * se_gwam
    ci_ipw_lo = mu_ipw - z_crit * se_ipw
    ci_ipw_hi = mu_ipw + z_crit * se_ipw
    ci_pet_lo = mu_pet - z_crit * se_pet
    ci_pet_hi = mu_pet + z_crit * se_pet
    ci_bayesian_gwam_lo = bayesian_result.cri_lo
    ci_bayesian_gwam_hi = bayesian_result.cri_hi

    return {
        "lambda_sim_pmid_only": lambda_pmid_only,
        "lambda_sim_non_ghost": lambda_non_ghost,
        "mu_re": mu_re,
        "se_re": se_re,
        "ci_re_lo": ci_re_lo,
        "ci_re_hi": ci_re_hi,
        "mu_gwam": mu_gwam,
        "se_gwam": se_gwam,
        "ci_gwam_lo": ci_gwam_lo,
        "ci_gwam_hi": ci_gwam_hi,
        "mu_oracle_ipw_re": mu_ipw,
        "se_oracle_ipw_re": se_ipw,
        "ci_oracle_ipw_re_lo": ci_ipw_lo,
        "ci_oracle_ipw_re_hi": ci_ipw_hi,
        "mu_pet_peese": mu_pet,
        "se_pet_peese": se_pet,
        "ci_pet_peese_lo": ci_pet_lo,
        "ci_pet_peese_hi": ci_pet_hi,
        "mu_bayesian_gwam": mu_bayesian_gwam,
        "se_bayesian_gwam": se_bayesian_gwam,
        "ci_bayesian_gwam_lo": ci_bayesian_gwam_lo,
        "ci_bayesian_gwam_hi": ci_bayesian_gwam_hi,
    }


def summarize_model(est: np.ndarray, ci_lo: np.ndarray, ci_hi: np.ndarray, mu_true: float) -> dict[str, float]:
    valid = np.isfinite(est) & np.isfinite(ci_lo) & np.isfinite(ci_hi)
    if not np.any(valid):
        return {
            "n_valid": 0.0,
            "mean_estimate": float("nan"),
            "bias": float("nan"),
            "rmse": float("nan"),
            "coverage_95": float("nan"),
            "mean_ci_width": float("nan"),
        }
    est_v = est[valid]
    ci_lo_v = ci_lo[valid]
    ci_hi_v = ci_hi[valid]
    err = est_v - mu_true
    return {
        "n_valid": float(np.sum(valid)),
        "mean_estimate": float(np.mean(est_v)),
        "bias": float(np.mean(err)),
        "rmse": float(math.sqrt(np.mean(err**2))),
        "coverage_95": float(np.mean((ci_lo_v <= mu_true) & (ci_hi_v >= mu_true))),
        "mean_ci_width": float(np.mean(ci_hi_v - ci_lo_v)),
    }


def summarize_model_with_multiplier(
    est: np.ndarray,
    se: np.ndarray,
    mu_true: float,
    multiplier: float,
) -> tuple[float, float]:
    valid = np.isfinite(est) & np.isfinite(se) & (se > 0)
    if not np.any(valid):
        return float("nan"), float("nan")
    est_v = est[valid]
    se_v = se[valid]
    ci_lo = est_v - (multiplier * se_v)
    ci_hi = est_v + (multiplier * se_v)
    coverage = float(np.mean((ci_lo <= mu_true) & (ci_hi >= mu_true)))
    mean_width = float(np.mean(ci_hi - ci_lo))
    return coverage, mean_width


def parse_mu_list(text: str) -> list[float]:
    out: list[float] = []
    for token in text.split(","):
        token = token.strip()
        if token:
            out.append(float(token))
    if not out:
        raise ValueError("mu-true-list is empty.")
    return out


def load_registry_weights(
    path: Path,
    weight_column: str,
    winsorize_pct: float,
) -> tuple[np.ndarray, np.ndarray, float, dict[str, float]]:
    with path.open("r", newline="", encoding="utf-8") as _fh:
        rows = list(csv.DictReader(_fh))
    if not rows:
        raise ValueError(f"Registry CSV is empty: {path}")
    if weight_column not in rows[0]:
        raise ValueError(f"Weight column '{weight_column}' not found in {path}.")

    weight_values: list[float] = []
    for row_idx, row in enumerate(rows, start=2):
        raw_val = row[weight_column]
        try:
            w = float(raw_val)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Non-numeric weight '{raw_val}' in row {row_idx}, column '{weight_column}'."
            ) from exc
        if not math.isfinite(w) or w <= 0:
            raise ValueError(
                f"Invalid weight {w} in row {row_idx}, column '{weight_column}'. Must be finite and > 0."
            )
        weight_values.append(w)
    raw_weights = np.asarray(weight_values, dtype=float)

    is_ghost = np.asarray([parse_bool(row.get("is_ghost_protocol", "")) for row in rows], dtype=bool)
    has_pmid = np.asarray([parse_bool(row.get("has_pmid", "")) for row in rows], dtype=bool)
    has_results = np.asarray([parse_bool(row.get("has_results", "")) for row in rows], dtype=bool)

    # Results-only/unknown observed stratum: non-ghost rows lacking PMID.
    # This mirrors model_gwam.py behavior, where non-ghost rows without PMID are
    # not treated as latent ghosts even if has_results is missing/false.
    is_results_only = (~is_ghost) & (~has_pmid)

    weights = raw_weights.copy()
    if winsorize_pct < 100.0:
        cap = float(np.quantile(weights, winsorize_pct / 100.0, method="linear"))
        weights = np.minimum(weights, cap)

    total_weight = float(np.sum(weights))
    if total_weight <= 0:
        raise ValueError("Total weight must be > 0.")

    observed = {
        "pmid_only_raw": float(np.sum(raw_weights[(~is_ghost) & has_pmid]) / np.sum(raw_weights)),
        "non_ghost_raw": float(np.sum(raw_weights[~is_ghost]) / np.sum(raw_weights)),
        "pmid_only_effective": float(np.sum(weights[(~is_ghost) & has_pmid]) / total_weight),
        "non_ghost_effective": float(np.sum(weights[~is_ghost]) / total_weight),
        "results_only_weight_raw": float(np.sum(raw_weights[is_results_only])),
        "results_only_weight_effective": float(np.sum(weights[is_results_only])),
        "selection_weight_effective": float(np.sum(weights[~is_results_only])),
        "n_results_only_trials": float(np.sum(is_results_only)),
        "n_results_only_with_results_trials": float(np.sum(is_results_only & has_results)),
    }
    weights_selection = weights[~is_results_only]
    results_only_weights = weights[is_results_only]
    return weights_selection, results_only_weights, total_weight, observed


def estimate_lambda_from_publication(
    *,
    weights_selection: np.ndarray,
    results_only_weights: np.ndarray,
    total_weight: float,
    lambda_target: str,
    mu_true: float,
    tau2: float,
    p_sig: float,
    p_nonsig: float,
    control_event_rate: float,
    control_event_sd_logit: float,
    allocation_min_frac_treat: float,
    allocation_max_frac_treat: float,
    ghost_sigma: float,
    results_only_mode: str,
    results_only_mu: float,
    results_only_sd: float,
    n_meta: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    lambdas: list[float] = []
    for _ in range(n_meta):
        result = simulate_one_meta(
            weights_selection=weights_selection,
            results_only_weights=results_only_weights,
            total_weight=total_weight,
            mu_true=mu_true,
            tau2=tau2,
            p_sig=p_sig,
            p_nonsig=p_nonsig,
            control_event_rate=control_event_rate,
            control_event_sd_logit=control_event_sd_logit,
            allocation_min_frac_treat=allocation_min_frac_treat,
            allocation_max_frac_treat=allocation_max_frac_treat,
            ghost_sigma=ghost_sigma,
            results_only_mode=results_only_mode,
            results_only_mu=results_only_mu,
            results_only_sd=results_only_sd,
            rng=rng,
        )
        if result:
            if lambda_target == "pmid_only":
                lambdas.append(result["lambda_sim_pmid_only"])
            else:
                lambdas.append(result["lambda_sim_non_ghost"])
    if not lambdas:
        return 0.0
    return float(np.mean(lambdas))


def calibrate_p_nonsig(
    *,
    weights_selection: np.ndarray,
    results_only_weights: np.ndarray,
    total_weight: float,
    lambda_target: str,
    mu_true: float,
    tau2: float,
    p_sig: float,
    control_event_rate: float,
    control_event_sd_logit: float,
    allocation_min_frac_treat: float,
    allocation_max_frac_treat: float,
    ghost_sigma: float,
    target_lambda: float,
    results_only_mode: str,
    results_only_mu: float,
    results_only_sd: float,
    n_probe: int,
    grid_size: int,
    refine_rounds: int,
    seed: int,
) -> dict[str, float | int | bool]:
    low = 0.01
    high = min(0.99, p_sig)
    grid_size = max(5, int(grid_size))
    refine_rounds = max(0, int(refine_rounds))

    cache: dict[float, float] = {}

    def eval_lambda(p_nonsig: float) -> float:
        p = float(min(high, max(low, p_nonsig)))
        key = round(p, 8)
        if key not in cache:
            cache[key] = estimate_lambda_from_publication(
                weights_selection=weights_selection,
                results_only_weights=results_only_weights,
                total_weight=total_weight,
                lambda_target=lambda_target,
                mu_true=mu_true,
                tau2=tau2,
                p_sig=p_sig,
                p_nonsig=p,
                control_event_rate=control_event_rate,
                control_event_sd_logit=control_event_sd_logit,
                allocation_min_frac_treat=allocation_min_frac_treat,
                allocation_max_frac_treat=allocation_max_frac_treat,
                ghost_sigma=ghost_sigma,
                results_only_mode=results_only_mode,
                results_only_mu=results_only_mu,
                results_only_sd=results_only_sd,
                n_meta=n_probe,
                seed=seed,
            )
        return float(cache[key])

    best_p = low
    best_lambda = eval_lambda(low)
    best_err = abs(best_lambda - target_lambda)

    interval_low = low
    interval_high = high
    step_width = (interval_high - interval_low) / (grid_size - 1)

    for _ in range(refine_rounds + 1):
        candidates = np.linspace(interval_low, interval_high, grid_size)
        for candidate in candidates:
            lam = eval_lambda(float(candidate))
            err = abs(lam - target_lambda)
            if err < best_err:
                best_err = err
                best_p = float(candidate)
                best_lambda = lam

        interval_low = max(low, best_p - step_width)
        interval_high = min(high, best_p + step_width)
        if interval_high <= interval_low:
            break
        step_width = (interval_high - interval_low) / (grid_size - 1)

    if cache:
        lambda_range_low = float(min(cache.values()))
        lambda_range_high = float(max(cache.values()))
    else:
        lambda_range_low = 0.0
        lambda_range_high = 0.0

    return {
        "p_nonsig_calibrated": best_p,
        "lambda_at_p_nonsig": best_lambda,
        "abs_error": best_err,
        "lambda_range_low": lambda_range_low,
        "lambda_range_high": lambda_range_high,
        "target_attainable": bool(lambda_range_low <= target_lambda <= lambda_range_high),
        "n_candidates_evaluated": int(len(cache)),
    }


def main() -> int:
    args = parse_args()
    if args.n_meta < 2:
        raise ValueError(f"--n-meta={args.n_meta} must be at least 2 for meaningful Monte Carlo.")
    if args.n_meta > 10_000_000:
        raise ValueError(f"--n-meta={args.n_meta} exceeds maximum of 10,000,000.")
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    if not (0.0 < args.control_event_rate < 1.0):
        raise ValueError("--control-event-rate must be in (0,1).")
    if not (0.0 < args.p_nonsig <= 1.0):
        raise ValueError("--p-nonsig must be in (0,1].")
    if not (0.0 < args.p_sig <= 1.0):
        raise ValueError("--p-sig must be in (0,1].")
    if args.results_only_sd < 0:
        raise ValueError("--results-only-sd must be >= 0.")
    if not (0.5 < args.ci_target_coverage < 1.0):
        raise ValueError("--ci-target-coverage must be in (0.5, 1.0).")
    if args.calibration_runs < 100:
        raise ValueError("--calibration-runs must be >= 100.")
    if args.calibration_grid_size < 5:
        raise ValueError("--calibration-grid-size must be >= 5.")
    if args.calibration_refine_rounds < 0:
        raise ValueError("--calibration-refine-rounds must be >= 0.")
    if args.calibration_max_abs_error < 0:
        raise ValueError("--calibration-max-abs-error must be >= 0.")
    if args.ci_calibration_runs < 100:
        raise ValueError("--ci-calibration-runs must be >= 100.")
    if args.control_event_sd_logit < 0:
        raise ValueError("--control-event-sd-logit must be >= 0.")
    if args.tau2 < 0:
        raise ValueError("--tau2 must be >= 0.")
    if args.ghost_sigma < 0:
        raise ValueError("--ghost-sigma must be >= 0.")
    if args.results_only_sd < 0:
        raise ValueError("--results-only-sd must be >= 0.")
    if not (0 < args.allocation_min_frac_treat < 1):
        raise ValueError("--allocation-min-frac-treat must be in (0,1).")
    if not (0 < args.allocation_max_frac_treat < 1):
        raise ValueError("--allocation-max-frac-treat must be in (0,1).")
    if args.allocation_min_frac_treat > args.allocation_max_frac_treat:
        raise ValueError("--allocation-min-frac-treat cannot exceed --allocation-max-frac-treat.")

    mu_true_values = parse_mu_list(args.mu_true_list)
    calibration_mu = mu_true_values[0] if args.calibration_mu is None else float(args.calibration_mu)
    (
        weights_selection,
        results_only_weights,
        total_weight,
        observed_lambda,
    ) = load_registry_weights(
        args.registry_csv, args.weight_column, args.winsorize_pct
    )
    if weights_selection.size < 2:
        raise ValueError("Need at least two selection-eligible trials after results-only split.")

    if args.lambda_target == "pmid_only":
        lambda_target_raw = observed_lambda["pmid_only_raw"]
        lambda_target_effective = observed_lambda["pmid_only_effective"]
    else:
        lambda_target_raw = observed_lambda["non_ghost_raw"]
        lambda_target_effective = observed_lambda["non_ghost_effective"]

    if args.calibration_objective == "auto":
        calibration_objective = (
            "pmid_only"
            if observed_lambda["n_results_only_trials"] > 0
            else args.lambda_target
        )
    else:
        calibration_objective = args.calibration_objective
    if calibration_objective == "pmid_only":
        calibration_target_effective = observed_lambda["pmid_only_effective"]
    else:
        calibration_target_effective = observed_lambda["non_ghost_effective"]

    p_nonsig_used = args.p_nonsig
    calibration_info: dict[str, float | int | bool] | None = None
    if args.calibrate_nonsig:
        calibration_info = calibrate_p_nonsig(
            weights_selection=weights_selection,
            results_only_weights=results_only_weights,
            total_weight=total_weight,
            lambda_target=calibration_objective,
            mu_true=calibration_mu,
            tau2=args.tau2,
            p_sig=args.p_sig,
            control_event_rate=args.control_event_rate,
            control_event_sd_logit=args.control_event_sd_logit,
            allocation_min_frac_treat=args.allocation_min_frac_treat,
            allocation_max_frac_treat=args.allocation_max_frac_treat,
            ghost_sigma=args.ghost_sigma,
            target_lambda=calibration_target_effective,
            results_only_mode=args.results_only_mode,
            results_only_mu=args.results_only_mu,
            results_only_sd=args.results_only_sd,
            n_probe=args.calibration_runs,
            grid_size=args.calibration_grid_size,
            refine_rounds=args.calibration_refine_rounds,
            seed=args.seed,
        )
        p_nonsig_used = float(calibration_info["p_nonsig_calibrated"])
        if args.calibration_enforce_tolerance and float(calibration_info["abs_error"]) > args.calibration_max_abs_error:
            raise RuntimeError(
                "Lambda calibration error exceeds tolerance: "
                f"{float(calibration_info['abs_error']):.4f} > {args.calibration_max_abs_error:.4f}. "
                "Adjust simulation assumptions or disable --calibration-enforce-tolerance."
            )

    scenario_results: list[dict[str, object]] = []
    scenario_buffers: list[dict[str, object]] = []
    for idx, mu_true in enumerate(mu_true_values):
        rng = np.random.default_rng(args.seed + 1000 + idx)
        sim_rows: list[dict[str, float]] = []
        for _ in range(args.n_meta):
            item = simulate_one_meta(
                weights_selection=weights_selection,
                results_only_weights=results_only_weights,
                total_weight=total_weight,
                mu_true=mu_true,
                tau2=args.tau2,
                p_sig=args.p_sig,
                p_nonsig=p_nonsig_used,
                control_event_rate=args.control_event_rate,
                control_event_sd_logit=args.control_event_sd_logit,
                allocation_min_frac_treat=args.allocation_min_frac_treat,
                allocation_max_frac_treat=args.allocation_max_frac_treat,
                ghost_sigma=args.ghost_sigma,
                results_only_mode=args.results_only_mode,
                results_only_mu=args.results_only_mu,
                results_only_sd=args.results_only_sd,
                rng=rng,
            )
            if item:
                sim_rows.append(item)

        if not sim_rows:
            raise RuntimeError("No successful simulation replicates; adjust parameters.")

        lambda_sim_pmid_only = np.asarray([row["lambda_sim_pmid_only"] for row in sim_rows], dtype=float)
        lambda_sim_non_ghost = np.asarray([row["lambda_sim_non_ghost"] for row in sim_rows], dtype=float)
        mu_re = np.asarray([row["mu_re"] for row in sim_rows], dtype=float)
        se_re = np.asarray([row["se_re"] for row in sim_rows], dtype=float)
        ci_re_lo = np.asarray([row["ci_re_lo"] for row in sim_rows], dtype=float)
        ci_re_hi = np.asarray([row["ci_re_hi"] for row in sim_rows], dtype=float)
        mu_g = np.asarray([row["mu_gwam"] for row in sim_rows], dtype=float)
        se_g = np.asarray([row["se_gwam"] for row in sim_rows], dtype=float)
        ci_g_lo = np.asarray([row["ci_gwam_lo"] for row in sim_rows], dtype=float)
        ci_g_hi = np.asarray([row["ci_gwam_hi"] for row in sim_rows], dtype=float)
        mu_ipw = np.asarray([row["mu_oracle_ipw_re"] for row in sim_rows], dtype=float)
        se_ipw = np.asarray([row["se_oracle_ipw_re"] for row in sim_rows], dtype=float)
        ci_ipw_lo = np.asarray([row["ci_oracle_ipw_re_lo"] for row in sim_rows], dtype=float)
        ci_ipw_hi = np.asarray([row["ci_oracle_ipw_re_hi"] for row in sim_rows], dtype=float)
        mu_pet = np.asarray([row["mu_pet_peese"] for row in sim_rows], dtype=float)
        se_pet = np.asarray([row["se_pet_peese"] for row in sim_rows], dtype=float)
        ci_pet_lo = np.asarray([row["ci_pet_peese_lo"] for row in sim_rows], dtype=float)
        ci_pet_hi = np.asarray([row["ci_pet_peese_hi"] for row in sim_rows], dtype=float)
        mu_bgwam = np.asarray([row["mu_bayesian_gwam"] for row in sim_rows], dtype=float)
        se_bgwam = np.asarray([row["se_bayesian_gwam"] for row in sim_rows], dtype=float)
        ci_bgwam_lo = np.asarray([row["ci_bayesian_gwam_lo"] for row in sim_rows], dtype=float)
        ci_bgwam_hi = np.asarray([row["ci_bayesian_gwam_hi"] for row in sim_rows], dtype=float)

        scenario_buffers.append(
            {
                "mu_true": mu_true,
                "mu_re": mu_re,
                "se_re": se_re,
                "mu_g": mu_g,
                "se_g": se_g,
                "mu_ipw": mu_ipw,
                "se_ipw": se_ipw,
                "mu_pet": mu_pet,
                "se_pet": se_pet,
                "mu_bgwam": mu_bgwam,
                "se_bgwam": se_bgwam,
            }
        )
        scenario_results.append(
            {
                "mu_true": mu_true,
                "n_meta_requested": args.n_meta,
                "n_meta_success": int(len(sim_rows)),
                "lambda_sim_mean_pmid_only": float(np.mean(lambda_sim_pmid_only)),
                "lambda_sim_mean_non_ghost": float(np.mean(lambda_sim_non_ghost)),
                "random_effects": summarize_model(mu_re, ci_re_lo, ci_re_hi, mu_true),
                "gwam_null": summarize_model(mu_g, ci_g_lo, ci_g_hi, mu_true),
                "oracle_ipw_random_effects": summarize_model(mu_ipw, ci_ipw_lo, ci_ipw_hi, mu_true),
                "pet_peese": summarize_model(mu_pet, ci_pet_lo, ci_pet_hi, mu_true),
                "bayesian_gwam": summarize_model(mu_bgwam, ci_bgwam_lo, ci_bgwam_hi, mu_true),
            }
        )

    ci_calibration_summary: dict[str, object] = {
        "mode": args.ci_calibration,
        "target_coverage": args.ci_target_coverage,
        "calibration_mu": calibration_mu,
        "out_of_sample": args.ci_calibration == "empirical",
        "calibration_runs_requested": args.ci_calibration_runs if args.ci_calibration == "empirical" else None,
        "calibration_runs_success": None,
        "scenario_mu_used": None,
        "multiplier_random_effects": None,
        "multiplier_gwam": None,
        "multiplier_oracle_ipw_random_effects": None,
        "multiplier_pet_peese": None,
        "multiplier_bayesian_gwam": None,
    }
    if args.ci_calibration == "empirical":
        cal_rng = np.random.default_rng(args.seed + 500000)
        cal_rows: list[dict[str, float]] = []
        for _ in range(args.ci_calibration_runs):
            item = simulate_one_meta(
                weights_selection=weights_selection,
                results_only_weights=results_only_weights,
                total_weight=total_weight,
                mu_true=calibration_mu,
                tau2=args.tau2,
                p_sig=args.p_sig,
                p_nonsig=p_nonsig_used,
                control_event_rate=args.control_event_rate,
                control_event_sd_logit=args.control_event_sd_logit,
                allocation_min_frac_treat=args.allocation_min_frac_treat,
                allocation_max_frac_treat=args.allocation_max_frac_treat,
                ghost_sigma=args.ghost_sigma,
                results_only_mode=args.results_only_mode,
                results_only_mu=args.results_only_mu,
                results_only_sd=args.results_only_sd,
                rng=cal_rng,
            )
            if item:
                cal_rows.append(item)
        if not cal_rows:
            raise RuntimeError("No successful simulation replicates for CI calibration batch.")

        cal_mu_true = float(calibration_mu)
        re_cal = np.asarray([row["mu_re"] for row in cal_rows], dtype=float)
        se_re_cal = np.asarray([row["se_re"] for row in cal_rows], dtype=float)
        g_cal = np.asarray([row["mu_gwam"] for row in cal_rows], dtype=float)
        se_g_cal = np.asarray([row["se_gwam"] for row in cal_rows], dtype=float)
        ipw_cal = np.asarray([row["mu_oracle_ipw_re"] for row in cal_rows], dtype=float)
        se_ipw_cal = np.asarray([row["se_oracle_ipw_re"] for row in cal_rows], dtype=float)
        pet_cal = np.asarray([row["mu_pet_peese"] for row in cal_rows], dtype=float)
        se_pet_cal = np.asarray([row["se_pet_peese"] for row in cal_rows], dtype=float)
        bgwam_cal = np.asarray([row["mu_bayesian_gwam"] for row in cal_rows], dtype=float)
        se_bgwam_cal = np.asarray([row["se_bayesian_gwam"] for row in cal_rows], dtype=float)
        eps = 1e-12
        re_std = np.abs(re_cal - cal_mu_true) / np.maximum(se_re_cal, eps)
        g_std = np.abs(g_cal - cal_mu_true) / np.maximum(se_g_cal, eps)
        ipw_std = np.abs(ipw_cal - cal_mu_true) / np.maximum(se_ipw_cal, eps)
        pet_std = np.abs(pet_cal - cal_mu_true) / np.maximum(se_pet_cal, eps)
        bgwam_std = np.abs(bgwam_cal - cal_mu_true) / np.maximum(se_bgwam_cal, eps)
        ipw_std = ipw_std[np.isfinite(ipw_std)]
        pet_std = pet_std[np.isfinite(pet_std)]
        bgwam_std = bgwam_std[np.isfinite(bgwam_std)]
        mult_re = float(np.quantile(re_std, args.ci_target_coverage, method="linear"))
        mult_g = float(np.quantile(g_std, args.ci_target_coverage, method="linear"))
        mult_ipw = float(np.quantile(ipw_std, args.ci_target_coverage, method="linear")) if ipw_std.size else float("nan")
        mult_pet = float(np.quantile(pet_std, args.ci_target_coverage, method="linear")) if pet_std.size else float("nan")
        mult_bgwam = float(np.quantile(bgwam_std, args.ci_target_coverage, method="linear")) if bgwam_std.size else float("nan")
        ci_calibration_summary["calibration_runs_success"] = int(len(cal_rows))
        ci_calibration_summary["scenario_mu_used"] = cal_mu_true
        ci_calibration_summary["multiplier_random_effects"] = mult_re
        ci_calibration_summary["multiplier_gwam"] = mult_g
        ci_calibration_summary["multiplier_oracle_ipw_random_effects"] = (
            mult_ipw if np.isfinite(mult_ipw) else None
        )
        ci_calibration_summary["multiplier_pet_peese"] = mult_pet if np.isfinite(mult_pet) else None
        ci_calibration_summary["multiplier_bayesian_gwam"] = mult_bgwam if np.isfinite(mult_bgwam) else None

        for idx, scenario in enumerate(scenario_results):
            buf = scenario_buffers[idx]
            mu_true = float(buf["mu_true"])  # type: ignore[index]
            re_cov, re_w = summarize_model_with_multiplier(
                buf["mu_re"],  # type: ignore[index]
                buf["se_re"],  # type: ignore[index]
                mu_true,
                mult_re,
            )
            g_cov, g_w = summarize_model_with_multiplier(
                buf["mu_g"],  # type: ignore[index]
                buf["se_g"],  # type: ignore[index]
                mu_true,
                mult_g,
            )
            scenario["random_effects"]["coverage_calibrated"] = re_cov  # type: ignore[index]
            scenario["random_effects"]["mean_ci_width_calibrated"] = re_w  # type: ignore[index]
            scenario["gwam_null"]["coverage_calibrated"] = g_cov  # type: ignore[index]
            scenario["gwam_null"]["mean_ci_width_calibrated"] = g_w  # type: ignore[index]
            if np.isfinite(mult_ipw):
                ipw_cov, ipw_w = summarize_model_with_multiplier(
                    buf["mu_ipw"],  # type: ignore[index]
                    buf["se_ipw"],  # type: ignore[index]
                    mu_true,
                    mult_ipw,
                )
                scenario["oracle_ipw_random_effects"]["coverage_calibrated"] = ipw_cov  # type: ignore[index]
                scenario["oracle_ipw_random_effects"]["mean_ci_width_calibrated"] = ipw_w  # type: ignore[index]
            if np.isfinite(mult_pet):
                pet_cov, pet_w = summarize_model_with_multiplier(
                    buf["mu_pet"],  # type: ignore[index]
                    buf["se_pet"],  # type: ignore[index]
                    mu_true,
                    mult_pet,
                )
                scenario["pet_peese"]["coverage_calibrated"] = pet_cov  # type: ignore[index]
                scenario["pet_peese"]["mean_ci_width_calibrated"] = pet_w  # type: ignore[index]
            if np.isfinite(mult_bgwam):
                bgwam_cov, bgwam_w = summarize_model_with_multiplier(
                    buf["mu_bgwam"],  # type: ignore[index]
                    buf["se_bgwam"],  # type: ignore[index]
                    mu_true,
                    mult_bgwam,
                )
                scenario["bayesian_gwam"]["coverage_calibrated"] = bgwam_cov  # type: ignore[index]
                scenario["bayesian_gwam"]["mean_ci_width_calibrated"] = bgwam_w  # type: ignore[index]

    summary = {
        "metadata": build_environment_metadata(),
        "registry_csv": str(args.registry_csv),
        "n_trials_selection_pool": int(weights_selection.size),
        "n_trials_results_only_pool": int(results_only_weights.size),
        "lambda_target": args.lambda_target,
        "calibration_objective": calibration_objective,
        "lambda_observed_registry_pmid_only_raw": observed_lambda["pmid_only_raw"],
        "lambda_observed_registry_non_ghost_raw": observed_lambda["non_ghost_raw"],
        "lambda_observed_registry_pmid_only_effective": observed_lambda["pmid_only_effective"],
        "lambda_observed_registry_non_ghost_effective": observed_lambda["non_ghost_effective"],
        "lambda_observed_registry_raw": lambda_target_raw,
        "lambda_observed_registry_effective": lambda_target_effective,
        "results_only_weight_raw": observed_lambda["results_only_weight_raw"],
        "results_only_weight_effective": observed_lambda["results_only_weight_effective"],
        "winsorize_pct": args.winsorize_pct,
        "simulation_settings": {
            "mu_true_list": mu_true_values,
            "calibration_mu": calibration_mu,
            "tau2": args.tau2,
            "n_meta": args.n_meta,
            "p_sig": args.p_sig,
            "p_nonsig_used": p_nonsig_used,
            "calibration_runs": args.calibration_runs,
            "calibration_grid_size": args.calibration_grid_size,
            "calibration_refine_rounds": args.calibration_refine_rounds,
            "control_event_rate": args.control_event_rate,
            "control_event_sd_logit": args.control_event_sd_logit,
            "allocation_min_frac_treat": args.allocation_min_frac_treat,
            "allocation_max_frac_treat": args.allocation_max_frac_treat,
            "results_only_mode": args.results_only_mode,
            "results_only_mu": args.results_only_mu,
            "results_only_sd": args.results_only_sd,
            "ci_calibration_runs": args.ci_calibration_runs,
            "seed": args.seed,
        },
        "calibration": {
            "enabled": bool(args.calibrate_nonsig),
            "target_lambda_reporting": lambda_target_effective if args.calibrate_nonsig else None,
            "target_lambda_objective": calibration_target_effective if args.calibrate_nonsig else None,
            "p_nonsig_calibrated": float(calibration_info["p_nonsig_calibrated"])
            if calibration_info
            else None,
            "lambda_calibrated_probe": float(calibration_info["lambda_at_p_nonsig"])
            if calibration_info
            else None,
            "lambda_calibration_abs_error": float(calibration_info["abs_error"]) if calibration_info else None,
            "lambda_calibration_max_abs_error": args.calibration_max_abs_error if args.calibrate_nonsig else None,
            "lambda_calibration_within_tolerance": (
                float(calibration_info["abs_error"]) <= args.calibration_max_abs_error
                if calibration_info
                else None
            ),
            "lambda_target_attainable_in_grid": bool(calibration_info["target_attainable"])
            if calibration_info
            else None,
            "lambda_calibration_range_grid_low": float(calibration_info["lambda_range_low"])
            if calibration_info
            else None,
            "lambda_calibration_range_grid_high": float(calibration_info["lambda_range_high"])
            if calibration_info
            else None,
            "n_calibration_candidates_evaluated": int(calibration_info["n_candidates_evaluated"])
            if calibration_info
            else None,
        },
        "ci_calibration": ci_calibration_summary,
        "scenario_results": scenario_results,
    }

    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Wrote {args.output_json}")
    print(
        f"lambda_target={args.lambda_target}, "
        f"calibration_objective={calibration_objective}, "
        f"observed_lambda_raw={lambda_target_raw:.3f}, "
        f"observed_lambda_effective={lambda_target_effective:.3f}, "
        f"p_sig={args.p_sig:.3f}, p_nonsig_used={p_nonsig_used:.3f}, "
        f"selection_trials={weights_selection.size}, results_only_trials={results_only_weights.size}"
    )
    if calibration_info:
        print(
            "lambda_calibration="
            f"objective_target={calibration_target_effective:.3f}, "
            f"probe={float(calibration_info['lambda_at_p_nonsig']):.3f}, "
            f"abs_error={float(calibration_info['abs_error']):.3f}, "
            f"within_tol={float(calibration_info['abs_error']) <= args.calibration_max_abs_error}, "
            f"attainable={bool(calibration_info['target_attainable'])}, "
            f"range=[{float(calibration_info['lambda_range_low']):.3f}, {float(calibration_info['lambda_range_high']):.3f}]"
        )
    if args.ci_calibration == "empirical":
        print(
            "ci_calibration="
            f"empirical(target={args.ci_target_coverage:.3f}, "
            f"mu={ci_calibration_summary['scenario_mu_used']}, "
            f"n={ci_calibration_summary['calibration_runs_success']}, "
            f"mult_re={ci_calibration_summary['multiplier_random_effects']:.3f}, "
            f"mult_gwam={ci_calibration_summary['multiplier_gwam']:.3f}, "
            f"mult_ipw={ci_calibration_summary['multiplier_oracle_ipw_random_effects']}, "
            f"mult_pet={ci_calibration_summary['multiplier_pet_peese']})"
        )
    for scenario in scenario_results:
        mu_true = scenario["mu_true"]
        n_succ = scenario["n_meta_success"]
        lambda_sim_mean_pmid = scenario["lambda_sim_mean_pmid_only"]
        lambda_sim_mean_non_ghost = scenario["lambda_sim_mean_non_ghost"]
        re_bias = scenario["random_effects"]["bias"]  # type: ignore[index]
        gwam_bias = scenario["gwam_null"]["bias"]  # type: ignore[index]
        ipw_bias = scenario["oracle_ipw_random_effects"]["bias"]  # type: ignore[index]
        pet_bias = scenario["pet_peese"]["bias"]  # type: ignore[index]
        bgwam_bias = scenario["bayesian_gwam"]["bias"]  # type: ignore[index]
        re_cov = scenario["random_effects"]["coverage_95"]  # type: ignore[index]
        gwam_cov = scenario["gwam_null"]["coverage_95"]  # type: ignore[index]
        bgwam_cov = scenario["bayesian_gwam"]["coverage_95"]  # type: ignore[index]
        re_cov_cal = scenario["random_effects"].get("coverage_calibrated", None)  # type: ignore[index]
        gwam_cov_cal = scenario["gwam_null"].get("coverage_calibrated", None)  # type: ignore[index]
        cov_tail = ""
        if re_cov_cal is not None and gwam_cov_cal is not None:
            cov_tail = (
                f", cov_re={re_cov:.3f}/{re_cov_cal:.3f}, "
                f"cov_gwam={gwam_cov:.3f}/{gwam_cov_cal:.3f}"
            )
        else:
            cov_tail = f", cov_re={re_cov:.3f}, cov_gwam={gwam_cov:.3f}"
        print(
            f"mu_true={mu_true:.3f}, n={n_succ}, "
            f"lambda_sim_pmid={lambda_sim_mean_pmid:.3f}, "
            f"lambda_sim_non_ghost={lambda_sim_mean_non_ghost:.3f}, "
            f"bias_re={re_bias:.3f}, bias_gwam={gwam_bias:.3f}, "
            f"bias_ipw={ipw_bias:.3f}, bias_pet={pet_bias:.3f}, "
            f"bias_bgwam={bgwam_bias:.3f}"
            f"{cov_tail}, cov_bgwam={bgwam_cov:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
