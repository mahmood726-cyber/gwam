#!/usr/bin/env python3
"""Run full GWAM pipeline for one intervention/condition pair."""

from __future__ import annotations

import argparse
import math
import platform
import subprocess
import sys
from pathlib import Path

from gwam_utils import sanitize_path_component

DEPRESSION_RESPONSE_KEYWORDS = [
    "response",
    "remission",
    "hamd",
    "ham d",
    "hamilton",
    "madrs",
    "montgomery asberg",
]


def run(cmd: list[str]) -> None:
    print(">", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intervention", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument(
        "--published-mu",
        type=float,
        required=True,
        help="Published pooled effect size (input to GWAM correction).",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed forwarded to model and simulation scripts.")
    parser.add_argument("--sim-n", type=int, default=5000)
    parser.add_argument("--ghost-mu", type=float, default=0.0)
    parser.add_argument("--ghost-sd", type=float, default=0.1)
    parser.add_argument(
        "--results-only-mode",
        choices=["as_unknown", "as_observed"],
        default="as_unknown",
    )
    parser.add_argument("--results-only-mu", type=float, default=0.0)
    parser.add_argument("--results-only-sd", type=float, default=0.1)
    parser.add_argument(
        "--lambda-target",
        choices=["pmid_only", "non_ghost"],
        default="non_ghost",
        help="Lambda definition used to calibrate simulation publication model.",
    )
    parser.add_argument("--tau2", type=float, default=0.1)
    parser.add_argument("--control-event-rate", type=float, default=0.35)
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
        help="Minimum treatment-arm fraction used to sample per-study allocation ratio in simulation.",
    )
    parser.add_argument(
        "--allocation-max-frac-treat",
        type=float,
        default=0.55,
        help="Maximum treatment-arm fraction used to sample per-study allocation ratio in simulation.",
    )
    parser.add_argument(
        "--calibration-runs",
        type=int,
        default=1200,
        help="Replicates used per lambda-calibration candidate evaluation.",
    )
    parser.add_argument("--calibration-grid-size", type=int, default=25)
    parser.add_argument("--calibration-refine-rounds", type=int, default=3)
    parser.add_argument(
        "--calibration-max-abs-error",
        type=float,
        default=0.03,
        help="Advisory maximum lambda-calibration absolute error.",
    )
    parser.add_argument(
        "--calibration-enforce-tolerance",
        action="store_true",
        help="Fail workflow if lambda calibration exceeds --calibration-max-abs-error.",
    )
    parser.add_argument(
        "--calibration-objective",
        choices=["auto", "pmid_only", "non_ghost"],
        default="auto",
        help="Objective lambda used to calibrate p-nonsig in simulation.",
    )
    parser.add_argument(
        "--ci-calibration",
        choices=["none", "empirical"],
        default="empirical",
    )
    parser.add_argument("--ci-target-coverage", type=float, default=0.95)
    parser.add_argument(
        "--ci-calibration-runs",
        type=int,
        default=2000,
        help="Independent replicates used to estimate empirical CI multipliers.",
    )
    parser.add_argument(
        "--ghost-definition",
        choices=["no_pmid", "no_pmid_no_results"],
        default="no_pmid_no_results",
    )
    parser.add_argument(
        "--ghost-require-no-any-pmid",
        action="store_true",
        default=True,
        help=(
            "Require no PMID at all when classifying ghost protocols. "
            "Prevents strict linkage modes from labeling non-linked PMID trials as latent ghosts."
        ),
    )
    parser.add_argument(
        "--allow-any-pmid-ghost",
        action="store_false",
        dest="ghost_require_no_any_pmid",
    )
    parser.add_argument(
        "--comparator-scope",
        choices=["any", "active", "placebo"],
        default="active",
    )
    parser.add_argument(
        "--publication-linkage",
        choices=["results_pmid_strict", "results_pmid", "any_pmid"],
        default="results_pmid_strict",
        help="Publication linkage rule used for has_pmid/ghost classification.",
    )
    parser.add_argument(
        "--missing-outcome-policy",
        choices=["exclude", "include_as_unknown"],
        default="include_as_unknown",
        help="How to handle trials with missing/non-parseable primary outcome metadata.",
    )
    parser.add_argument(
        "--estimand-profile",
        choices=["auto", "none", "depression_response_acute"],
        default="auto",
        help="Optional default estimand profile; explicit outcome/timeframe args override.",
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
    parser.add_argument("--outcome-timeframe-min-weeks", type=float, default=None)
    parser.add_argument("--outcome-timeframe-max-weeks", type=float, default=None)
    parser.add_argument(
        "--require-randomized",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no-require-randomized",
        action="store_false",
        dest="require_randomized",
    )
    parser.add_argument(
        "--require-treatment-purpose",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no-require-treatment-purpose",
        action="store_false",
        dest="require_treatment_purpose",
    )
    parser.add_argument(
        "--calibration-mu",
        type=float,
        default=None,
        help="Mu value used to calibrate publication model; defaults to --published-mu.",
    )
    parser.add_argument(
        "--published-mu-comparator-scope",
        choices=["active", "placebo", "any"],
        required=True,
        help="Comparator scope represented by the published pooled effect source.",
    )
    parser.add_argument(
        "--skip-pubmed-enrichment",
        action="store_true",
        help="Skip PubMed metadata enrichment step.",
    )
    parser.add_argument(
        "--allow-pubmed-enrichment-failure",
        action="store_true",
        default=False,
        help="Continue workflow if PubMed enrichment fails (not recommended for final reproducible runs).",
    )
    parser.add_argument(
        "--strict-pubmed-enrichment",
        action="store_false",
        dest="allow_pubmed_enrichment_failure",
        help="Fail workflow if PubMed enrichment fails.",
    )
    parser.add_argument("--mu-true-list", default="0.0,0.1,0.2")

    # --- Published SE (for HVP model) ---
    parser.add_argument(
        "--published-se",
        type=float,
        default=None,
        help="SE of the published pooled effect. Required for HVP model unless CI bounds are provided.",
    )
    parser.add_argument(
        "--published-ci-lower",
        type=float,
        default=None,
        help="Lower bound of published pooled effect CI (used to compute SE if --published-se not given).",
    )
    parser.add_argument(
        "--published-ci-upper",
        type=float,
        default=None,
        help="Upper bound of published pooled effect CI (used to compute SE if --published-se not given).",
    )
    parser.add_argument(
        "--published-effect-scale",
        choices=["log_ratio", "difference"],
        default="log_ratio",
        help="Scale of published CI bounds: 'log_ratio' (OR/RR/HR) or 'difference' (MD/SMD/RD). Default: log_ratio.",
    )

    # --- CT.gov results extraction ---
    parser.add_argument(
        "--extract-ctgov-results",
        action="store_true",
        help="Run CT.gov results extraction for results-only studies before GWAM modelling.",
    )
    parser.add_argument(
        "--outcome-keywords-results",
        default="",
        help="Comma-separated outcome keywords for CT.gov results extraction.",
    )

    # --- HVP model args ---
    parser.add_argument(
        "--skip-bayesian",
        action="store_true",
        help="Skip the HVP (hierarchical variance-propagation) GWAM model step.",
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
        help="Comma-separated within-ghost heterogeneity values for HVP sensitivity grid.",
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
        help="Comma-separated within-results-only heterogeneity values for HVP sensitivity grid.",
    )

    return parser.parse_args()


def safe_stem(text: str) -> str:
    """Sanitize text for use as a path component (no traversal characters)."""
    text = text.strip().lower().replace(" ", "_")
    return sanitize_path_component(text)


def resolve_estimand_filters(args: argparse.Namespace) -> tuple[list[str], list[str], float | None, float | None, str]:
    allowed_phase = list(args.allowed_phase)
    keywords = [k for k in args.required_outcome_keyword if str(k).strip()]
    min_weeks = args.outcome_timeframe_min_weeks
    max_weeks = args.outcome_timeframe_max_weeks

    profile = args.estimand_profile
    if profile == "auto":
        profile = "depression_response_acute" if "depression" in args.condition.lower() else "none"

    if profile == "depression_response_acute":
        if not keywords:
            keywords = list(DEPRESSION_RESPONSE_KEYWORDS)
        if min_weeks is None:
            min_weeks = 4.0
        if max_weeks is None:
            max_weeks = 16.0
    return allowed_phase, keywords, min_weeks, max_weeks, profile


def compute_se_from_ci(
    ci_lower: float, ci_upper: float, effect_scale: str, conf_level: float = 0.95,
) -> float:
    """Derive SE from confidence interval bounds.

    For log_ratio scale (OR/RR/HR): SE = (log(ci_upper) - log(ci_lower)) / (2*z)
    For difference scale (MD/SMD/RD): SE = (ci_upper - ci_lower) / (2*z)
    """
    from gwam_utils import normal_quantile  # local import to avoid top-level dependency issues

    z = normal_quantile(1.0 - (1.0 - conf_level) / 2.0)
    z_divisor = 2.0 * z

    if effect_scale == "log_ratio":
        if ci_lower <= 0 or ci_upper <= 0:
            raise ValueError(f"CI bounds must be > 0 for log_ratio scale, got ({ci_lower}, {ci_upper}).")
        se = (math.log(ci_upper) - math.log(ci_lower)) / z_divisor
    else:
        se = (ci_upper - ci_lower) / z_divisor

    if se <= 0:
        raise ValueError(f"Computed SE={se:.6f} is not positive.")
    return se


def main() -> int:
    args = parse_args()
    if args.published_mu_comparator_scope != "any" and args.published_mu_comparator_scope != args.comparator_scope:
        raise ValueError(
            "--published-mu-comparator-scope must match --comparator-scope "
            "unless set to 'any'."
        )
    if args.outcome_timeframe_min_weeks is not None and args.outcome_timeframe_max_weeks is not None:
        if args.outcome_timeframe_min_weeks > args.outcome_timeframe_max_weeks:
            raise ValueError("--outcome-timeframe-min-weeks cannot exceed --outcome-timeframe-max-weeks.")
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
    if args.ci_calibration_runs < 200:
        raise ValueError("--ci-calibration-runs must be >= 200.")
    if args.control_event_sd_logit < 0:
        raise ValueError("--control-event-sd-logit must be >= 0.")
    if not (0 < args.allocation_min_frac_treat < 1):
        raise ValueError("--allocation-min-frac-treat must be in (0,1).")
    if not (0 < args.allocation_max_frac_treat < 1):
        raise ValueError("--allocation-max-frac-treat must be in (0,1).")
    if args.allocation_min_frac_treat > args.allocation_max_frac_treat:
        raise ValueError("--allocation-min-frac-treat cannot exceed --allocation-max-frac-treat.")

    # --- Resolve published SE ---
    published_se: float | None = args.published_se
    if published_se is None and args.published_ci_lower is not None and args.published_ci_upper is not None:
        published_se = compute_se_from_ci(
            args.published_ci_lower, args.published_ci_upper, args.published_effect_scale,
        )
        print(f"Computed published SE from CI ({args.published_ci_lower}, {args.published_ci_upper}) "
              f"on {args.published_effect_scale} scale: SE = {published_se:.6f}")

    if not args.skip_bayesian and published_se is None:
        raise ValueError(
            "HVP model requires --published-se or (--published-ci-lower + --published-ci-upper). "
            "Use --skip-bayesian to skip the HVP step."
        )
    if published_se is not None and published_se <= 0:
        raise ValueError("--published-se must be > 0.")

    project_root = Path(__file__).resolve().parents[1]
    scripts = project_root / "scripts"

    intr = safe_stem(args.intervention)
    cond = safe_stem(args.condition)
    if not intr or not cond:
        raise ValueError("Intervention and condition must not be empty after sanitization.")
    registry_csv = project_root / "data" / "raw" / f"{intr}__{cond}.csv"
    registry_csv_enriched = project_root / "data" / "raw" / f"{intr}__{cond}_with_results.csv"
    candidate_summary_global = project_root / "data" / "registry_candidate_summary.csv"
    candidate_summary_targeted = project_root / "data" / "registry_candidate_summary_targeted.csv"
    gwam_json = project_root / "data" / "analysis" / f"{intr}_{cond}_gwam.json"
    bayesian_json = project_root / "data" / "analysis" / f"{intr}_{cond}_gwam_bayesian.json"
    sim_json = project_root / "data" / "analysis" / f"simulation_{intr}_{cond}.json"
    calibration_mu = args.published_mu if args.calibration_mu is None else args.calibration_mu
    allowed_phase, required_keywords, min_weeks, max_weeks, active_profile = resolve_estimand_filters(args)

    print(
        "Estimand settings:",
        f"profile={active_profile}",
        f"publication_linkage={args.publication_linkage}",
        f"ghost_require_no_any_pmid={args.ghost_require_no_any_pmid}",
        f"allowed_phase={allowed_phase or 'ANY'}",
        f"required_outcome_keyword={required_keywords or 'ANY'}",
        f"outcome_timeframe_min_weeks={min_weeks}",
        f"outcome_timeframe_max_weeks={max_weeks}",
        f"missing_outcome_policy={args.missing_outcome_policy}",
    )

    # Global scan across default candidates for cross-domain context.
    run(
        [
            sys.executable,
            str(scripts / "run_candidate_scan.py"),
            "--output",
            str(candidate_summary_global),
            "--ghost-definition",
            args.ghost_definition,
            "--ghost-require-no-any-pmid"
            if args.ghost_require_no_any_pmid
            else "--allow-any-pmid-ghost",
            "--publication-linkage",
            args.publication_linkage,
            "--comparator-scope",
            args.comparator_scope,
            "--require-randomized" if args.require_randomized else "--no-require-randomized",
            "--require-treatment-purpose"
            if args.require_treatment_purpose
            else "--no-require-treatment-purpose",
            "--missing-outcome-policy",
            args.missing_outcome_policy,
        ]
        + sum([["--allowed-phase", phase] for phase in allowed_phase], [])
        + sum([["--required-outcome-keyword", keyword] for keyword in required_keywords], [])
        + (
            ["--outcome-timeframe-min-weeks", str(min_weeks)]
            if min_weeks is not None
            else []
        )
        + (
            ["--outcome-timeframe-max-weeks", str(max_weeks)]
            if max_weeks is not None
            else []
        )
    )
    # Targeted scan aligned exactly with selected intervention/condition.
    run(
        [
            sys.executable,
            str(scripts / "run_candidate_scan.py"),
            "--output",
            str(candidate_summary_targeted),
            "--candidate",
            f"{args.intervention}::{args.condition}",
            "--ghost-definition",
            args.ghost_definition,
            "--ghost-require-no-any-pmid"
            if args.ghost_require_no_any_pmid
            else "--allow-any-pmid-ghost",
            "--publication-linkage",
            args.publication_linkage,
            "--comparator-scope",
            args.comparator_scope,
            "--require-randomized" if args.require_randomized else "--no-require-randomized",
            "--require-treatment-purpose"
            if args.require_treatment_purpose
            else "--no-require-treatment-purpose",
            "--missing-outcome-policy",
            args.missing_outcome_policy,
        ]
        + sum([["--allowed-phase", phase] for phase in allowed_phase], [])
        + sum([["--required-outcome-keyword", keyword] for keyword in required_keywords], [])
        + (
            ["--outcome-timeframe-min-weeks", str(min_weeks)]
            if min_weeks is not None
            else []
        )
        + (
            ["--outcome-timeframe-max-weeks", str(max_weeks)]
            if max_weeks is not None
            else []
        )
    )
    run(
        [
            sys.executable,
            str(scripts / "fetch_ctgov_registry.py"),
            "--intervention",
            args.intervention,
            "--condition",
            args.condition,
            "--ghost-definition",
            args.ghost_definition,
            "--ghost-require-no-any-pmid"
            if args.ghost_require_no_any_pmid
            else "--allow-any-pmid-ghost",
            "--publication-linkage",
            args.publication_linkage,
            "--comparator-scope",
            args.comparator_scope,
            "--require-randomized" if args.require_randomized else "--no-require-randomized",
            "--require-treatment-purpose"
            if args.require_treatment_purpose
            else "--no-require-treatment-purpose",
            "--missing-outcome-policy",
            args.missing_outcome_policy,
        ]
        + sum([["--allowed-phase", phase] for phase in allowed_phase], [])
        + sum([["--required-outcome-keyword", keyword] for keyword in required_keywords], [])
        + (
            ["--outcome-timeframe-min-weeks", str(min_weeks)]
            if min_weeks is not None
            else []
        )
        + (
            ["--outcome-timeframe-max-weeks", str(max_weeks)]
            if max_weeks is not None
            else []
        )
    )
    pmid_column = "pmids_results" if args.publication_linkage in {"results_pmid_strict", "results_pmid"} else "pmids_any"
    if args.skip_pubmed_enrichment:
        print("Skipping PubMed enrichment by request.")
    else:
        try:
            run(
                [
                    sys.executable,
                    str(scripts / "enrich_pubmed_links.py"),
                    "--registry-csv",
                    str(registry_csv),
                    "--pmid-column",
                    pmid_column,
                ]
            )
        except subprocess.CalledProcessError:
            if not args.allow_pubmed_enrichment_failure:
                raise
            print("WARNING: PubMed enrichment failed; continuing without metadata enrichment.")

    # --- CT.gov results extraction ---
    if args.extract_ctgov_results:
        ctgov_results_cmd = [
            sys.executable,
            str(scripts / "extract_ctgov_results.py"),
            "--registry-csv",
            str(registry_csv),
            "--output-csv",
            str(registry_csv_enriched),
        ]
        if args.outcome_keywords_results:
            ctgov_results_cmd += ["--outcome-keywords", args.outcome_keywords_results]
        run(ctgov_results_cmd)
        # Use enriched CSV for downstream steps
        model_csv = registry_csv_enriched
    else:
        model_csv = registry_csv

    run(
        [
            sys.executable,
            str(scripts / "model_gwam.py"),
            "--registry-csv",
            str(model_csv),
            "--published-mu",
            str(args.published_mu),
            "--sim-n",
            str(args.sim_n),
            "--ghost-mu",
            str(args.ghost_mu),
            "--ghost-sd",
            str(args.ghost_sd),
            "--results-only-mode",
            args.results_only_mode,
            "--results-only-mu",
            str(args.results_only_mu),
            "--results-only-sd",
            str(args.results_only_sd),
            "--seed",
            str(args.seed),
            "--output-json",
            str(gwam_json),
        ]
    )

    # --- HVP (hierarchical variance-propagation) GWAM ---
    if not args.skip_bayesian:
        bayesian_cmd = [
            sys.executable,
            str(scripts / "model_gwam_bayesian.py"),
            "--registry-csv",
            str(model_csv),
            "--published-mu",
            str(args.published_mu),
            "--published-se",
            str(published_se),
            "--prior-ghost-mu-mean",
            str(args.prior_ghost_mu_mean),
            "--prior-ghost-mu-sd",
            str(args.prior_ghost_mu_sd),
            "--prior-ghost-sigma-grid",
            args.prior_ghost_sigma_grid,
            "--prior-ro-mu-mean",
            str(args.prior_ro_mu_mean),
            "--prior-ro-mu-sd",
            str(args.prior_ro_mu_sd),
            "--prior-ro-sigma-grid",
            args.prior_ro_sigma_grid,
            "--output-json",
            str(bayesian_json),
        ]
        if args.extract_ctgov_results:
            bayesian_cmd.append("--use-observed-results")
        run(bayesian_cmd)
    else:
        print("Skipping HVP-GWAM model by request.")

    run(
        [
            sys.executable,
            str(scripts / "simulate_gwam_vs_re.py"),
            "--registry-csv",
            str(model_csv),
            "--calibrate-nonsig",
            "--n-meta",
            str(args.sim_n),
            "--calibration-runs",
            str(args.calibration_runs),
            "--calibration-grid-size",
            str(args.calibration_grid_size),
            "--calibration-refine-rounds",
            str(args.calibration_refine_rounds),
            "--calibration-max-abs-error",
            str(args.calibration_max_abs_error),
            *(["--calibration-enforce-tolerance"] if args.calibration_enforce_tolerance else []),
            "--calibration-objective",
            args.calibration_objective,
            "--tau2",
            str(args.tau2),
            "--mu-true-list",
            args.mu_true_list,
            "--calibration-mu",
            str(calibration_mu),
            "--control-event-rate",
            str(args.control_event_rate),
            "--control-event-sd-logit",
            str(args.control_event_sd_logit),
            "--allocation-min-frac-treat",
            str(args.allocation_min_frac_treat),
            "--allocation-max-frac-treat",
            str(args.allocation_max_frac_treat),
            "--results-only-mode",
            args.results_only_mode,
            "--results-only-mu",
            str(args.results_only_mu),
            "--results-only-sd",
            str(args.results_only_sd),
            "--lambda-target",
            args.lambda_target,
            "--ci-calibration",
            args.ci_calibration,
            "--ci-target-coverage",
            str(args.ci_target_coverage),
            "--ci-calibration-runs",
            str(args.ci_calibration_runs),
            "--seed",
            str(args.seed),
            "--output-json",
            str(sim_json),
        ]
    )

    print("\nWorkflow complete.")
    print(f"Environment: Python {sys.version.split()[0]}, {platform.system()} {platform.release()}")
    try:
        import numpy as np
        print(f"  numpy={np.__version__}", end="")
        import requests as req
        print(f", requests={req.__version__}", end="")
    except ImportError:
        pass
    print()
    print(f"- Candidate summary (global): {candidate_summary_global}")
    print(f"- Candidate summary (targeted): {candidate_summary_targeted}")
    print(f"- Registry: {registry_csv}")
    if args.extract_ctgov_results:
        print(f"- Registry (with results): {registry_csv_enriched}")
    print(f"- GWAM summary: {gwam_json}")
    if not args.skip_bayesian:
        print(f"- HVP-GWAM: {bayesian_json}")
    print(f"- Simulation summary: {sim_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
