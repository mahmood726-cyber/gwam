# GWAM Real-Data Workflow

**Requires Python >= 3.10** (uses `X | Y` union type syntax).

This project provides a reproducible GWAM pipeline on ClinicalTrials.gov registry data:

1. Pull completed trial records for selected intervention/condition pairs.
2. Classify publication linkage from CT.gov reference types.
3. Identify strict ghost protocols (no linked PMID and no posted results).
4. Compute integrity ratios (`lambda_pmid_only`, `lambda_non_ghost`).
5. Run GWAM correction and RE-vs-GWAM Monte Carlo simulation with optional CI calibration.
6. Compare against standard methods (IV-FE, DL-RE, PM-RE, HK, Huber, Student-t) and GRMA (`grey_meta_v8.py`).

See `LICENSE` for terms.

## Setup

```bash
python -m pip install -r requirements.txt
```

## Running Tests

```bash
python -m pytest tests/ -v
```

## 1) Candidate scan

Global scan:

```bash
python scripts/run_candidate_scan.py \
  --ghost-definition no_pmid_no_results \
  --ghost-require-no-any-pmid \
  --publication-linkage results_pmid_strict \
  --comparator-scope active \
  --require-randomized \
  --require-treatment-purpose
```

Targeted scan:

```bash
python scripts/run_candidate_scan.py \
  --candidate "escitalopram::depression" \
  --publication-linkage results_pmid_strict
```

## 2) Trial-level extraction

Example aligned to acute depression-response estimand:

```bash
python scripts/fetch_ctgov_registry.py \
  --intervention escitalopram \
  --condition depression \
  --ghost-definition no_pmid_no_results \
  --ghost-require-no-any-pmid \
  --publication-linkage results_pmid_strict \
  --comparator-scope active \
  --require-randomized \
  --require-treatment-purpose \
  --missing-outcome-policy include_as_unknown \
  --required-outcome-keyword response \
  --required-outcome-keyword remission \
  --required-outcome-keyword hamd \
  --required-outcome-keyword "ham d" \
  --required-outcome-keyword hamilton \
  --required-outcome-keyword madrs \
  --required-outcome-keyword "montgomery asberg" \
  --outcome-timeframe-min-weeks 4 \
  --outcome-timeframe-max-weeks 16
```

## 3) GWAM model

```bash
python scripts/model_gwam.py \
  --registry-csv data/raw/escitalopram__depression.csv \
  --published-mu 0.25464221837358075 \
  --sim-n 5000 \
  --results-only-mode as_unknown \
  --results-only-mu 0.0 \
  --results-only-sd 0.1 \
  --ghost-mu 0.0 \
  --ghost-sd 0.1 \
  --output-json data/analysis/escitalopram_depression_gwam.json
```

## 4) PubMed enrichment

```bash
python scripts/enrich_pubmed_links.py \
  --registry-csv data/raw/escitalopram__depression.csv \
  --pmid-column pmids_results
```

## 5) One-command workflow

```bash
python scripts/run_full_workflow.py \
  --intervention escitalopram \
  --condition depression \
  --published-mu 0.25464221837358075 \
  --published-mu-comparator-scope active \
  --sim-n 5000 \
  --calibration-runs 1200 \
  --calibration-grid-size 25 \
  --calibration-refine-rounds 3 \
  --calibration-max-abs-error 0.03 \
  --calibration-objective auto \
  --results-only-mode as_unknown \
  --results-only-mu 0.0 \
  --results-only-sd 0.1 \
  --lambda-target non_ghost \
  --ci-calibration empirical \
  --ci-target-coverage 0.95 \
  --ci-calibration-runs 2000 \
  --mu-true-list 0.0,0.1,0.25464221837358075 \
  --control-event-rate 0.35 \
  --control-event-sd-logit 0.35 \
  --allocation-min-frac-treat 0.45 \
  --allocation-max-frac-treat 0.55 \
  --comparator-scope active \
  --require-randomized \
  --require-treatment-purpose \
  --ghost-definition no_pmid_no_results \
  --ghost-require-no-any-pmid \
  --skip-pubmed-enrichment \
  --calibration-mu 0.25464221837358075
```

Outputs:
- `data/registry_candidate_summary.csv` (global scan)
- `data/registry_candidate_summary_targeted.csv` (selected pair)
- `data/raw/<intervention>__<condition>.csv`
- `data/raw/pubmed_links_<intervention>__<condition>.csv`
- `data/analysis/<intervention>_<condition>_gwam.json`
- `data/analysis/simulation_<intervention>_<condition>.json`

## 6) Pairwise70 real-data benchmark

Run all binary Cochrane analyses from Pairwise70 with RE, selection-IPW RE,
PET-PEESE, and GWAM shrinkage sensitivity:

1) Build CT.gov linkage summary (query-level completed status + exact NCT augmentation when available):

```bash
python scripts/build_pairwise70_ctgov_linkage_summary.py \
  --query-csv "<PATH_TO_PAIRWISE70>\analysis\transportability\ctgov_query_terms.csv" \
  --output-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --cache-dir pairwise70_benchmark_sensitivity/ctgov_linkage_raw \
  --max-pages 3 \
  --page-size 100 \
  --sleep-ms 60 \
  --refresh-cache
```

2) Run benchmark with linkage-preferred lambda source:

```bash
python scripts/run_pairwise70_benchmark.py \
  --pairwise-data-dir "<PATH_TO_PAIRWISE70>\data" \
  --ctgov-covariates-csv "<PATH_TO_PAIRWISE70>\analysis\transportability\ctgov_target_covariates.csv" \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --output-dir pairwise70_benchmark \
  --lambda-grid 0.2,0.4,0.6,0.8 \
  --default-lambda-proxy 0.7 \
  --lambda-proxy-lower 0.2 \
  --lambda-proxy-upper 0.95 \
  --p-sig 0.95
```

Outputs:
- `pairwise70_benchmark/analysis_results.csv` (analysis-level estimates)
- `pairwise70_benchmark/top_shifted_analyses.csv` (largest GWAM proxy shifts)
- `pairwise70_benchmark/summary.json` (global benchmark summary)
- `pairwise70_benchmark/stratified_by_lambda_source.csv` (required transparency stratification)
- `pairwise70_benchmark/stratified_by_k_band.csv` (required sparse-k stratification)
- `pairwise70_benchmark/clustered_by_review.csv` (review-level clustered attenuation summaries)
- `pairwise70_benchmark/robust_sensitivity_excluding_abs_mu_re_gt_4.json` (outlier sensitivity check)

Interpretation notes:
- Use attenuation metrics (`crossed below 0.10/0.20`) as primary correction outputs.
- Treat significance-drop counts as diagnostic only (GWAM proxy scales both estimate and SE).
- `lambda-source-mode=auto` now prefers exact Pairwise NCT linkage, then CT.gov linkage-derived completed-study lambda, then transport proxy, then default fallback.

Lambda sensitivity reruns and side-by-side table:

```bash
python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/base \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto

python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/strict \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --lambda-proxy-lower 0.10 \
  --lambda-proxy-upper 0.70 \
  --default-lambda-proxy 0.50 \
  --lambda-grid 0.1,0.2,0.4,0.6

python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/loose \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --lambda-proxy-lower 0.40 \
  --lambda-proxy-upper 1.00 \
  --default-lambda-proxy 0.90 \
  --lambda-grid 0.4,0.6,0.8,1.0

# PET-PEESE strict sensitivity (k >= 5)
python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/base_petk5 \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --min-k-pet 5

python scripts/compare_pairwise70_sensitivity.py \
  --scenario base=pairwise70_benchmark_sensitivity/base/summary.json \
  --scenario strict=pairwise70_benchmark_sensitivity/strict/summary.json \
  --scenario loose=pairwise70_benchmark_sensitivity/loose/summary.json \
  --output-csv pairwise70_benchmark_sensitivity/side_by_side_summary.csv \
  --output-md reports/pairwise70_lambda_sensitivity_comparison.md
```

Methods-grade uncertainty and clustering mode:

```bash
# Bayesian lambda uncertainty + review-cluster bootstrap
python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/methods_base \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --lambda-uncertainty-model beta \
  --lambda-posterior-draws 400 \
  --review-bootstrap-iters 400

# Bound-free lambda sensitivity (disable hard clipping bounds)
python scripts/run_pairwise70_benchmark.py \
  --output-dir pairwise70_benchmark_sensitivity/methods_base_unclipped \
  --ctgov-linkage-csv pairwise70_benchmark_sensitivity/ctgov_linkage_summary.csv \
  --lambda-source-mode auto \
  --lambda-uncertainty-model beta \
  --lambda-posterior-draws 400 \
  --review-bootstrap-iters 400 \
  --disable-lambda-clipping
```

## Notes

- `results_pmid_strict` links publication only via CT.gov reference types `RESULT`/`PRIMARY`.
- `results_pmid` is a broader legacy mode that also includes `DERIVED`.
- `ghost-require-no-any-pmid` keeps trials with any PMID metadata out of the fully latent ghost stratum by default.
- `missing-outcome-policy=include_as_unknown` avoids automatic exclusion of metadata-sparse trials.
- Simulation reports nominal coverage (`coverage_95`) and out-of-sample empirically calibrated coverage (`coverage_calibrated`) when enabled.
- Lambda calibration now reports absolute error and attainable range diagnostics in the output JSON.
- Simulation includes comparator models beyond RE/GWAM: selection-IPW random effects and PET-PEESE.
- Pairwise70 benchmark now supports `--lambda-uncertainty-model beta` (Beta posterior propagation), review-level clustered bootstrap summaries, and `--disable-lambda-clipping` bound-free sensitivity.
- PubMed enrichment failure is strict by default; use `--allow-pubmed-enrichment-failure` only for exploratory/offline runs.
- Use `--skip-pubmed-enrichment` if NCBI E-utilities are unavailable in your environment.
- Current GWAM implementation is a weighted correction/sensitivity framework, not a full hierarchical selection model.

## Linkage sensitivity (recommended)

Run the broad linkage sensitivity case:

```bash
python scripts/run_full_workflow.py \
  --intervention escitalopram \
  --condition depression \
  --published-mu 0.25464221837358075 \
  --published-mu-comparator-scope active \
  --publication-linkage results_pmid \
  --sim-n 5000 \
  --calibration-runs 1200 \
  --calibration-grid-size 25 \
  --calibration-refine-rounds 3 \
  --calibration-max-abs-error 0.03 \
  --calibration-objective auto \
  --results-only-mode as_unknown \
  --results-only-mu 0.0 \
  --results-only-sd 0.1 \
  --lambda-target non_ghost \
  --ci-calibration empirical \
  --ci-target-coverage 0.95 \
  --ci-calibration-runs 2000 \
  --mu-true-list 0.0,0.1,0.25464221837358075 \
  --control-event-rate 0.35 \
  --control-event-sd-logit 0.35 \
  --allocation-min-frac-treat 0.45 \
  --allocation-max-frac-treat 0.55 \
  --comparator-scope active \
  --require-randomized \
  --require-treatment-purpose \
  --ghost-definition no_pmid_no_results \
  --ghost-require-no-any-pmid \
  --skip-pubmed-enrichment \
  --calibration-mu 0.25464221837358075
```

Archive strict and broad outputs under `data/analysis/sensitivity/` and `data/raw/sensitivity/` for side-by-side reporting.
