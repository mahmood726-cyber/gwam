# GWAM: reviewer rerun manifest

This manifest is the shortest reviewer-facing rerun path for the local software package. It lists the files that should be sufficient to recreate one worked example, inspect saved outputs, and verify that the manuscript claims remain bounded to what the repository actually demonstrates.

## Reviewer Entry Points
- Project directory: `C:\Models\GWAM`.
- Preferred documentation start points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected public repository root: `https://github.com/mahmood726-cyber/gwam`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/gwam/tree/8b9c3f6495175dde32525705332a0ad1d8cb46cf`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture files: `requirements.txt`.
- Validation/test artifacts: `f1000_artifacts/validation_summary.md`, `tests/conftest.py`, `tests/test_bayesian_gwam.py`, `tests/test_ctgov_linkage_builder.py`, `tests/test_ctgov_results_extraction.py`, `tests/test_estimators.py`, `tests/test_gwam_logic.py`, `tests/test_pairwise70_benchmark.py`.

## Worked Example Inputs
- Manuscript-named example paths: `scripts/run_full_workflow.py` for end-to-end registry-to-model execution; `paper_review.md` for the strongest current reviewer critique of assumptions and reporting; `pairwise70_benchmark/` and `pairwise70_benchmark_sensitivity/` for large-scale sensitivity outputs; f1000_artifacts/example_dataset.csv.
- Auto-detected sample/example files: `f1000_artifacts/example_dataset.csv`.

## Expected Outputs To Inspect
- Registry candidate scans and CT.gov linkage summaries.
- GWAM-adjusted analyses and simulation summaries.
- Large-scale Pairwise70 attenuation and sensitivity tables.

## Minimal Reviewer Rerun Sequence
- Start with the README/tutorial files listed below and keep the manuscript paths synchronized with the public archive.
- Create the local runtime from the detected environment capture files if available: `requirements.txt`.
- Run at least one named example path from the manuscript and confirm that the generated outputs match the saved validation materials.
- Quote one concrete numeric result from the local validation snippets below when preparing the final software paper.

## Local Numeric Evidence Available
- `f1000_artifacts/validation_summary.md` reports Example dataset used for walkthrough: pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json.
- `paper/gwam_manuscript.md` reports The published random-effects pooled estimate was log-OR = 0.255 (SE_pub = 0.096; 95% CI: 0.066, 0.443; OR = 1.29), derived from the PMID-linked active-comparator trials in the escitalopram registry subset using DerSimonian-Laird random effects.
