# GWAM: concrete submission fixes

This file converts the multi-persona review into repository-side actions that should be checked before external submission of the F1000 software paper for `GWAM`.

## Detectable Local State
- Documentation files detected: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Environment lock or container files detected: `requirements.txt`.
- Package manifests detected: none detected.
- Example data files detected: `f1000_artifacts/example_dataset.csv`.
- Validation artifacts detected: `f1000_artifacts/validation_summary.md`, `tests/conftest.py`, `tests/test_bayesian_gwam.py`, `tests/test_ctgov_linkage_builder.py`, `tests/test_ctgov_results_extraction.py`, `tests/test_estimators.py`, `tests/test_gwam_logic.py`, `tests/test_pairwise70_benchmark.py`.
- Detected public repository root: `https://github.com/mahmood726-cyber/gwam`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/gwam/tree/8b9c3f6495175dde32525705332a0ad1d8cb46cf`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## High-Priority Fixes
- Check that the manuscript's named example paths exist in the public archive and can be run without repository archaeology.
- Confirm that the cited repository root (`https://github.com/mahmood726-cyber/gwam`) resolves to the same fixed public source snapshot used for submission.
- Archive the tagged release and insert the Zenodo DOI or record URL once it has been minted; no project-specific archive DOI was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text matches the shipped artifacts.

## Numeric Evidence Available To Quote
- `f1000_artifacts/validation_summary.md` reports Example dataset used for walkthrough: pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json.
- `paper/gwam_manuscript.md` reports The published random-effects pooled estimate was log-OR = 0.255 (SE_pub = 0.096; 95% CI: 0.066, 0.443; OR = 1.29), derived from the PMID-linked active-comparator trials in the escitalopram registry subset using DerSimonian-Laird random effects.

## Manuscript Files To Keep In Sync
- `F1000_Software_Tool_Article.md`
- `F1000_Reviewer_Rerun_Manifest.md`
- `F1000_MultiPersona_Review.md`
- `F1000_Submission_Checklist_RealReview.md` where present
- README/tutorial files and the public repository release metadata
