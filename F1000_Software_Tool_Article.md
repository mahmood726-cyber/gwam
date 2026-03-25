# GWAM: a software tool for reviewer-auditable evidence synthesis

## Authors
- Mahmood Ahmad [1,2]
- Niraj Kumar [1]
- Bilaal Dar [3]
- Laiba Khan [1]
- Andrew Woo [4]
- Corresponding author: Andrew Woo (andy2709w@gmail.com)

## Affiliations
1. Royal Free Hospital
2. Tahir Heart Institute Rabwah
3. King's College Medical School
4. St George's Medical School

## Abstract
**Background:** Standard publication-bias diagnostics infer missing evidence from the published literature alone. GWAM instead treats completed-but-unpublished protocols in trial registries as an observable missing-evidence layer that can be used to quantify attenuation and transparency loss.

**Methods:** The Python project implements a registry-first workflow: candidate scanning on ClinicalTrials.gov, trial-level extraction, publication-linkage classification, ghost-protocol identification, GWAM correction, and Monte Carlo sensitivity analyses with optional Pairwise70 benchmarking.

**Results:** The repository includes full workflow scripts, benchmark outputs, linkage summaries, sensitivity reruns, and clear methodological review notes documenting where the current implementation behaves as a conservative correction framework rather than a full hierarchical selection model.

**Conclusions:** GWAM is reported as a registry-informed sensitivity framework for publication bias, with explicit caveats about calibration, prior assumptions, and comparator breadth.

## Keywords
publication bias; trial registries; ghost trials; ClinicalTrials.gov; meta-analysis sensitivity analysis; software tool

## Introduction
The software contribution lies in making registry linkage auditable at scale. It preserves raw registry pulls, query-level linkage summaries, and scenario-wise benchmark outputs so that reviewers can inspect how a ghost-trial assumption changes a pooled conclusion.

The package should be positioned alongside conventional random-effects models, IPW sensitivity analyses, PET-PEESE, trim-and-fill, and formal selection models. Reviewer critique in the local project files makes it clear that claim discipline and comparator transparency are essential.

The manuscript structure below is deliberately aligned to common open-software review requests: the rationale is stated explicitly, at least one runnable example path is named, local validation artifacts are listed, and conclusions are bounded to the functions and outputs documented in the repository.

## Methods
### Software architecture and workflow
Scripts are arranged into candidate scanning, CT.gov extraction, PubMed enrichment, GWAM modeling, one-command workflow execution, and Pairwise70 benchmarking. Output directories preserve analysis JSON, registry tables, and scenario summaries for side-by-side comparison.

### Installation, runtime, and reviewer reruns
The local implementation is packaged under `C:\Models\GWAM`. The manuscript identifies the local entry points, dependency manifest, fixed example input, and expected saved outputs so that reviewers can rerun the documented workflow without reconstructing it from scratch.

- Entry directory: `C:\Models\GWAM`.
- Detected documentation entry points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `requirements.txt`.
- Named worked-example paths in this draft: `scripts/run_full_workflow.py` for end-to-end registry-to-model execution; `paper_review.md` for the strongest current reviewer critique of assumptions and reporting; `pairwise70_benchmark/` and `pairwise70_benchmark_sensitivity/` for large-scale sensitivity outputs.
- Detected validation or regression artifacts: `f1000_artifacts/validation_summary.md`, `tests/conftest.py`, `tests/test_bayesian_gwam.py`, `tests/test_ctgov_linkage_builder.py`, `tests/test_ctgov_results_extraction.py`, `tests/test_estimators.py`, `tests/test_gwam_logic.py`, `tests/test_pairwise70_benchmark.py`.
- Detected example or sample data files: `f1000_artifacts/example_dataset.csv`.

### Worked examples and validation materials
**Example or fixed demonstration paths**
- `scripts/run_full_workflow.py` for end-to-end registry-to-model execution.
- `paper_review.md` for the strongest current reviewer critique of assumptions and reporting.
- `pairwise70_benchmark/` and `pairwise70_benchmark_sensitivity/` for large-scale sensitivity outputs.

**Validation and reporting artifacts**
- `tests/` for code-level regression checks.
- Simulation coverage and calibration outputs embedded in workflow JSON products.
- Pairwise70 benchmark summaries stratified by lambda source and sparsity bands.

### Typical outputs and user-facing deliverables
- Registry candidate scans and CT.gov linkage summaries.
- GWAM-adjusted analyses and simulation summaries.
- Large-scale Pairwise70 attenuation and sensitivity tables.

### Reviewer-informed safeguards
- Provides a named example workflow or fixed demonstration path.
- Documents local validation artifacts rather than relying on unsupported claims.
- Positions the software against existing tools without claiming blanket superiority.
- States limitations and interpretation boundaries in the manuscript itself.
- Requires explicit environment capture and public example accessibility in the released archive.

## Review-Driven Revisions
This draft has been tightened against recurring open peer-review objections taken from the supplied reviewer reports.
- Reproducibility: the draft names a reviewer rerun path and points readers to validation artifacts instead of assuming interface availability is proof of correctness.
- Validation: claims are anchored to local tests, validation summaries, simulations, or consistency checks rather than to unsupported assertions of performance.
- Comparators and niche: the manuscript now names the relevant comparison class and keeps the claimed niche bounded instead of implying universal superiority.
- Documentation and interpretation: the text expects a worked example, input transparency, and reviewer-verifiable outputs rather than a high-level feature list alone.
- Claims discipline: conclusions are moderated to the documented scope of GWAM and paired with explicit limitations.

## Use Cases and Results
The software outputs should be described in terms of concrete reviewer-verifiable workflows: running the packaged example, inspecting the generated results, and checking that the reported interpretation matches the saved local artifacts. In this project, the most important result layer is the availability of a transparent execution path from input to analysis output.

Representative local result: `f1000_artifacts/validation_summary.md` reports Example dataset used for walkthrough: pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json.

### Concrete local quantitative evidence
- `f1000_artifacts/validation_summary.md` reports Example dataset used for walkthrough: pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json.
- `paper/gwam_manuscript.md` reports The published random-effects pooled estimate was log-OR = 0.255 (SE_pub = 0.096; 95% CI: 0.066, 0.443; OR = 1.29), derived from the PMID-linked active-comparator trials in the escitalopram registry subset using DerSimonian-Laird random effects.

## Discussion
Representative local result: `f1000_artifacts/validation_summary.md` reports Example dataset used for walkthrough: pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json.

The real reviewer report in this directory is useful because it identifies exactly what the F1000 paper must do: describe the model transparently, moderate claims, state that fixed-zero or conservative assumptions are sensitivity choices rather than truth, and present comparator methods beyond a plain random-effects baseline.

### Limitations
- Current interval calibration is imperfect and must be described honestly.
- The present implementation is a conservative sensitivity framework, not a full Bayesian selection model.
- Comparator coverage is broader than before but still requires careful framing relative to established publication-bias methods.

## Software Availability
- Local source package: `GWAM` under `C:\Models`.
- Public repository: `https://github.com/mahmood726-cyber/gwam`.
- Public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/gwam/tree/8b9c3f6495175dde32525705332a0ad1d8cb46cf`.
- DOI/archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture detected locally: `requirements.txt`.
- Reviewer-facing documentation detected locally: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Reproducibility walkthrough: `f1000_artifacts/tutorial_walkthrough.md` where present.
- Validation summary: `f1000_artifacts/validation_summary.md` where present.
- Reviewer rerun manifest: `F1000_Reviewer_Rerun_Manifest.md`.
- Multi-persona review memo: `F1000_MultiPersona_Review.md`.
- Concrete submission-fix note: `F1000_Concrete_Submission_Fixes.md`.
- License: see the local `LICENSE` file.

## Data Availability
Registry queries, raw linkage tables, and benchmark outputs are stored locally in the project tree. External release materials should freeze the exact registry snapshot and document any PubMed-enrichment failures.

## Reporting Checklist
Real-peer-review-aligned checklist: `F1000_Submission_Checklist_RealReview.md`.
Reviewer rerun companion: `F1000_Reviewer_Rerun_Manifest.md`.
Companion reviewer-response artifact: `F1000_MultiPersona_Review.md`.
Project-level concrete fix list: `F1000_Concrete_Submission_Fixes.md`.

## Declarations
### Competing interests
The authors declare that no competing interests were disclosed.

### Grant information
No specific grant was declared for this manuscript draft.

### Author contributions (CRediT)
| Author | CRediT roles |
|---|---|
| Mahmood Ahmad | Conceptualization; Software; Validation; Data curation; Writing - original draft; Writing - review and editing |
| Niraj Kumar | Conceptualization |
| Bilaal Dar | Conceptualization |
| Laiba Khan | Conceptualization |
| Andrew Woo | Conceptualization |

### Acknowledgements
The authors acknowledge contributors to open statistical methods, reproducible research software, and reviewer-led software quality improvement.

## References
1. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188.
2. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Statistics in Medicine. 2002;21(11):1539-1558.
3. Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of Statistical Software. 2010;36(3):1-48.
4. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
5. Fay C, Rochette S, Guyader V, Girard C. Engineering Production-Grade Shiny Apps. Chapman and Hall/CRC. 2022.
