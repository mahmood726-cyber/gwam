# GWAM: a software tool for reproducible evidence-synthesis and modeling workflows

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
**Background:** Reproducible evidence-synthesis model development often requires repeated fitting, heterogeneity checks, and structured reporting outputs. GWAM was prepared to support these workflows in a traceable software package.

**Methods:** The project is implemented using Python. The documented workflow covers data import, parameter selection, model execution, diagnostics, and export. The manuscript presents a reviewer-facing reproducibility path with explicit artifact pointers and bounded claims.

**Results:** The package supports end-to-end local workflows from input through model execution and exportable outputs, with validation evidence linked through project artifacts where available.

**Conclusions:** GWAM provides a practical workflow for transparent modeling and evidence synthesis; external submission readiness depends on final public repository metadata and DOI-archived release records.

## Keywords
evidence synthesis; reproducibility; model validation; software tool

## Visual Abstract
### Workflow overview
| Panel | Key message | What the software does | Reviewer-check evidence |
|---|---|---|---|
| Clinical problem | Evidence-synthesis workflows are often fragmented and hard to audit. | Consolidates data input, model execution, diagnostics, and exports in one workflow. | Manuscript Methods + reproducibility checklist. |
| Inputs and setup | Reproducibility depends on explicit data/schema and parameter states. | Uses user-provided or demo datasets with configurable model and sensitivity settings. | Use-case walkthrough + saved run configuration. |
| Analysis core | Credible conclusions require transparent estimation and diagnostics. | Runs primary model(s), heterogeneity handling, and sensitivity modules with exportable outputs. | Results/diagnostic tables + validation artifact paths. |
| Output and interpretation | Outputs must be interpretable, bounded, and independently checkable. | Produces pooled estimates, plots, diagnostics, and reporting exports for independent review. | Validation evidence table + checklist-linked artifacts. |
| Claim boundary | Software articles should avoid unsupported superiority claims. | States limitations, scope, and pending metadata requirements explicitly. | Discussion limitations + submission blockers checklist. |

## Introduction
Complex modeling projects can become difficult to audit when data handling, estimation, diagnostics, and reporting are separated across multiple ad hoc steps. This manuscript documents a consolidated workflow with reproducibility-first structure.

### Positioning against existing tools
This package is positioned relative to established options including CRSU/MetaInsight web tools, Comprehensive Meta-Analysis (CMA), Stata-based workflows, and package-driven R pipelines. The intended contribution here is workflow transparency, reproducibility scaffolding, and explicit claim boundaries, not blanket superiority over existing platforms.

### Table 1. Positioning matrix
| Dimension | This package | Established alternatives | Claim boundary |
|---|---|---|---|
| Primary goal | Transparent, reproducible end-to-end workflow | Mature GUIs/statistical packages with broad legacy adoption | Scope limited to demonstrated workflows |
| User profile | Clinicians/researchers needing guided reproducibility | Advanced analysts and mixed-skill teams | Complementary use is recommended |
| Strength emphasis | Auditability, artifact linkage, structured outputs | Feature breadth and ecosystem maturity | Interpret strengths relative to use case |
| Reproducibility support | Walkthrough + validation summary + checklist | Varies by tool/package and setup | Claims remain artifact-bounded |

## Methods
### Implementation
The software package is organized for local execution with explicit analysis states and reproducible outputs. Components cover data handling, model settings, estimation routines, diagnostics, and export.

### Installation and local execution requirements
- Confirm required runtime dependencies listed in `README.md` and project environment files.
- Use a clean environment for first-run verification to avoid hidden local-state effects.
- Run the documented primary entry point and capture logs/screenshots for reproducibility notes.
- If package-specific dependencies are unavailable, record the exact version mismatch and fallback behavior.

### Operation
- Open and run one primary project entry point (e.g., `README.md`).
- Load demonstration or test data (example reference: ``pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json``).
- Configure default model and sensitivity settings, then execute analysis.
- Review outputs and export artifacts for independent verification.

### Table 2. Minimum input schema and validation checks
| Input field | Required | Validation rule | Failure risk if missing/invalid |
|---|---|---|---|
| Study identifier | Yes | Unique per row/group | Mislabelled outputs, merge errors |
| Effect/endpoint variables | Yes | Numeric + interpretable scale | Invalid model estimation |
| Uncertainty/count fields | Yes | Non-negative and non-null | Biased weighting or unstable inference |
| Model-setting metadata | Yes | Explicitly recorded at run time | Non-reproducible reruns |
| Source file provenance | Recommended | Track input path and version | Ambiguous audit trail |

### Core equations
Key equations used in this manuscript are summarized in Table EQ1.

### Equation summary table
| Eq. ID | Model component | Expression | Interpretation role |
|---|---|---|---|
| E1 | General model specification | `y_i = g(x_i; \theta) + u_i + \varepsilon_i` | Defines core model structure with parameterized signal and noise terms. |
| E2 | Random-effects variance structure | `u_i \sim \mathcal{N}(0,\tau^2),\; \varepsilon_i \sim \mathcal{N}(0,v_i)` | Encodes heterogeneity and sampling variability assumptions. |
| E3 | Precision-weighted estimator | `\hat{\theta} = \frac{\sum_i w_i \theta_i}{\sum_i w_i},\; w_i = \frac{1}{v_i + \tau^2}` | Summarizes pooled effects under random-effects weighting. |

### Reproducibility and validation
The package is reported with explicit artifact references so claims can be independently checked.

#### Validation evidence table
| Validation dimension | Evidence summary | Artifact source |
|---|---|---|
| Local execution context | Project execution and output generation are documented for local reruns. | `README.md` |
| Validation scope | Validation, test, or benchmark artifacts were identified for package-level checking. | `scripts\__pycache__\run_pairwise70_benchmark.cpython-313.pyc` |
| Evidence policy | Claims are bounded to artifact-supported local outputs; no unsupported superiority claims are made. | `F1000_Submission_Checklist_RealReview.md` |
| Release requirements | Public repository URL and DOI archive are required before external submission. | `F1000_Submission_Checklist_RealReview.md` |

### Core functionality exposed in the package
- Data ingestion and preprocessing hooks
- Configurable modeling and sensitivity parameters
- Diagnostic and interpretive output pathways
- Exportable artifacts for reporting and reproducibility

### Reviewer-informed reproducibility safeguards
- End-to-end walkthrough included for independent rerun.
- Evidence statements tied to local artifact paths.
- Limitations and claim boundaries stated explicitly.
- Submission blockers clearly listed in checklist form.

### Output interpretation guidance
Interpret outputs jointly across effect estimates, uncertainty intervals, heterogeneity diagnostics, and sensitivity results. For small study counts, rare-event settings, or model-mismatch scenarios, treat asymmetry tests and pooled estimates cautiously. When assumptions are only approximately met (e.g., large-sample approximations), results should be reported with explicit caveats.

### Table 4. Output-to-decision interpretation guide
| Output type | What it tells you | What it does not guarantee | Reporting recommendation |
|---|---|---|---|
| Primary pooled/model estimate | Central tendency under stated assumptions | Universal validity across all settings | Report with assumptions and uncertainty |
| Heterogeneity metrics | Between-study variability signal | Definitive cause of heterogeneity | Pair with subgroup/sensitivity rationale |
| Bias/asymmetry checks | Potential small-study/publication-bias signal | Definitive proof of bias mechanism | Report small-k limitations explicitly |
| Sensitivity analyses | Robustness under alternate assumptions | Immunity to all model misspecification | Present scenario-wise evidence table |

## Use cases
### Demonstration dataset used for manuscript walkthrough
- Dataset reference: ``pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json``
- Rationale: fixed demonstration artifacts improve reviewer reproducibility.

### Use case 1: Primary model execution
Workflow:
- Load demonstration data and run default model configuration.
- Review primary outputs and uncertainty diagnostics.
- Export results for reproducibility archive.

Expected outputs for manuscript:
- Primary estimate summary
- Diagnostics and model-status outputs
- Exportable reporting artifacts

### Use case 2: Sensitivity and robustness check
Workflow:
- Re-run with alternative settings (e.g., heterogeneity/sensitivity options).
- Compare directional and magnitude stability of key outputs.
- Document any materially different conclusions.

Expected outputs for manuscript:
- Sensitivity-aware interpretation
- Traceable robustness evidence

### Reviewer-facing walkthrough (replicable package check)
1. Run one documented entry point (e.g., `README.md`).
2. Load a demonstration/test dataset (e.g., ``pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json``).
3. Execute default analysis and record outputs.
4. Execute at least one sensitivity rerun and compare conclusions.
5. Cross-check outputs against listed validation evidence paths.

### User tutorial and onboarding
- Start with packaged demonstration data before using external data.
- Verify exported results against on-screen summaries.
- Use checklist items before submission or release tagging.
- Report warnings and convergence caveats with analysis outputs.

### Table 5. Assumptions, diagnostics, and caution flags
| Component | Assumption | Recommended diagnostic | Caution flag |
|---|---|---|---|
| Effect model | Chosen form reflects study design and outcome scale | Residual pattern + sensitivity reruns | Large directional shifts across settings |
| Heterogeneity handling | Random-effects assumptions are plausible | Tau/I2/Q-related diagnostics | Small-k instability or extreme heterogeneity |
| Approximation regime | Large-sample approximations adequate for data context | Rare-event and small-k checks | Sparse events / unstable variance |
| Sensitivity module | Alternative settings should not reverse core interpretation without explanation | Structured scenario comparison table | Inference changes without transparent rationale |

## Discussion
The primary strength is operational transparency: workflow steps, outputs, and evidence references are all documented in one package-level path. This improves reviewer auditability and reduces undocumented handoffs.

### Limitations and claim boundaries (review-informed)
- The package is not framed as a universal replacement for all existing methods/tools.
- Claims are restricted to artifact-supported local workflows.
- Interpretation quality remains dependent on methodological fit and data quality.
- Repository URL and DOI metadata are still required before submission.

## Conclusions
GWAM is structured for reproducible modeling workflows with explicit evidence mapping and bounded claims. With final repository/DOI metadata completed, the package is prepared for software-tool reporting.

<!-- FLOWCHART_BLOCK_START -->
## Workflow Figure Blueprint

### Figure FA1. End-to-end analytical flowchart
Recommended node sequence:
1. Data input and schema checks
2. Model setup and assumptions
3. Primary estimation
4. Diagnostics and heterogeneity review
5. Sensitivity and robustness analysis
6. Export, reporting, and reproducibility checks

Design specifications:
- Use clean vector geometry, no decorative backgrounds.
- Keep labels short, method-focused, and assumption-aware.
- Ensure grayscale legibility for print workflows.
- Keep scientific claims in manuscript text/tables; flowchart is explanatory only.

Proposed figure files:
- `figures/figure00_workflow_flowchart.png` (working draft)
- `figures/figure00_workflow_flowchart.tiff` (submission-ready raster)
- `figures/figure00_workflow_flowchart.eps` (submission-ready vector)

### Infographic-style quick panel
If needed, add a one-panel “study at a glance” infographic summarizing:
- Problem context
- Workflow contribution
- Validation evidence anchor
- Claim boundary statement
<!-- FLOWCHART_BLOCK_END -->

## Figures and visual walkthrough
Figure 1. Full-page interface or primary run overview.

Figure 2. Model configuration and parameter workflow.

Figure 3. Primary results and diagnostics view.

Figure 4. Sensitivity/robustness view.

Figure 5. Export/reporting workflow view.

Figure assets should be generated from representative full-page runs and mapped to Figure 1-5 during final submission.

### Submission figure files (separate, 300 DPI; screenshots only)
- `figures/figure01_overview_fullpage.tiff` and `figures/figure01_overview_fullpage.eps`
- `figures/figure02_model_fullpage.tiff` and `figures/figure02_model_fullpage.eps`
- `figures/figure03_results_fullpage.tiff` and `figures/figure03_results_fullpage.eps`
- `figures/figure04_bias_fullpage.tiff` and `figures/figure04_bias_fullpage.eps`
- `figures/figure05_report_fullpage.tiff` and `figures/figure05_report_fullpage.eps`

### Table 3. Reproducibility and submission readiness map
| Item | Local artifact | Current status | Action before external submission |
|---|---|---|---|
| Example walkthrough dataset | `pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json` | Present | Verify rerun on clean machine |
| Validation summary | `f1000_artifacts/validation_summary.md` | Present | Confirm numbers and paths |
| User walkthrough | `f1000_artifacts/tutorial_walkthrough.md` | Present | Align screenshots/captions |
| Repository metadata | `[TO_BE_ADDED_GITHUB_OR_GITLAB_URL]` | Placeholder | Replace after final tagging |
| DOI metadata | `[TO_BE_ADDED_ZENODO_DOI]` | Placeholder | Replace after Zenodo archive creation |

## Software availability
- Local source package: `GWAM` under `C:\Models`.
- Public repository (placeholder): `[TO_BE_ADDED_GITHUB_OR_GITLAB_URL]`
- Zenodo DOI (placeholder): `[TO_BE_ADDED_ZENODO_DOI]`
- Version: development build (version label to be finalized)
- Reproducibility walkthrough: `f1000_artifacts/tutorial_walkthrough.md`
- Validation summary: `f1000_artifacts/validation_summary.md`
- License: see package `LICENSE` file.
- Note: repository and DOI placeholders are intentionally retained until release archival is finalized.
## Data availability
No new participant-level clinical data were generated for this software article package. Example dataset for reviewer walkthrough: `pairwise70_benchmark_sensitivity\ctgov_linkage_raw\CD014089_data.json`. Additional project data assets, where present, remain available within the local package tree.
## Reporting guidelines
Real-peer-review-aligned checklist included: `F1000_Submission_Checklist_RealReview.md`.

## Declarations
### Competing interests
The author declares that no competing interests were disclosed.

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
The author acknowledges contributors to open statistical methods and reproducible software engineering practices.

## References
1. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188.
2. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Statistics in Medicine. 2002;21(11):1539-1558.
3. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
4. Guyatt GH, Oxman AD, Vist GE, et al. GRADE: an emerging consensus on rating quality of evidence and strength of recommendations. BMJ. 2008;336(7650):924-926.
5. Fay C, Rochette S, Guyader V, Girard C. Engineering Production-Grade Shiny Apps. Chapman and Hall/CRC. 2022.
