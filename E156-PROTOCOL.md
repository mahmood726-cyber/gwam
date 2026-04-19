# E156 Protocol — `GWAM`

This repository is the source code and dashboard backing an E156 micro-paper on the [E156 Student Board](https://mahmood726-cyber.github.io/e156/students.html).

---

## `[63]` GWAM: Ghost-Weighted Aggregate Meta-Analysis for Registry-Based Publication Bias Sensitivity

**Type:** methods  |  ESTIMAND: Integrity ratio (lambda)  
**Data:** ClinicalTrials.gov registry + 398 Cochrane reviews (Pairwise70)

### 156-word body

Can trial registry data provide a quantitative sensitivity framework for publication bias in meta-analysis? GWAM classifies ClinicalTrials.gov records as published, results-posted, or ghost protocols across 7,467 meta-analyses from 398 Cochrane reviews, computing enrollment-weighted integrity ratios for each drug-condition pair. The method scales published pooled estimates by the integrity ratio lambda, with simulation and hierarchical variance-propagation extensions quantifying uncertainty about unobserved effects. For escitalopram in depression, the adjusted OR on the log scale fell to 0.025 (95% CI negative 0.041 to 0.087) versus the published random-effects estimate of 0.255, with lambda of 0.096. Simulation at 2,000 replications confirmed minimal bias at the null and near-nominal hierarchical coverage, while across the Pairwise70 benchmark 7.9 percent of meta-analyses with effects exceeding 0.10 were attenuated below that threshold. GWAM provides transparent registry-anchored sensitivity analysis that complements existing statistical bias corrections. A limitation is that incomplete ClinicalTrials.gov publication linkage inflates ghost classification, making lambda a lower bound on true publication fraction.

### Submission metadata

```
Corresponding author: Mahmood Ahmad <mahmood.ahmad2@nhs.net>
ORCID: 0000-0001-9107-3704
Affiliation: Tahir Heart Institute, Rabwah, Pakistan

Links:
  Code:      https://github.com/mahmood726-cyber/GWAM
  Protocol:  https://github.com/mahmood726-cyber/GWAM/blob/main/E156-PROTOCOL.md
  Dashboard: https://mahmood726-cyber.github.io/gwam/

References (topic pack: publication bias / selection):
  1. Egger M, Davey Smith G, Schneider M, Minder C. 1997. Bias in meta-analysis detected by a simple, graphical test. BMJ. 315(7109):629-634. doi:10.1136/bmj.315.7109.629
  2. Duval S, Tweedie R. 2000. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 56(2):455-463. doi:10.1111/j.0006-341X.2000.00455.x

Data availability: No patient-level data used. Analysis derived exclusively
  from publicly available aggregate records. All source identifiers are in
  the protocol document linked above.

Ethics: Not required. Study uses only publicly available aggregate data; no
  human participants; no patient-identifiable information; no individual-
  participant data. No institutional review board approval sought or required
  under standard research-ethics guidelines for secondary methodological
  research on published literature.

Funding: None.

Competing interests: MA serves on the editorial board of Synthēsis (the
  target journal); MA had no role in editorial decisions on this
  manuscript, which was handled by an independent editor of the journal.

Author contributions (CRediT):
  [STUDENT REWRITER, first author] — Writing – original draft, Writing –
    review & editing, Validation.
  [SUPERVISING FACULTY, last/senior author] — Supervision, Validation,
    Writing – review & editing.
  Mahmood Ahmad (middle author, NOT first or last) — Conceptualization,
    Methodology, Software, Data curation, Formal analysis, Resources.

AI disclosure: Computational tooling (including AI-assisted coding via
  Claude Code [Anthropic]) was used to develop analysis scripts and assist
  with data extraction. The final manuscript was human-written, reviewed,
  and approved by the author; the submitted text is not AI-generated. All
  quantitative claims were verified against source data; cross-validation
  was performed where applicable. The author retains full responsibility for
  the final content.

Preprint: Not preprinted.

Reporting checklist: PRISMA 2020 (methods-paper variant — reports on review corpus).

Target journal: ◆ Synthēsis (https://www.synthesis-medicine.org/index.php/journal)
  Section: Methods Note — submit the 156-word E156 body verbatim as the main text.
  The journal caps main text at ≤400 words; E156's 156-word, 7-sentence
  contract sits well inside that ceiling. Do NOT pad to 400 — the
  micro-paper length is the point of the format.

Manuscript license: CC-BY-4.0.
Code license: MIT.

SUBMITTED: [ ]
```


---

_Auto-generated from the workbook by `C:/E156/scripts/create_missing_protocols.py`. If something is wrong, edit `rewrite-workbook.txt` and re-run the script — it will overwrite this file via the GitHub API._