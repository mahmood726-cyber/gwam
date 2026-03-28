Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

GWAM: Ghost-Weighted Aggregate Meta-Analysis for Registry-Based Publication Bias Sensitivity

Can trial registry data provide a quantitative sensitivity framework for publication bias in meta-analysis? GWAM classifies ClinicalTrials.gov records as published, results-posted, or ghost protocols across 7,467 meta-analyses from 398 Cochrane reviews, computing enrollment-weighted integrity ratios for each drug-condition pair. The method scales published pooled estimates by the integrity ratio lambda, with simulation and hierarchical variance-propagation extensions quantifying uncertainty about unobserved effects. For escitalopram in depression, the adjusted OR on the log scale fell to 0.025 (95% CI negative 0.041 to 0.087) versus the published random-effects estimate of 0.255, with lambda of 0.096. Simulation at 2,000 replications confirmed minimal bias at the null and near-nominal hierarchical coverage, while across the Pairwise70 benchmark 7.9 percent of meta-analyses with effects exceeding 0.10 were attenuated below that threshold. GWAM provides transparent registry-anchored sensitivity analysis that complements existing statistical bias corrections. A limitation is that incomplete ClinicalTrials.gov publication linkage inflates ghost classification, making lambda a lower bound on true publication fraction.

Outside Notes

Type: methods
Primary estimand: Integrity ratio (lambda)
App: GWAM v1.0.5
Data: ClinicalTrials.gov registry + 398 Cochrane reviews (Pairwise70)
Code: https://github.com/mahmood726-cyber/gwam
Version: 1.0.5
Validation: DRAFT

References

1. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. BMJ. 1997;315(7109):629-634.
2. Duval S, Tweedie R. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
