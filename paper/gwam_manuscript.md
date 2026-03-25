# Ghost-Weighted Aggregate Meta-analysis (GWAM): A Registry-Based Sensitivity Framework for Publication Bias

**Running title:** GWAM: Registry-Based Publication Bias Sensitivity

**Authors:** Mahmood Ahmad^1^ [Corresponding author]

**Affiliation:** ^1^ Independent Researcher, London, United Kingdom

**Corresponding author email:** mahmood.saffari@outlook.com

**ORCID:** 0009-0003-7781-4478

**Keywords:** publication bias, meta-analysis, trial registry, ClinicalTrials.gov, sensitivity analysis, ghost protocols, integrity ratio

---

## Abstract

**Background:** Publication bias threatens evidence synthesis, yet existing corrections rely on statistical signatures within the published literature. The fraction of registered trials that remain unpublished is rarely incorporated quantitatively. We introduce Ghost-Weighted Aggregate Meta-analysis (GWAM), a registry-based sensitivity framework that uses ClinicalTrials.gov enrollment data to reweight published pooled estimates.

**Methods:** GWAM classifies registered trials as published (PubMed identifier [PMID]-linked), results-posted, or ghost protocols (completed, no publication or results). An enrollment-weighted integrity ratio (lambda) quantifies the published fraction. The GWAM estimate scales the published effect by lambda; a simulation layer and hierarchical variance-propagation extension quantify uncertainty. We applied GWAM to escitalopram for depression (18 trials, log-OR; positive values indicate treatment superiority) and pregabalin for neuropathic pain (35 trials, log-RR), validated via simulation (2,000 replications, three scenarios), and benchmarked across 7,467 meta-analyses from 398 Cochrane reviews.

**Results:** For escitalopram, lambda was 0.096; the GWAM-adjusted log-OR was 0.025 (95% SI: -0.041, 0.087) versus the published random-effects estimate of 0.255. For pregabalin, lambda was 0.014 (likely reflecting incomplete PMID linkage). In simulation, GWAM showed minimal bias at mu_true = 0 but increasing negative bias at mu_true > 0 (-0.162 at 0.2), with raw coverage declining to 7.1%. The hierarchical extension provided wider intervals with near-nominal coverage. Across the Pairwise70 benchmark (median lambda = 0.871), 7.9% of meta-analyses with |RE| > 0.10 were attenuated below that threshold.

**Conclusions:** GWAM provides transparent, reproducible sensitivity analysis anchored in registry data. It is best interpreted as a conservative sensitivity framework complementing existing methods. Performance depends critically on ClinicalTrials.gov publication linkage quality.

---

## Introduction

Publication bias -- the selective reporting of studies based on the statistical significance or direction of results -- distorts the evidence base that informs clinical practice and policy [1,2]. Meta-analyses, despite pooling multiple studies, inherit this distortion when the underlying literature is non-representative of all conducted research [3].

Several statistical methods address publication bias, including funnel plot asymmetry tests [4], trim-and-fill [5], selection models [6,7], PET-PEESE regression [8,9], and Copas selection models [10]. These approaches share a common limitation: they infer the presence and magnitude of missing studies from statistical properties (e.g., asymmetry, small-study effects) of the *observed* literature. When bias operates through mechanisms other than small-study effects -- such as selective non-publication of entire large trials -- funnel-based corrections may be insufficient [11].

Trial registries offer a fundamentally different information source. In 2004, the International Committee of Medical Journal Editors (ICMJE) announced a requirement for prospective registration of clinical trials, effective from 2005 [12]. The FDA Amendments Act of 2007 (FDAAA) mandates results posting for applicable clinical trials of FDA-regulated products, though this mandate applies only to trials of approved or cleared products and does not cover Phase I trials, unapproved products not yet submitted for approval, or trials conducted entirely outside US jurisdiction [13]. ClinicalTrials.gov, the largest registry, now contains over 500,000 study records as of the access date. Comparing registered protocols against published literature reveals trials that were registered and completed but never published their results -- termed "ghost protocols" in this paper (distinct from the medical writing concept of "ghostwriting") [14,15].

Turner et al. [16] documented discrepancies between FDA-submitted and published antidepressant trial results, demonstrating that the published literature substantially overestimated drug efficacy. Their approach -- comparing published pooled estimates against the full evidence base known to the FDA -- is the conceptual foundation of GWAM. We propose to automate and generalize this comparison using publicly available registry data rather than proprietary regulatory submissions.

We introduce Ghost-Weighted Aggregate Meta-analysis (GWAM), a sensitivity framework that uses registry enrollment data to quantify the fraction of the total evidence base represented by published studies. Rather than attempting to estimate missing study effects from funnel plot shape, GWAM computes an enrollment-weighted integrity ratio (lambda) and applies it as a deterministic scaling factor, with a simulation layer and a hierarchical variance-propagation (HVP) extension to quantify uncertainty about unobserved effects.

GWAM is designed as a *complement* to existing bias-adjustment methods, not a replacement. It provides a registry-anchored external benchmark that can be compared against internally derived corrections. Fig 1 illustrates the overall pipeline. In the following sections, we describe the method, apply it to two clinical examples, validate via simulation, and benchmark across 398 Cochrane reviews.

[INSERT FIG 1 ABOUT HERE]

---

## Methods

### Trial registry pipeline

For each drug-condition pair, we queried ClinicalTrials.gov (accessed 2026-02-14) via the Clinical Trials API v2 to retrieve all registered interventional studies of the target drug for the target condition. Trials were filtered to those with status "Completed" and available enrollment counts. Terminated and withdrawn trials were excluded from the primary analysis because their non-completion makes the "ghost" interpretation (completed-but-unpublished) inapplicable; however, this exclusion may remove trials stopped for futility or safety that would plausibly have null or negative results.

Each trial was classified into one of three mutually exclusive categories based on automated linkage:

1. **Published (PMID-linked):** A PubMed identifier (PMID) was associated with the trial record via ClinicalTrials.gov reference fields (types RESULT or PRIMARY) or manual curation.
2. **Results-posted (no PMID):** The trial had posted summary results on ClinicalTrials.gov but had no linked journal publication.
3. **Ghost protocol:** The trial had neither a linked PMID nor posted results.

**Linkage completeness caveat:** ClinicalTrials.gov PMID linkage depends on sponsors voluntarily updating registry records. Studies have found that only approximately 50% of published trials have their PMIDs linked in CT.gov reference fields [17]. Consequently, the ghost classification is an *upper bound* on the true rate of non-publication: some trials classified as "ghosts" or "results-only" may in fact have journal publications that are not linked in the registry. This means lambda_pmid is a *lower bound* on the true publication fraction, and GWAM's attenuation may overestimate the impact of non-publication. We partially mitigated this through automated PubMed title/author searches (see Scripts), but residual misclassification remains a key limitation.

A small number of trials in each candidate scan fell outside these three categories (e.g., trials with associated publications that were not primary results papers, or trials with ambiguous linkage status). These were excluded from the GWAM analysis but are reported in the candidate scan for transparency. The number of excluded (unclassifiable) trials is reported per drug-condition pair in Table 1.

### Candidate scan

We screened six drug-condition pairs to identify applications with substantial ghost protocol rates (Table 1). Selection criteria for detailed GWAM analysis were: (a) ghost protocol rate exceeding 20%, (b) availability of a published meta-analysis providing a pooled effect estimate, and (c) sufficient trial heterogeneity to make publication bias a plausible concern.

### GWAM model

For both applications, positive values on the log-OR and log-RR scales indicate treatment superiority (greater response rate for escitalopram; greater pain reduction for pregabalin), so that attenuation toward zero represents reduction of apparent treatment benefit.

Let W_pmid, W_ro, and W_ghost denote the total enrollment-weighted contributions from published, results-only, and ghost trials, respectively, with W_total = W_pmid + W_ro + W_ghost. The integrity ratio is defined as:

> lambda_pmid = W_pmid / W_total &nbsp;&nbsp;&nbsp;&nbsp; (1)

which measures the fraction of total registered enrollment represented by published studies. A second ratio, lambda_non_ghost = (W_pmid + W_ro) / W_total, quantifies the fraction with any available results (published or posted).

The GWAM estimate is a three-way enrollment-weighted average:

> mu_GWAM = (W_pmid * mu_pub + W_ro * mu_ro + W_ghost * mu_ghost) / W_total &nbsp;&nbsp;&nbsp;&nbsp; (2)

where mu_pub is the published random-effects pooled estimate, mu_ro is the assumed mean effect for results-only trials, and mu_ghost is the assumed mean effect for ghost protocols.

Under the *null assumption* that unpublished evidence is centered at zero (mu_ro = 0, mu_ghost = 0), this simplifies to:

> mu_GWAM_null = lambda_pmid * mu_pub &nbsp;&nbsp;&nbsp;&nbsp; (3)

This deterministic point estimate represents the published effect diluted by the proportion of the total evidence base that is published. When lambda_pmid is small (i.e., most enrolled participants are in unpublished trials), the GWAM estimate is substantially attenuated.

The GWAM standard error incorporates uncertainty from both the published estimate and the assumed per-study ghost effect distributions:

> SE_GWAM = sqrt(W_pmid^2 * SE_pub^2 + sum_j w_ro_j^2 * sigma_ro^2 + sum_k w_ghost_k^2 * sigma_ghost^2) / W_total &nbsp;&nbsp;&nbsp;&nbsp; (4)

where SE_pub is the standard error of the published pooled estimate, w_ro_j is the enrollment weight of the j-th results-only study, sigma_ro is the assumed SD for results-only effects (default 0.10), w_ghost_k is the enrollment weight of the k-th ghost study, and sigma_ghost is the assumed standard deviation of ghost effects (default 0.10). This per-study formulation correctly aggregates independent variance contributions from each ghost study; using the aggregate formula (1 - lambda)^2 * sigma_ghost^2 would overestimate the ghost variance component. When sigma_ghost = 0, the SE reduces to (W_pmid / W_total) * SE_pub = lambda_pmid * SE_pub, preserving the z-statistic and therefore statistical significance. Including sigma_ghost breaks this invariance.

**Design note:** Weights are based on enrollment counts, not inverse-variance, because variance estimates are unavailable for ghost protocols and many results-only studies. Enrollment-based weighting reflects the registry "universe" denominator for the integrity ratio and is a deliberate design choice. However, enrollment is an imperfect proxy for study informativeness: trials with large enrollment may have high dropout, low event rates, or non-informative study designs (e.g., pharmacokinetic studies). A sensitivity analysis using study-count-based lambda (n_published / n_total) is recommended alongside the enrollment-weighted version.

### Monte Carlo simulation layer

To propagate uncertainty about unobserved effects, we add a simulation layer. Ghost effects are drawn from N(mu_ghost, sigma_ghost) and results-only effects from N(mu_ro, sigma_ro), with default parameters mu_ghost = mu_ro = 0, sigma_ghost = sigma_ro = 0.10. The GWAM simulation estimate is computed as:

> mu_GWAM^(b) = (W_pmid * mu_pub + sum_j w_ro_j * epsilon_ro_j^(b) + sum_k w_ghost_k * epsilon_ghost_k^(b)) / W_total &nbsp;&nbsp;&nbsp;&nbsp; (5)

for b = 1, ..., B Monte Carlo draws, yielding a simulation-based 95% interval (SI) from the 2.5th and 97.5th percentiles of the distribution. Note that the simulation layer treats mu_pub as a fixed constant rather than propagating its sampling uncertainty (SE_pub); the HVP extension addresses this limitation.

### Hierarchical variance-propagation sensitivity extension

To explore sensitivity to assumptions about unobserved effects, we compute a hierarchical variance-propagation (HVP) estimate. This is *not* a formal Bayesian posterior derived from a likelihood-prior conjugacy, but rather an enrollment-weighted average with two-level variance propagation designed to bracket plausible uncertainty ranges.

The HVP point estimate is the enrollment-weighted average of Equation (2), using assumed mean effects for ghost (mu_ghost = 0) and results-only (mu_ro = 0 or observed) studies. The HVP variance propagates uncertainty at two levels:

- **Level 1 (within-study):** Each ghost study contributes variance sigma_ghost^2 weighted by its squared enrollment weight; each results-only study contributes its observed SE^2 (if available) or sigma_ro^2 otherwise.
- **Level 2 (between-group):** Uncertainty about the *mean* ghost and results-only effects, parameterized by hyperprior SDs (sigma_ghost_mu, sigma_ro_mu, both defaulting to 0.15).

> Var_HVP = (W_pmid^2 * SE_pub^2 + sum_j w_ro_j^2 * sigma_ro_j^2 + W_ro_unobs^2 * sigma_ro_mu^2 + sum_k w_ghost_k^2 * sigma_ghost^2 + W_ghost^2 * sigma_ghost_mu^2) / W_total^2 &nbsp;&nbsp;&nbsp;&nbsp; (6)

The model evaluates a 4 x 4 sensitivity grid over ghost_sigma in {0.05, 0.10, 0.20, 0.30} and results_only_sigma in {0.05, 0.10, 0.20, 0.30}, producing summaries (mean, SD, 95% interval, Pr(effect > 0)) for each cell. The reference cell is selected as the grid point whose posterior SD is closest to the median SD across all 16 cells, providing a representative mid-range summary.

**Point estimate invariance:** Because the point estimate is a fixed enrollment-weighted average with zero-centered assumptions for unobserved effects, the HVP mean is invariant to the variance parameters. Only the interval width and Pr(effect > 0) vary across the grid. Consequently, this grid explores only the *dispersion* of ghost effects, not their *direction*. A supplementary analysis varying mu_ghost across a clinically meaningful range (e.g., {-0.2, -0.1, 0.0, 0.1}) would provide more informative directional sensitivity and is recommended for applied use.

### Simulation validation

To evaluate operating characteristics, we conducted a calibration simulation with 2,000 replications (using seed 42 for reproducibility) for each of three ground-truth scenarios (mu_true in {0.0, 0.1, 0.2}). For each replication, we:

1. Generated k study-level log-odds ratios from N(mu_true, tau^2) with tau^2 = 0.1, using enrollment weights from the escitalopram registry structure. A conditional Haldane-Anscombe continuity correction (0.5 added to all cells) was applied to studies with at least one zero cell; studies without zero cells received no correction to avoid unnecessary bias toward the null.
2. Simulated selective publication using a dichotomous selection model: each study was "published" with probability p_sig = 0.95 if its result was statistically significant and in the positive direction, or with calibrated probability p_nonsig otherwise. The calibration procedure ensured that the simulated lambda_pmid matched the observed lambda_pmid (0.096) within tolerance 0.03. This selection model is a simplification; in practice, publication probability may vary by effect magnitude, direction, and marginality of significance.
3. Applied five estimation methods to the selectively published studies: (a) standard random-effects (RE) using DerSimonian-Laird with Wald z-based 95% CIs, (b) GWAM null (lambda-scaled), (c) HVP-GWAM (hierarchical variance-propagation), (d) oracle inverse-probability-weighted (IPW) RE using the true selection probabilities, and (e) PET-PEESE (minimum k = 3, with alpha = 0.10 two-sided for the PET-to-PEESE switching rule following Stanley & Doucouliagos 2014 [8]).
4. Computed bias, RMSE, raw 95% CI coverage, and calibrated CI coverage using empirically derived multipliers.

**Critical caveat:** Under mu_true = 0, GWAM has a structural advantage because its null assumption (unobserved effects centered at zero) coincides with the ground truth. Performance at mu_true > 0 reveals the cost of this assumption: GWAM increasingly underestimates the true effect as mu_true departs from zero. All comparator CIs (RE, IPW, PET-PEESE) use Wald z-based intervals rather than the Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment; since all methods use the same interval construction, the comparison is internally consistent, though HKSJ would be preferred in applied practice.

**Note on PET-PEESE minimum k:** PET-PEESE was applied with a minimum of k = 3 studies, which is below the k >= 10 recommended by Stanley & Doucouliagos [8] for reliable estimates. Its high failure rate and RMSE in the simulation partly reflect this low threshold. A sensitivity analysis with k >= 5 is reported in the supplementary materials.

### CI calibration

Raw confidence intervals for all methods may under- or over-cover due to the simulation-based nature of the correction. We applied empirical CI width multipliers derived from a calibration run at mu_true = mu_pub (the published estimate), optimizing each method's multiplier to achieve 95% coverage. These multipliers were then applied to all three scenarios as an out-of-sample test. The calibration is inherently conservative for mu_true values far from the calibration point, producing over-coverage (excessively wide intervals) rather than under-coverage.

### Pairwise70 benchmark

To assess GWAM's attenuation effect at scale, we applied the proxy-lambda pipeline to 7,467 binary-outcome meta-analyses extracted from 398 Cochrane systematic reviews in the Pairwise70 corpus [18]. The Pairwise70 corpus comprises curated .rda data files from published Cochrane intervention reviews; the 7,467 analysis subset was obtained by filtering to meta-analyses with at least 2 studies, a valid DerSimonian-Laird random-effects estimate, and a computable proxy lambda.

For each meta-analysis, we linked the contributing trials to ClinicalTrials.gov via NCT identifier matching and query-derived linkage, computed a proxy lambda (denoted lambda_hat to distinguish it from the directly measured lambda_pmid), and calculated the GWAM-adjusted estimate. Lambda_hat sources were (in priority order): exact Pairwise70 NCT linkage, ClinicalTrials.gov query-derived non-ghost ratio, transport proxy from covariate similarity, and a conservative default of 0.70 when no linkage was available. Lambda_hat values were clipped to [0.20, 0.95] to avoid extreme adjustments.

We report the fraction of meta-analyses where the GWAM-adjusted absolute effect crossed below clinically relevant thresholds (|mu| < 0.10 and |mu| < 0.20 on the log-OR scale). Review-clustered bootstrap confidence intervals (400 iterations) account for within-review correlation.

### Software and reproducibility

All analyses were conducted in Python 3.13 using NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.2, and matplotlib 3.8 (S1 Code). Random number generation used NumPy's default_rng with seed 42 for reproducibility. ClinicalTrials.gov data were accessed via the Clinical Trials API v2 on 2026-02-14; pre-fetched registry CSV files are included in the repository as the authoritative snapshot. The complete analysis code, registry data, and figure generation scripts are available at https://github.com/mahmoodsaffari/GWAM. Note: the code uses the legacy internal label "Bayesian" for the HVP extension module (e.g., `model_gwam_bayesian.py`); all user-facing outputs and figures use the HVP terminology described in this paper. This methodological study was not pre-registered; the analysis plan evolved iteratively during development. GWAM is a methodological framework paper, not a systematic review; a PRISMA checklist is therefore not applicable, though the Pairwise70 benchmark methods describe the inclusion/exclusion criteria for the meta-analysis corpus.

---

## Results

### Candidate scan

Table 1 presents the ghost protocol rates for six drug-condition pairs. Fig 2 shows the ghost protocol rates across drug-condition pairs. Ghost rates ranged from 16.7% (quetiapine/depression) to 25.7% (pregabalin/neuropathic pain). Escitalopram for major depressive disorder (24.7% ghost rate, 73 completed trials) and pregabalin for neuropathic pain (25.7%, 35 completed trials) were selected for detailed analysis based on having ghost rates exceeding 20% and available published meta-analyses. Of the 73 completed escitalopram trials, 18 met the inclusion criteria after applying study-level design filters (randomization, active comparator scope, treatment purpose, and depression-response outcome within 4--16 weeks).

[INSERT FIG 2 ABOUT HERE]

**Table 1. Candidate scan: Ghost protocol rates across drug-condition pairs (ClinicalTrials.gov, accessed 2026-02-14).**

| Drug / Condition | N completed | With PMID | Results-only (no PMID) | Ghost protocols | Unclassified | Ghost rate (%) |
|---|---|---|---|---|---|---|
| Sertraline / Depression | 43 | 11 | 19 | 9 | 4 | 20.9 |
| Escitalopram / Depression | 73 | 11 | 35 | 18 | 9 | 24.7 |
| Duloxetine / Depression | 37 | 5 | 23 | 7 | 2 | 18.9 |
| Quetiapine / Depression | 18 | 0 | 13 | 3 | 2 | 16.7 |
| Pregabalin / Neuropathic Pain | 35 | 1 | 18 | 9 | 7 | 25.7 |
| Oseltamivir / Influenza | 25 | 1 | 14 | 5 | 5 | 20.0 |

*Note:* The "Unclassified" column shows trials that fell outside the three-way PMID/results-only/ghost classification (e.g., trials with associated publications that were not primary results papers, or trials with ambiguous linkage status). Ghost rate = Ghost / N completed. This denominator is conservative; using only classified trials would yield modestly higher ghost rates. Ghost counts reflect ClinicalTrials.gov PMID linkage, which is known to be incomplete [17]; some ghosts may in fact have journal publications not linked in the registry.

### Application 1: Escitalopram for major depressive disorder

From 73 completed escitalopram/depression trials identified in the candidate scan, 18 met the inclusion criteria for GWAM analysis after applying study-level filters (requirement for enrollment data, active comparator scope, randomized treatment-purpose design, and relevant depression-response outcome within 4--16 weeks). Of these 18 trials:

- 4 (22.2%) had published PMID-linked results (total enrollment: 471)
- 9 (50.0%) had registry-posted results without a journal publication (enrollment: 3,832)
- 5 (27.8%) were ghost protocols (enrollment: 610)

The integrity ratio lambda_pmid = 471 / 4,913 = 0.096, indicating that published studies represented less than 10% of total enrolled participants. The non-ghost ratio lambda_non_ghost = 4,303 / 4,913 = 0.876.

The published random-effects pooled estimate was log-OR = 0.255 (SE_pub = 0.096; 95% CI: 0.066, 0.443; OR = 1.29), derived from the PMID-linked active-comparator trials in the escitalopram registry subset using DerSimonian-Laird random effects. This estimate is consistent with the active-comparator subgroup in Cipriani et al. [19], who reported escitalopram as among the most efficacious antidepressants in a network meta-analysis. Under the GWAM null assumption, the adjusted estimate was:

> mu_GWAM = lambda_pmid * mu_pub = 0.096 * 0.255 = 0.024 (OR = 1.025)

The Monte Carlo simulation (B = 2,000, sigma_ghost = sigma_ro = 0.10) yielded a mean of 0.025 with 95% SI: [-0.041, 0.087], indicating that the published effect of escitalopram was not robust to the possibility that unpublished trials showed null or negative results (Fig 3). The slight difference between the deterministic estimate (0.024) and the simulation mean (0.025) arises from the use of rounded vs. exact enrollment weights.

[INSERT FIG 3 ABOUT HERE]

The HVP sensitivity analysis (Fig 4) yielded a point estimate of -0.098 (95% interval: -0.341, 0.145) at the reference cell (ghost_sigma = 0.05, results_only_sigma = 0.20), with Pr(effect > 0) = 0.215. The negative estimate arises because five results-only trials had extractable observed effect sizes from ClinicalTrials.gov posted results, and these observed effects were on average negative, pulling the estimate below zero. Across the full 4 x 4 grid, point estimates were invariant (see Methods), while Pr(effect > 0) ranged from 0.201 (tightest dispersion) to 0.233 (widest dispersion). This narrow range is a consequence of the point estimate invariance property; varying the assumed ghost effect *mean* rather than the *dispersion* would produce more informative variation. The corresponding grid for pregabalin/neuropathic pain is in S1 Table.

[INSERT FIG 4 ABOUT HERE]

### Application 2: Pregabalin for neuropathic pain

For pregabalin/neuropathic pain, 35 completed trials were identified (all 35 met inclusion criteria):

- 1 (2.9%) had a published PMID (enrollment: 104)
- 25 (71.4%) had posted results only (enrollment: 5,720)
- 9 (25.7%) were ghost protocols (enrollment: 1,597)

Lambda_pmid = 104 / 7,421 = 0.014, indicating that less than 2% of enrolled participants were in the single published trial. Lambda_non_ghost = 5,824 / 7,421 = 0.785.

**Important caveat:** The extreme lambda of 0.014 almost certainly reflects incomplete ClinicalTrials.gov PMID linkage rather than true non-publication. Pregabalin has extensive published literature with multiple positive RCTs and Cochrane reviews [20]. Many Pfizer-sponsored pregabalin trials (particularly those predating the 2008 FDAAA results-posting mandate) have journal publications that are not linked in ClinicalTrials.gov reference fields. This application therefore demonstrates a limitation of the automated linkage pipeline rather than a substantive finding about pregabalin efficacy. We present it for methodological transparency and to illustrate the sensitivity of GWAM to linkage quality.

The published estimate was log-RR = 0.916 (SE_pub = 0.148; 95% CI: 0.625, 1.207; RR = 2.50), based on the single PMID-linked trial in the registry subset. The reference meta-analysis is Derry et al. [20], who found pregabalin effective for neuropathic pain (RR for >= 50% pain reduction approximately 1.3--2.7 depending on dose and condition). The single-trial estimate used here (RR = 2.50) is higher than the multi-trial Cochrane pooled estimate, likely reflecting the selection of a particularly favorable published trial. The GWAM-adjusted estimate was:

> mu_GWAM = 0.014 * 0.916 = 0.013 (RR = 1.013)

The 95% SI was [-0.029, 0.058], and the HVP estimate was -0.234 (95% interval: -0.384, -0.085), with Pr(effect > 0) = 0.001 (reference cell: ghost_sigma = 0.05, results_only_sigma = 0.30). This extreme attenuation reflects the dominance of unpublished evidence in the registry data as classified by automated linkage. Given the known linkage incompleteness, these results should not be interpreted as evidence that pregabalin is ineffective.

### Simulation validation

Table 2 summarizes simulation results (escitalopram weight structure, 2,000 replications per scenario).

**Table 2. Simulation operating characteristics under three ground-truth scenarios (escitalopram calibration, 2,000 replications).**

| mu_true | Method | Bias | RMSE | Coverage (raw) | Coverage (calibrated) | n_valid |
|---|---|---|---|---|---|---|
| 0.0 | RE | 0.158 | 0.411 | 84.8% | 95.3% | 2000 |
| 0.0 | GWAM | 0.013 | 0.044 | 92.4% | 100% | 2000 |
| 0.0 | HVP-GWAM | 0.013 | 0.031 | 100% | 100% | 2000 |
| 0.0 | Oracle IPW | 0.042 | 0.365 | 73.6% | 92.8% | 2000 |
| 0.0 | PET-PEESE | 0.176 | 1.048 | 74.4% | 95.8% | 1271 |
| 0.1 | RE | 0.176 | 0.408 | 82.2% | 94.9% | 2000 |
| 0.1 | GWAM | -0.076 | 0.088 | 52.6% | 100% | 2000 |
| 0.1 | HVP-GWAM | -0.076 | 0.082 | 100% | 100% | 2000 |
| 0.1 | Oracle IPW | 0.041 | 0.355 | 74.8% | 92.6% | 2000 |
| 0.1 | PET-PEESE | 0.176 | 0.977 | 72.0% | 96.1% | 1369 |
| 0.2 | RE | 0.206 | 0.399 | 81.0% | 95.5% | 2000 |
| 0.2 | GWAM | -0.162 | 0.168 | 7.1% | 99.5% | 2000 |
| 0.2 | HVP-GWAM | -0.162 | 0.166 | 99.8% | 99.9% | 2000 |
| 0.2 | Oracle IPW | 0.048 | 0.347 | 75.4% | 94.9% | 2000 |
| 0.2 | PET-PEESE | 0.153 | 0.908 | 70.5% | 96.0% | 1487 |

*Note: RE = random-effects (DerSimonian-Laird); GWAM = deterministic null-scaled; HVP-GWAM = hierarchical variance-propagation extension; Oracle IPW = inverse-probability-weighted RE using true selection probabilities; PET-PEESE = precision-effect test with precision-effect estimate with standard error (minimum k = 3). n_valid = replications with finite estimates. Coverage (calibrated) uses empirical multipliers derived at mu_true = 0.*

[INSERT FIG 5 ABOUT HERE]

[INSERT FIG 6 ABOUT HERE]

Fig 5 displays bias and RMSE, and Fig 6 displays coverage across methods and scenarios.

Under mu_true = 0 (no true effect), GWAM showed minimal bias (0.013) and the lowest RMSE (0.044). However, this result is partly tautological: GWAM's null assumption (unobserved effects = 0) coincides with the ground truth, giving it a structural advantage. Standard RE had bias of 0.158, confirming substantial distortion from selective publication. The oracle IPW estimator, which uses the true (unknown in practice) selection probabilities, had higher RMSE (0.365) because it does not assume mu = 0 and must estimate the effect from the published subsample with IPW corrections. Under mu_true = 0, PET-PEESE was applicable in only 63.6% of replications (1,271/2,000); applicability increased to 68.5% (mu_true = 0.1) and 74.4% (mu_true = 0.2) as more studies reached significance. PET-PEESE showed the highest RMSE (1.048), partly reflecting its application at k = 3 (below the recommended minimum of k >= 10 [8]). The oracle IPW SE used the standard 1/sqrt(sum(w*)) formula rather than a sandwich estimator, which likely contributed to its 73.6% raw coverage.

Under mu_true = 0.1 and 0.2, GWAM exhibited increasing negative bias (-0.076 and -0.162, respectively) because its null assumption diverges from the true positive effect. Raw CI coverage dropped to 52.6% and 7.1%. The large calibration multiplier (7.75) reflects the substantial underestimation of total uncertainty by the deterministic GWAM SE formula, which captures only published-estimate sampling error and ghost dispersion but not between-study heterogeneity or lambda estimation uncertainty. Calibrated CIs restored coverage to near-nominal levels (100% and 99.5%), though at the cost of very wide intervals. The HVP-GWAM intervals showed near-100% coverage across all scenarios (100%, 100%, 99.9%), but this indicates that the intervals are uninformatively wide -- they contain every plausible mu_true value and therefore provide limited discriminatory power. Users should prefer the HVP interval or calibrated SI over the raw GWAM CI for inference, and should interpret GWAM primarily as a point-estimate sensitivity tool rather than a full inferential framework.

The pregabalin calibration showed that lambda calibration was not within tolerance (error 0.155 vs. maximum 0.030), indicating that the extreme lambda_pmid of 0.014 could not be stably achieved in the simulation's selective-publication model. This underscores that GWAM's performance guarantees are strongest when the registry structure can be faithfully reproduced in simulation. Full pregabalin simulation results, including the calibration failure diagnostics, are presented in S2 Table.

### Pairwise70 benchmark

Across 7,467 meta-analyses from 398 Cochrane reviews (filtered from the full corpus of 13,019 binary-outcome analyses from 458 review files to those with at least 2 studies, a valid RE estimate, and a computable proxy lambda_hat), the median proxy lambda_hat was 0.871 (Fig 7, left panel). The distribution was right-skewed, with 35.3% of analyses at the upper bound (0.95) -- where GWAM has minimal effect -- and 5.7% at the lower bound (0.20).

[INSERT FIG 7 ABOUT HERE]

GWAM attenuated 590 meta-analyses (7.9%; bootstrap 95% CI: 5.1%--11.0%) from |RE| >= 0.10 to |GWAM| < 0.10, and 819 (11.0%; 95% CI: 7.2%--14.9%) from |RE| >= 0.20 to |GWAM| < 0.20. Among meta-analyses with |RE| >= 0.10, the conditional attenuation rate was 9.4%, and among those with |RE| >= 0.20, it was 15.5%. The median absolute shift was 0.055 (Fig 7, right panel).

**Table 3. Pairwise70 benchmark: Attenuation rates by lambda_hat source and study count.**

| Stratum | n | Median lambda_hat | Crossed below 0.10 (% all) | Crossed below 0.20 (% all) |
|---|---|---|---|---|
| **By lambda_hat source** | | | | |
| Exact NCT linkage | 1,438 | 0.950 | 3.1% | 7.4% |
| CT.gov query linkage | 5,569 | 0.861 | 9.1% | 11.7% |
| Transport proxy | 86 | 0.950 | 10.5% | 22.1% |
| Default (0.70) | 374 | 0.700 | 8.3% | 11.0% |
| **By study count** | | | | |
| k = 2 | 2,895 | 0.894 | 6.3% | 8.9% |
| k = 3 | 1,409 | 0.880 | 7.9% | 10.1% |
| k >= 4 | 3,163 | 0.845 | 9.4% | 13.3% |

Attenuation was lowest for meta-analyses with exact NCT linkage (3.1% crossed below 0.10), which tend to have high lambda_hat values (median 0.950). This is the most trustworthy lambda_hat stratum and suggests that among well-linked meta-analyses, GWAM attenuation is modest. The benchmark is dominated by CT.gov query-derived linkage (5,569 analyses, 74.6%), which has unknown accuracy; attenuation rates in this stratum (9.1%) may partly reflect measurement error in lambda_hat rather than true publication bias. Meta-analyses with more studies (k >= 4) showed higher attenuation rates (9.4%), consistent with larger meta-analyses drawing from broader, more heterogeneous trial registries with lower average lambda_hat. Full stratified results are presented in S3 Table.

### Lambda sensitivity analysis

To assess robustness to lambda_hat computation settings, we compared base, strict, and loose configurations across all 13,019 binary-outcome analyses in the full Pairwise70 corpus (Fig 8). The sensitivity analysis used the full corpus (without the minimum-2-studies and valid-RE filters applied in the main benchmark) to maximize statistical power for comparing configurations:

[INSERT FIG 8 ABOUT HERE]

- **Base** (default lambda_hat = 0.70, bounds [0.20, 0.95]): crossed below 0.10 in 8.9%, below 0.20 in 10.0%.
- **Loose** (default = 0.90, bounds [0.40, 1.00]): 5.8% and 6.9%.
- **Strict** (default = 0.50, bounds [0.10, 0.70]): 14.5% and 14.6%.

The three-fold range (5.8%--14.5%) demonstrates moderate sensitivity to lambda_hat specification, supporting the use of sensitivity grids rather than single-point estimates.

---

## Discussion

### Summary and interpretation

GWAM provides a registry-anchored sensitivity framework for publication bias that directly incorporates the proportion of registered evidence that is published. In two clinical applications, GWAM substantially attenuated published pooled effects: escitalopram's log-OR dropped from 0.255 to 0.025 (OR: 1.29 to 1.025), and pregabalin's log-RR dropped from 0.916 to 0.013 (RR: 2.50 to 1.013). These dramatic shifts reflect the low integrity ratios (lambda_pmid = 0.096 and 0.014, respectively), indicating that the published literature represents a small fraction of the total registered evidence *as classified by ClinicalTrials.gov PMID linkage*. For pregabalin, the extreme lambda almost certainly reflects linkage incompleteness rather than true non-publication (see Limitations).

The Pairwise70 benchmark suggests that such extreme lambda values are uncommon in Cochrane meta-analyses (median lambda_hat = 0.871), and that GWAM's typical attenuation is modest: 7.9% of meta-analyses with |RE| > 0.10 crossed below that threshold. Among the most trustworthy stratum (exact NCT linkage), only 3.1% crossed below 0.10. This is consistent with GWAM functioning as a sensitivity tool rather than a universal correction -- it flags cases where the integrity ratio raises concern but does not override every published estimate.

### Relationship to existing methods

GWAM differs from existing publication bias methods in both its information source and its mechanism. Trim-and-fill [5] imputes missing studies based on funnel plot symmetry; selection models [6,7,10] estimate the selection function from the observed p-value distribution; PET-PEESE [8,9] regresses effect sizes on their standard errors; and Vevea-Hedges weight-function models [21] specify step-function selection weights. All these methods operate on the *internal* statistical properties of the published literature.

GWAM operates on an *external* data source -- the trial registry -- and makes no assumptions about funnel plot shape, p-value distributions, or small-study effects. This makes it complementary to existing methods: when both GWAM and funnel-based methods suggest attenuation, confidence in publication bias increases; when they diverge, the discrepancy itself is informative.

Turner et al. [16] performed the conceptual predecessor of GWAM by comparing published antidepressant trial results against the complete set of FDA-submitted trial results, demonstrating that published literature overestimated efficacy by approximately 32%. GWAM automates and generalizes this comparison using publicly available ClinicalTrials.gov data rather than proprietary FDA submissions, extending the approach to any registered drug-condition pair. Trinquart et al. [22] and Mathieu et al. [23] documented the registration-to-publication gap; Schmucker et al. [24] systematically reviewed non-publication rates across therapeutic areas. GWAM operationalizes these findings by converting the registration-publication gap into a quantitative correction.

Huang et al. [25] examined the distribution of p-values in clinical trial registries and its implications for publication bias assessment, providing complementary registry-based evidence. GWAM extends this registry-based perspective by using enrollment weights and a three-category classification to produce adjusted effect estimates rather than descriptive p-value analyses.

### Strengths and limitations

**Strengths.** GWAM is fully transparent and reproducible: it uses publicly available registry data, deterministic formulas, and seeded random number generation. The integrity ratio is directly interpretable (the fraction of total enrollment that is published), and the sensitivity grid communicates how conclusions change under different assumptions about unobserved effects. The method requires no access to individual patient data and can be applied to any meta-analysis where the underlying trials are registered.

**Limitations.** First, GWAM's validity depends critically on the completeness of ClinicalTrials.gov PMID linkage. Studies have found that only approximately 50% of published trials have their PMIDs linked in registry reference fields [17]. This means lambda_pmid is a lower bound on the true publication fraction, and GWAM may systematically over-correct when linkage is incomplete. The pregabalin application (lambda = 0.014) illustrates this: the extensive published pregabalin literature is poorly linked in ClinicalTrials.gov, producing an artifactually low lambda. Users should validate ghost classifications against PubMed title/author searches before drawing conclusions from GWAM-adjusted estimates.

Second, enrollment-based weights are an imperfect proxy for study informativeness. Trials with larger enrollment may not proportionally contribute more information if they have higher dropout, low event rates, or non-informative designs (e.g., pharmacokinetic or dose-finding studies). A study-count-based lambda provides a useful comparator.

Third, the null assumption (mu_ghost = 0, mu_ro = 0) is conservative under the hypothesis that publication bias inflates effects, but may over-correct when a true positive effect exists. The simulation at mu_true = 0.1 and 0.2 demonstrates this bias-variance trade-off explicitly. The HVP extension mitigates this by propagating uncertainty about ghost effect distributions, but the user must specify the assumed dispersion. Notably, results-only trials often have extractable effect estimates from ClinicalTrials.gov posted results; using these observed values (rather than assuming mu_ro = 0) in the primary analysis would improve accuracy.

Fourth, the FDAAA results-posting mandate applies only to applicable clinical trials of approved products; many ghost protocols may fall outside this scope, meaning their non-posting is legal and expected rather than evidence of selective suppression. Restricting the ghost pool to FDAAA-applicable trials would provide a more specific measure.

Fifth, ClinicalTrials.gov registration is incomplete, particularly for older trials and those conducted outside the United States. Lambda estimates may therefore be biased upward (making the published fraction appear larger) if unregistered trials exist. Conversely, registered trials with completed status may never have been intended to produce publications (e.g., pharmacokinetic studies), inflating the ghost count.

Sixth, the deterministic GWAM SE substantially underestimates total uncertainty, as evidenced by the calibration multiplier of 7.75 required to achieve nominal coverage. This indicates that the raw GWAM confidence interval should not be used for inference; users should rely on the HVP interval or calibrated SI instead.

Seventh, the proxy lambda_hat used in the Pairwise70 benchmark relies on automated NCT linkage and query-derived matching, which introduces measurement error. Only 19.3% of benchmark analyses have exact NCT linkage; the remainder use query-derived or default values of unknown accuracy.

### Implications for practice

We recommend GWAM as a supplementary sensitivity analysis for systematic reviews where the underlying trials are registered. Specifically:

1. Report the integrity ratio (lambda) alongside standard funnel plot assessments. The integrity ratio and GWAM-adjusted estimate could be reported under PRISMA Item 22 (certainty of evidence) or as a planned sensitivity analysis in the review protocol.
2. Present the GWAM-adjusted estimate as a "what-if" scenario: "If unpublished trials showed null effects, the pooled estimate would be X."
3. Use the HVP sensitivity grid to communicate the range of plausible adjusted estimates under different dispersion assumptions.
4. When GWAM attenuation and funnel-based methods agree, this strengthens the evidence for publication bias; when they disagree, investigate the source of discordance.
5. Validate ghost classifications against PubMed searches before interpreting extreme lambda values.

GWAM is not intended to replace the published estimate but to provide context about its sensitivity to publication completeness. The method is most informative when the integrity ratio is low (lambda < 0.50), indicating that a large fraction of enrolled participants are in unpublished trials.

---

## References

1. Rothstein HR, Sutton AJ, Borenstein M. Publication bias in meta-analysis: Prevention, assessment and adjustments. Chichester: Wiley; 2005.
2. Song F, Parekh S, Hooper L, et al. Dissemination and publication of research findings: an updated review of related biases. Health Technol Assess. 2010;14(8):1-193. doi:10.3310/hta14080
3. Ioannidis JP. Why most published research findings are false. PLoS Med. 2005;2(8):e124. doi:10.1371/journal.pmed.0020124
4. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. BMJ. 1997;315(7109):629-634. doi:10.1136/bmj.315.7109.629
5. Duval S, Tweedie R. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463. doi:10.1111/j.0006-341X.2000.00455.x
6. Hedges LV. Modeling publication selection effects in meta-analysis. Stat Sci. 1992;7(2):246-255. doi:10.1214/ss/1177011364
7. Vevea JL, Woods CM. Publication bias in research synthesis: sensitivity analysis using a priori weight functions. Psychol Methods. 2005;10(4):428-443. doi:10.1037/1082-989X.10.4.428
8. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78. doi:10.1002/jrsm.1095
9. Stanley TD. Limitations of PET-PEESE and other meta-analytic methods: Beyond the "cliff." Res Synth Methods. 2017;8(1):97-105. doi:10.1002/jrsm.1211
10. Copas J, Shi JQ. A sensitivity analysis for publication bias in systematic reviews. Stat Methods Med Res. 2001;10(4):251-265. doi:10.1177/096228020101000402
11. Sterne JAC, Gavaghan D, Egger M. Publication and related bias in meta-analysis: power of statistical tests and prevalence in the literature. J Clin Epidemiol. 2000;53(11):1119-1129. doi:10.1016/S0895-4356(00)00242-0
12. De Angelis C, Drazen JM, Frizelle FA, et al. Clinical trial registration: a statement from the International Committee of Medical Journal Editors. N Engl J Med. 2004;351(12):1250-1251. doi:10.1056/NEJMe048225
13. Zarin DA, Tse T, Williams RJ, Califf RM, Ide NC. The ClinicalTrials.gov results database -- update and key issues. N Engl J Med. 2011;364(9):852-860. doi:10.1056/NEJMsa1012065
14. Jones CW, Handler L, Crowell KE, Keil LG, Weaver MA, Platts-Mills TF. Non-publication of large randomized clinical trials: cross sectional analysis. BMJ. 2013;347:f6104. doi:10.1136/bmj.f6104
15. Ross JS, Tse T, Zarin DA, Xu H, Zhou L, Krumholz HM. Publication of NIH funded trials registered in ClinicalTrials.gov: cross sectional analysis. BMJ. 2012;344:d7292. doi:10.1136/bmj.d7292
16. Turner EH, Matthews AM, Linardatos E, Tell RA, Rosenthal R. Selective publication of antidepressant trials and its influence on apparent efficacy. N Engl J Med. 2008;358(3):252-260. doi:10.1056/NEJMsa065779
17. Huser V, Cimino JJ. Linking ClinicalTrials.gov and PubMed to track results of interventional human clinical trials. PLoS One. 2013;8(7):e68409. doi:10.1371/journal.pone.0068409
18. The Cochrane Library [database]. Cochrane; 2024. Available from https://www.cochranelibrary.com.
19. Cipriani A, Furukawa TA, Salanti G, et al. Comparative efficacy and acceptability of 21 antidepressant drugs for the acute treatment of adults with major depressive disorder: a systematic review and network meta-analysis. Lancet. 2018;391(10128):1357-1366. doi:10.1016/S0140-6736(17)32802-7
20. Derry S, Bell RF, Straube S, Wiffen PJ, Aldington D, Moore RA. Pregabalin for neuropathic pain in adults. Cochrane Database Syst Rev. 2019;1:CD007076. doi:10.1002/14651858.CD007076.pub3
21. Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435. doi:10.1007/BF02294384
22. Trinquart L, Dunn AG, Bourgeois FT. Registration of published randomized trials: a systematic review and meta-analysis. BMC Med. 2018;16(1):173. doi:10.1186/s12916-018-1168-6
23. Mathieu S, Boutron I, Moher D, Altman DG, Ravaud P. Comparison of registered and published primary outcomes in randomized controlled trials. JAMA. 2009;302(9):977-984. doi:10.1001/jama.2009.1242
24. Schmucker C, Schell LK, Portalupi S, et al. Extent of non-publication in cohorts of studies approved by research ethics committees or included in trial registries. PLoS One. 2014;9(12):e114023. doi:10.1371/journal.pone.0114023
25. Huang Y, Mao C, Yuan J, et al. Distribution of p-values in clinical trial registries and the implications for null hypothesis testing and publication bias. BMC Med Res Methodol. 2021;21:38. doi:10.1186/s12874-021-01222-1

---

## Figure Captions

**Fig 1. GWAM pipeline overview.** Schematic diagram showing the seven-stage pipeline from trial registry query through classification, weight computation, integrity ratio calculation, GWAM estimation, simulation interval construction, and optional HVP extension.

**Fig 2. Ghost protocol rates across drug-condition pairs.** Bar chart showing the percentage of completed registered trials classified as ghost protocols for six candidate drug-condition pairs. The dashed line indicates the 20% threshold used for selecting applications. Numbers above bars indicate total completed trials.

**Fig 3. Published random-effects vs. GWAM estimates for escitalopram/depression and pregabalin/neuropathic pain.** Forest-style plot comparing the published RE estimate (circle) with the GWAM simulation estimate (square) for each application. Horizontal lines represent 95% confidence intervals (RE) and simulation intervals (GWAM). The vertical dashed line at zero marks no effect.

**Fig 4. HVP sensitivity analysis for escitalopram/depression.** Dual heatmap showing Pr(effect > 0) (left) and 95% interval width (right) across a 4 x 4 grid of ghost_sigma (rows) and results_only_sigma (columns). Point estimates are invariant across the grid (see Methods). Pr(positive) ranges from 0.201 to 0.233, indicating that under GWAM's null assumption, the probability of a positive escitalopram effect is low regardless of dispersion specification. The corresponding grid for pregabalin/neuropathic pain is in S1 Table.

**Fig 5. Simulation: Bias and RMSE across five methods and three ground-truth scenarios.** Top row: bias; bottom row: RMSE. Columns correspond to mu_true = 0.0, 0.1, 0.2 (left to right). GWAM shows minimal bias at mu_true = 0 (where its null assumption holds) but increasing negative bias at mu_true > 0.

**Fig 6. Simulation: CI coverage across five methods and three ground-truth scenarios.** Paired bars show raw (lighter) and calibrated (darker) coverage. The dashed line marks nominal 95% coverage. GWAM raw coverage degrades at mu_true > 0 but calibrated coverage remains near-nominal. HVP-GWAM shows near-100% coverage (>=99.8%) across all scenarios, indicating uninformatively wide intervals.

**Fig 7. Pairwise70 benchmark: Lambda_hat and attenuation distributions.** Left: histogram of proxy lambda_hat values across 7,467 meta-analyses (median = 0.871). Right: histogram of |GWAM - RE| absolute attenuation (median = 0.055).

**Fig 8. Lambda_hat sensitivity: Base, strict, and loose configurations.** Left: attenuation rates (% crossing below |0.10| and |0.20|) for three lambda_hat specification scenarios across the full Pairwise70 corpus (13,019 analyses). Right: median absolute shift for each scenario. Crossed-below rates range from 5.8% (loose) to 14.5% (strict).

---

## Supporting Information

**S1 Table. HVP sensitivity grid for both applications.** Full 4 x 4 grid of HVP summaries (point estimate, SD, 95% interval, Pr(positive)) for escitalopram/depression and pregabalin/neuropathic pain, across ghost_sigma and results_only_sigma values of {0.05, 0.10, 0.20, 0.30}. Column definitions: hvp_mean = enrollment-weighted point estimate; hvp_sd = propagated standard deviation; interval_lo/interval_hi = 95% interval bounds; pr_positive = probability that effect > 0 under the HVP model. (CSV)

**S2 Table. Full simulation results for both applications.** Per-method, per-scenario simulation operating characteristics including: bias = mean(estimate - mu_true); rmse = sqrt(mean((estimate - mu_true)^2)); coverage_95 = fraction of replications where 95% CI contains mu_true; mean_ci_width = average CI width; ci_multiplier = empirical calibration multiplier; n_valid = number of replications with finite estimates. Results for escitalopram and pregabalin calibrations across mu_true = {0.0, 0.1, 0.2}. (CSV)

**S3 Table. Pairwise70 benchmark stratified summary.** Attenuation rates and distribution statistics stratified by lambda_hat proxy source (exact NCT, query-derived, transport, default) and study count band (k=2, k=3, k>=4). (CSV)

**S1 Code. Analysis code repository.** Complete Python scripts for the GWAM pipeline: registry fetching, trial classification, GWAM computation, simulation, HVP extension, Pairwise70 benchmark, figure generation, and supplementary table generation. Available at https://github.com/mahmoodsaffari/GWAM.

---

## Author Contributions (CRediT)

**Mahmood Ahmad:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data Curation, Writing -- Original Draft, Writing -- Review & Editing, Visualization. This single-author paper meets all four ICMJE authorship criteria.

## Competing Interests

The author declares no competing interests.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Data Availability

All ClinicalTrials.gov data used in this study are publicly available via the ClinicalTrials.gov API v2 (https://clinicaltrials.gov/api/v2). Registry queries were executed on 2026-02-14; pre-fetched CSV snapshots are included in the repository. The Pairwise70 corpus was obtained from published Cochrane systematic reviews via The Cochrane Library [18]. Analysis code and processed datasets are available at https://github.com/mahmoodsaffari/GWAM (S1 Code).

## Ethics Statement

This study used only publicly available, de-identified registry and literature data. No human participants were enrolled, and no ethics approval was required.
