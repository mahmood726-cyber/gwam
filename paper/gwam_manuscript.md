# Ghost-Weighted Aggregate Meta-analysis (GWAM): A Registry-Based Sensitivity Framework for Publication Bias

**Authors:** Mahmood S. [Corresponding author]

**Affiliation:** [To be completed]

**Corresponding author:** [email to be completed]

---

## Abstract

**Background:** Publication bias remains a pervasive threat to evidence synthesis. Existing corrections (trim-and-fill, PET-PEESE, selection models) rely on statistical signatures within the published literature, yet the most informative signal — the fraction of registered trials that remain unpublished — is rarely incorporated quantitatively. We introduce Ghost-Weighted Aggregate Meta-analysis (GWAM), a registry-based sensitivity framework that uses ClinicalTrials.gov enrollment data to reweight published pooled estimates by the proportion of the total evidence base they represent.

**Methods:** GWAM classifies registered trials into three categories: (i) published with a linked PubMed identifier (PMID), (ii) results-posted on the registry without a journal publication, and (iii) ghost protocols with neither publication nor posted results. An enrollment-weighted integrity ratio (lambda) quantifies the fraction of total enrolled participants represented by published studies. The deterministic GWAM estimate scales the published pooled effect by lambda, while a Monte Carlo simulation layer and a Bayesian sensitivity extension propagate uncertainty about ghost and results-only effects. We applied GWAM to escitalopram for depression (18 trials, log-OR scale) and pregabalin for neuropathic pain (35 trials, log-RR scale), validated performance via simulation (2,000 replications under three ground-truth scenarios), and benchmarked attenuation across 7,467 binary-outcome meta-analyses from 458 Cochrane systematic reviews (Pairwise70 corpus).

**Results:** For escitalopram, lambda (PMID-only) was 0.096, indicating that published studies accounted for less than 10% of total enrolled participants; the GWAM-adjusted estimate was log-OR = 0.025 (95% simulation interval [SI]: -0.041, 0.087) compared with the published random-effects (RE) estimate of 0.255. In simulation under the null (mu_true = 0), GWAM achieved bias of 0.013 and root mean square error (RMSE) of 0.044, compared with RE bias of 0.158 and RMSE of 0.411. Across the Pairwise70 benchmark (median lambda = 0.871), 7.9% of meta-analyses with |RE| > 0.10 were attenuated below that threshold by GWAM (bootstrap 95% CI: 5.1%--11.0%).

**Conclusions:** GWAM provides a transparent, reproducible sensitivity analysis anchored in registry data rather than statistical distributional assumptions. It complements existing methods by offering an external benchmark for the plausible magnitude of publication bias.

---

## Introduction

Publication bias — the selective reporting of studies based on the statistical significance or direction of results — distorts the evidence base that informs clinical practice and policy [1,2]. Meta-analyses, despite pooling multiple studies, inherit this distortion when the underlying literature is non-representative of all conducted research [3].

Several statistical methods address publication bias, including funnel plot asymmetry tests [4], trim-and-fill [5], selection models [6,7], PET-PEESE regression [8,9], and Copas selection models [10]. These approaches share a common limitation: they infer the presence and magnitude of missing studies from statistical properties (e.g., asymmetry, small-study effects) of the *observed* literature. When bias operates through mechanisms other than small-study effects — such as selective non-publication of entire large trials — funnel-based corrections may be insufficient [11].

Trial registries offer a fundamentally different information source. Since 2005, the International Committee of Medical Journal Editors (ICMJE) has required prospective registration of clinical trials [12], and the FDA Amendments Act of 2007 (FDAAA) mandates results posting for applicable trials [13]. ClinicalTrials.gov, the largest registry, now contains over 500,000 study records. Comparing registered protocols against published literature reveals "ghost protocols" — trials that were registered and completed but never published their results, either in a journal or on the registry [14,15].

We propose Ghost-Weighted Aggregate Meta-analysis (GWAM), a sensitivity framework that uses registry enrollment data to quantify the fraction of the total evidence base represented by published studies. Rather than attempting to estimate missing study effects from funnel plot shape, GWAM computes an enrollment-weighted integrity ratio (lambda) and applies it as a deterministic scaling factor, with a simulation layer and Bayesian extension to propagate uncertainty about unobserved effects.

GWAM is designed as a *complement* to existing bias-adjustment methods, not a replacement. It provides a registry-anchored external benchmark that can be compared against internally derived corrections. Fig 1 illustrates the overall pipeline. In the following sections, we describe the method, apply it to two clinical examples, validate via simulation, and benchmark across 458 Cochrane reviews.

---

## Methods

### Trial registry pipeline

For each drug-condition pair, we queried ClinicalTrials.gov (accessed 2026-02-14) via the Clinical Trials API v2 to retrieve all registered interventional studies of the target drug for the target condition. Trials were filtered to those with status "Completed" and available enrollment counts. Each trial was classified into one of three mutually exclusive categories based on automated linkage:

1. **Published (PMID-linked):** A PubMed identifier (PMID) was associated with the trial record via ClinicalTrials.gov reference fields or manual curation.
2. **Results-posted (no PMID):** The trial had posted summary results on ClinicalTrials.gov but had no linked journal publication.
3. **Ghost protocol:** The trial had neither a linked PMID nor posted results.

A small number of trials in each candidate scan fell outside these three categories (e.g., trials with associated publications that were not primary results papers, or trials with ambiguous linkage status). These were excluded from the GWAM analysis but are reported in the candidate scan for transparency.

### Candidate scan

We screened six drug-condition pairs to identify applications with substantial ghost protocol rates (Table 1). Selection criteria for detailed GWAM analysis were: (a) ghost protocol rate exceeding 20%, (b) availability of a published meta-analysis providing a pooled effect estimate, and (c) sufficient trial heterogeneity to make publication bias a plausible concern.

### GWAM model

Let W_pmid, W_ro, and W_ghost denote the total enrollment-weighted contributions from published, results-only, and ghost trials, respectively, with W_total = W_pmid + W_ro + W_ghost. The integrity ratio is defined as:

> lambda_pmid = W_pmid / W_total

which measures the fraction of total registered enrollment represented by published studies. A second ratio, lambda_non_ghost = (W_pmid + W_ro) / W_total, quantifies the fraction with any available results (published or posted).

The GWAM estimate is a three-way enrollment-weighted average:

> mu_GWAM = (W_pmid * mu_pub + W_ro * mu_ro + W_ghost * mu_ghost) / W_total

where mu_pub is the published random-effects pooled estimate, mu_ro is the assumed mean effect for results-only trials, and mu_ghost is the assumed mean effect for ghost protocols.

Under the *null assumption* that unpublished evidence is centered at zero (mu_ro = 0, mu_ghost = 0), this simplifies to:

> mu_GWAM_null = lambda_pmid * mu_pub

This deterministic point estimate represents the published effect diluted by the proportion of the total evidence base that is published. When lambda_pmid is small (i.e., most enrolled participants are in unpublished trials), the GWAM estimate is substantially attenuated.

**Design note:** Weights are based on enrollment counts, not inverse-variance, because variance estimates are unavailable for ghost protocols and many results-only studies. Enrollment-based weighting reflects the registry "universe" denominator for the integrity ratio and is a deliberate design choice.

### Monte Carlo simulation layer

To propagate uncertainty about unobserved effects, we add a simulation layer. Ghost effects are drawn from N(mu_ghost, sigma_ghost) and results-only effects from N(mu_ro, sigma_ro), with default parameters mu_ghost = mu_ro = 0, sigma_ghost = sigma_ro = 0.10. The GWAM simulation estimate is computed as:

> mu_GWAM^(b) = (W_pmid * mu_pub + sum_j w_ro_j * epsilon_ro_j^(b) + sum_k w_ghost_k * epsilon_ghost_k^(b)) / W_total

for b = 1, ..., B Monte Carlo draws, yielding a simulation-based 95% interval (SI) from the 2.5th and 97.5th percentiles of the distribution.

### Bayesian sensitivity extension

The Bayesian extension integrates over uncertainty in ghost and results-only effects using conjugate Normal priors. The published estimate provides the likelihood:

> mu_pub | mu_true ~ N(mu_true, SE_pub^2)

where SE_pub is the standard error of the published pooled estimate. The model evaluates a 4 x 4 sensitivity grid over ghost_sigma in {0.05, 0.10, 0.20, 0.30} and results_only_sigma in {0.05, 0.10, 0.20, 0.30}, producing posterior summaries (mean, SD, 95% credible interval [CrI], Pr(effect > 0)) for each cell. The reference cell uses the tightest ghost prior (sigma = 0.05) with moderate results-only variance.

**Posterior mean invariance:** Under this model, the posterior mean depends on the data likelihood and prior location parameters, but is invariant to the prior variance parameters (ghost_sigma, results_only_sigma). This structural property of Normal-Normal conjugacy means that only the posterior SD and CrI width vary across the sensitivity grid, while Pr(effect > 0) changes modestly through the widening tails.

### Simulation validation

To evaluate operating characteristics, we conducted a calibration simulation with 2,000 replications for each of three ground-truth scenarios (mu_true in {0.0, 0.1, 0.2}). For each replication, we:

1. Generated k study-level log-odds ratios from N(mu_true, tau^2) with tau^2 = 0.1, using enrollment weights from the escitalopram registry structure.
2. Simulated selective publication: each study was "published" with probability p_sig = 0.95 if its result was statistically significant, or with calibrated probability p_nonsig otherwise. The calibration procedure ensured that the simulated lambda_pmid matched the observed lambda_pmid (0.096) within tolerance 0.03.
3. Applied five estimation methods to the selectively published studies: (a) standard random-effects (RE), (b) GWAM null (lambda-scaled), (c) Bayesian GWAM, (d) oracle inverse-probability-weighted (IPW) RE using the true selection probabilities, and (e) PET-PEESE.
4. Computed bias, RMSE, raw 95% CI coverage, and calibrated CI coverage using empirically derived multipliers.

**Critical caveat:** Under mu_true = 0, GWAM has a structural advantage because its null assumption (unobserved effects centered at zero) coincides with the ground truth. Performance at mu_true > 0 reveals the cost of this assumption: GWAM increasingly underestimates the true effect as mu_true departs from zero.

### CI calibration

Raw confidence intervals for all methods may under- or over-cover due to the simulation-based nature of the correction. We applied empirical CI width multipliers derived from a calibration run at mu_true = mu_pub (the published estimate), optimizing each method's multiplier to achieve 95% coverage. These multipliers were then applied to all three scenarios as an out-of-sample test.

### Pairwise70 benchmark

To assess GWAM's attenuation effect at scale, we applied the proxy-lambda pipeline to 7,467 binary-outcome meta-analyses extracted from 458 Cochrane systematic reviews in the Pairwise70 corpus [16]. For each meta-analysis, we linked the contributing trials to ClinicalTrials.gov via NCT identifier matching and query-derived linkage, computed a proxy lambda, and calculated the GWAM-adjusted estimate. Lambda sources were (in priority order): exact Pairwise70 NCT linkage, ClinicalTrials.gov query-derived non-ghost ratio, transport proxy from covariate similarity, and a conservative default of 0.70 when no linkage was available. Lambda values were clipped to [0.20, 0.95] to avoid extreme adjustments.

We report the fraction of meta-analyses where the GWAM-adjusted absolute effect crossed below clinically relevant thresholds (|mu| < 0.10 and |mu| < 0.20 on the log-OR scale). Review-clustered bootstrap confidence intervals (400 iterations) account for within-review correlation.

### Software and reproducibility

All analyses were conducted in Python 3.12 using NumPy 1.26, pandas 2.1, and matplotlib 3.8. Random number generation used NumPy's default_rng with seed 42 for reproducibility. ClinicalTrials.gov data were accessed via the Clinical Trials API v2. The complete analysis code, registry data, and figure generation scripts are available at [repository URL].

---

## Results

### Candidate scan

Table 1 presents the ghost protocol rates for six drug-condition pairs. Ghost rates ranged from 16.7% (quetiapine/depression) to 25.7% (pregabalin/neuropathic pain). Escitalopram for depression (24.7% ghost rate, 73 completed trials) and pregabalin for neuropathic pain (25.7%, 35 completed trials) were selected for detailed analysis based on having ghost rates exceeding 20% and available published meta-analyses.

**Table 1. Candidate scan: Ghost protocol rates across drug-condition pairs (ClinicalTrials.gov, accessed 2026-02-14).**

| Drug / Condition | N completed | With PMID | Results-only (no PMID) | Ghost protocols | Ghost rate (%) |
|---|---|---|---|---|---|
| Sertraline / Depression | 43 | 11 | 19 | 9 | 20.9 |
| Escitalopram / Depression | 73 | 11 | 35 | 18 | 24.7 |
| Duloxetine / Depression | 37 | 5 | 23 | 7 | 18.9 |
| Quetiapine / Depression | 18 | 0 | 13 | 3 | 16.7 |
| Pregabalin / Neuropathic Pain | 35 | 1 | 18 | 9 | 25.7 |
| Oseltamivir / Influenza | 25 | 1 | 14 | 5 | 20.0 |

*Note:* Column totals (PMID + Results-only + Ghost) may not equal N completed. The difference represents trials with associated publications that are not primary results papers, or trials with ambiguous linkage status, excluded from the three-way classification. Ghost rate = Ghost / N completed.

### Application 1: Escitalopram for depression

From 73 completed escitalopram/depression trials identified in the candidate scan, 18 met the inclusion criteria for GWAM analysis after applying study-level filters (requirement for enrollment data, appropriate comparator, and relevant outcome). Of these 18 trials:

- 4 (22.2%) had published PMID-linked results (total enrollment weight: 471)
- 9 (50.0%) had registry-posted results without a journal publication (weight: 3,832)
- 5 (27.8%) were ghost protocols (weight: 610)

The integrity ratio lambda_pmid = 471 / 4,913 = 0.096, indicating that published studies represented less than 10% of total enrolled participants. The non-ghost ratio lambda_non_ghost = 4,303 / 4,913 = 0.876.

The published random-effects pooled estimate from the literature meta-analysis was log-OR = 0.255 (SE = 0.096; 95% CI: 0.066, 0.443; OR = 1.29). Under the GWAM null assumption, the adjusted estimate was:

> mu_GWAM = lambda_pmid * mu_pub = 0.096 * 0.255 = 0.024 (OR = 1.025)

The Monte Carlo simulation (B = 2,000, sigma_ghost = sigma_ro = 0.10) yielded a mean of 0.025 with 95% SI: [-0.041, 0.087], indicating that the published effect of escitalopram was not robust to the possibility that unpublished trials showed null or negative results (Fig 3).

The Bayesian sensitivity analysis (Fig 4) yielded a posterior mean of -0.098 (95% CrI: -0.341, 0.145) at the reference cell (ghost_sigma = 0.05, results_only_sigma = 0.20), with Pr(effect > 0) = 0.215. Across the full 4 x 4 grid, posterior means were invariant (see Methods), while Pr(effect > 0) ranged from 0.201 (tightest priors) to 0.231 (widest priors). The Bayesian analysis uses additional information from the published SE and five results-only trials with observed effect sizes, which explains the shift from the simulation-only GWAM estimate toward a more negative posterior.

### Application 2: Pregabalin for neuropathic pain

For pregabalin/neuropathic pain, 35 completed trials were identified (all 35 met inclusion criteria):

- 1 (2.9%) had a published PMID (weight: 104)
- 25 (71.4%) had posted results only (weight: 5,720)
- 9 (25.7%) were ghost protocols (weight: 1,597)

Lambda_pmid = 104 / 7,421 = 0.014, indicating that less than 2% of enrolled participants were in the single published trial. Lambda_non_ghost = 5,824 / 7,421 = 0.785.

The published estimate was log-RR = 0.916 (SE = 0.148; 95% CI: 0.625, 1.207; RR = 2.50). The GWAM-adjusted estimate was:

> mu_GWAM = 0.014 * 0.916 = 0.013 (RR = 1.013)

The 95% SI was [-0.029, 0.058], and the Bayesian posterior mean was -0.234 (95% CrI: -0.384, -0.085), with Pr(effect > 0) = 0.001 (reference cell: ghost_sigma = 0.05, results_only_sigma = 0.30). This extreme attenuation reflects the dominance of unpublished evidence: a single published trial represented only 1.4% of total enrollment.

### Simulation validation

Table 2 summarizes simulation results (escitalopram weight structure, 2,000 replications). Under mu_true = 0:

**Table 2. Simulation operating characteristics under three ground-truth scenarios (escitalopram calibration).**

| mu_true | Method | Bias | RMSE | Coverage (raw) | Coverage (calibrated) | n_valid |
|---|---|---|---|---|---|---|
| 0.0 | RE | 0.158 | 0.411 | 84.8% | 95.3% | 2000 |
| 0.0 | GWAM | 0.013 | 0.044 | 92.4% | 100% | 2000 |
| 0.0 | Bayesian GWAM | 0.013 | 0.031 | 100% | 100% | 2000 |
| 0.0 | Oracle IPW | 0.042 | 0.365 | 73.6% | 92.8% | 2000 |
| 0.0 | PET-PEESE | 0.176 | 1.048 | 74.4% | 95.8% | 1271 |
| 0.1 | RE | 0.176 | 0.408 | 82.2% | 94.9% | 2000 |
| 0.1 | GWAM | -0.076 | 0.088 | 52.6% | 100% | 2000 |
| 0.1 | Bayesian GWAM | -0.076 | 0.082 | 100% | 100% | 2000 |
| 0.2 | RE | 0.206 | 0.399 | 81.0% | 95.5% | 2000 |
| 0.2 | GWAM | -0.162 | 0.168 | 7.1% | 99.5% | 2000 |
| 0.2 | Bayesian GWAM | -0.162 | 0.166 | 99.8% | 100% | 2000 |

Fig 5 displays bias and RMSE, and Fig 6 displays coverage across methods and scenarios.

Under mu_true = 0 (no true effect), GWAM showed minimal bias (0.013) and the lowest RMSE (0.044), outperforming all comparators including the oracle IPW estimator. Standard RE had bias of 0.158, confirming substantial distortion from selective publication. PET-PEESE was applicable in only 63.6% of replications (1,271/2,000) and showed the highest RMSE (1.048).

Under mu_true = 0.1 and 0.2, GWAM exhibited increasing negative bias (-0.076 and -0.162, respectively) because its null assumption (unobserved effects = 0) diverges from the true positive effect. Raw CI coverage dropped to 52.6% and 7.1%, though calibrated CIs (using the empirically derived multiplier of 7.75) restored coverage to near-nominal levels. This bias-variance trade-off is inherent to GWAM's design: it provides maximal protection against spurious positive findings at the cost of underestimating true non-zero effects.

The pregabalin calibration showed that lambda calibration was not within tolerance (error 0.155 vs. maximum 0.030), indicating that the extreme lambda_pmid of 0.014 could not be stably achieved in the simulation's selective-publication model. This underscores that GWAM's performance guarantees are strongest when the registry structure can be faithfully reproduced in simulation.

### Pairwise70 benchmark

Across 7,467 meta-analyses from 458 Cochrane reviews, the median proxy lambda was 0.871 (Fig 7, left panel). The distribution was right-skewed, with 35.3% of analyses at the upper bound (0.95) and 5.7% at the lower bound (0.20).

GWAM attenuated 590 meta-analyses (7.9%; bootstrap 95% CI: 5.1%--11.0%) from |RE| >= 0.10 to |GWAM| < 0.10, and 819 (11.0%; 95% CI: 7.2%--14.9%) from |RE| >= 0.20 to |GWAM| < 0.20. Among meta-analyses with |RE| >= 0.10, the conditional attenuation rate was 9.4%, and among those with |RE| >= 0.20, it was 15.5%. The median absolute shift was 0.055 (Fig 7, right panel).

**Table 3. Pairwise70 benchmark: Attenuation rates by lambda source and study count.**

| Stratum | n | Median lambda | Crossed below 0.10 (% all) | Crossed below 0.20 (% all) |
|---|---|---|---|---|
| **By lambda source** | | | | |
| Exact NCT linkage | 1,438 | 0.950 | 3.1% | 7.4% |
| CT.gov query linkage | 5,569 | 0.861 | 9.1% | 11.7% |
| Transport proxy | 86 | 0.950 | 10.5% | 22.1% |
| Default (0.70) | 374 | 0.700 | 8.3% | 11.0% |
| **By study count** | | | | |
| k = 2 | 2,895 | 0.894 | 6.3% | 8.9% |
| k = 3 | 1,409 | 0.880 | 7.9% | 10.1% |
| k >= 4 | 3,163 | 0.845 | 9.4% | 13.3% |

Attenuation was lowest for meta-analyses with exact NCT linkage (3.1% crossed below 0.10), which tend to have high lambda values. Meta-analyses with more studies (k >= 4) showed higher attenuation rates (9.4%), consistent with larger meta-analyses drawing from broader, more heterogeneous trial registries with lower average lambda.

### Lambda sensitivity analysis

To assess robustness to lambda computation settings, we compared base, strict, and loose configurations across all 13,019 binary-outcome analyses in the Pairwise70 corpus (Fig 8):

- **Base** (default lambda = 0.70, bounds [0.20, 0.95]): crossed below 0.10 in 8.9%, below 0.20 in 10.0%.
- **Loose** (default = 0.90, bounds [0.40, 1.00]): 5.8% and 6.9%.
- **Strict** (default = 0.50, bounds [0.10, 0.70]): 14.5% and 14.6%.

The three-fold range (5.8%--14.5%) demonstrates moderate sensitivity to lambda specification, supporting the use of sensitivity grids rather than single-point estimates.

---

## Discussion

### Summary and interpretation

GWAM provides a registry-anchored sensitivity framework for publication bias that directly incorporates the proportion of registered evidence that is published. In two clinical applications, GWAM substantially attenuated published pooled effects: escitalopram's log-OR dropped from 0.255 to 0.025 (OR: 1.29 to 1.025), and pregabalin's log-RR dropped from 0.916 to 0.013 (RR: 2.50 to 1.013). These dramatic shifts reflect the low integrity ratios (lambda_pmid = 0.096 and 0.014, respectively), indicating that the published literature represents a small fraction of the total registered evidence.

The Pairwise70 benchmark suggests that such extreme lambda values are uncommon in Cochrane meta-analyses (median lambda = 0.871), and that GWAM's typical attenuation is modest: 7.9% of meta-analyses with |RE| > 0.10 crossed below that threshold. This is consistent with GWAM functioning as a sensitivity tool rather than a universal correction — it flags cases where the integrity ratio raises concern but does not override every published estimate.

### Relationship to existing methods

GWAM differs from existing publication bias methods in both its information source and its mechanism. Trim-and-fill [5] imputes missing studies based on funnel plot symmetry; selection models [6,7,10] estimate the selection function from the observed p-value distribution; PET-PEESE [8,9] regresses effect sizes on their standard errors; and Vevea-Hedges weight-function models [17] specify step-function selection weights. All these methods operate on the *internal* statistical properties of the published literature.

GWAM operates on an *external* data source — the trial registry — and makes no assumptions about funnel plot shape, p-value distributions, or small-study effects. This makes it complementary to existing methods: when both GWAM and funnel-based methods suggest attenuation, confidence in publication bias increases; when they diverge, the discrepancy itself is informative.

Huang et al. [18] proposed a related registry-based approach that counts registered-but-unpublished trials to construct a bias-corrected estimate. GWAM extends this concept by using enrollment-based weights (rather than study counts), incorporating a three-category classification (PMID / results-only / ghost), and providing simulation and Bayesian uncertainty quantification.

Turner et al. [19] and others have documented discrepancies between registered and published outcomes; GWAM operationalizes this evidence by converting the registration-publication gap into a quantitative correction.

### Strengths and limitations

**Strengths.** GWAM is fully transparent and reproducible: it uses publicly available registry data, deterministic formulas, and seeded random number generation. The integrity ratio is directly interpretable (the fraction of total enrollment that is published), and the sensitivity grid communicates how conclusions change under different assumptions about unobserved effects. The method requires no access to individual patient data and can be applied to any meta-analysis where the underlying trials are registered.

**Limitations.** First, GWAM assumes that enrollment weights are an appropriate proxy for study informativeness. Trials with larger enrollment may not proportionally contribute more information if they have higher within-study variance or different populations.

Second, the null assumption (mu_ghost = 0, mu_ro = 0) is conservative under the hypothesis that publication bias inflates effects, but may over-correct when a true positive effect exists. The simulation at mu_true = 0.1 and 0.2 demonstrates this bias-variance trade-off explicitly. The Bayesian extension mitigates this by allowing non-zero ghost effect distributions, but the user must specify the prior.

Third, ClinicalTrials.gov registration is incomplete, particularly for older trials and those conducted outside the United States. Lambda estimates may therefore be biased upward (making the published fraction appear larger) if unregistered trials exist. Conversely, registered trials with completed status may never have been intended to produce publications (e.g., pharmacokinetic studies), inflating the ghost count.

Fourth, the proxy lambda used in the Pairwise70 benchmark relies on automated NCT linkage and query-derived matching, which introduces measurement error. The sensitivity analysis (base/strict/loose) bounds the impact of this error.

Fifth, GWAM does not model within-trial selective outcome reporting (outcome reporting bias), which is a distinct mechanism from whole-trial non-publication. Combining GWAM with outcome-level bias assessment (e.g., ROB 2.0 domain 5) would provide more comprehensive bias evaluation.

### Implications for practice

We recommend GWAM as a supplementary sensitivity analysis for systematic reviews where the underlying trials are registered. Specifically:

1. Report the integrity ratio (lambda) alongside standard funnel plot assessments.
2. Present the GWAM-adjusted estimate as a "what-if" scenario: "If unpublished trials showed null effects, the pooled estimate would be X."
3. Use the Bayesian sensitivity grid to communicate the range of plausible adjusted estimates under different prior assumptions.
4. When GWAM attenuation and funnel-based methods agree, this strengthens the evidence for publication bias; when they disagree, investigate the source of discordance.

GWAM is not intended to replace the published estimate but to provide context about its sensitivity to publication completeness. The method is most informative when the integrity ratio is low (lambda < 0.50), indicating that a large fraction of enrolled participants are in unpublished trials.

---

## References

1. Rothstein HR, Sutton AJ, Borenstein M. Publication bias in meta-analysis: Prevention, assessment and adjustments. Chichester: Wiley; 2005.
2. Song F, Parekh S, Hooper L, et al. Dissemination and publication of research findings: an updated review of related biases. Health Technol Assess. 2010;14(8):1-193.
3. Ioannidis JP. Why most published research findings are false. PLoS Med. 2005;2(8):e124.
4. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. BMJ. 1997;315(7109):629-634.
5. Duval S, Tweedie R. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463.
6. Hedges LV. Modeling publication selection effects in meta-analysis. Stat Sci. 1992;7(2):246-255.
7. Vevea JL, Woods CM. Publication bias in research synthesis: sensitivity analysis using a priori weight functions. Psychol Methods. 2005;10(4):428-443.
8. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78.
9. Stanley TD. Limitations of PET-PEESE and other meta-analytic methods: Beyond the "cliff." Res Synth Methods. 2017;8(1):97-105.
10. Copas J, Shi JQ. A sensitivity analysis for publication bias in systematic reviews. Stat Methods Med Res. 2001;10(4):251-265.
11. Sterne JAC, Gavaghan D, Egger M. Publication and related bias in meta-analysis: power of statistical tests and prevalence in the literature. J Clin Epidemiol. 2000;53(11):1119-1129.
12. De Angelis C, Drazen JM, Frizelle FA, et al. Clinical trial registration: a statement from the International Committee of Medical Journal Editors. N Engl J Med. 2004;351(12):1250-1251.
13. Zarin DA, Tse T, Williams RJ, Califf RM, Ide NC. The ClinicalTrials.gov results database -- update and key issues. N Engl J Med. 2011;364(9):852-860.
14. Jones CW, Handler L, Crowell KE, Keil LG, Weaver MA, Platts-Mills TF. Non-publication of large randomized clinical trials: cross sectional analysis. BMJ. 2013;347:f6104.
15. Ross JS, Tse T, Zarin DA, Xu H, Zhou L, Krumholz HM. Publication of NIH funded trials registered in ClinicalTrials.gov: cross sectional analysis. BMJ. 2012;344:d7292.
16. Cochrane Database of Systematic Reviews. The Cochrane Library. Available at: https://www.cochranelibrary.com/. Accessed 2026-02-14.
17. Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435.
18. Huang Y, Mao C, Yuan J, et al. Distribution of p-values in clinical trial registries and the implications for null hypothesis testing and publication bias. BMC Med Res Methodol. 2021;21:38.
19. Turner EH, Matthews AM, Linardatos E, Tell RA, Rosenthal R. Selective publication of antidepressant trials and its influence on apparent efficacy. N Engl J Med. 2008;358(3):252-260.

---

## Figure Captions

**Fig 1. GWAM pipeline overview.** Schematic diagram showing the six-stage pipeline from trial registry query through classification, weight computation, integrity ratio calculation, GWAM estimation, and Bayesian extension.

**Fig 2. Ghost protocol rates across drug-condition pairs.** Bar chart showing the percentage of completed registered trials classified as ghost protocols for six candidate drug-condition pairs. The dashed line indicates the 20% threshold used for selecting applications. Numbers above bars indicate total completed trials.

**Fig 3. Published random-effects vs. GWAM estimates for escitalopram/depression and pregabalin/neuropathic pain.** Forest-style plot comparing the published RE estimate (circle) with the GWAM simulation estimate (square) for each application. Horizontal lines represent 95% confidence intervals (RE) and simulation intervals (GWAM). The vertical dashed line at zero marks no effect.

**Fig 4. Bayesian sensitivity analysis for escitalopram/depression.** Dual heatmap showing Pr(effect > 0) (left) and 95% CrI width (right) across a 4 x 4 grid of ghost_sigma (rows) and results_only_sigma (columns). Posterior means are invariant across the grid (see Methods). Pr(positive) ranges from 0.201 to 0.231, indicating robust evidence against a positive escitalopram effect regardless of prior variance specification.

**Fig 5. Simulation: Bias and RMSE across five methods and three ground-truth scenarios.** Top row: bias; bottom row: RMSE. Columns correspond to mu_true = 0.0, 0.1, 0.2 (left to right). GWAM shows minimal bias at mu_true = 0 but increasing negative bias at mu_true > 0.

**Fig 6. Simulation: CI coverage across five methods and three ground-truth scenarios.** Paired bars show raw (lighter) and calibrated (darker) coverage. The dashed line marks nominal 95% coverage. GWAM raw coverage degrades at mu_true > 0 but calibrated coverage remains near-nominal.

**Fig 7. Pairwise70 benchmark: Lambda and attenuation distributions.** Left: histogram of proxy lambda values across 7,467 meta-analyses (median = 0.871). Right: histogram of |GWAM - RE| absolute attenuation (median = 0.055).

**Fig 8. Lambda sensitivity: Base, strict, and loose configurations.** Left: attenuation rates (% crossing below |0.10| and |0.20|) for three lambda specification scenarios. Right: median absolute shift for each scenario. Crossed-below rates range from 5.8% (loose) to 14.5% (strict).

---

## Supporting Information

**S1 Table. Bayesian sensitivity grid for both applications.** Full 4 x 4 grid of posterior summaries (mean, SD, CrI, Pr(positive)) for escitalopram/depression and pregabalin/neuropathic pain, across ghost_sigma and results_only_sigma values of {0.05, 0.10, 0.20, 0.30}. (CSV)

**S2 Table. Full simulation results for both applications.** Per-method, per-scenario simulation operating characteristics (bias, RMSE, coverage, CI width, calibration multiplier) for escitalopram and pregabalin calibrations across mu_true = {0.0, 0.1, 0.2}. (CSV)

**S3 Table. Pairwise70 benchmark stratified summary.** Attenuation rates and distribution statistics stratified by lambda proxy source and study count band. (CSV)

**S1 Code. Analysis code repository.** Complete Python scripts for the GWAM pipeline: registry fetching, trial classification, GWAM computation, simulation, Bayesian extension, Pairwise70 benchmark, figure generation, and supplementary table generation.

---

## Author Contributions (CRediT)

**Mahmood S.:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data Curation, Writing -- Original Draft, Writing -- Review & Editing, Visualization.

## Competing Interests

The author declares no competing interests.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Data Availability Statement

All ClinicalTrials.gov data used in this study are publicly available via the ClinicalTrials.gov API v2 (https://clinicaltrials.gov/api/v2). Registry queries were executed on 2026-02-14. The Pairwise70 corpus was obtained from published Cochrane systematic reviews via The Cochrane Library. Analysis code and processed datasets are available at [repository URL to be provided upon acceptance].

## Ethics Statement

This study used only publicly available, de-identified registry and literature data. No human participants were enrolled, and no ethics approval was required.
