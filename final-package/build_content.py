#!/usr/bin/env python3
"""GWAM final-package content + Table 1 + visual abstract generator.

All quantitative values are taken verbatim from the verified GWAM source pack
(F:\\Models\\GWAM): the raw registry CSV, the analysis JSON artifacts, and the
canonical manuscript paper/gwam_manuscript.md (whose numbers were independently
reproduced in-session: lambda=0.0959, mu_GWAM_null=0.024408, sim_mean=0.024619,
95% SI [-0.041357, 0.087373], to 6 decimals).
"""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- identity
GNOSIS = dict(
    title="Gnosis",
    subtitle="Global Health Evidence",
    tag="Diamond Open Access",
    accent="#182650",   # midnight blue
    accent2="#384C82",
    bg="#F2F4FA",
    ink="#0A1020",
)

# ---------------------------------------------------------------- article meta (final, no "preliminary")
META = dict(
    title="GWAM: Ghost-Weighted Aggregate Meta-Analysis for registry-based publication-bias sensitivity",
    author="Mahmood Ahmad",
    orcid="0009-0003-7781-4478",
    affiliation="Royal Free Hospital, London, United Kingdom",  # FLAG: confirm vs OJS author record
    email="mahmood.ahmad2@nhs.net",
    journal="Gnosis",                # per author framing + PDF identity; OJS install currently labels it Synthesis (flag at routing)
    journal_section="Methods Note",
    volume="2", issue="4", year="2026",
    published="2026-06-16",
    license="CC BY 4.0",
    doi="",        # suppressed site-wide / pending -- DO NOT fabricate
    issn="",       # pending
    code="https://github.com/mahmood726-cyber/GWAM",
)

# ---------------------------------------------------------------- final structured abstract (preliminary notice REMOVED)
ABSTRACT = [
    ("Background", "Publication bias threatens the validity of meta-analysis, yet routine detection tools "
        "(funnel-plot asymmetry, trim-and-fill) rest on assumptions that are hard to verify. Trial registries "
        "offer an external anchor for how much of the evidence base is actually visible."),
    ("Methods", "GWAM (Ghost-Weighted Aggregate Meta-Analysis) is a registry-anchored sensitivity method. It "
        "classifies ClinicalTrials.gov records as ghost (completed, with no publication and no posted results), "
        "results-posted, or published (PMID-linked); computes an enrollment-weighted integrity ratio λ (the "
        "published-participant fraction); and scales the published pooled effect by λ. A Monte-Carlo layer and "
        "a hierarchical variance-propagation extension carry the uncertainty about unobserved trials."),
    ("Results", "For escitalopram in depression (18 trials) λ was 0.096; the GWAM-adjusted log-odds-ratio was "
        "0.025 (95% SI −0.041 to 0.087; SI denotes a simulation interval, not a confidence interval) versus a "
        "published random-effects log-odds-ratio of 0.255 (OR 1.29). Across a benchmark of 7,467 meta-analyses from "
        "398 Cochrane reviews the median λ was 0.871, and 7.9% of meta-analyses with a random-effects estimate "
        "of magnitude at least 0.10 were attenuated below 0.10. For pregabalin (35 trials, log-risk-ratio) λ was "
        "0.014, flagged as likely incomplete PMID linkage. In simulation (2,000 replications) bias was minimal at "
        "μ=0 but grew negative for μ>0 (−0.162 at μ=0.2) and raw coverage fell to 7.1%; the "
        "hierarchical extension restored near-nominal coverage."),
    ("Conclusions", "GWAM is a transparent, reproducible registry-anchored sensitivity analysis, best interpreted "
        "as a conservative framework that complements existing publication-bias methods; its performance depends on "
        "the quality of ClinicalTrials.gov publication linkage."),
]

# ---------------------------------------------------------------- 3 references (final)
REFERENCES = [
    dict(authors="Egger M, Davey Smith G, Schneider M, Minder C.",
         title="Bias in meta-analysis detected by a simple, graphical test.",
         source="BMJ.", year="1997", vol="315", issue="7109", fpage="629", lpage="634",
         doi="10.1136/bmj.315.7109.629", pmid="9310563"),
    dict(authors="Duval S, Tweedie R.",
         title="Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis.",
         source="Biometrics.", year="2000", vol="56", issue="2", fpage="455", lpage="463",
         doi="10.1111/j.0006-341X.2000.00455.x", pmid="10877304"),
    dict(authors="Cipriani A, Furukawa TA, Salanti G, et al.",
         title="Comparative efficacy and acceptability of 21 antidepressant drugs for the acute treatment of adults with major depressive disorder: a systematic review and network meta-analysis.",
         source="Lancet.", year="2018", vol="391", issue="10128", fpage="1357", lpage="1366",
         doi="10.1016/S0140-6736(17)32802-7", pmid="29477251"),
]

# ---------------------------------------------------------------- Table 1 (10 data rows): lambda / effect / benchmark
# Columns: panel, stratum, n, lambda, published_effect, gwam_result
# All values verbatim from the verified manuscript (Tables 1-3 + Results text) and analysis JSONs.
TABLE1_COLS = ["Panel", "Stratum / drug–condition", "n", "λ",
               "Published effect (log scale)", "GWAM result"]
TABLE1_ROWS = [
    # Panel A -- worked applications
    ["Applications", "Escitalopram / depression", "18 trials", "0.096",
        "0.255 (OR 1.29)", "0.025 (95% SI −0.041, 0.087)"],
    ["Applications", "Pregabalin / neuropathic pain", "35 trials", "0.014",
        "0.916 (RR 2.50)", "0.013 (95% SI −0.029, 0.058)"],
    # Panel B -- Pairwise70 benchmark, overall
    ["Benchmark (overall)", "All Cochrane meta-analyses", "7,467", "0.871 (median)",
        "median |Δ| 0.055", "7.9% attenuated <0.10 (95% CI 5.1–11.0)"],
    # Panel C -- benchmark by lambda-hat source
    ["Benchmark by λ̂ source", "Exact NCT linkage", "1,438", "0.950", "—", "3.1% attenuated <0.10"],
    ["Benchmark by λ̂ source", "CT.gov query linkage", "5,569", "0.861", "—", "9.1% attenuated <0.10"],
    ["Benchmark by λ̂ source", "Transport proxy", "86", "0.950", "—", "10.5% attenuated <0.10"],
    ["Benchmark by λ̂ source", "Default (0.70)", "374", "0.700", "—", "8.3% attenuated <0.10"],
    # Panel D -- benchmark by study count
    ["Benchmark by study count", "k = 2", "2,895", "0.894", "—", "6.3% attenuated <0.10"],
    ["Benchmark by study count", "k = 3", "1,409", "0.880", "—", "7.9% attenuated <0.10"],
    ["Benchmark by study count", "k ≥ 4", "3,163", "0.845", "—", "9.4% attenuated <0.10"],
]
TABLE1_TITLE = ("Table 1. Integrity ratio (λ) and GWAM attenuation across the two worked applications and the "
                "Pairwise70 benchmark of 7,467 meta-analyses from 398 Cochrane reviews.")
TABLE1_FOOTNOTE = ("λ = enrollment-weighted integrity ratio (published-participant fraction); for applications "
                   "λ is the directly measured λ_pmid, for the benchmark it is the proxy λ̂. "
                   "Effects are on the log-OR (escitalopram) or log-RR (pregabalin) scale; positive values indicate "
                   "treatment superiority. SI = simulation interval (not a confidence interval). "
                   "“% attenuated <0.10” = share of meta-analyses with |random-effects estimate| ≥ 0.10 "
                   "whose GWAM-adjusted |estimate| fell below 0.10. Δ = GWAM − random-effects shift. "
                   "Pregabalin's λ = 0.014 reflects incomplete ClinicalTrials.gov PMID linkage, not true "
                   "non-publication.")

FIG1_CAPTION = ("Figure 1. GWAM pipeline overview. The seven-stage workflow runs from the ClinicalTrials.gov registry "
                "query through trial classification (published / results-posted / ghost), enrollment-weight "
                "computation, integrity-ratio (λ) calculation, GWAM point estimation, Monte-Carlo simulation-"
                "interval construction, and the optional hierarchical variance-propagation extension.")

def write_outputs():
    # Table 1 CSV
    with open(os.path.join(HERE, "table1.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(TABLE1_COLS)
        w.writerows(TABLE1_ROWS)
    # content bundle JSON
    bundle = dict(identity=GNOSIS, meta=META, abstract=ABSTRACT, references=REFERENCES,
                  table1=dict(title=TABLE1_TITLE, cols=TABLE1_COLS, rows=TABLE1_ROWS, footnote=TABLE1_FOOTNOTE),
                  fig1_caption=FIG1_CAPTION)
    with open(os.path.join(HERE, "content_bundle.json"), "w", encoding="utf-8") as f:
        json.dump(bundle, f, indent=2, ensure_ascii=False)
    print("wrote table1.csv, content_bundle.json  (rows:", len(TABLE1_ROWS), ")")

if __name__ == "__main__":
    write_outputs()
