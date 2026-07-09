# GWAM — Verification report (finalization to FINAL form)

**Date:** 2026-06-20 · **Verifier:** Claude (Opus 4.8), in-session
**Source pack:** `F:\Models\GWAM` (located by search of F:\ and C:\Users\mahmo for `gwam` / `ghost-weighted`)
**Live preliminary article:** Synthēsis OJS, slug `gwam-ghost-weighted-meta-analysis`, Vol. 2 No. 4 (June 2026),
Methods Note, published 2026-06-16 (captured HTML: `C:\Users\mahmo\AppData\Local\Temp\gw.html`).

---

## 1. Headline — independently reproduced (not merely read from the stored artifact)

The escitalopram/depression headline was re-derived **from the raw registry CSV**
(`data/raw/escitalopram__depression_with_results.csv`, 18 trials) using an independent classification +
weighting + Monte-Carlo reimplementation, then compared to the stored pipeline artifact
(`data/analysis/escitalopram_depression_gwam.json`). **Agreement to 6 decimals.**

| Quantity | Independent recompute | Stored artifact | Published text | Match |
|---|---|---|---|---|
| Trial classification (pmid / results-only / ghost) | 4 / 9 / 5 | 4 / 9 / 5 | 4 / 9 / 5 | ✅ |
| W_published / W_total | 471 / 4913 | 471 / 4913 | 471 / 4913 | ✅ |
| λ_pmid | 0.095868 → **0.096** | 0.095868 | 0.096 | ✅ |
| μ_GWAM (null point) | 0.024408 | 0.024408 | 0.024 | ✅ |
| μ_GWAM (sim mean) | 0.024619 → **0.025** | 0.024619 | 0.025 | ✅ |
| sim SD | 0.032699 | 0.032699 | — | ✅ |
| 95% **SI** q025 | −0.041357 → **−0.041** | −0.041357 | −0.041 | ✅ |
| 95% **SI** q975 | 0.087373 → **0.087** | 0.087373 | 0.087 | ✅ |
| Published RE log-OR | 0.2546 → **0.255** | 0.2546 | 0.255 | ✅ |
| Published OR = exp(0.2546) | 1.2899 → **1.29** | — | 1.29 | ✅ |

**Conclusion:** the published headline — *GWAM-adjusted log-OR 0.025 (95% SI −0.041 to 0.087) at λ=0.096,
vs published RE 0.255 (OR 1.29)* — is **correct and reproducible**. No correction needed.

The simulation reproduction used `numpy.random.default_rng(42)`, B=2000, σ_ghost=σ_ro=0.10, drawing results-only
effects then ghost effects in the same order as `scripts/model_gwam.py` — matching `mean/sd/q025/q975` exactly.
(The repo's own pipeline + tests are a further check; an `agy`/Codex pass was not required because the in-session
reproduction already matches the stored artifact to 6 decimals on an independent code path.)

## 2. Remaining abstract / table claims vs verified source (manuscript `paper/gwam_manuscript.md` + JSONs)

| Claim (final text) | Source value | Match |
|---|---|---|
| escitalopram 18 trials, λ 0.096 | raw CSV + JSON | ✅ |
| benchmark 7,467 MAs from 398 Cochrane reviews | manuscript §Pairwise70 | ✅ |
| median λ 0.871 | manuscript (Fig 7 / Results) | ✅ |
| 7.9% of |RE|≥0.10 attenuated below 0.10 (590/7467) | manuscript Results | ✅ |
| pregabalin 35 trials, λ 0.014 (104/7421) | manuscript App. 2 | ✅ |
| pregabalin GWAM 0.013, SI −0.029 to 0.058; log-RR 0.916, RR 2.50 | manuscript App. 2 | ✅ |
| simulation 2,000 reps; bias −0.162 at μ=0.2; raw coverage 7.1% | manuscript Table 2 | ✅ |
| Table 1 benchmark strata (NCT/query/transport/default; k=2/3/≥4) | manuscript Table 3 | ✅ |

All numbers in the final text, **Table 1**, and **visual abstract** trace to the verified source pack.

## 3. Honest distinction preserved (truth-first)

- The interval is labelled **SI = simulation interval, not a confidence interval** in the abstract, Table 1
  footnote, visual abstract, and a dedicated "Honest-reporting note" box in the PDF — **not relabelled as a CI**.
- λ is stated as a **lower bound** on the true publication fraction (incomplete CT.gov PMID linkage).
- GWAM framed as a **conservative sensitivity** tool that complements, not replaces, existing methods.

## 4. Finalization changes applied

1. **Removed all "preliminary" markers** — title suffix `[preliminary version]` dropped; the italic
   "Preliminary version… not peer reviewed" abstract notice deleted. (0 occurrences of "preliminary" remain in the
   PDF, JATS, or HTML — verified by grep.)
2. **JATS 1.3** now contains the previously-missing structured content: `<ref-list>` (3 refs: Egger 1997,
   Duval/Tweedie 2000, Cipriani 2018, with DOIs + PMIDs), `<table-wrap>` (Table 1, 10 data rows), and `<fig>`
   (Figure 1). Well-formedness validated; structured `<abstract>`, `<contrib>`+ORCID, `<permissions>` CC BY present.
3. **Regenerated PDF galley** in the **Gnosis** visual identity (midnight-blue #182650) with preliminary notices
   gone, structured abstract, Figure 1, Table 1, references, and a **visual abstract**.

## 5. Flags for you to resolve at routing time (NOT fabricated here)

- **DOI:** suppressed site-wide / pending — intentionally **not** emitted in JATS and shown as "pending" in the PDF.
  Do not mint a placeholder.
- **ISSN:** pending — shown as "pending"; insert into JATS `<issn>` when assigned.
- **Submission / publication IDs:** the live page's citation links show `submissionId=129 / publicationId=281`
  (galleys: PDF=395, JATS=396). Your brief said **sub 131 / pub 283** — please confirm which IDs the OJS update
  should target.
- **Journal name:** the PDF + JATS use **Gnosis** (per your framing and the Gnosis PDF identity), but the live OJS
  install currently records the journal title as **Synthēsis**. A note to this effect is in the JATS as an XML
  comment. Decide whether the final galley should carry Gnosis or Synthēsis branding before republish.
- **Affiliation:** set to "Royal Free Hospital, London, United Kingdom" (from the canonical manuscript). The live
  OJS page displayed no affiliation; the E156 protocol says "Tahir Heart Institute". Confirm the correct one.
- **Live galley files (395/396) were not on disk** and the live download URLs returned 404 at finalization time;
  this package was therefore rebuilt from the verified source pack + captured live HTML, not patched from the
  original galley binaries.
