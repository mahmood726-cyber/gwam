# GWAM — finalized publication package

Final-form package for **"GWAM: Ghost-Weighted Aggregate Meta-Analysis for registry-based publication-bias
sensitivity"** (Gnosis / Synthēsis OJS, Vol. 2 No. 4, Methods Note). Brings the preliminary article to FINAL
form: preliminary markers removed, structured JATS (ref-list + table + fig), regenerated Gnosis PDF galley,
and a visual abstract. **Not yet published** — routing left to the author.

## Deliverables

| File | What it is |
|---|---|
| `gwam_final_galley.pdf` | **PDF galley** (Gnosis identity) — 3 pp; structured abstract, visual abstract, Fig 1, Table 1, refs, honest-note. Upload as the PDF galley. |
| `gwam_final.jats.xml` | **JATS 1.3** — structured `<abstract>` + `<ref-list>` (3 refs) + `<table-wrap>` (Table 1, 10 rows) + `<fig>` (Fig 1). Upload as the JATS galley. |
| `abstract_final.txt` | Plain structured abstract (no "preliminary") for the OJS abstract field. |
| `paper_final.md` | Final E156-standard text (title, author, abstract, Table 1, Fig 1, refs, article info). |
| `figures/visual_abstract.png` | Visual abstract (2460×1340). |
| `figures/figure1_pipeline.png` | Figure 1 (pipeline) — dependent file for JATS. |
| `figures/gnosis_hero.png`, `gnosis_cover.png` | Gnosis identity assets. |
| `table1.csv` | Table 1 data (10 rows). |
| `content_bundle.json` | Single source of all final content (drives JATS + HTML). |
| `VERIFICATION_REPORT.md` | Independent reproduction of the headline + claim-by-claim source check + flags. |
| `build_content.py`, `build_visual_abstract.py`, `build_jats_and_html.py` | Regenerators (deterministic). |
| `gwam_galley.html` | HTML source the PDF is rendered from. |

## Regenerate
```
python build_content.py            # table1.csv + content_bundle.json
python build_visual_abstract.py    # figures/visual_abstract.png
python build_jats_and_html.py      # gwam_final.jats.xml + gwam_galley.html
# PDF: render gwam_galley.html with Playwright/Chromium A4, no header/footer
```

## Before republish — confirm (see VERIFICATION_REPORT.md §5)
DOI pending (do not fabricate) · ISSN pending · sub/pub IDs (live=129/281 vs brief=131/283) ·
journal name Gnosis vs Synthēsis · author affiliation · attach `figures/*.png` as JATS dependent files.

**Headline (verified, reproduced to 6 dp):** escitalopram λ=0.096 → GWAM-adjusted log-OR **0.025
(95% simulation interval −0.041 to 0.087)** vs published RE **0.255 (OR 1.29)**.
