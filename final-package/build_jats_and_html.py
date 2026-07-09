#!/usr/bin/env python3
"""Emit final JATS 1.3 XML (with ref-list, table-wrap, fig) and the Gnosis-branded
HTML galley (rendered to PDF separately via Chrome headless).
Reads content_bundle.json produced by build_content.py."""
import json, os, html

HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, "content_bundle.json"), encoding="utf-8"))
M, ID, ABS, REFS, T1, FIG1 = B["meta"], B["identity"], B["abstract"], B["references"], B["table1"], B["fig1_caption"]

def esc(s): return html.escape(str(s), quote=True)

# =====================================================================  JATS 1.3
def jats():
    o = []
    o.append('<?xml version="1.0" encoding="UTF-8"?>')
    o.append('<!DOCTYPE article PUBLIC "-//NLM//DTD JATS (Z39.96) Journal Publishing DTD v1.3 20210610//EN" '
             '"JATS-journalpublishing-1-3.dtd">')
    o.append('<article xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:mml="http://www.w3.org/1998/Math/MathML" '
             'dtd-version="1.3" article-type="other" xml:lang="en">')
    # ---- FRONT
    o.append("  <front>")
    o.append("    <journal-meta>")
    o.append('      <journal-id journal-id-type="publisher">gnosis</journal-id>')
    o.append("      <journal-title-group>")
    o.append(f"        <journal-title>{esc(M['journal'])}</journal-title>")
    o.append("        <abbrev-journal-title>Gnosis</abbrev-journal-title>")
    o.append("      </journal-title-group>")
    o.append("      <!-- ISSN pending; insert <issn pub-type='epub'> when assigned -->")
    o.append("      <publisher><publisher-name>Perfervid Consultancy Services</publisher-name></publisher>")
    o.append("      <!-- NOTE: the live OJS install currently records the journal title as 'Synthesis'; "
             "reconcile journal-title with the OJS journal record at routing time. -->")
    o.append("    </journal-meta>")
    o.append("    <article-meta>")
    o.append("      <!-- DOI suppressed site-wide / pending: no <article-id pub-id-type='doi'> emitted by design -->")
    o.append('      <article-categories><subj-group subj-group-type="heading">'
             f'<subject>{esc(M["journal_section"])}</subject></subj-group></article-categories>')
    o.append("      <title-group>")
    o.append(f"        <article-title>{esc(M['title'])}</article-title>")
    o.append("      </title-group>")
    o.append("      <contrib-group>")
    o.append('        <contrib contrib-type="author">')
    o.append(f'          <contrib-id contrib-id-type="orcid">https://orcid.org/{esc(M["orcid"])}</contrib-id>')
    o.append("          <name><surname>Ahmad</surname><given-names>Mahmood</given-names></name>")
    o.append('          <xref ref-type="aff" rid="aff1"/>')
    o.append(f"          <email>{esc(M['email'])}</email>")
    o.append("        </contrib>")
    o.append(f'        <aff id="aff1">{esc(M["affiliation"])}</aff>')
    o.append("      </contrib-group>")
    y, mth, d = M["published"].split("-")
    o.append(f'      <pub-date publication-format="electronic" date-type="pub"><day>{d}</day><month>{mth}</month><year>{y}</year></pub-date>')
    o.append(f"      <volume>{esc(M['volume'])}</volume>")
    o.append(f"      <issue>{esc(M['issue'])}</issue>")
    o.append("      <permissions>")
    o.append(f"        <copyright-statement>Copyright (c) {y} {esc(M['author'])}</copyright-statement>")
    o.append(f"        <copyright-year>{y}</copyright-year>")
    o.append('        <license license-type="open-access" xlink:href="https://creativecommons.org/licenses/by/4.0/">')
    o.append("          <license-p>This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0).</license-p>")
    o.append("        </license>")
    o.append("      </permissions>")
    o.append("      <self-uri content-type='code' xlink:href='%s'/>" % esc(M['code']))
    # structured abstract
    o.append("      <abstract>")
    for head, text in ABS:
        o.append(f"        <sec><title>{esc(head)}</title><p>{esc(text)}</p></sec>")
    o.append("      </abstract>")
    o.append("      <kwd-group>")
    for kw in ["publication bias", "meta-analysis", "trial registry", "ClinicalTrials.gov",
               "sensitivity analysis", "ghost protocols", "integrity ratio"]:
        o.append(f"        <kwd>{esc(kw)}</kwd>")
    o.append("      </kwd-group>")
    o.append("    </article-meta>")
    o.append("  </front>")
    # ---- BODY
    o.append("  <body>")
    o.append("    <sec id=\"s1\"><title>Summary</title>")
    o.append("      <p>GWAM (Ghost-Weighted Aggregate Meta-Analysis) is a registry-anchored sensitivity method for "
             "publication bias. It classifies ClinicalTrials.gov records as published (PMID-linked), results-posted, "
             "or ghost protocols, computes an enrollment-weighted integrity ratio "
             "<inline-formula><mml:math><mml:mi>&#955;</mml:mi></mml:math></inline-formula> (the published-participant "
             "fraction), and scales the published pooled effect by "
             "<inline-formula><mml:math><mml:mi>&#955;</mml:mi></mml:math></inline-formula>. A Monte-Carlo layer "
             "returns a simulation interval (SI, not a confidence interval) and a hierarchical variance-propagation "
             "extension carries the uncertainty about unobserved trials. The pipeline is shown in "
             "<xref ref-type=\"fig\" rid=\"fig1\">Fig 1</xref>; integrity ratios and attenuation across the two worked "
             "applications and the Pairwise70 benchmark are reported in <xref ref-type=\"table\" rid=\"tab1\">Table 1</xref>.</p>")
    o.append("    </sec>")
    # FIG 1
    o.append('    <fig id="fig1">')
    o.append("      <label>Fig 1</label>")
    o.append(f"      <caption><p>{esc(FIG1.replace('Figure 1. ',''))}</p></caption>")
    o.append('      <graphic xlink:href="figures/figure1_pipeline.png" mimetype="image" mime-subtype="png"/>')
    o.append("    </fig>")
    # TABLE 1 (table-wrap)
    o.append('    <table-wrap id="tab1">')
    o.append("      <label>Table 1</label>")
    o.append(f"      <caption><p>{esc(T1['title'].replace('Table 1. ',''))}</p></caption>")
    o.append('      <table frame="hsides" rules="groups">')
    o.append("        <thead><tr>" + "".join(f"<th>{esc(c)}</th>" for c in T1["cols"]) + "</tr></thead>")
    o.append("        <tbody>")
    for row in T1["rows"]:
        o.append("          <tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>")
    o.append("        </tbody>")
    o.append("      </table>")
    o.append(f"      <table-wrap-foot><fn><p>{esc(T1['footnote'])}</p></fn></table-wrap-foot>")
    o.append("    </table-wrap>")
    o.append("  </body>")
    # ---- BACK / ref-list
    o.append("  <back>")
    o.append('    <ref-list><title>References</title>')
    for i, r in enumerate(REFS, 1):
        o.append(f'      <ref id="R{i}">')
        o.append(f'        <element-citation publication-type="journal">')
        o.append(f"          <person-group person-group-type='author'><collab>{esc(r['authors'])}</collab></person-group>")
        o.append(f"          <article-title>{esc(r['title'])}</article-title>")
        o.append(f"          <source>{esc(r['source'].rstrip('.'))}</source>")
        o.append(f"          <year>{esc(r['year'])}</year>")
        o.append(f"          <volume>{esc(r['vol'])}</volume>")
        o.append(f"          <issue>{esc(r['issue'])}</issue>")
        o.append(f"          <fpage>{esc(r['fpage'])}</fpage><lpage>{esc(r['lpage'])}</lpage>")
        o.append(f"          <pub-id pub-id-type='doi'>{esc(r['doi'])}</pub-id>")
        o.append(f"          <pub-id pub-id-type='pmid'>{esc(r['pmid'])}</pub-id>")
        o.append("        </element-citation>")
        o.append("      </ref>")
    o.append("    </ref-list>")
    o.append("  </back>")
    o.append("</article>")
    return "\n".join(o)

# =====================================================================  HTML galley (Gnosis identity -> PDF)
def html_galley():
    accent, accent2, ink = ID["accent"], ID["accent2"], ID["ink"]
    refs_html = ""
    for i, r in enumerate(REFS, 1):
        refs_html += (f'<li><span class="rn">{i}.</span> {esc(r["authors"])} {esc(r["title"])} '
                      f'<i>{esc(r["source"])}</i> {esc(r["year"])};{esc(r["vol"])}({esc(r["issue"])}):'
                      f'{esc(r["fpage"])}&#8211;{esc(r["lpage"])}. doi:{esc(r["doi"])}. PMID: {esc(r["pmid"])}.</li>')
    abs_html = "".join(f'<p><span class="ahead">{esc(h)}.</span> {esc(t)}</p>' for h, t in ABS)
    # table
    thead = "".join(f"<th>{esc(c)}</th>" for c in T1["cols"])
    trows = ""
    prev_panel = None
    for row in T1["rows"]:
        panel = row[0]
        cls = ' class="panelstart"' if panel != prev_panel else ""
        prev_panel = panel
        cells = "".join(f"<td>{esc(c)}</td>" for c in row)
        trows += f"<tr{cls}>{cells}</tr>"
    doi_line = (f'DOI {esc(M["doi"])}' if M["doi"] else 'DOI: pending (suppressed site-wide)')
    issn_line = (esc(M["issn"]) if M["issn"] else "ISSN: pending")
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<style>
@page {{ size: A4; margin: 15mm 16mm 16mm 16mm; }}
* {{ box-sizing: border-box; }}
body {{ font-family: Georgia,'Times New Roman',serif; color:{ink}; line-height:1.5; font-size:10.3pt; margin:0; }}
.hero {{ background:{accent}; color:#fff; padding:14px 18px; display:flex; justify-content:space-between; align-items:flex-start; }}
.hero .j {{ font-family:Georgia,serif; font-size:20px; font-weight:bold; letter-spacing:.06em; }}
.hero .sub {{ font-size:8.5px; color:#C9D2EA; letter-spacing:.12em; text-transform:uppercase; margin-top:3px; }}
.hero .meta {{ text-align:right; font-size:8.5px; color:#C9D2EA; font-family:system-ui,sans-serif; line-height:1.5; }}
.titlewrap {{ border-bottom:2.5px solid {accent}; padding:14px 0 12px; margin-bottom:12px; }}
h1 {{ font-size:16.5pt; line-height:1.25; margin:0 0 8px; color:{ink}; }}
.byline {{ font-family:system-ui,sans-serif; font-size:9px; color:#444; }}
.byline a {{ color:{accent2}; text-decoration:none; }}
.section-label {{ font-family:system-ui,sans-serif; font-size:8px; font-weight:700; letter-spacing:.2em;
  text-transform:uppercase; color:{accent}; margin:16px 0 6px; border-bottom:1px solid #D8DEEC; padding-bottom:3px; }}
.abstract p {{ margin:0 0 6px; text-align:justify; }}
.ahead {{ font-weight:bold; color:{accent}; }}
.vabox {{ text-align:center; margin:8px 0 2px; }}
.vabox img {{ width:100%; border:1px solid #D8DEEC; border-radius:4px; }}
.cap {{ font-family:system-ui,sans-serif; font-size:8px; color:#555; margin:4px 2px 0; line-height:1.4; }}
figure {{ margin:8px 0; text-align:center; }}
figure img {{ width:86%; }}
table {{ width:100%; border-collapse:collapse; font-family:system-ui,sans-serif; font-size:7.6px; margin-top:4px; }}
th {{ background:{accent}; color:#fff; text-align:left; padding:4px 5px; font-weight:600; }}
td {{ padding:3px 5px; border-bottom:.5px solid #E2E6F0; vertical-align:top; }}
tr.panelstart td {{ border-top:1.2px solid {accent2}; }}
tbody tr:nth-child(even) td {{ background:#F6F8FC; }}
.refs ol {{ list-style:none; padding:0; margin:0; font-family:system-ui,sans-serif; font-size:8px; color:#333; }}
.refs li {{ margin-bottom:5px; padding-left:16px; text-indent:-16px; line-height:1.45; }}
.rn {{ color:{accent}; font-weight:700; }}
.meta-grid {{ font-family:system-ui,sans-serif; font-size:8px; color:#444; display:grid;
  grid-template-columns:auto 1fr; gap:2px 12px; margin-top:4px; }}
.meta-grid dt {{ font-weight:600; color:{ink}; }}
.flag {{ background:#FBF7EC; border-left:3px solid #B98A2E; padding:6px 9px; font-family:system-ui,sans-serif;
  font-size:8px; color:#5b4410; margin-top:8px; border-radius:0 3px 3px 0; }}
.footer {{ margin-top:14px; border-top:1px solid #D8DEEC; padding-top:6px; font-family:system-ui,sans-serif;
  font-size:7.5px; color:#777; display:flex; justify-content:space-between; }}
</style></head><body>
<div class="hero"><div><div class="j">{esc(ID['title'])}</div><div class="sub">{esc(ID['subtitle'])} &middot; {esc(ID['tag'])}</div></div>
<div class="meta">{esc(M['journal_section'])}<br>Vol. {esc(M['volume'])} No. {esc(M['issue'])} &middot; {esc(M['year'])}<br>Published {esc(M['published'])}</div></div>

<div class="titlewrap">
<h1>{esc(M['title'])}</h1>
<div class="byline"><b>{esc(M['author'])}</b> &nbsp;
<a href="https://orcid.org/{esc(M['orcid'])}">ORCID {esc(M['orcid'])}</a> &middot; {esc(M['affiliation'])} &middot;
<a href="mailto:{esc(M['email'])}">{esc(M['email'])}</a></div>
</div>

<div class="vabox"><img src="figures/visual_abstract.png" alt="Visual abstract">
<div class="cap">Visual abstract. Registry classification &rarr; integrity ratio &lambda; &rarr; &lambda;-scaled effect, with the escitalopram application and the Pairwise70 benchmark. All values are verbatim from the verified GWAM source pack.</div></div>

<div class="section-label">Abstract</div>
<div class="abstract">{abs_html}</div>

<div class="section-label">Pipeline</div>
<figure><img src="figures/figure1_pipeline.png" alt="GWAM pipeline">
<figcaption class="cap">{esc(FIG1)}</figcaption></figure>

<div class="section-label">Table 1 &mdash; Integrity ratio (&lambda;) and GWAM attenuation</div>
<table><thead><tr>{thead}</tr></thead><tbody>{trows}</tbody></table>
<div class="cap">{esc(T1['footnote'])}</div>

<div class="section-label">References</div>
<div class="refs"><ol>{refs_html}</ol></div>

<div class="section-label">Article information</div>
<dl class="meta-grid">
<dt>Journal</dt><dd>{esc(M['journal'])} &mdash; {esc(ID['subtitle'])} ({issn_line})</dd>
<dt>Section</dt><dd>{esc(M['journal_section'])}</dd>
<dt>Identifier</dt><dd>{doi_line}</dd>
<dt>Licence</dt><dd>{esc(M['license'])} &mdash; Copyright (c) {esc(M['year'])} {esc(M['author'])}</dd>
<dt>Code &amp; data</dt><dd>{esc(M['code'])}</dd>
</dl>
<div class="flag"><b>Honest-reporting note.</b> The 95% interval on the GWAM-adjusted estimate is a
<b>simulation interval (SI)</b>, not a confidence interval. The integrity ratio &lambda; is a <b>lower bound</b> on
the true publication fraction because ClinicalTrials.gov PMID linkage is incomplete; GWAM is a conservative
<b>sensitivity</b> framework that complements, and does not replace, existing publication-bias methods.</div>

<div class="footer"><span>{esc(ID['title'])} &middot; {esc(ID['tag'])}</span><span>Numbers verbatim from the verified GWAM source pack</span></div>
</body></html>"""

with open(os.path.join(HERE, "gwam_final.jats.xml"), "w", encoding="utf-8") as f:
    f.write(jats())
with open(os.path.join(HERE, "gwam_galley.html"), "w", encoding="utf-8") as f:
    f.write(html_galley())
print("wrote gwam_final.jats.xml and gwam_galley.html")
