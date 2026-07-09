#!/usr/bin/env python3
"""GWAM visual abstract -- Gnosis identity (midnight blue).
All numbers verbatim from the verified source pack."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
INK   = "#0A1020"; ACCENT = "#182650"; ACCENT2 = "#384C82"
BG    = "#FFFFFF"; PANEL = "#F2F4FA"; GHOST = "#9AA6C4"; PUB = "#182650"; RES = "#6E80B0"
GOLD  = "#B98A2E"

fig = plt.figure(figsize=(12, 6.4), dpi=200)
fig.patch.set_facecolor(BG)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 120); ax.set_ylim(0, 64); ax.axis("off")

def box(x, y, w, h, fc, ec=ACCENT, lw=1.2, r=0.04):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.2,rounding_size={r*100}",
                 fc=fc, ec=ec, lw=lw, mutation_aspect=0.5))
def arrow(x1, y1, x2, y2, c=ACCENT2, lw=2.2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
                 color=c, lw=lw, shrinkA=2, shrinkB=2))

# ---- header band (Gnosis identity)
ax.add_patch(plt.Rectangle((0, 57), 120, 7, color=ACCENT))
ax.text(2.5, 60.5, "GNOSIS", color="white", fontsize=15, fontweight="bold", va="center",
        family="DejaVu Serif")
ax.text(15.5, 61.4, "Global Health Evidence", color="#C9D2EA", fontsize=8, va="center")
ax.text(15.5, 59.2, "Diamond Open Access · Methods Note", color="#9AA6C4", fontsize=7, va="center")
ax.text(117.5, 60.5, "Vol. 2 (4) · 2026", color="#C9D2EA", fontsize=8, va="center", ha="right")

ax.text(2.5, 53.2, "GWAM: Ghost-Weighted Aggregate Meta-Analysis",
        color=INK, fontsize=16.5, fontweight="bold", va="center")
ax.text(2.5, 49.6, "A registry-anchored sensitivity analysis for publication bias",
        color=ACCENT2, fontsize=10.5, va="center", style="italic")

# ---- Stage 1: registry classification
sx = 2.5
ax.text(sx+9, 45.0, "1 · Classify registered trials", color=ACCENT, fontsize=9.5, fontweight="bold", ha="center")
labels = [("Published (PMID)", PUB, "white"), ("Results-posted", RES, "white"), ("Ghost protocols", GHOST, INK)]
yy = 41.5
for txt, fc, tc in labels:
    box(sx, yy, 18, 3.0, fc, ec=fc, r=0.05)
    ax.text(sx+9, yy+1.5, txt, color=tc, fontsize=8.2, ha="center", va="center", fontweight="bold")
    yy -= 4.0
ax.text(sx+9, 28.5, "weight by enrollment", color=ACCENT2, fontsize=7.6, ha="center", style="italic")

# ---- Stage 2: integrity ratio
ix = 30
arrow(21, 36.5, ix-0.5, 36.5)
ax.text(ix+11, 45.0, "2 · Integrity ratio λ", color=ACCENT, fontsize=9.5, fontweight="bold", ha="center")
box(ix, 33.0, 22, 9.5, PANEL, ec=ACCENT, lw=1.4, r=0.05)
ax.text(ix+11, 39.7, "λ = W$_{published}$ / W$_{total}$", color=INK, fontsize=10.5, ha="center", va="center")
ax.text(ix+11, 36.2, "escitalopram:", color=ACCENT2, fontsize=8, ha="center")
ax.text(ix+11, 34.3, "λ = 471 / 4913 = 0.096", color=GOLD, fontsize=9.5, ha="center", fontweight="bold")

# ---- Stage 3: attenuation result (escitalopram)
rx = 56
arrow(52.5, 37.7, rx-0.5, 37.7)
ax.text(rx+24, 45.0, "3 · Scale the published effect by λ", color=ACCENT, fontsize=9.5, fontweight="bold", ha="center")
box(rx, 30.5, 48, 12.0, "white", ec=ACCENT, lw=1.4, r=0.04)
# mini forest
fx0, fx1 = rx+3, rx+45
ax.plot([fx0, fx1], [40.0, 40.0], color="#C9D2EA", lw=0.8)  # published row
ax.plot([fx0, fx1], [34.0, 34.0], color="#C9D2EA", lw=0.8)  # gwam row
# zero line
import numpy as np
def xmap(v):  # map log-OR [-0.10, 0.50] to [fx0, fx1]
    return fx0 + (v - (-0.10)) / (0.50 - (-0.10)) * (fx1 - fx0)
ax.plot([xmap(0), xmap(0)], [31.5, 41.5], color=GHOST, lw=1.0, ls="--")
ax.text(xmap(0), 31.0, "no effect", color=GHOST, fontsize=6.6, ha="center", va="top")
# published point + CI 0.255 (0.066, 0.443)
ax.plot([xmap(0.066), xmap(0.443)], [40.0, 40.0], color=ACCENT, lw=2.0)
ax.plot(xmap(0.255), 40.0, "o", color=ACCENT, ms=8)
ax.text(fx0-0.5, 40.0, "Published RE", color=INK, fontsize=7.6, ha="right", va="center")
ax.text(xmap(0.255), 41.6, "log-OR 0.255 (OR 1.29)", color=ACCENT, fontsize=7.4, ha="center")
# gwam point + SI (-0.041, 0.087)
ax.plot([xmap(-0.041), xmap(0.087)], [34.0, 34.0], color=GOLD, lw=2.0)
ax.plot(xmap(0.025), 34.0, "s", color=GOLD, ms=8)
ax.text(fx0-0.5, 34.0, "GWAM-adjusted", color=INK, fontsize=7.6, ha="right", va="center")
ax.text(xmap(0.025), 32.4, "0.025 (95% SI −0.041, 0.087)", color=GOLD, fontsize=7.4, ha="center")

# ---- bottom strip: benchmark + honesty note
box(2.5, 14.0, 53, 12.5, PANEL, ec=ACCENT2, lw=1.2, r=0.03)
ax.text(29, 24.0, "Benchmark — 7,467 Cochrane meta-analyses", color=ACCENT, fontsize=9.5, fontweight="bold", ha="center")
ax.text(29, 20.6, "median λ = 0.871", color=INK, fontsize=10.5, ha="center")
ax.text(29, 17.4, "7.9%  of effects ≥ 0.10 attenuated below 0.10", color=GOLD, fontsize=10, ha="center", fontweight="bold")

box(58.5, 14.0, 59, 12.5, "#FBF7EC", ec=GOLD, lw=1.2, r=0.03)
ax.text(60.5, 24.0, "Reads as a sensitivity tool, not a correction", color="#7A5A10", fontsize=9.2,
        fontweight="bold", ha="left")
ax.text(60.5, 20.8, "• Interval is a simulation interval (SI), not a CI.", color=INK, fontsize=8.0, ha="left")
ax.text(60.5, 18.3, "• λ is a lower bound — CT.gov PMID linkage is incomplete.", color=INK, fontsize=8.0, ha="left")
ax.text(60.5, 15.8, "• Conservative null assumption can over-correct true effects.", color=INK, fontsize=8.0, ha="left")

# ---- footer
ax.add_patch(plt.Rectangle((0, 0), 120, 9.0, color=PANEL))
ax.text(2.5, 5.6, "Mahmood Ahmad", color=INK, fontsize=8.6, fontweight="bold", va="center")
ax.text(2.5, 3.0, "ORCID 0009-0003-7781-4478 · CC BY 4.0 · github.com/mahmood726-cyber/GWAM",
        color=ACCENT2, fontsize=7.2, va="center")
ax.text(117.5, 4.3, "Numbers verbatim from the verified GWAM source pack", color=GHOST, fontsize=7.0,
        va="center", ha="right", style="italic")

out_png = os.path.join(HERE, "figures", "visual_abstract.png")
fig.savefig(out_png, dpi=200, facecolor=BG, bbox_inches="tight", pad_inches=0.15)
print("wrote", out_png)
