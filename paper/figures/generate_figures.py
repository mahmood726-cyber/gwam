#!/usr/bin/env python
"""
GWAM Paper — Publication-Quality Figure Generator
Reads analysis JSONs/CSVs and produces 8 figures (PDF + PNG).

Usage:
    python paper/figures/generate_figures.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent.parent
DATA = ROOT / "data"
ANALYSIS = DATA / "analysis"
P70_GRMA = ROOT / "pairwise70_benchmark_grma"
P70_SENS = ROOT / "pairwise70_benchmark_sensitivity"
OUTDIR = Path(__file__).resolve().parent

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

COLORS = {
    'pmid': '#2166ac',
    'results_only': '#67a9cf',
    'ghost': '#ef8a62',
    're': '#d73027',
    'gwam': '#4575b4',
    'bayesian': '#91bfdb',
    'oracle': '#fc8d59',
    'pet_peese': '#fee090',
    'base': '#4575b4',
    'strict': '#d73027',
    'loose': '#1a9850',
}


def save_fig(fig, name):
    """Save figure as both PDF and PNG."""
    fig.savefig(OUTDIR / f"{name}.pdf", format='pdf')
    fig.savefig(OUTDIR / f"{name}.png", format='png')
    plt.close(fig)
    print(f"  Saved {name}.pdf + {name}.png")


def find_scenario(scenario_results, mu_true):
    """Find scenario by mu_true value instead of fragile index access."""
    for s in scenario_results:
        if abs(s['mu_true'] - mu_true) < 1e-9:
            return s
    raise ValueError(f"No scenario with mu_true={mu_true}")


# ── Fig 1: Pipeline Flow Diagram ──────────────────────────────────────────
def fig1_pipeline():
    """Programmatic pipeline flow diagram."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis('off')
    ax.set_title("Fig 1. GWAM Pipeline Overview", fontweight='bold', pad=15)

    boxes = [
        (1.0, 5.5, "Trial Registry\n(ClinicalTrials.gov)", '#e0e0e0'),
        (4.0, 5.5, "Classify Trials\n(PMID / Results / Ghost)", '#b3cde3'),
        (7.5, 5.5, "Compute Weights\n(Enrollment-based)", '#ccebc5'),
        (1.0, 3.0, "Integrity Ratio\n(\u03bb = W_obs / W_total)", '#decbe4'),
        (4.0, 3.0, "GWAM Estimate\n(Weighted average)", '#fed9a6'),
        (7.5, 3.0, "Simulation CI\n(Calibrated)", '#fddaec'),
        (4.0, 0.8, "Bayesian Extension\n(Prior on ghost effects)", '#f2f0f7'),
    ]

    for x, y, text, color in boxes:
        box = FancyBboxPatch((x - 1.1, y - 0.55), 2.2, 1.1,
                             boxstyle="round,pad=0.1",
                             facecolor=color, edgecolor='#333333', linewidth=1.2)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=8.5,
                fontweight='bold', wrap=True)

    arrows = [
        (2.1, 5.5, 2.9, 5.5),
        (5.1, 5.5, 6.4, 5.5),
        (1.0, 4.95, 1.0, 3.55),
        (7.5, 4.95, 7.5, 3.55),
        (2.1, 3.0, 2.9, 3.0),
        (5.1, 3.0, 6.4, 3.0),
        (4.0, 2.45, 4.0, 1.35),
    ]
    for x1, y1, x2, y2 in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))

    save_fig(fig, 'fig1_pipeline')


# ── Fig 2: Candidate Scan Ghost Ratios ────────────────────────────────────
def fig2_candidate_scan():
    """Bar chart of ghost protocol ratios across drug-condition pairs."""
    df = pd.read_csv(DATA / "registry_candidate_summary.csv")

    labels = [f"{row['intervention'].title()}\n({row['condition'].title()})"
              for _, row in df.iterrows()]
    ghost_pct = df['ghost_ratio'].values * 100
    total = df['total_completed'].values

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(labels))
    width = 0.55

    bars = ax.bar(x, ghost_pct, width, color=COLORS['ghost'], edgecolor='#333',
                  linewidth=0.8, zorder=3)

    for i, (bar, n) in enumerate(zip(bars, total)):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"n={n}", ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Ghost Protocol Rate (%)")
    ax.set_title("Fig 2. Ghost Protocol Rates Across Drug-Condition Pairs",
                 fontweight='bold')
    ax.set_ylim(0, max(ghost_pct) + 5)
    ax.grid(axis='y', alpha=0.3, zorder=0)
    ax.axhline(y=20, color='grey', linestyle='--', alpha=0.5, label='20% reference')
    ax.legend(fontsize=8)

    save_fig(fig, 'fig2_candidate_scan')


# ── Fig 3: Forest-Style Published vs GWAM ─────────────────────────────────
def fig3_forest():
    """Forest-style comparison: Published RE vs GWAM for both applications."""
    esc = json.loads((ANALYSIS / "escitalopram_depression_gwam.json").read_text())
    preg = json.loads((ANALYSIS / "pregabalin_neuropathic_pain_gwam.json").read_text())
    esc_bay = json.loads((ANALYSIS / "escitalopram_depression_gwam_bayesian.json").read_text())
    preg_bay = json.loads((ANALYSIS / "pregabalin_neuropathic_pain_gwam_bayesian.json").read_text())

    # Published RE CIs from published_se (z-based)
    esc_mu_re = esc['estimates']['mu_random_effects_published_only']
    esc_se = esc_bay['inputs']['published_se']
    preg_mu_re = preg['estimates']['mu_random_effects_published_only']
    preg_se = preg_bay['inputs']['published_se']

    items = [
        ("Escitalopram/Depression", "Published RE",
         esc_mu_re, esc_mu_re - 1.96 * esc_se, esc_mu_re + 1.96 * esc_se),
        ("Escitalopram/Depression", "GWAM",
         esc['estimates']['mu_gwam_sim_mean'],
         esc['estimates']['mu_gwam_sim_q025'],
         esc['estimates']['mu_gwam_sim_q975']),
        ("Pregabalin/Neuropathic Pain", "Published RE",
         preg_mu_re, preg_mu_re - 1.96 * preg_se, preg_mu_re + 1.96 * preg_se),
        ("Pregabalin/Neuropathic Pain", "GWAM",
         preg['estimates']['mu_gwam_sim_mean'],
         preg['estimates']['mu_gwam_sim_q025'],
         preg['estimates']['mu_gwam_sim_q975']),
    ]

    fig, ax = plt.subplots(figsize=(8, 4))
    y_positions = [3.5, 3.0, 1.5, 1.0]
    colors_list = [COLORS['re'], COLORS['gwam'], COLORS['re'], COLORS['gwam']]

    for i, (app, method, mu, lo, hi) in enumerate(items):
        y = y_positions[i]
        color = colors_list[i]
        marker = 'o' if method == "Published RE" else 's'
        ax.plot(mu, y, marker, color=color, markersize=10, zorder=5)
        ax.hlines(y, lo, hi, color=color, linewidth=2.5, zorder=4)
        label_text = f"{method}: {mu:.3f} [{lo:.3f}, {hi:.3f}]"
        ax.text(max(mu, hi) + 0.02, y, label_text, va='center', fontsize=8.5)

    ax.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
    ax.set_yticks([3.25, 1.25])
    ax.set_yticklabels(["Escitalopram /\nDepression", "Pregabalin /\nNeuropathic Pain"],
                       fontsize=10)
    ax.set_xlabel("Effect Size (log-OR / log-RR)")
    ax.set_title("Fig 3. Published vs GWAM Estimates", fontweight='bold')
    ax.set_ylim(0.3, 4.2)

    re_patch = mpatches.Patch(color=COLORS['re'], label='Published RE')
    gwam_patch = mpatches.Patch(color=COLORS['gwam'], label='GWAM')
    ax.legend(handles=[re_patch, gwam_patch], loc='upper right', fontsize=9)
    ax.grid(axis='x', alpha=0.2)

    save_fig(fig, 'fig3_forest')


# ── Fig 4: Bayesian Sensitivity Heatmap ───────────────────────────────────
def fig4_bayesian_heatmap():
    """Bayesian sensitivity heatmap for escitalopram: Pr(positive)."""
    bayes = json.loads((ANALYSIS / "escitalopram_depression_gwam_bayesian.json").read_text())

    grid_sigmas = [0.05, 0.1, 0.2, 0.3]
    n = len(grid_sigmas)
    pr_positive = np.zeros((n, n))
    cri_width = np.zeros((n, n))

    for entry in bayes['sensitivity_table']:
        i = grid_sigmas.index(entry['ghost_sigma'])
        j = grid_sigmas.index(entry['results_only_sigma'])
        pr_positive[i, j] = entry['pr_positive']
        cri_width[i, j] = entry['cri_width']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: Pr(positive) — data-driven color range
    pr_min, pr_max = pr_positive.min(), pr_positive.max()
    im1 = ax1.imshow(pr_positive, cmap='RdYlBu_r', vmin=pr_min - 0.005,
                     vmax=pr_max + 0.005, origin='lower', aspect='auto')
    ax1.set_xticks(range(n))
    ax1.set_xticklabels(grid_sigmas)
    ax1.set_yticks(range(n))
    ax1.set_yticklabels(grid_sigmas)
    ax1.set_xlabel("Results-only \u03c3")
    ax1.set_ylabel("Ghost \u03c3")
    ax1.set_title("Pr(effect > 0)")
    mid_val = (pr_min + pr_max) / 2
    for i in range(n):
        for j in range(n):
            ax1.text(j, i, f"{pr_positive[i,j]:.3f}", ha='center', va='center',
                     fontsize=8, color='white' if pr_positive[i,j] > mid_val else 'black')
    fig.colorbar(im1, ax=ax1, shrink=0.8)

    # Right: CrI width — data-driven color range
    im2 = ax2.imshow(cri_width, cmap='YlOrRd', origin='lower', aspect='auto')
    ax2.set_xticks(range(n))
    ax2.set_xticklabels(grid_sigmas)
    ax2.set_yticks(range(n))
    ax2.set_yticklabels(grid_sigmas)
    ax2.set_xlabel("Results-only \u03c3")
    ax2.set_ylabel("Ghost \u03c3")
    ax2.set_title("95% CrI Width")
    for i in range(n):
        for j in range(n):
            ax2.text(j, i, f"{cri_width[i,j]:.3f}", ha='center', va='center',
                     fontsize=8)
    fig.colorbar(im2, ax=ax2, shrink=0.8)

    fig.suptitle("Fig 4. Bayesian Sensitivity Analysis (Escitalopram/Depression)",
                 fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'fig4_bayesian_sensitivity')


# ── Fig 5: Simulation Bias + RMSE Panel ───────────────────────────────────
def fig5_simulation_bias_rmse():
    """Simulation: bias + RMSE panel (3 mu_true x 5 methods) for escitalopram."""
    sim = json.loads((ANALYSIS / "simulation_escitalopram_depression.json").read_text())

    methods = ['random_effects', 'gwam_null', 'bayesian_gwam',
               'oracle_ipw_random_effects', 'pet_peese']
    method_labels = ['RE', 'GWAM', 'Bayesian\nGWAM', 'Oracle\nIPW', 'PET-\nPEESE']
    method_colors = [COLORS['re'], COLORS['gwam'], COLORS['bayesian'],
                     COLORS['oracle'], COLORS['pet_peese']]
    mu_trues = [0.0, 0.1, 0.2]

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharey='row')

    for col, mu_true in enumerate(mu_trues):
        scenario = find_scenario(sim['scenario_results'], mu_true)
        biases = []
        rmses = []
        for m in methods:
            d = scenario.get(m, {})
            biases.append(d.get('bias', np.nan))
            rmses.append(d.get('rmse', np.nan))

        x = np.arange(len(methods))

        # Bias
        ax_b = axes[0, col]
        ax_b.bar(x, biases, color=method_colors, edgecolor='#333',
                 linewidth=0.6, zorder=3)
        ax_b.axhline(0, color='grey', linestyle='--', alpha=0.5)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(method_labels, fontsize=7)
        ax_b.set_title(f"\u03bc_true = {mu_true}", fontweight='bold')
        if col == 0:
            ax_b.set_ylabel("Bias")
        ax_b.grid(axis='y', alpha=0.2, zorder=0)

        # RMSE
        ax_r = axes[1, col]
        ax_r.bar(x, rmses, color=method_colors, edgecolor='#333',
                 linewidth=0.6, zorder=3)
        ax_r.set_xticks(x)
        ax_r.set_xticklabels(method_labels, fontsize=7)
        if col == 0:
            ax_r.set_ylabel("RMSE")
        ax_r.grid(axis='y', alpha=0.2, zorder=0)

    fig.suptitle("Fig 5. Simulation: Bias and RMSE (Escitalopram Calibration, n=2000)",
                 fontweight='bold', y=1.01)
    fig.tight_layout()
    save_fig(fig, 'fig5_simulation_bias_rmse')


# ── Fig 6: Simulation Coverage Comparison ─────────────────────────────────
def fig6_simulation_coverage():
    """Simulation: calibrated coverage comparison panel."""
    sim = json.loads((ANALYSIS / "simulation_escitalopram_depression.json").read_text())

    methods = ['random_effects', 'gwam_null', 'bayesian_gwam',
               'oracle_ipw_random_effects', 'pet_peese']
    method_labels = ['RE', 'GWAM', 'Bayesian GWAM', 'Oracle IPW', 'PET-PEESE']
    method_colors = [COLORS['re'], COLORS['gwam'], COLORS['bayesian'],
                     COLORS['oracle'], COLORS['pet_peese']]
    mu_trues = [0.0, 0.1, 0.2]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)

    for col, mu_true in enumerate(mu_trues):
        scenario = find_scenario(sim['scenario_results'], mu_true)
        ax = axes[col]

        uncalibrated = []
        calibrated = []
        for m in methods:
            d = scenario.get(m, {})
            uncalibrated.append(d.get('coverage_95', np.nan))
            calibrated.append(d.get('coverage_calibrated', np.nan))

        x = np.arange(len(methods))
        w = 0.35
        ax.bar(x - w/2, uncalibrated, w,
               color=[c for c in method_colors], alpha=0.4,
               edgecolor='#333', linewidth=0.6, label='Raw 95% CI', zorder=3)
        ax.bar(x + w/2, calibrated, w, color=method_colors,
               edgecolor='#333', linewidth=0.6, label='Calibrated CI', zorder=3)

        ax.axhline(0.95, color='black', linestyle='--', linewidth=1.2,
                   label='Nominal 95%', zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels(method_labels, fontsize=7, rotation=30, ha='right')
        ax.set_title(f"\u03bc_true = {mu_true}", fontweight='bold')
        ax.set_ylim(0, 1.08)
        if col == 0:
            ax.set_ylabel("Coverage")
            ax.legend(fontsize=7, loc='lower left')
        ax.grid(axis='y', alpha=0.2, zorder=0)

    fig.suptitle("Fig 6. Simulation: CI Coverage (Escitalopram Calibration, n=2000)",
                 fontweight='bold', y=1.01)
    fig.tight_layout()
    save_fig(fig, 'fig6_simulation_coverage')


# ── Fig 7: Pairwise70 Attenuation Distribution ────────────────────────────
def fig7_attenuation_distribution():
    """Histogram of lambda_proxy and |shift| distribution across Pairwise70."""
    df = pd.read_csv(P70_GRMA / "analysis_results.csv")
    n_analyses = len(df)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: Lambda proxy distribution
    lam = df['lambda_proxy'].dropna()
    ax1.hist(lam, bins=50, color=COLORS['gwam'], edgecolor='#333', linewidth=0.4,
             alpha=0.85, zorder=3)
    ax1.axvline(lam.median(), color=COLORS['re'], linestyle='--', linewidth=1.5,
                label=f"Median = {lam.median():.3f}")
    ax1.set_xlabel("\u03bb (Integrity Ratio)")
    ax1.set_ylabel("Count")
    ax1.set_title("Distribution of \u03bb_proxy")
    ax1.legend(fontsize=9)
    ax1.grid(axis='y', alpha=0.2, zorder=0)

    # Right: Absolute shift distribution
    shift = (df['mu_gwam_proxy'] - df['mu_re']).abs().dropna()
    ax2.hist(shift[shift < 1.0], bins=60, color=COLORS['ghost'], edgecolor='#333',
             linewidth=0.4, alpha=0.85, zorder=3)
    ax2.axvline(shift.median(), color=COLORS['re'], linestyle='--', linewidth=1.5,
                label=f"Median = {shift.median():.3f}")
    ax2.set_xlabel("|GWAM - RE| (absolute shift)")
    ax2.set_ylabel("Count")
    ax2.set_title("Distribution of Attenuation")
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', alpha=0.2, zorder=0)

    fig.suptitle(f"Fig 7. Pairwise70 Benchmark: \u03bb and Attenuation Distributions "
                 f"(n={n_analyses:,})", fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'fig7_attenuation_distribution')


# ── Fig 8: Lambda Sensitivity (Base/Strict/Loose) ─────────────────────────
def fig8_lambda_sensitivity():
    """Side-by-side comparison of base/strict/loose lambda settings."""
    df = pd.read_csv(P70_SENS / "side_by_side_summary.csv")

    scenarios = df['scenario'].tolist()
    cross_010 = df['crossed_below_0p10_pct'].values * 100
    cross_020 = df['crossed_below_0p20_pct'].values * 100
    med_shift = df['median_abs_shift_proxy'].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    x = np.arange(len(scenarios))
    w = 0.3
    scenario_colors = [COLORS.get(s, '#999999') for s in scenarios]

    # Left: Crossing rates
    ax1.bar(x - w/2, cross_010, w, color=scenario_colors, edgecolor='#333',
            linewidth=0.8, alpha=0.7, zorder=3)
    ax1.bar(x + w/2, cross_020, w, color=scenario_colors, edgecolor='#333',
            linewidth=0.8, alpha=1.0, zorder=3)

    for i in range(len(scenarios)):
        ax1.text(x[i] - w/2, cross_010[i] + 0.3, f"{cross_010[i]:.1f}%",
                 ha='center', fontsize=7.5, fontweight='bold')
        ax1.text(x[i] + w/2, cross_020[i] + 0.3, f"{cross_020[i]:.1f}%",
                 ha='center', fontsize=7.5, fontweight='bold')

    ax1.set_xticks(x)
    ax1.set_xticklabels([s.title() for s in scenarios])
    ax1.set_ylabel("Attenuation Rate (%)")
    ax1.set_title("Crossed-Below Rates")
    ax1.grid(axis='y', alpha=0.2, zorder=0)
    ax1.legend(handles=[
        mpatches.Patch(facecolor='grey', alpha=0.7, label='Crossed |0.10|'),
        mpatches.Patch(facecolor='grey', alpha=1.0, label='Crossed |0.20|'),
    ], fontsize=8)

    # Right: Median absolute shift
    ax2.bar(x, med_shift, 0.5, color=scenario_colors, edgecolor='#333',
            linewidth=0.8, zorder=3)
    for i in range(len(scenarios)):
        ax2.text(x[i], med_shift[i] + 0.002, f"{med_shift[i]:.4f}",
                 ha='center', fontsize=8, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels([s.title() for s in scenarios])
    ax2.set_ylabel("Median |GWAM - RE|")
    ax2.set_title("Median Absolute Shift")
    ax2.grid(axis='y', alpha=0.2, zorder=0)

    fig.suptitle("Fig 8. Lambda Sensitivity: Base / Strict / Loose Configurations",
                 fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'fig8_lambda_sensitivity')


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("GWAM Figure Generator")
    print(f"  Root: {ROOT}")
    print(f"  Output: {OUTDIR}")
    print()

    fig1_pipeline()
    fig2_candidate_scan()
    fig3_forest()
    fig4_bayesian_heatmap()
    fig5_simulation_bias_rmse()
    fig6_simulation_coverage()
    fig7_attenuation_distribution()
    fig8_lambda_sensitivity()

    print("\nAll 8 figures generated successfully.")


if __name__ == "__main__":
    main()
