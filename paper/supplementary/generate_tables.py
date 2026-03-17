#!/usr/bin/env python
"""Generate supplementary CSV tables from GWAM analysis outputs."""

import json
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
ANALYSIS = ROOT / "data" / "analysis"
P70_GRMA = ROOT / "pairwise70_benchmark_grma"
OUTDIR = Path(__file__).resolve().parent

# Mapping from method names used in loops to JSON ci_calibration keys
MULTIPLIER_KEY_MAP = {
    'random_effects': 'multiplier_random_effects',
    'gwam_null': 'multiplier_gwam',
    'bayesian_gwam': 'multiplier_bayesian_gwam',
    'oracle_ipw_random_effects': 'multiplier_oracle_ipw_random_effects',
    'pet_peese': 'multiplier_pet_peese',
}


def s1_hvp_sensitivity():
    """S1: HVP sensitivity grid for both applications."""
    rows = []
    for app_name, fname in [
        ("Escitalopram/Depression", "escitalopram_depression_gwam_bayesian.json"),
        ("Pregabalin/Neuropathic Pain", "pregabalin_neuropathic_pain_gwam_bayesian.json"),
    ]:
        data = json.loads((ANALYSIS / fname).read_text())
        for entry in data['sensitivity_table']:
            rows.append({
                'application': app_name,
                'ghost_sigma': entry['ghost_sigma'],
                'results_only_sigma': entry['results_only_sigma'],
                'posterior_mean': round(entry['posterior_mean'], 6),
                'posterior_sd': round(entry['posterior_sd'], 6),
                'cri_lo': round(entry['cri_lo'], 6),
                'cri_hi': round(entry['cri_hi'], 6),
                'cri_width': round(entry['cri_width'], 6),
                'pr_positive': round(entry['pr_positive'], 6),
            })

    outpath = OUTDIR / "S1_hvp_sensitivity.csv"
    with open(outpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  S1: {len(rows)} rows -> {outpath.name}")


def s2_simulation_results():
    """S2: Full simulation results for both applications."""
    rows = []
    methods = ['random_effects', 'gwam_null', 'bayesian_gwam',
               'oracle_ipw_random_effects', 'pet_peese']

    for app_name, fname in [
        ("Escitalopram/Depression", "simulation_escitalopram_depression.json"),
        ("Pregabalin/Neuropathic Pain", "simulation_pregabalin_neuropathic_pain.json"),
    ]:
        data = json.loads((ANALYSIS / fname).read_text())
        calib = data.get('calibration', {})
        ci_calib = data.get('ci_calibration', {})

        for scenario in data['scenario_results']:
            mu_true = scenario['mu_true']
            for method in methods:
                d = scenario.get(method, {})
                if not d:
                    continue
                multiplier_key = MULTIPLIER_KEY_MAP.get(method, f'multiplier_{method}')
                multiplier_val = ci_calib.get(multiplier_key)
                # Map internal key to user-facing label
                method_label = 'hvp_gwam' if method == 'bayesian_gwam' else method
                rows.append({
                    'application': app_name,
                    'mu_true': mu_true,
                    'method': method_label,
                    'n_valid': int(d.get('n_valid', 0)),
                    'mean_estimate': round(d.get('mean_estimate', 0), 6),
                    'bias': round(d.get('bias', 0), 6),
                    'rmse': round(d.get('rmse', 0), 6),
                    'coverage_raw': round(d.get('coverage_95', 0), 4),
                    'coverage_calibrated': round(d.get('coverage_calibrated', 0), 4),
                    'mean_ci_width_raw': round(d.get('mean_ci_width', 0), 4),
                    'mean_ci_width_calibrated': round(d.get('mean_ci_width_calibrated', 0), 4),
                    'lambda_calibrated': calib.get('lambda_calibration_within_tolerance', ''),
                    'ci_multiplier': round(multiplier_val, 4) if multiplier_val is not None else '',
                })

    outpath = OUTDIR / "S2_simulation_results.csv"
    with open(outpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  S2: {len(rows)} rows -> {outpath.name}")


def s3_pairwise70_stratified():
    """S3: Pairwise70 benchmark stratified summary."""
    summary = json.loads((P70_GRMA / "summary.json").read_text())
    rows = []

    # By lambda source
    for source, d in summary['stratified']['by_lambda_proxy_source'].items():
        dist = d['distribution']
        att = d['attenuation']
        rows.append({
            'stratification': 'lambda_source',
            'subset': source,
            'n_analyses': d['n_analyses'],
            'median_k': dist['median_k'],
            'median_lambda_proxy': round(dist['median_lambda_proxy'], 4),
            'median_abs_shift': round(dist['median_abs_shift_proxy'], 4),
            'crossed_below_010_pct_all': round(att['re_to_gwam_proxy_crossed_below_0p10_pct_all'], 4),
            'crossed_below_020_pct_all': round(att['re_to_gwam_proxy_crossed_below_0p20_pct_all'], 4),
            'crossed_below_010_pct_within_re': round(att['re_to_gwam_proxy_crossed_below_0p10_pct_within_re_ge'], 4),
            'crossed_below_020_pct_within_re': round(att['re_to_gwam_proxy_crossed_below_0p20_pct_within_re_ge'], 4),
        })

    # By k band
    for band, d in summary['stratified']['by_k_band'].items():
        dist = d['distribution']
        att = d['attenuation']
        rows.append({
            'stratification': 'k_band',
            'subset': band,
            'n_analyses': d['n_analyses'],
            'median_k': dist['median_k'],
            'median_lambda_proxy': round(dist['median_lambda_proxy'], 4),
            'median_abs_shift': round(dist['median_abs_shift_proxy'], 4),
            'crossed_below_010_pct_all': round(att['re_to_gwam_proxy_crossed_below_0p10_pct_all'], 4),
            'crossed_below_020_pct_all': round(att['re_to_gwam_proxy_crossed_below_0p20_pct_all'], 4),
            'crossed_below_010_pct_within_re': round(att['re_to_gwam_proxy_crossed_below_0p10_pct_within_re_ge'], 4),
            'crossed_below_020_pct_within_re': round(att['re_to_gwam_proxy_crossed_below_0p20_pct_within_re_ge'], 4),
        })

    outpath = OUTDIR / "S3_pairwise70_stratified.csv"
    with open(outpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  S3: {len(rows)} rows -> {outpath.name}")


if __name__ == "__main__":
    print("Generating supplementary tables...")
    s1_hvp_sensitivity()
    s2_simulation_results()
    s3_pairwise70_stratified()
    print("Done.")
