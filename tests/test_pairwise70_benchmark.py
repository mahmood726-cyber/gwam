#!/usr/bin/env python3
"""Unit tests for Pairwise70 benchmark helpers."""

from __future__ import annotations

import unittest

import pandas as pd

# sys.path setup handled by conftest.py
from run_pairwise70_benchmark import (
    apply_lambda_bounds,
    build_lambda_posterior_params,
    build_review_cluster_summary,
    calibrate_p_nonsig_from_lambda,
    extract_review_id,
    make_analysis_key,
    parse_lambda_grid,
)


class TestPairwise70BenchmarkHelpers(unittest.TestCase):
    def test_parse_lambda_grid(self) -> None:
        vals = parse_lambda_grid("0.2,0.4,0.8")
        self.assertEqual(vals, [0.2, 0.4, 0.8])

    def test_extract_review_id(self) -> None:
        self.assertEqual(extract_review_id("CD002042_pub6_data.rda"), "CD002042")
        self.assertEqual(extract_review_id("CD014089_data.rda"), "CD014089")
        self.assertEqual(extract_review_id("custom_file.rda"), "custom_file")

    def test_make_analysis_key(self) -> None:
        frame = pd.DataFrame(
            {
                "Analysis.group": [1, 1],
                "Analysis.number": [2, 2],
                "Analysis.name": ["All-cause mortality", "All-cause mortality"],
                "Subgroup": ["overall", "overall"],
            }
        )
        out = make_analysis_key(frame)
        self.assertEqual(out.iloc[0], "1|2|All-cause mortality|overall")

    def test_calibrate_p_nonsig_from_lambda(self) -> None:
        # With no significant-positive studies, p_nonsig tracks lambda.
        self.assertAlmostEqual(calibrate_p_nonsig_from_lambda(0.7, q_sig=0.0, p_sig=0.95), 0.7, places=6)
        # Clipping behavior at lower bound.
        self.assertGreaterEqual(calibrate_p_nonsig_from_lambda(0.001, q_sig=0.2, p_sig=0.95), 0.01)
        # Never exceeds p_sig.
        self.assertLessEqual(calibrate_p_nonsig_from_lambda(0.99, q_sig=0.1, p_sig=0.95), 0.95)

    def test_apply_lambda_bounds(self) -> None:
        self.assertAlmostEqual(
            apply_lambda_bounds(0.05, lower=0.2, upper=0.95, clipping_enabled=True),
            0.2,
            places=6,
        )
        self.assertAlmostEqual(
            apply_lambda_bounds(1.2, lower=0.2, upper=0.95, clipping_enabled=False),
            1.2,
            places=6,
        )
        # lambda=1.0 (no ghost trials): should clamp to upper bound
        self.assertAlmostEqual(
            apply_lambda_bounds(1.0, lower=0.2, upper=0.95, clipping_enabled=True),
            0.95,
            places=6,
        )

    def test_build_lambda_posterior_params_count_beta(self) -> None:
        out = build_lambda_posterior_params(
            lambda_source="ctgov_linkage_non_ghost",
            lambda_center=0.7,
            review_link={"n_completed_non_ghost": 8, "n_completed_ghost": 2},
            prior_alpha=1.0,
            prior_beta=1.0,
            pseudo_concentration=20.0,
        )
        self.assertEqual(out["kind"], "count_beta")
        self.assertAlmostEqual(float(out["mean"]), 9.0 / 12.0, places=6)

    def test_build_review_cluster_summary(self) -> None:
        frame = pd.DataFrame(
            {
                "review_id": ["A", "A", "B", "B", "B"],
                "cross_prob_010": [1.0, 0.0, 0.5, 0.5, 0.0],
                "cross_prob_020": [1.0, 0.0, 1.0, 0.0, 0.0],
            }
        )
        summary, table = build_review_cluster_summary(frame, bootstrap_iters=20, bootstrap_seed=1)
        self.assertEqual(summary["n_reviews"], 2)
        self.assertEqual(int(summary["n_analyses"]), 5)
        self.assertEqual(len(table), 2)
        self.assertGreaterEqual(summary["crossed_below_0p10_pct_all_weighted"], 0.0)


if __name__ == "__main__":
    unittest.main()
