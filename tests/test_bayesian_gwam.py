#!/usr/bin/env python3
"""Unit tests for Bayesian hierarchical GWAM model."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


# sys.path setup handled by conftest.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

from gwam_utils import normal_cdf  # noqa: E402
from model_gwam_bayesian import (  # noqa: E402
    BayesianResult,
    PriorSpec,
    bayesian_gwam_posterior,
    run_sensitivity_grid,
)


class TestBayesianGWAMPosterior(unittest.TestCase):
    """Test the analytical Bayesian GWAM posterior computation."""

    def test_converges_to_deterministic_when_priors_vanish(self) -> None:
        """When sigma_G -> 0 and s_mG -> 0, posterior should match deterministic GWAM."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=1e-12,  # near-zero hyperprior uncertainty
            ghost_sigma_grid=[1e-12],  # near-zero within-ghost heterogeneity
            results_only_mu_mean=0.0,
            results_only_mu_sd=1e-12,
            results_only_sigma_grid=[1e-12],
        )
        w_pub = 100.0
        mu_pub = 0.5
        se_pub = 0.1
        w_ro = np.array([50.0])
        w_g = np.array([50.0])
        w_total = 200.0

        result = bayesian_gwam_posterior(
            w_pub=w_pub,
            mu_pub=mu_pub,
            se_pub=se_pub,
            w_results_only_individual=w_ro,
            w_ghost_individual=w_g,
            w_total=w_total,
            prior=prior,
            ghost_sigma=1e-12,
            results_only_sigma=1e-12,
        )

        # Deterministic: mu_gwam = (100 * 0.5 + 50 * 0 + 50 * 0) / 200 = 0.25
        expected_mean = (w_pub * mu_pub) / w_total
        self.assertAlmostEqual(result.posterior_mean, expected_mean, places=6)

        # With near-zero ghost/results-only uncertainty, posterior SD should
        # be dominated by published SE: w_pub * se_pub / w_total = 0.05
        expected_sd = (w_pub * se_pub) / w_total
        self.assertAlmostEqual(result.posterior_sd, expected_sd, places=4)

    def test_posterior_sd_increases_with_hyperprior_uncertainty(self) -> None:
        """Larger hyperprior SD (s_mG) should widen the posterior."""
        w_pub = 100.0
        mu_pub = 0.5
        se_pub = 0.1
        w_ro = np.array([50.0])
        w_g = np.array([50.0])
        w_total = 200.0

        sds: list[float] = []
        for s_mg in [0.0, 0.1, 0.2, 0.5]:
            prior = PriorSpec(
                ghost_mu_mean=0.0,
                ghost_mu_sd=s_mg,
                ghost_sigma_grid=[0.1],
                results_only_mu_mean=0.0,
                results_only_mu_sd=0.1,
                results_only_sigma_grid=[0.1],
            )
            r = bayesian_gwam_posterior(
                w_pub=w_pub,
                mu_pub=mu_pub,
                se_pub=se_pub,
                w_results_only_individual=w_ro,
                w_ghost_individual=w_g,
                w_total=w_total,
                prior=prior,
                ghost_sigma=0.1,
                results_only_sigma=0.1,
            )
            sds.append(r.posterior_sd)

        # SD should be monotonically non-decreasing with hyperprior uncertainty
        for i in range(1, len(sds)):
            self.assertGreaterEqual(sds[i], sds[i - 1] - 1e-12)

    def test_posterior_sd_increases_with_ghost_sigma(self) -> None:
        """Larger within-ghost sigma should widen the posterior."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
        )
        w_pub = 100.0
        mu_pub = 0.5
        se_pub = 0.1
        w_ro = np.array([50.0])
        w_g = np.array([50.0, 30.0])
        w_total = 230.0

        sds: list[float] = []
        for gs in [0.0, 0.1, 0.2, 0.5]:
            r = bayesian_gwam_posterior(
                w_pub=w_pub,
                mu_pub=mu_pub,
                se_pub=se_pub,
                w_results_only_individual=w_ro,
                w_ghost_individual=w_g,
                w_total=w_total,
                prior=prior,
                ghost_sigma=gs,
                results_only_sigma=0.1,
            )
            sds.append(r.posterior_sd)

        for i in range(1, len(sds)):
            self.assertGreaterEqual(sds[i], sds[i - 1] - 1e-12)

    def test_cri_width_sensitivity_to_prior(self) -> None:
        """CrI width should vary meaningfully across sensitivity grid."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.05, 0.3],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.05, 0.3],
        )
        w_pub = 100.0
        mu_pub = 0.3
        se_pub = 0.08
        w_ro = np.array([40.0, 60.0])
        w_g = np.array([30.0, 20.0, 50.0])
        w_total = 300.0

        grid = run_sensitivity_grid(
            w_pub=w_pub,
            mu_pub=mu_pub,
            se_pub=se_pub,
            w_results_only_individual=w_ro,
            w_ghost_individual=w_g,
            w_total=w_total,
            prior=prior,
        )

        widths = [r.cri_hi - r.cri_lo for r in grid]
        self.assertEqual(len(grid), 4)  # 2 ghost * 2 ro
        self.assertGreater(max(widths) - min(widths), 0.01)

    def test_no_ghost_studies(self) -> None:
        """With no ghost studies, posterior should be close to published with prior on RO."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.1],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.1],
        )
        w_pub = 200.0
        mu_pub = 0.4
        se_pub = 0.08
        w_ro = np.array([50.0])
        w_g = np.array([], dtype=float)  # no ghosts
        w_total = 250.0

        r = bayesian_gwam_posterior(
            w_pub=w_pub,
            mu_pub=mu_pub,
            se_pub=se_pub,
            w_results_only_individual=w_ro,
            w_ghost_individual=w_g,
            w_total=w_total,
            prior=prior,
            ghost_sigma=0.1,
            results_only_sigma=0.1,
        )

        # Mean should be between published_mu and 0
        self.assertGreater(r.posterior_mean, 0.0)
        self.assertLessEqual(r.posterior_mean, mu_pub)
        # CrI should contain the mean
        self.assertLess(r.cri_lo, r.posterior_mean)
        self.assertGreater(r.cri_hi, r.posterior_mean)

    def test_observed_results_reduce_uncertainty(self) -> None:
        """When results-only studies have observed effects, posterior should tighten."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.1],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.1],
        )
        w_pub = 100.0
        mu_pub = 0.5
        se_pub = 0.1
        w_ro = np.array([50.0, 40.0, 30.0])
        w_g = np.array([50.0])
        w_total = 270.0

        # Without observed effects
        r_unobs = bayesian_gwam_posterior(
            w_pub=w_pub,
            mu_pub=mu_pub,
            se_pub=se_pub,
            w_results_only_individual=w_ro,
            w_ghost_individual=w_g,
            w_total=w_total,
            prior=prior,
            ghost_sigma=0.1,
            results_only_sigma=0.1,
        )

        # With observed effects for 2 of 3 results-only
        obs = np.array([0.2, 0.1, float("nan")])
        r_obs = bayesian_gwam_posterior(
            w_pub=w_pub,
            mu_pub=mu_pub,
            se_pub=se_pub,
            w_results_only_individual=w_ro,
            w_ghost_individual=w_g,
            w_total=w_total,
            prior=prior,
            ghost_sigma=0.1,
            results_only_sigma=0.1,
            results_only_observed=obs,
        )

        # Observed should have tighter posterior (less uncertainty)
        self.assertLess(r_obs.posterior_sd, r_unobs.posterior_sd)

    def test_pr_positive_correct(self) -> None:
        """Pr(positive) should be consistent with mean and SD."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.1],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.1],
        )
        r = bayesian_gwam_posterior(
            w_pub=100.0,
            mu_pub=0.5,
            se_pub=0.1,
            w_results_only_individual=np.array([50.0]),
            w_ghost_individual=np.array([50.0]),
            w_total=200.0,
            prior=prior,
            ghost_sigma=0.1,
            results_only_sigma=0.1,
        )

        expected_pr = normal_cdf(r.posterior_mean / r.posterior_sd)
        self.assertAlmostEqual(r.pr_positive, expected_pr, places=8)

    def test_observed_length_mismatch_raises(self) -> None:
        """Mismatched results_only_observed length should raise ValueError."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.1],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.1],
        )
        with self.assertRaises(ValueError):
            bayesian_gwam_posterior(
                w_pub=100.0,
                mu_pub=0.5,
                se_pub=0.1,
                w_results_only_individual=np.array([50.0, 40.0]),
                w_ghost_individual=np.array([50.0]),
                w_total=240.0,
                prior=prior,
                ghost_sigma=0.1,
                results_only_sigma=0.1,
                results_only_observed=np.array([0.2]),  # length 1 vs 2
            )

    def test_posterior_mean_is_zero_when_all_zero(self) -> None:
        """When published_mu = 0 and prior means = 0, posterior mean should be 0."""
        prior = PriorSpec(
            ghost_mu_mean=0.0,
            ghost_mu_sd=0.15,
            ghost_sigma_grid=[0.1],
            results_only_mu_mean=0.0,
            results_only_mu_sd=0.15,
            results_only_sigma_grid=[0.1],
        )
        r = bayesian_gwam_posterior(
            w_pub=100.0,
            mu_pub=0.0,
            se_pub=0.1,
            w_results_only_individual=np.array([50.0]),
            w_ghost_individual=np.array([50.0]),
            w_total=200.0,
            prior=prior,
            ghost_sigma=0.1,
            results_only_sigma=0.1,
        )

        self.assertAlmostEqual(r.posterior_mean, 0.0, places=10)
        self.assertAlmostEqual(r.pr_positive, 0.5, places=3)


class TestBayesianGWAMScript(unittest.TestCase):
    """Integration tests for the CLI script."""

    def _make_registry_csv(self, tmpdir: Path) -> Path:
        csv_path = tmpdir / "registry.csv"
        rows = [
            {"is_ghost_protocol": "false", "has_pmid": "true", "has_results": "true", "weight_n": "100"},
            {"is_ghost_protocol": "false", "has_pmid": "false", "has_results": "true", "weight_n": "80"},
            {"is_ghost_protocol": "true", "has_pmid": "false", "has_results": "false", "weight_n": "60"},
            {"is_ghost_protocol": "true", "has_pmid": "false", "has_results": "false", "weight_n": "40"},
        ]
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"])
            writer.writeheader()
            writer.writerows(rows)
        return csv_path

    def test_script_produces_valid_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = self._make_registry_csv(tmp)
            out_path = tmp / "bayesian.json"

            cmd = [
                sys.executable,
                str(SCRIPTS_DIR / "model_gwam_bayesian.py"),
                "--registry-csv", str(csv_path),
                "--published-mu", "0.25",
                "--published-se", "0.10",
                "--prior-ghost-sigma-grid", "0.1,0.2",
                "--prior-ro-sigma-grid", "0.1,0.2",
                "--output-json", str(out_path),
            ]
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)

            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertIn("bayesian_gwam", payload)
            self.assertIn("sensitivity_table", payload)
            self.assertEqual(len(payload["sensitivity_table"]), 4)  # 2 x 2 grid
            self.assertEqual(payload["n_trials_total"], 4)
            self.assertEqual(payload["n_trials_with_pmid"], 1)
            self.assertEqual(payload["n_trials_results_only"], 1)
            self.assertEqual(payload["n_trials_ghost"], 2)

            # Bayesian GWAM mean should be between 0 and published_mu
            ref = payload["bayesian_gwam"]
            self.assertGreater(ref["reference_posterior_mean"], -1.0)
            self.assertLess(ref["reference_posterior_mean"], 0.5)
            self.assertGreater(ref["reference_posterior_sd"], 0.0)
            self.assertLess(ref["reference_cri_lo"], ref["reference_cri_hi"])

    def test_lambda_matches_model_gwam(self) -> None:
        """Lambda values should match between bayesian and deterministic models."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = self._make_registry_csv(tmp)
            out_path = tmp / "bayesian.json"

            cmd = [
                sys.executable,
                str(SCRIPTS_DIR / "model_gwam_bayesian.py"),
                "--registry-csv", str(csv_path),
                "--published-mu", "0.25",
                "--published-se", "0.10",
                "--prior-ghost-sigma-grid", "0.1",
                "--prior-ro-sigma-grid", "0.1",
                "--output-json", str(out_path),
            ]
            subprocess.run(cmd, check=True)

            payload = json.loads(out_path.read_text(encoding="utf-8"))
            # w_pub=100, w_ro=80, w_g=100, w_total=280
            expected_lambda_pmid = 100.0 / 280.0
            expected_lambda_non_ghost = 180.0 / 280.0
            self.assertAlmostEqual(
                payload["weights"]["integrity_ratio_lambda_pmid_only"],
                expected_lambda_pmid,
                places=6,
            )
            self.assertAlmostEqual(
                payload["weights"]["integrity_ratio_lambda_non_ghost"],
                expected_lambda_non_ghost,
                places=6,
            )


if __name__ == "__main__":
    unittest.main()
