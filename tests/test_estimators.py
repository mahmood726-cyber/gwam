#!/usr/bin/env python3
"""Unit tests for random_effects_dl() and pet_peese_estimate()."""

from __future__ import annotations

import math
import unittest

import numpy as np

from simulate_gwam_vs_re import pet_peese_estimate, random_effects_dl


class TestRandomEffectsDL(unittest.TestCase):
    """Tests for DerSimonian-Laird random-effects estimator."""

    def test_homogeneous_studies(self) -> None:
        """Identical effects → mu equals the common value, tau2 = 0."""
        y = np.array([0.5, 0.5, 0.5, 0.5])
        v = np.array([0.04, 0.04, 0.04, 0.04])
        result = random_effects_dl(y, v)
        self.assertIsNotNone(result)
        mu, se, tau2 = result  # type: ignore[misc]
        self.assertAlmostEqual(mu, 0.5, places=6)
        self.assertAlmostEqual(tau2, 0.0, places=6)
        self.assertGreater(se, 0.0)

    def test_known_heterogeneous(self) -> None:
        """Effects with spread → tau2 > 0."""
        y = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
        v = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        result = random_effects_dl(y, v)
        self.assertIsNotNone(result)
        mu, se, tau2 = result  # type: ignore[misc]
        self.assertGreater(tau2, 0.0)
        # Mean should be near 1.0 (midpoint)
        self.assertAlmostEqual(mu, 1.0, places=2)

    def test_single_study_returns_none(self) -> None:
        """k=1 → None (cannot estimate heterogeneity)."""
        result = random_effects_dl(np.array([0.3]), np.array([0.05]))
        self.assertIsNone(result)

    def test_two_studies(self) -> None:
        """k=2 → valid result (minimum for heterogeneity)."""
        y = np.array([0.2, 0.8])
        v = np.array([0.05, 0.05])
        result = random_effects_dl(y, v)
        self.assertIsNotNone(result)
        mu, se, tau2 = result  # type: ignore[misc]
        self.assertTrue(math.isfinite(mu))
        self.assertGreater(se, 0.0)

    def test_weight_multiplier(self) -> None:
        """IPW multiplier changes the estimate."""
        y = np.array([0.1, 0.3, 0.5, 0.7])
        v = np.array([0.04, 0.04, 0.04, 0.04])
        result_plain = random_effects_dl(y, v)
        # Up-weight the last two studies
        mult = np.array([1.0, 1.0, 3.0, 3.0])
        result_ipw = random_effects_dl(y, v, weight_multiplier=mult)
        self.assertIsNotNone(result_plain)
        self.assertIsNotNone(result_ipw)
        # IPW estimate should be pulled toward 0.5/0.7 (heavier weights)
        self.assertGreater(result_ipw[0], result_plain[0])  # type: ignore[index]

    def test_zero_variance_returns_none(self) -> None:
        """Zero variance → None (invalid input)."""
        result = random_effects_dl(np.array([0.1, 0.2, 0.3]), np.array([0.0, 0.1, 0.1]))
        self.assertIsNone(result)


class TestPetPeeseEstimate(unittest.TestCase):
    """Tests for the PET-PEESE publication-bias estimator."""

    def test_fewer_than_three_returns_none(self) -> None:
        """k < 3 → None (insufficient for WLS regression)."""
        self.assertIsNone(pet_peese_estimate(np.array([0.1, 0.2]), np.array([0.01, 0.01])))

    def test_zero_variance_returns_none(self) -> None:
        """Any v <= 0 → None."""
        self.assertIsNone(
            pet_peese_estimate(np.array([0.1, 0.2, 0.3]), np.array([0.01, 0.0, 0.01]))
        )

    def test_returns_mu_and_se(self) -> None:
        """Valid input → (mu, se) tuple."""
        rng = np.random.default_rng(42)
        k = 20
        true_mu = 0.3
        y = rng.normal(true_mu, 0.2, size=k)
        v = rng.uniform(0.01, 0.1, size=k)
        result = pet_peese_estimate(y, v)
        self.assertIsNotNone(result)
        mu, se = result  # type: ignore[misc]
        self.assertTrue(math.isfinite(mu))
        self.assertGreater(se, 0.0)

    def test_no_bias_returns_near_true(self) -> None:
        """No funnel asymmetry → estimate near true value."""
        rng = np.random.default_rng(123)
        k = 50
        true_mu = 0.5
        # Varying variances needed so SE regression is non-degenerate
        v = rng.uniform(0.01, 0.10, size=k)
        y = rng.normal(true_mu, np.sqrt(v))
        result = pet_peese_estimate(y, v)
        self.assertIsNotNone(result)
        mu, _ = result  # type: ignore[misc]
        # Should be within ~0.4 of true value
        self.assertAlmostEqual(mu, true_mu, delta=0.4)


if __name__ == "__main__":
    unittest.main()
