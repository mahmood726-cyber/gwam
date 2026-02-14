#!/usr/bin/env python3
"""Unit tests for CT.gov results extraction (extract_ctgov_results.py)."""

from __future__ import annotations

import math
import unittest
from unittest.mock import patch

# sys.path setup handled by conftest.py
from extract_ctgov_results import (
    extract_binary_from_measurements,
    extract_continuous_from_measurements,
    extract_from_analyses,
    extract_outcome_effect,
    safe_float,
    _match_outcome,
)
from gwam_utils import normal_quantile


class TestSafeFloat(unittest.TestCase):
    def test_valid_float(self) -> None:
        self.assertEqual(safe_float(3.14), 3.14)
        self.assertEqual(safe_float("2.5"), 2.5)
        self.assertEqual(safe_float(0), 0.0)

    def test_none_returns_nan(self) -> None:
        self.assertTrue(math.isnan(safe_float(None)))

    def test_invalid_string_returns_nan(self) -> None:
        self.assertTrue(math.isnan(safe_float("abc")))

    def test_inf_returns_nan(self) -> None:
        self.assertTrue(math.isnan(safe_float(float("inf"))))
        self.assertTrue(math.isnan(safe_float(float("-inf"))))


class TestMatchOutcome(unittest.TestCase):
    def test_empty_keywords_matches_all(self) -> None:
        self.assertTrue(_match_outcome("Anything", []))

    def test_keyword_match_case_insensitive(self) -> None:
        self.assertTrue(_match_outcome("Change in MADRS Score", ["madrs"]))
        self.assertTrue(_match_outcome("HAM-D Response Rate", ["ham-d", "response"]))

    def test_keyword_no_match(self) -> None:
        self.assertFalse(_match_outcome("Blood Pressure", ["madrs", "depression"]))


class TestExtractBinaryFromMeasurements(unittest.TestCase):
    """Test binary outcome extraction from 2x2 table."""

    def _make_binary_outcome(
        self,
        e1: float, n1: int, e2: float, n2: int,
    ) -> dict:
        """Create a mock outcome dict with binary measurements."""
        return {
            "groups": [
                {"id": "OG000", "title": "Treatment"},
                {"id": "OG001", "title": "Placebo"},
            ],
            "denoms": [
                {
                    "units": "Participants",
                    "counts": [
                        {"groupId": "OG000", "value": str(n1)},
                        {"groupId": "OG001", "value": str(n2)},
                    ],
                }
            ],
            "classes": [
                {
                    "categories": [
                        {
                            "measurements": [
                                {"groupId": "OG000", "value": str(e1)},
                                {"groupId": "OG001", "value": str(e2)},
                            ]
                        }
                    ]
                }
            ],
        }

    def test_basic_binary_extraction(self) -> None:
        """Standard 2x2: 30/100 vs 20/100 -> log(OR) without correction (no zero cells)."""
        outcome = self._make_binary_outcome(30, 100, 20, 100)
        result = extract_binary_from_measurements(outcome, [])
        self.assertIsNotNone(result)
        log_or, se = result

        # No zero cells -> no Haldane correction applied
        a, b, c, d = 30, 70, 20, 80
        expected_log_or = math.log((a * d) / (b * c))
        expected_se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)

        self.assertAlmostEqual(log_or, expected_log_or, places=6)
        self.assertAlmostEqual(se, expected_se, places=6)
        # Should be positive (higher events in treatment)
        self.assertGreater(log_or, 0)

    def test_zero_cell_correction(self) -> None:
        """Zero events in one arm -> Haldane correction prevents log(0)."""
        outcome = self._make_binary_outcome(10, 50, 0, 50)
        result = extract_binary_from_measurements(outcome, [])
        self.assertIsNotNone(result)
        log_or, se = result
        self.assertTrue(math.isfinite(log_or))
        self.assertTrue(math.isfinite(se))
        self.assertGreater(se, 0)

    def test_too_few_groups_returns_none(self) -> None:
        outcome = {
            "groups": [{"id": "OG000", "title": "Single arm"}],
            "denoms": [],
            "classes": [],
        }
        self.assertIsNone(extract_binary_from_measurements(outcome, []))

    def test_no_classes_returns_none(self) -> None:
        outcome = {
            "groups": [
                {"id": "OG000", "title": "A"},
                {"id": "OG001", "title": "B"},
            ],
            "denoms": [],
            "classes": [],
        }
        self.assertIsNone(extract_binary_from_measurements(outcome, []))

    def test_events_exceeding_denom_returns_none(self) -> None:
        """Events > denominator should be rejected."""
        outcome = self._make_binary_outcome(120, 100, 20, 100)
        self.assertIsNone(extract_binary_from_measurements(outcome, []))

    def test_double_all_events_returns_none(self) -> None:
        """Both arms 100% events (e1==n1, e2==n2) -> double-all-events exclusion."""
        outcome = self._make_binary_outcome(100, 100, 80, 80)
        self.assertIsNone(extract_binary_from_measurements(outcome, []))


class TestExtractContinuousFromMeasurements(unittest.TestCase):
    """Test continuous outcome extraction for SMD (Hedges' g)."""

    def _make_continuous_outcome(
        self,
        mean1: float, sd1: float, n1: int,
        mean2: float, sd2: float, n2: int,
        dispersion_type: str = "Standard Deviation",
    ) -> dict:
        """Create a mock outcome dict with continuous measurements.

        CT.gov API v2: dispersionType is at the outcome level,
        spread (numeric value) is at the measurement level.
        """
        return {
            "dispersionType": dispersion_type,
            "groups": [
                {"id": "OG000", "title": "Treatment"},
                {"id": "OG001", "title": "Control"},
            ],
            "denoms": [
                {
                    "units": "Participants",
                    "counts": [
                        {"groupId": "OG000", "value": str(n1)},
                        {"groupId": "OG001", "value": str(n2)},
                    ],
                }
            ],
            "classes": [
                {
                    "categories": [
                        {
                            "measurements": [
                                {
                                    "groupId": "OG000",
                                    "value": str(mean1),
                                    "spread": str(sd1),
                                },
                                {
                                    "groupId": "OG001",
                                    "value": str(mean2),
                                    "spread": str(sd2),
                                },
                            ]
                        }
                    ]
                }
            ],
        }

    def test_basic_continuous_extraction(self) -> None:
        """mean1=10, sd1=3, n1=50 vs mean2=8, sd2=3, n2=50 -> Hedges' g."""
        outcome = self._make_continuous_outcome(10.0, 3.0, 50, 8.0, 3.0, 50)
        result = extract_continuous_from_measurements(outcome)
        self.assertIsNotNone(result)
        smd, se = result

        # Manual: pooled SD = 3.0, raw Cohen's d = (10-8)/3 = 0.667
        # Hedges' g correction: j = 1 - 3/(4*(50+50-2)-1) = 1 - 3/391 ≈ 0.99233
        n1, n2 = 50, 50
        s_pooled = 3.0
        raw_d = (10.0 - 8.0) / s_pooled
        df = n1 + n2 - 2
        j = 1.0 - 3.0 / (4.0 * df - 1.0)
        expected_smd = raw_d * j
        # SE uses df=n1+n2-2 in denominator (not n1+n2)
        expected_se = math.sqrt(
            (n1 + n2) / (n1 * n2) + expected_smd ** 2 / (2.0 * df)
        )

        self.assertAlmostEqual(smd, expected_smd, places=4)
        self.assertAlmostEqual(se, expected_se, places=4)
        self.assertGreater(smd, 0)
        self.assertGreater(se, 0)

    def test_equal_means_gives_zero_smd(self) -> None:
        """Same means -> SMD should be ~0."""
        outcome = self._make_continuous_outcome(5.0, 2.0, 30, 5.0, 2.0, 30)
        result = extract_continuous_from_measurements(outcome)
        self.assertIsNotNone(result)
        smd, se = result
        self.assertAlmostEqual(smd, 0.0, places=6)
        self.assertGreater(se, 0)

    def test_small_n_returns_none(self) -> None:
        """n=1 in either arm should be rejected (n<=1)."""
        outcome = self._make_continuous_outcome(10.0, 3.0, 1, 8.0, 3.0, 50)
        self.assertIsNone(extract_continuous_from_measurements(outcome))

    def test_zero_sd_returns_none(self) -> None:
        """SD=0 should be rejected."""
        outcome = self._make_continuous_outcome(10.0, 0.0, 50, 8.0, 3.0, 50)
        self.assertIsNone(extract_continuous_from_measurements(outcome))


    def test_standard_error_dispersion_converts_to_sd(self) -> None:
        """When dispersionType is Standard Error, spread values should be converted to SD."""
        # SE=1.0 with n=25 -> SD = 1.0 * sqrt(25) = 5.0
        outcome = self._make_continuous_outcome(
            10.0, 1.0, 25, 8.0, 1.0, 25,
            dispersion_type="Standard Error of the Mean",
        )
        result = extract_continuous_from_measurements(outcome)
        self.assertIsNotNone(result)
        smd, se = result
        # After SE->SD conversion: sd1 = 1.0*sqrt(25)=5.0, sd2 = 1.0*sqrt(25)=5.0
        # s_pooled = 5.0, raw d = 2.0/5.0 = 0.4
        n1, n2 = 25, 25
        sd_converted = 1.0 * math.sqrt(25)
        df = n1 + n2 - 2
        raw_d = (10.0 - 8.0) / sd_converted
        j = 1.0 - 3.0 / (4.0 * df - 1.0)
        expected_smd = raw_d * j
        self.assertAlmostEqual(smd, expected_smd, places=4)

    def test_missing_dispersion_type_returns_none(self) -> None:
        """Outcome without dispersionType should return None."""
        outcome = {
            "groups": [
                {"id": "OG000", "title": "Treatment"},
                {"id": "OG001", "title": "Control"},
            ],
            "denoms": [
                {
                    "counts": [
                        {"groupId": "OG000", "value": "50"},
                        {"groupId": "OG001", "value": "50"},
                    ],
                }
            ],
            "classes": [
                {
                    "categories": [
                        {
                            "measurements": [
                                {"groupId": "OG000", "value": "10.0", "spread": "3.0"},
                                {"groupId": "OG001", "value": "8.0", "spread": "3.0"},
                            ]
                        }
                    ]
                }
            ],
        }
        self.assertIsNone(extract_continuous_from_measurements(outcome))


class TestExtractFromAnalyses(unittest.TestCase):
    """Test extraction from the analyses section."""

    def test_odds_ratio_extraction(self) -> None:
        """OR=2.0, CI=[1.2, 3.33] -> log scale."""
        outcome = {
            "analyses": [
                {
                    "paramType": "Odds Ratio (Non-Inferiority or Superiority)",
                    "paramValue": "2.0",
                    "ciLowerLimit": "1.2",
                    "ciUpperLimit": "3.33",
                    "ciPctValue": "95",
                }
            ]
        }
        result = extract_from_analyses(outcome)
        self.assertIsNotNone(result)
        effect, se, effect_type = result

        self.assertAlmostEqual(effect, math.log(2.0), places=6)
        z = normal_quantile(0.975)
        expected_se = (math.log(3.33) - math.log(1.2)) / (2 * z)
        self.assertAlmostEqual(se, expected_se, places=4)
        self.assertIn("odds_ratio", effect_type)

    def test_risk_ratio_extraction(self) -> None:
        outcome = {
            "analyses": [
                {
                    "paramType": "Risk Ratio (RR)",
                    "paramValue": "0.75",
                    "ciLowerLimit": "0.55",
                    "ciUpperLimit": "1.02",
                }
            ]
        }
        result = extract_from_analyses(outcome)
        self.assertIsNotNone(result)
        effect, se, effect_type = result
        self.assertAlmostEqual(effect, math.log(0.75), places=6)
        self.assertIn("risk_ratio", effect_type)

    def test_hazard_ratio_extraction(self) -> None:
        outcome = {
            "analyses": [
                {
                    "paramType": "Hazard Ratio (HR)",
                    "paramValue": "0.80",
                    "ciLowerLimit": "0.60",
                    "ciUpperLimit": "1.10",
                }
            ]
        }
        result = extract_from_analyses(outcome)
        self.assertIsNotNone(result)
        effect, se, effect_type = result
        self.assertAlmostEqual(effect, math.log(0.80), places=6)
        self.assertIn("hazard_ratio", effect_type)

    def test_mean_difference_extraction(self) -> None:
        outcome = {
            "analyses": [
                {
                    "paramType": "Mean Difference (Net)",
                    "paramValue": "-2.5",
                    "ciLowerLimit": "-4.0",
                    "ciUpperLimit": "-1.0",
                }
            ]
        }
        result = extract_from_analyses(outcome)
        self.assertIsNotNone(result)
        effect, se, effect_type = result
        self.assertAlmostEqual(effect, -2.5, places=6)
        z = normal_quantile(0.975)
        expected_se = (-1.0 - (-4.0)) / (2 * z)
        self.assertAlmostEqual(se, expected_se, places=4)
        self.assertIn("difference", effect_type)

    def test_missing_ci_returns_none(self) -> None:
        outcome = {
            "analyses": [
                {
                    "paramType": "Odds Ratio",
                    "paramValue": "2.0",
                    # No CI limits
                }
            ]
        }
        self.assertIsNone(extract_from_analyses(outcome))

    def test_zero_or_negative_ratio_returns_none(self) -> None:
        outcome = {
            "analyses": [
                {
                    "paramType": "Odds Ratio",
                    "paramValue": "-1.0",
                    "ciLowerLimit": "0.5",
                    "ciUpperLimit": "2.0",
                }
            ]
        }
        self.assertIsNone(extract_from_analyses(outcome))

    def test_90_percent_ci_uses_correct_z(self) -> None:
        """90% CI should use z=1.6449, not 1.96."""
        outcome = {
            "analyses": [
                {
                    "paramType": "Mean Difference",
                    "paramValue": "-3.0",
                    "ciLowerLimit": "-5.0",
                    "ciUpperLimit": "-1.0",
                    "ciPctValue": "90",
                }
            ]
        }
        result = extract_from_analyses(outcome)
        self.assertIsNotNone(result)
        effect, se, _ = result
        # z for 90% two-sided: alpha/2 = 0.05, z = normal_quantile(0.95)
        z90 = normal_quantile(0.95)
        expected_se = (-1.0 - (-5.0)) / (2 * z90)
        self.assertAlmostEqual(se, expected_se, places=3)

    def test_ci_bracket_violation_skipped(self) -> None:
        """If param_value is outside [ci_lower, ci_upper], should skip."""
        outcome = {
            "analyses": [
                {
                    "paramType": "Odds Ratio",
                    "paramValue": "5.0",  # outside CI [1.2, 3.0]
                    "ciLowerLimit": "1.2",
                    "ciUpperLimit": "3.0",
                }
            ]
        }
        self.assertIsNone(extract_from_analyses(outcome))

    def test_empty_analyses_returns_none(self) -> None:
        self.assertIsNone(extract_from_analyses({"analyses": []}))
        self.assertIsNone(extract_from_analyses({}))


class TestExtractOutcomeEffect(unittest.TestCase):
    """Test the top-level extraction priority logic."""

    def test_analyses_preferred_over_measurements(self) -> None:
        """Analyses section should be tried first (priority 1)."""
        payload = {
            "resultsSection": {
                "outcomeMeasuresModule": {
                    "outcomeMeasures": [
                        {
                            "type": "PRIMARY",
                            "title": "Response Rate",
                            "analyses": [
                                {
                                    "paramType": "Odds Ratio",
                                    "paramValue": "1.5",
                                    "ciLowerLimit": "0.9",
                                    "ciUpperLimit": "2.5",
                                }
                            ],
                            "groups": [
                                {"id": "OG000", "title": "Treatment"},
                                {"id": "OG001", "title": "Placebo"},
                            ],
                            "denoms": [
                                {
                                    "counts": [
                                        {"groupId": "OG000", "value": "100"},
                                        {"groupId": "OG001", "value": "100"},
                                    ]
                                }
                            ],
                            "classes": [
                                {
                                    "categories": [
                                        {
                                            "measurements": [
                                                {"groupId": "OG000", "value": "40"},
                                                {"groupId": "OG001", "value": "30"},
                                            ]
                                        }
                                    ]
                                }
                            ],
                        }
                    ]
                }
            },
            "protocolSection": {},
        }
        result = extract_outcome_effect(payload, [])
        self.assertIsNotNone(result)
        _, _, effect_type = result
        # Should come from analyses, not from 2x2 table
        self.assertIn("odds_ratio", effect_type)
        self.assertNotEqual(effect_type, "log_or_2x2")

    def test_binary_fallback_when_no_analyses(self) -> None:
        """Without analyses, should fall back to binary 2x2 extraction."""
        payload = {
            "resultsSection": {
                "outcomeMeasuresModule": {
                    "outcomeMeasures": [
                        {
                            "type": "PRIMARY",
                            "title": "Responders",
                            "groups": [
                                {"id": "OG000", "title": "Drug"},
                                {"id": "OG001", "title": "Placebo"},
                            ],
                            "denoms": [
                                {
                                    "counts": [
                                        {"groupId": "OG000", "value": "80"},
                                        {"groupId": "OG001", "value": "80"},
                                    ]
                                }
                            ],
                            "classes": [
                                {
                                    "categories": [
                                        {
                                            "measurements": [
                                                {"groupId": "OG000", "value": "25"},
                                                {"groupId": "OG001", "value": "15"},
                                            ]
                                        }
                                    ]
                                }
                            ],
                        }
                    ]
                }
            },
            "protocolSection": {},
        }
        result = extract_outcome_effect(payload, [])
        self.assertIsNotNone(result)
        _, _, effect_type = result
        self.assertEqual(effect_type, "log_or_2x2")

    def test_continuous_fallback(self) -> None:
        """When no binary data available, should fall back to SMD.

        Use negative change scores so binary extraction fails (values < 0).
        dispersionType is at the outcome level per CT.gov API v2.
        """
        payload = {
            "resultsSection": {
                "outcomeMeasuresModule": {
                    "outcomeMeasures": [
                        {
                            "type": "PRIMARY",
                            "title": "Change in Score",
                            "dispersionType": "Standard Deviation",
                            "groups": [
                                {"id": "OG000", "title": "Drug"},
                                {"id": "OG001", "title": "Placebo"},
                            ],
                            "denoms": [
                                {
                                    "counts": [
                                        {"groupId": "OG000", "value": "50"},
                                        {"groupId": "OG001", "value": "50"},
                                    ]
                                }
                            ],
                            "classes": [
                                {
                                    "categories": [
                                        {
                                            "measurements": [
                                                {
                                                    "groupId": "OG000",
                                                    "value": "-8.5",
                                                    "spread": "4.0",
                                                },
                                                {
                                                    "groupId": "OG001",
                                                    "value": "-5.2",
                                                    "spread": "4.0",
                                                },
                                            ]
                                        }
                                    ]
                                }
                            ],
                        }
                    ]
                }
            },
            "protocolSection": {},
        }
        result = extract_outcome_effect(payload, [])
        self.assertIsNotNone(result)
        _, _, effect_type = result
        self.assertEqual(effect_type, "smd_hedges_g")

    def test_keyword_filtering(self) -> None:
        """With keywords, should only match relevant outcomes."""
        payload = {
            "resultsSection": {
                "outcomeMeasuresModule": {
                    "outcomeMeasures": [
                        {
                            "type": "PRIMARY",
                            "title": "Blood Pressure",
                            "analyses": [
                                {
                                    "paramType": "Mean Difference",
                                    "paramValue": "-5.0",
                                    "ciLowerLimit": "-8.0",
                                    "ciUpperLimit": "-2.0",
                                }
                            ],
                        },
                        {
                            "type": "PRIMARY",
                            "title": "MADRS Score Change",
                            "analyses": [
                                {
                                    "paramType": "Mean Difference",
                                    "paramValue": "-3.0",
                                    "ciLowerLimit": "-5.0",
                                    "ciUpperLimit": "-1.0",
                                }
                            ],
                        },
                    ]
                }
            },
            "protocolSection": {},
        }
        # Without keywords: should get first outcome
        result_any = extract_outcome_effect(payload, [])
        self.assertIsNotNone(result_any)
        effect_any, _, _ = result_any
        self.assertAlmostEqual(effect_any, -5.0, places=4)

        # With keyword "MADRS": should get second outcome
        result_madrs = extract_outcome_effect(payload, ["madrs"])
        self.assertIsNotNone(result_madrs)
        effect_madrs, _, _ = result_madrs
        self.assertAlmostEqual(effect_madrs, -3.0, places=4)

    def test_no_outcomes_returns_none(self) -> None:
        payload = {
            "resultsSection": {
                "outcomeMeasuresModule": {"outcomeMeasures": []}
            },
            "protocolSection": {},
        }
        self.assertIsNone(extract_outcome_effect(payload, []))

    def test_no_results_section_returns_none(self) -> None:
        self.assertIsNone(extract_outcome_effect({}, []))
        self.assertIsNone(extract_outcome_effect({"resultsSection": {}}, []))


if __name__ == "__main__":
    unittest.main()
