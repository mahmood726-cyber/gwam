#!/usr/bin/env python3
"""Unit tests for GWAM core logic."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


# sys.path setup handled by conftest.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

from gwam_utils import (  # noqa: E402
    classify_publication_status,
    extract_pmids,
    intervention_matches,
    trial_passes_design_filters,
)


class TestGWAMUtils(unittest.TestCase):
    def test_intervention_matches_alias(self) -> None:
        names = ["Lexapro", "Placebo"]
        self.assertTrue(intervention_matches(names, intervention="escitalopram"))
        self.assertFalse(intervention_matches(["citalopram"], intervention="escitalopram"))

    def test_publication_classification(self) -> None:
        self.assertEqual(
            classify_publication_status(
                has_pmid=False, has_results=True, ghost_definition="no_pmid_no_results"
            ),
            (True, False),
        )
        self.assertEqual(
            classify_publication_status(
                has_pmid=False, has_results=False, ghost_definition="no_pmid_no_results"
            ),
            (False, True),
        )
        self.assertEqual(
            classify_publication_status(
                has_pmid=False,
                has_results=True,
                ghost_definition="no_pmid",
            ),
            (True, True),
        )

    def test_design_filters_active(self) -> None:
        protocol = {
            "designModule": {
                "studyType": "INTERVENTIONAL",
                "designInfo": {"allocation": "RANDOMIZED", "primaryPurpose": "TREATMENT"},
            },
            "armsInterventionsModule": {
                "armGroups": [
                    {"label": "Escitalopram", "type": "EXPERIMENTAL"},
                    {"label": "Sertraline", "type": "ACTIVE_COMPARATOR"},
                ],
                "interventions": [
                    {"name": "Escitalopram", "type": "DRUG"},
                    {"name": "Sertraline", "type": "DRUG"},
                ],
            },
        }
        self.assertTrue(
            trial_passes_design_filters(
                protocol,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
            )
        )

        protocol_placebo = {
            **protocol,
            "armsInterventionsModule": {
                "armGroups": [{"label": "Escitalopram"}, {"label": "Placebo"}],
                "interventions": [
                    {"name": "Escitalopram", "type": "DRUG"},
                    {"name": "Placebo", "type": "PLACEBO"},
                ],
            },
        }
        self.assertFalse(
            trial_passes_design_filters(
                protocol_placebo,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
            )
        )

    def test_extract_pmids_respects_linkage_rule(self) -> None:
        refs = [
            {"pmid": "100", "type": "BACKGROUND"},
            {"pmid": "200", "type": "RESULT"},
            {"pmid": "300", "type": "DERIVED"},
        ]
        _, pm_res_strict, has_results_link_strict = extract_pmids(refs, publication_linkage="results_pmid_strict")
        self.assertEqual(pm_res_strict, ["200"])
        self.assertTrue(has_results_link_strict)

        pm_any, pm_res, has_results_link = extract_pmids(refs, publication_linkage="results_pmid")
        self.assertEqual(pm_any, ["100", "200", "300"])
        self.assertEqual(pm_res, ["200", "300"])
        self.assertTrue(has_results_link)

        _, _, has_any_link = extract_pmids(refs, publication_linkage="any_pmid")
        self.assertTrue(has_any_link)

        _, _, has_results_only = extract_pmids([{"pmid": "100", "type": "BACKGROUND"}], publication_linkage="results_pmid")
        self.assertFalse(has_results_only)

    def test_publication_classification_any_pmid_guard(self) -> None:
        self.assertEqual(
            classify_publication_status(
                has_pmid=False,
                has_results=False,
                has_any_pmid=True,
                ghost_require_no_any_pmid=True,
                ghost_definition="no_pmid_no_results",
            ),
            (False, False),
        )
        self.assertEqual(
            classify_publication_status(
                has_pmid=False,
                has_results=False,
                has_any_pmid=True,
                ghost_require_no_any_pmid=False,
                ghost_definition="no_pmid_no_results",
            ),
            (False, True),
        )

    def test_missing_outcome_policy_behavior(self) -> None:
        protocol = {
            "designModule": {
                "studyType": "INTERVENTIONAL",
                "designInfo": {"allocation": "RANDOMIZED", "primaryPurpose": "TREATMENT"},
            },
            "armsInterventionsModule": {
                "armGroups": [
                    {"label": "Escitalopram", "type": "EXPERIMENTAL"},
                    {"label": "Sertraline", "type": "ACTIVE_COMPARATOR"},
                ],
                "interventions": [
                    {"name": "Escitalopram", "type": "DRUG"},
                    {"name": "Sertraline", "type": "DRUG"},
                ],
            },
            # outcomesModule intentionally missing
        }
        self.assertFalse(
            trial_passes_design_filters(
                protocol,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
                required_outcome_keywords=["response"],
                missing_outcome_policy="exclude",
            )
        )
        self.assertTrue(
            trial_passes_design_filters(
                protocol,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
                required_outcome_keywords=["response"],
                missing_outcome_policy="include_as_unknown",
            )
        )

    def test_design_filters_estimand_keyword_and_timeframe(self) -> None:
        protocol = {
            "designModule": {
                "studyType": "INTERVENTIONAL",
                "designInfo": {"allocation": "RANDOMIZED", "primaryPurpose": "TREATMENT"},
            },
            "armsInterventionsModule": {
                "armGroups": [
                    {"label": "Escitalopram", "type": "EXPERIMENTAL"},
                    {"label": "Sertraline", "type": "ACTIVE_COMPARATOR"},
                ],
                "interventions": [
                    {"name": "Escitalopram", "type": "DRUG"},
                    {"name": "Sertraline", "type": "DRUG"},
                ],
            },
            "outcomesModule": {
                "primaryOutcomes": [
                    {
                        "measure": "Response rate on MADRS",
                        "description": "Clinical response",
                        "timeFrame": "Baseline, Week 8 and Week 24",
                    }
                ]
            },
        }
        self.assertTrue(
            trial_passes_design_filters(
                protocol,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
                required_outcome_keywords=["response", "madrs"],
                outcome_timeframe_min_weeks=4.0,
                outcome_timeframe_max_weeks=16.0,
            )
        )
        self.assertFalse(
            trial_passes_design_filters(
                protocol,
                require_randomized=True,
                require_treatment_purpose=True,
                comparator_scope="active",
                required_outcome_keywords=["response", "madrs"],
                outcome_timeframe_min_weeks=9.0,
                outcome_timeframe_max_weeks=12.0,
            )
        )


class TestModelGWAMScript(unittest.TestCase):
    def test_results_only_not_forced_as_observed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "registry.csv"
            out_path = tmp / "summary.json"

            rows = [
                # Published proxy
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "true",
                    "has_results": "true",
                    "weight_n": "100",
                },
                # Results-only (no PMID)
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "false",
                    "has_results": "true",
                    "weight_n": "100",
                },
                # Strict ghost
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "100",
                },
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"]
                )
                writer.writeheader()
                writer.writerows(rows)

            cmd = [
                sys.executable,
                str(SCRIPTS_DIR / "model_gwam.py"),
                "--registry-csv",
                str(csv_path),
                "--published-mu",
                "1.0",
                "--sim-n",
                "100",
                "--results-only-mode",
                "as_unknown",
                "--results-only-mu",
                "0.0",
                "--results-only-sd",
                "0.0",
                "--ghost-mu",
                "0.0",
                "--ghost-sd",
                "0.0",
                "--output-json",
                str(out_path),
            ]
            subprocess.run(cmd, check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))

            # PMID-only lambda should be 1/3 (only first row contributes to published evidence).
            self.assertAlmostEqual(
                payload["weights"]["integrity_ratio_lambda_pmid_only"],
                1.0 / 3.0,
                places=6,
            )
            # With conservative unknown assumptions, deterministic GWAM point should also be 1/3.
            self.assertAlmostEqual(payload["estimates"]["mu_gwam_null_point"], 1.0 / 3.0, places=6)

    def test_model_respects_ghost_flag_for_no_pmid(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "registry.csv"
            out_path = tmp / "summary.json"

            rows = [
                # PMID published
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "true",
                    "has_results": "true",
                    "weight_n": "100",
                },
                # no PMID but has results, flagged as ghost (no_pmid mode)
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "true",
                    "weight_n": "100",
                },
                # strict ghost
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "100",
                },
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"]
                )
                writer.writeheader()
                writer.writerows(rows)

            cmd = [
                sys.executable,
                str(SCRIPTS_DIR / "model_gwam.py"),
                "--registry-csv",
                str(csv_path),
                "--published-mu",
                "1.0",
                "--sim-n",
                "30",
                "--output-json",
                str(out_path),
            ]
            subprocess.run(cmd, check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))

            self.assertEqual(payload["n_trials_results_only"], 0)
            self.assertEqual(payload["n_trials_ghost_proxy"], 2)
            self.assertEqual(payload["n_trials_ghost_with_results_no_pmid"], 1)
            self.assertAlmostEqual(
                payload["weights"]["integrity_ratio_lambda_pmid_only"],
                1.0 / 3.0,
                places=6,
            )

    def test_model_and_sim_lambda_consistent_no_pmid(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "registry.csv"
            model_out = tmp / "model.json"
            sim_out = tmp / "sim.json"

            rows = [
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "true",
                    "has_results": "true",
                    "weight_n": "200",
                },
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "true",
                    "weight_n": "100",
                },
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "100",
                },
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"]
                )
                writer.writeheader()
                writer.writerows(rows)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "model_gwam.py"),
                    "--registry-csv",
                    str(csv_path),
                    "--published-mu",
                    "0.2",
                    "--sim-n",
                    "30",
                    "--output-json",
                    str(model_out),
                ],
                check=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "simulate_gwam_vs_re.py"),
                    "--registry-csv",
                    str(csv_path),
                    "--n-meta",
                    "30",
                    "--mu-true-list",
                    "0.2",
                    "--output-json",
                    str(sim_out),
                ],
                check=True,
            )

            model_payload = json.loads(model_out.read_text(encoding="utf-8"))
            sim_payload = json.loads(sim_out.read_text(encoding="utf-8"))
            self.assertAlmostEqual(
                model_payload["weights"]["integrity_ratio_lambda_pmid_only"],
                sim_payload["lambda_observed_registry_raw"],
                places=6,
            )

    def test_sim_treats_non_ghost_no_pmid_as_results_only_pool(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "registry.csv"
            sim_out = tmp / "sim.json"

            rows = [
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "true",
                    "has_results": "true",
                    "weight_n": "200",
                },
                {
                    # Non-ghost with missing publication flags should be in results-only/unknown pool.
                    "is_ghost_protocol": "false",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "100",
                },
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "100",
                },
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"]
                )
                writer.writeheader()
                writer.writerows(rows)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "simulate_gwam_vs_re.py"),
                    "--registry-csv",
                    str(csv_path),
                    "--n-meta",
                    "30",
                    "--mu-true-list",
                    "0.2",
                    "--output-json",
                    str(sim_out),
                ],
                check=True,
            )
            sim_payload = json.loads(sim_out.read_text(encoding="utf-8"))
            self.assertEqual(sim_payload["n_trials_results_only_pool"], 1)

    def test_sim_calibration_reports_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "registry.csv"
            sim_out = tmp / "sim.json"

            rows = [
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "true",
                    "has_results": "true",
                    "weight_n": "220",
                },
                {
                    "is_ghost_protocol": "false",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "120",
                },
                {
                    "is_ghost_protocol": "true",
                    "has_pmid": "false",
                    "has_results": "false",
                    "weight_n": "80",
                },
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["is_ghost_protocol", "has_pmid", "has_results", "weight_n"]
                )
                writer.writeheader()
                writer.writerows(rows)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "simulate_gwam_vs_re.py"),
                    "--registry-csv",
                    str(csv_path),
                    "--n-meta",
                    "40",
                    "--mu-true-list",
                    "0.0,0.2",
                    "--calibrate-nonsig",
                    "--calibration-runs",
                    "120",
                    "--calibration-grid-size",
                    "9",
                    "--calibration-refine-rounds",
                    "1",
                    "--ci-calibration-runs",
                    "120",
                    "--output-json",
                    str(sim_out),
                ],
                check=True,
            )
            sim_payload = json.loads(sim_out.read_text(encoding="utf-8"))
            self.assertTrue(sim_payload["calibration"]["enabled"])
            self.assertIn("lambda_calibration_abs_error", sim_payload["calibration"])
            self.assertIsNotNone(sim_payload["calibration"]["p_nonsig_calibrated"])
            self.assertEqual(sim_payload["calibration_objective"], "pmid_only")
            self.assertEqual(sim_payload["ci_calibration"]["mode"], "empirical")
            self.assertTrue(sim_payload["ci_calibration"]["out_of_sample"])
            self.assertGreaterEqual(sim_payload["ci_calibration"]["calibration_runs_success"], 1)
            scen = sim_payload["scenario_results"][0]
            self.assertIn("oracle_ipw_random_effects", scen)
            self.assertIn("pet_peese", scen)
            self.assertIn("n_valid", scen["oracle_ipw_random_effects"])
            self.assertIn("n_valid", scen["pet_peese"])


    def test_model_gwam_all_ghosts_lambda_zero(self) -> None:
        """model_gwam.py handles lambda=0 (all trials are ghosts) without crashing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            csv_path = tmp / "all_ghost.csv"
            model_out = tmp / "gwam.json"

            rows = [
                {"is_ghost_protocol": "true", "has_pmid": "false", "has_results": "false", "weight_n": "100"},
                {"is_ghost_protocol": "true", "has_pmid": "false", "has_results": "false", "weight_n": "200"},
            ]
            fieldnames = list(rows[0].keys())
            with csv_path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

            proc = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "model_gwam.py"),
                    "--registry-csv", str(csv_path),
                    "--published-mu", "0.3",
                    "--sim-n", "100",
                    "--output-json", str(model_out),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            payload = json.loads(model_out.read_text(encoding="utf-8"))
            self.assertAlmostEqual(payload["weights"]["integrity_ratio_lambda_pmid_only"], 0.0)
            self.assertAlmostEqual(payload["weights"]["integrity_ratio_lambda_non_ghost"], 0.0)


class TestRunWorkflowGuards(unittest.TestCase):
    def test_published_mu_scope_mismatch_fails_fast(self) -> None:
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "run_full_workflow.py"),
            "--intervention",
            "escitalopram",
            "--condition",
            "depression",
            "--published-mu",
            "0.2",
            "--published-mu-comparator-scope",
            "placebo",
            "--comparator-scope",
            "active",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("must match --comparator-scope", proc.stderr)


if __name__ == "__main__":
    unittest.main()
