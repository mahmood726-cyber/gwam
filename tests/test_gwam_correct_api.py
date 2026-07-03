#!/usr/bin/env python3
"""Unit tests for the reusable GWAM correction API in gwam_utils.

Covers:
- Deterministic correction math (lambda ratios + weighted null-point mean).
- Input validation on the public entry points.
- A parity check that gwam_correct() reproduces the numbers that the
  model_gwam.py CLI writes to its JSON summary on bundled registry data.
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# sys.path setup handled by conftest.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

from gwam_utils import (  # noqa: E402
    GwamResult,
    gwam_correct,
    partition_registry_weights,
)


def _row(*, ghost: bool, pmid: bool, weight: float, results: bool = False) -> dict:
    return {
        "is_ghost_protocol": "True" if ghost else "False",
        "has_pmid": "True" if pmid else "False",
        "has_results": "True" if results else "False",
        "weight_n": str(weight),
    }


class TestPartition(unittest.TestCase):
    def test_three_strata(self) -> None:
        rows = [
            _row(ghost=False, pmid=True, weight=10),
            _row(ghost=False, pmid=False, weight=20),  # results-only
            _row(ghost=True, pmid=False, weight=30),
        ]
        pmid_w, ro_w, ghost_w = partition_registry_weights(rows)
        self.assertEqual(pmid_w, [10.0])
        self.assertEqual(ro_w, [20.0])
        self.assertEqual(ghost_w, [30.0])

    def test_non_ghost_no_pmid_goes_to_results_only(self) -> None:
        # A non-ghost row without a PMID lands in results-only regardless of
        # the has_results flag (matches model_gwam.py stratification).
        rows = [
            _row(ghost=False, pmid=False, weight=5, results=True),
            _row(ghost=False, pmid=False, weight=7, results=False),
        ]
        pmid_w, ro_w, ghost_w = partition_registry_weights(rows)
        self.assertEqual(pmid_w, [])
        self.assertEqual(ro_w, [5.0, 7.0])
        self.assertEqual(ghost_w, [])

    def test_empty_rows_raises(self) -> None:
        with self.assertRaises(ValueError):
            partition_registry_weights([])

    def test_missing_weight_column_raises(self) -> None:
        with self.assertRaises(ValueError):
            partition_registry_weights([{"is_ghost_protocol": "False"}])

    def test_non_numeric_weight_raises(self) -> None:
        with self.assertRaises(ValueError):
            partition_registry_weights([_row(ghost=False, pmid=True, weight="abc")])  # type: ignore[arg-type]

    def test_nonpositive_weight_raises(self) -> None:
        with self.assertRaises(ValueError):
            partition_registry_weights([_row(ghost=False, pmid=True, weight=0)])
        with self.assertRaises(ValueError):
            partition_registry_weights([_row(ghost=False, pmid=True, weight=-3)])


class TestGwamCorrect(unittest.TestCase):
    def test_known_values(self) -> None:
        # 100 pmid, 100 results-only, 100 ghost. lambda_pmid = 1/3,
        # lambda_non_ghost = 2/3. With ghost_mu = results_only_mu = 0 and
        # published_mu = 0.6 -> mu = (100*0.6)/300 = 0.2.
        rows = [
            _row(ghost=False, pmid=True, weight=100),
            _row(ghost=False, pmid=False, weight=100),
            _row(ghost=True, pmid=False, weight=100),
        ]
        r = gwam_correct(rows, published_mu=0.6)
        self.assertIsInstance(r, GwamResult)
        self.assertAlmostEqual(r.lambda_pmid_only, 1.0 / 3.0, places=12)
        self.assertAlmostEqual(r.lambda_non_ghost, 2.0 / 3.0, places=12)
        self.assertAlmostEqual(r.mu_gwam_null_point, 0.2, places=12)
        self.assertEqual(r.n_pmid_only, 1)
        self.assertEqual(r.n_results_only, 1)
        self.assertEqual(r.n_ghost, 1)

    def test_all_published_recovers_full_effect(self) -> None:
        # No ghosts, no results-only: lambda == 1 and correction == published.
        rows = [_row(ghost=False, pmid=True, weight=50) for _ in range(4)]
        r = gwam_correct(rows, published_mu=0.42)
        self.assertAlmostEqual(r.lambda_pmid_only, 1.0, places=12)
        self.assertAlmostEqual(r.lambda_non_ghost, 1.0, places=12)
        self.assertAlmostEqual(r.mu_gwam_null_point, 0.42, places=12)

    def test_as_observed_mode(self) -> None:
        # Under as_observed, results-only studies are held at published_mu.
        rows = [
            _row(ghost=False, pmid=True, weight=100),
            _row(ghost=False, pmid=False, weight=100),
            _row(ghost=True, pmid=False, weight=100),
        ]
        r = gwam_correct(rows, published_mu=0.6, results_only_mode="as_observed")
        # (100*0.6 + 100*0.6 + 100*0) / 300 = 0.4
        self.assertAlmostEqual(r.mu_gwam_null_point, 0.4, places=12)

    def test_nonzero_ghost_mu(self) -> None:
        rows = [
            _row(ghost=False, pmid=True, weight=100),
            _row(ghost=True, pmid=False, weight=100),
        ]
        r = gwam_correct(rows, published_mu=0.5, ghost_mu=0.5)
        # (100*0.5 + 100*0.5) / 200 = 0.5
        self.assertAlmostEqual(r.mu_gwam_null_point, 0.5, places=12)

    def test_nonfinite_published_mu_raises(self) -> None:
        rows = [_row(ghost=False, pmid=True, weight=10)]
        with self.assertRaises(ValueError):
            gwam_correct(rows, published_mu=math.inf)
        with self.assertRaises(ValueError):
            gwam_correct(rows, published_mu=math.nan)

    def test_bad_results_only_mode_raises(self) -> None:
        rows = [_row(ghost=False, pmid=True, weight=10)]
        with self.assertRaises(ValueError):
            gwam_correct(rows, published_mu=0.1, results_only_mode="bogus")

    def test_to_dict_roundtrip(self) -> None:
        rows = [_row(ghost=False, pmid=True, weight=10)]
        d = gwam_correct(rows, published_mu=0.1).to_dict()
        self.assertIn("lambda_pmid_only", d)
        self.assertIn("mu_gwam_null_point", d)
        # JSON-serialisable.
        json.dumps(d)


class TestParityWithCLI(unittest.TestCase):
    """gwam_correct() must reproduce the CLI's JSON numbers on bundled data."""

    def test_matches_model_gwam_cli(self) -> None:
        registry = PROJECT_ROOT / "data" / "raw" / "escitalopram__depression.csv"
        if not registry.exists():
            self.skipTest("bundled registry CSV not present")

        published_mu = 0.25464221837358075
        with registry.open("r", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        api = gwam_correct(rows, published_mu=published_mu)

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "cli.json"
            proc = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS_DIR / "model_gwam.py"),
                    "--registry-csv",
                    str(registry),
                    "--published-mu",
                    repr(published_mu),
                    "--sim-n",
                    "2",
                    "--output-json",
                    str(out),
                ],
                capture_output=True,
                text=True,
                cwd=str(SCRIPTS_DIR),
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            cli = json.loads(out.read_text(encoding="utf-8"))

        w = cli["weights"]
        self.assertEqual(api.lambda_pmid_only, w["integrity_ratio_lambda_pmid_only"])
        self.assertEqual(api.lambda_non_ghost, w["integrity_ratio_lambda_non_ghost"])
        self.assertEqual(api.mu_gwam_null_point, cli["estimates"]["mu_gwam_null_point"])
        self.assertEqual(api.weight_total, w["total"])
        self.assertEqual(api.n_ghost, cli["n_trials_ghost_proxy"])


if __name__ == "__main__":
    unittest.main()
