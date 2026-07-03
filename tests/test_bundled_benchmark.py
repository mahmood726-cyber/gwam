#!/usr/bin/env python3
"""Tests for the self-contained bundled GWAM benchmark."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

# sys.path setup handled by conftest.py
PROJECT_ROOT = Path(__file__).resolve().parents[1]

from run_bundled_benchmark import (  # noqa: E402
    benchmark_dataset,
    discover_registry_csvs,
    format_table,
    main,
)


class TestBundledBenchmark(unittest.TestCase):
    def test_discovers_bundled_csvs(self) -> None:
        csvs = discover_registry_csvs()
        self.assertGreater(len(csvs), 0)
        # PubMed link tables must be excluded.
        self.assertFalse(any("pubmed" in p.name.lower() for p in csvs))

    def test_benchmark_record_shape(self) -> None:
        csvs = discover_registry_csvs()
        rec = benchmark_dataset(csvs[0])
        for key in (
            "dataset",
            "n_trials",
            "lambda_pmid_only",
            "lambda_non_ghost",
            "attenuation",
            "published_mu_source",
            "mu_gwam_null_point",
        ):
            self.assertIn(key, rec)
        self.assertGreaterEqual(rec["lambda_pmid_only"], 0.0)
        self.assertLessEqual(rec["lambda_pmid_only"], 1.0)
        self.assertLessEqual(rec["lambda_pmid_only"], rec["lambda_non_ghost"] + 1e-12)

    def test_deterministic(self) -> None:
        csvs = discover_registry_csvs()
        first = [benchmark_dataset(p) for p in csvs]
        second = [benchmark_dataset(p) for p in csvs]
        self.assertEqual(first, second)

    def test_documented_dataset_attenuation_equals_lambda(self) -> None:
        # Under the null-ghost assumption the corrected mean is a pure
        # multiplicative attenuation, so attenuation == lambda_pmid_only.
        csvs = {p.stem: p for p in discover_registry_csvs()}
        self.assertIn("escitalopram__depression", csvs)
        rec = benchmark_dataset(csvs["escitalopram__depression"])
        self.assertEqual(rec["published_mu_source"], "documented")
        self.assertAlmostEqual(rec["attenuation"], rec["lambda_pmid_only"], places=12)
        # mu_gwam = lambda * published_mu.
        self.assertAlmostEqual(
            rec["mu_gwam_null_point"],
            rec["lambda_pmid_only"] * rec["published_mu"],
            places=12,
        )

    def test_format_table_runs(self) -> None:
        csvs = discover_registry_csvs()
        recs = [benchmark_dataset(p) for p in csvs]
        table = format_table(recs)
        self.assertIn("dataset", table)
        self.assertIn("lam_pmid", table)

    def test_main_writes_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_csv = Path(tmp) / "bench.csv"
            out_json = Path(tmp) / "bench.json"
            rc = main(["--output-csv", str(out_csv), "--output-json", str(out_json)])
            self.assertEqual(rc, 0)
            self.assertTrue(out_csv.exists())
            self.assertTrue(out_json.exists())
            with out_csv.open("r", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertGreater(len(rows), 0)
            payload = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertIn("metadata", payload)
            self.assertIn("results", payload)
            self.assertEqual(len(payload["results"]), len(rows))


if __name__ == "__main__":
    unittest.main()
