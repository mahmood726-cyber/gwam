#!/usr/bin/env python3
"""Unit tests for Pairwise70 CT.gov linkage summary builder."""

from __future__ import annotations

import unittest

# sys.path setup handled by conftest.py
from build_pairwise70_ctgov_linkage_summary import (
    extract_nct_ids_from_text,
    normalize_studies_shape,
    parse_review_id_from_filename,
    summarize_payload,
)


class TestCTGovLinkageBuilder(unittest.TestCase):
    def test_extract_nct_ids_from_text(self) -> None:
        text = "Protocol IDs NCT12345678 and nct00001111 appear here."
        out = extract_nct_ids_from_text(text)
        self.assertEqual(out, {"NCT12345678", "NCT00001111"})

    def test_parse_review_id_from_filename(self) -> None:
        self.assertEqual(parse_review_id_from_filename("CD001155_pub3_data.rda"), "CD001155")
        self.assertEqual(parse_review_id_from_filename("bad_name.rda"), "")

    def test_normalize_studies_shape_list(self) -> None:
        studies = [{"protocolSection": {"statusModule": {"overallStatus": "COMPLETED"}}}]
        out = normalize_studies_shape(studies)
        self.assertEqual(len(out), 1)
        self.assertIn("protocolSection", out[0])

    def test_normalize_studies_shape_columnar(self) -> None:
        studies = {
            "protocolSection": [{"statusModule": {"overallStatus": "COMPLETED"}}],
            "protocolSection.1": [{"statusModule": {"overallStatus": "RECRUITING"}}],
        }
        out = normalize_studies_shape(studies)
        self.assertEqual(len(out), 2)
        self.assertIn("protocolSection", out[0])

    def test_summarize_payload_basic(self) -> None:
        payload = {
            "query": "demo",
            "studies": [
                {
                    "hasResults": False,
                    "protocolSection": {
                        "statusModule": {"overallStatus": "COMPLETED"},
                        "designModule": {"enrollmentInfo": {"count": 100}},
                        "referencesModule": {"references": []},
                    },
                },
                {
                    "hasResults": True,
                    "protocolSection": {
                        "statusModule": {"overallStatus": "COMPLETED"},
                        "designModule": {"enrollmentInfo": {"count": 200}},
                        "referencesModule": {"references": []},
                    },
                },
                {
                    "hasResults": False,
                    "protocolSection": {
                        "statusModule": {"overallStatus": "COMPLETED"},
                        "designModule": {"enrollmentInfo": {"count": 300}},
                        "referencesModule": {
                            "references": [{"pmid": "12345", "type": "RESULT", "citation": "x"}]
                        },
                    },
                },
                {
                    "hasResults": False,
                    "protocolSection": {
                        "statusModule": {"overallStatus": "RECRUITING"},
                        "designModule": {"enrollmentInfo": {"count": 999}},
                        "referencesModule": {"references": []},
                    },
                },
            ],
        }
        out = summarize_payload("CDTEST", payload)
        self.assertEqual(out["review_id"], "CDTEST")
        self.assertEqual(out["n_completed"], 3)
        self.assertEqual(out["n_completed_ghost"], 1)
        self.assertEqual(out["n_completed_with_results"], 1)
        self.assertEqual(out["n_completed_with_pmid_any"], 1)
        self.assertAlmostEqual(out["w_completed_total"], 600.0, places=6)
        self.assertAlmostEqual(out["lambda_completed_non_ghost_weighted"], 500.0 / 600.0, places=6)
        self.assertAlmostEqual(out["lambda_completed_pmid_any_weighted"], 300.0 / 600.0, places=6)


if __name__ == "__main__":
    unittest.main()
