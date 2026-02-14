#!/usr/bin/env python3
"""Unit tests for PubMed bridge and fuzzy matching in linkage builder."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from build_pairwise70_ctgov_linkage_summary import (
    fuzzy_match_ctgov,
    pubmed_to_nct,
)


# --- Mock PubMed XML for testing ---
MOCK_PUBMED_XML_TWO_ARTICLES = """\
<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>12345678</PMID>
      <Article>
        <ArticleTitle>Study of Drug X for Depression</ArticleTitle>
        <Abstract>
          <AbstractText>
            This trial was registered at ClinicalTrials.gov (NCT00112233).
            Results showed significant improvement.
          </AbstractText>
        </Abstract>
        <DataBankList>
          <DataBank>
            <DataBankName>ClinicalTrials.gov</DataBankName>
            <AccessionNumberList>
              <AccessionNumber>NCT00112233</AccessionNumber>
            </AccessionNumberList>
          </DataBank>
        </DataBankList>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>87654321</PMID>
      <Article>
        <ArticleTitle>Another Study</ArticleTitle>
        <Abstract>
          <AbstractText>
            Registered under NCT99887766 and NCT11223344.
          </AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>
"""

MOCK_PUBMED_XML_NO_NCT = """\
<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>11111111</PMID>
      <Article>
        <ArticleTitle>Study without NCT ID</ArticleTitle>
        <Abstract>
          <AbstractText>No registration mentioned.</AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>
"""


class TestPubmedToNct(unittest.TestCase):
    """Test NCT ID extraction from PubMed records."""

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    @patch("build_pairwise70_ctgov_linkage_summary.time.sleep")
    def test_extracts_nct_from_databank_and_abstract(
        self, mock_sleep: MagicMock, mock_get: MagicMock
    ) -> None:
        """Should find NCT IDs from both DataBankList and abstract text."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.text = MOCK_PUBMED_XML_TWO_ARTICLES
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = pubmed_to_nct(["12345678", "87654321"])

        # PMID 12345678: NCT00112233 from DataBankList + abstract
        self.assertIn("12345678", result)
        self.assertIn("NCT00112233", result["12345678"])

        # PMID 87654321: NCT99887766 and NCT11223344 from abstract
        self.assertIn("87654321", result)
        self.assertIn("NCT99887766", result["87654321"])
        self.assertIn("NCT11223344", result["87654321"])

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    @patch("build_pairwise70_ctgov_linkage_summary.time.sleep")
    def test_no_nct_found_returns_empty(
        self, mock_sleep: MagicMock, mock_get: MagicMock
    ) -> None:
        """PMIDs without NCT references should not appear in results."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.text = MOCK_PUBMED_XML_NO_NCT
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = pubmed_to_nct(["11111111"])
        self.assertEqual(result, {})

    def test_empty_pmid_list(self) -> None:
        """Empty input should return empty dict without API call."""
        result = pubmed_to_nct([])
        self.assertEqual(result, {})

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    @patch("build_pairwise70_ctgov_linkage_summary.time.sleep")
    def test_api_failure_returns_empty(
        self, mock_sleep: MagicMock, mock_get: MagicMock
    ) -> None:
        """API errors should be handled gracefully."""
        mock_get.side_effect = Exception("Network error")
        result = pubmed_to_nct(["12345678"])
        self.assertEqual(result, {})


class TestFuzzyMatchCtgov(unittest.TestCase):
    """Test fuzzy matching against CT.gov."""

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    def test_unique_match_returns_nct_id(self, mock_get: MagicMock) -> None:
        """Single match within enrollment tolerance should return NCT ID."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {"nctId": "NCT55667788"},
                        "designModule": {
                            "enrollmentInfo": {"count": 105}
                        },
                        "statusModule": {
                            "completionDateStruct": {"date": "2020-06-15"}
                        },
                    }
                }
            ]
        }
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = fuzzy_match_ctgov(
            condition="depression",
            intervention="escitalopram",
            enrollment_count=100,  # 105 is within 20%
            tolerance=0.2,
        )
        self.assertEqual(result, "NCT55667788")

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    def test_multiple_matches_returns_none(self, mock_get: MagicMock) -> None:
        """Multiple matches should return None (conservative)."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {"nctId": "NCT11111111"},
                        "designModule": {"enrollmentInfo": {"count": 100}},
                        "statusModule": {},
                    }
                },
                {
                    "protocolSection": {
                        "identificationModule": {"nctId": "NCT22222222"},
                        "designModule": {"enrollmentInfo": {"count": 98}},
                        "statusModule": {},
                    }
                },
            ]
        }
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = fuzzy_match_ctgov(
            condition="diabetes",
            intervention="metformin",
            enrollment_count=100,
            tolerance=0.2,
        )
        self.assertIsNone(result)

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    def test_enrollment_outside_tolerance_returns_none(self, mock_get: MagicMock) -> None:
        """Match outside enrollment tolerance should be excluded."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {"nctId": "NCT33333333"},
                        "designModule": {"enrollmentInfo": {"count": 200}},
                        "statusModule": {},
                    }
                }
            ]
        }
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        result = fuzzy_match_ctgov(
            condition="cancer",
            intervention="pembrolizumab",
            enrollment_count=100,  # 200 is outside 20% tolerance
            tolerance=0.2,
        )
        self.assertIsNone(result)

    def test_missing_inputs_returns_none(self) -> None:
        """Missing condition/intervention/enrollment should return None."""
        self.assertIsNone(fuzzy_match_ctgov(
            condition="", intervention="drug", enrollment_count=100
        ))
        self.assertIsNone(fuzzy_match_ctgov(
            condition="disease", intervention="", enrollment_count=100
        ))
        self.assertIsNone(fuzzy_match_ctgov(
            condition="disease", intervention="drug", enrollment_count=0
        ))

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    def test_api_failure_returns_none(self, mock_get: MagicMock) -> None:
        """API errors should return None gracefully."""
        mock_get.side_effect = Exception("Timeout")
        result = fuzzy_match_ctgov(
            condition="depression",
            intervention="escitalopram",
            enrollment_count=100,
        )
        self.assertIsNone(result)

    @patch("build_pairwise70_ctgov_linkage_summary.requests.get")
    def test_completion_year_filter(self, mock_get: MagicMock) -> None:
        """When completion_year is specified, only matching year should be returned."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {"nctId": "NCT44444444"},
                        "designModule": {"enrollmentInfo": {"count": 100}},
                        "statusModule": {
                            "completionDateStruct": {"date": "2019-03-01"}
                        },
                    }
                }
            ]
        }
        mock_resp.raise_for_status.return_value = None
        mock_get.return_value = mock_resp

        # Wrong year -> None
        result = fuzzy_match_ctgov(
            condition="depression",
            intervention="drug",
            enrollment_count=100,
            completion_year=2020,
        )
        self.assertIsNone(result)

        # Correct year -> match
        result = fuzzy_match_ctgov(
            condition="depression",
            intervention="drug",
            enrollment_count=100,
            completion_year=2019,
        )
        self.assertEqual(result, "NCT44444444")


if __name__ == "__main__":
    unittest.main()
