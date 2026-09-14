from __future__ import annotations

import csv
import json
import struct
from collections import Counter
from pathlib import Path
from unittest import TestCase


ROOT = Path(__file__).parents[1]
RESULT = ROOT / "results" / "authentic_comparison"


def rows(relative_path: str) -> list[dict[str, str]]:
    with (RESULT / relative_path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class AuthenticComparisonResultTests(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(
            (RESULT / "run_manifest.json").read_text(encoding="utf-8")
        )
        cls.burdens = rows("tables/gene_consequence_burden.csv")
        cls.comparisons = rows("tables/baseline_ortholog_comparison.csv")
        cls.review = rows("tables/review_candidate_orthology.csv")

    def test_manifest_has_exact_cross_artifact_accounting(self) -> None:
        counts = self.manifest["counts"]
        self.assertEqual(self.manifest["status"], "pass")
        self.assertEqual(counts["references"], 6)
        self.assertEqual(counts["accepted_variants"], 357_581)
        self.assertEqual(counts["annotation_effects"], 357_959)
        self.assertEqual(counts["genes"], 28_859)
        self.assertEqual(counts["variant_gene_assignments"], 327_542)
        self.assertEqual(counts["effects_without_exact_gene_key"], 30_382)
        self.assertEqual(counts["effects_collapsed_within_variant_gene"], 35)
        self.assertEqual(
            counts["variant_gene_assignments"]
            + counts["effects_without_exact_gene_key"]
            + counts["effects_collapsed_within_variant_gene"],
            counts["annotation_effects"],
        )
        self.assertEqual(counts["comparable_gene_pairs"], 6_356)
        self.assertTrue(self.manifest["source_inputs_unchanged"])

    def test_all_six_callsets_reproduce_project13_counts(self) -> None:
        expected = {
            "GCF_000005845.2": (30_973, 31_059),
            "GCF_003697165.2": (92_665, 92_706),
            "GCF_036503815.1": (125_327, 125_378),
            "GCF_005843885.1": (49_201, 49_230),
            "GCF_000026225.1": (15_264, 15_266),
            "GCF_002900365.1": (44_151, 44_320),
        }
        observed = {
            row["accession"]: (
                int(row["accepted_variants"]),
                int(row["annotation_effects"]),
            )
            for row in self.manifest["per_reference_annotations"]
        }
        self.assertEqual(observed, expected)
        self.assertTrue(self.manifest["project14_effects_exactly_reproduced"])
        reproduction = rows("summaries/project14_reproduction.csv")
        self.assertEqual(reproduction[0]["project14_effects"], "31059")
        self.assertEqual(reproduction[0]["project15_effects"], "31059")
        self.assertEqual(reproduction[0]["passed"], "True")

    def test_gene_burden_has_one_valid_row_per_gene(self) -> None:
        keys = [(row["accession"], row["gene_id"]) for row in self.burdens]
        self.assertEqual(len(keys), 28_859)
        self.assertEqual(len(keys), len(set(keys)))
        self.assertTrue(
            all(
                int(row["callable_cds_bases"]) <= int(row["cds_bases"])
                for row in self.burdens
            )
        )
        self.assertEqual(
            sum(int(row["distinct_variant_count"]) for row in self.burdens),
            327_542,
        )
        for row in self.burdens:
            impact_total = sum(
                int(row[name])
                for name in (
                    "high_impact_variants",
                    "moderate_impact_variants",
                    "low_impact_variants",
                    "modifier_impact_variants",
                )
            )
            self.assertEqual(impact_total, int(row["distinct_variant_count"]))

    def test_comparisons_are_unique_and_comparability_is_strict(self) -> None:
        keys = [
            (row["baseline_gene_id"], row["alternative_accession"])
            for row in self.comparisons
        ]
        self.assertEqual(len(keys), 23_255)
        self.assertEqual(len(keys), len(set(keys)))
        per_alternative = Counter(
            row["alternative_accession"] for row in self.comparisons
        )
        self.assertEqual(set(per_alternative.values()), {4_651})
        comparable = [row for row in self.comparisons if row["comparable"] == "True"]
        self.assertEqual(len(comparable), 6_356)
        for row in comparable:
            self.assertEqual(row["relationship"], "one_to_one")
            self.assertGreaterEqual(float(row["baseline_callable_fraction"]), 0.90)
            self.assertGreaterEqual(float(row["alternative_callable_fraction"]), 0.90)
            self.assertEqual(row["project13_alternative_eligible"], "False")

    def test_comparable_density_summary_matches_detailed_rows(self) -> None:
        summary = rows("summaries/comparable_density_summary.csv")
        expected_counts = {
            "GCF_003697165.2": 3_206,
            "GCF_036503815.1": 3_137,
            "GCF_005843885.1": 9,
            "GCF_000026225.1": 1,
            "GCF_002900365.1": 3,
        }
        self.assertEqual(
            {
                row["alternative_accession"]: int(row["comparable_gene_pairs"])
                for row in summary
            },
            expected_counts,
        )
        for row in summary:
            accession = row["alternative_accession"]
            detailed = [
                item
                for item in self.comparisons
                if item["alternative_accession"] == accession
                and item["comparable"] == "True"
            ]
            self.assertEqual(len(detailed), int(row["comparable_gene_pairs"]))
            self.assertEqual(
                sum(int(item["baseline_distinct_variant_count"]) for item in detailed),
                int(row["baseline_variant_gene_assignments"]),
            )
            self.assertEqual(
                sum(
                    int(item["alternative_distinct_variant_count"])
                    for item in detailed
                ),
                int(row["alternative_variant_gene_assignments"]),
            )

    def test_review_candidates_have_five_gene_level_contexts(self) -> None:
        self.assertEqual(len(self.review), 400)
        per_variant = Counter(row["baseline_variant_key"] for row in self.review)
        self.assertEqual(len(per_variant), 80)
        self.assertEqual(set(per_variant.values()), {5})
        self.assertTrue(
            all(
                row["comparison_scope"] == "gene_level_only_no_variant_equivalence"
                for row in self.review
            )
        )

    def test_dashboard_is_a_full_resolution_png(self) -> None:
        path = RESULT / "dashboard" / "comparative_gene_dashboard.png"
        with path.open("rb") as handle:
            header = handle.read(24)
        self.assertEqual(header[:8], b"\x89PNG\r\n\x1a\n")
        width, height = struct.unpack(">II", header[16:24])
        self.assertEqual((width, height), (2880, 1980))
        self.assertGreater(path.stat().st_size, 100_000)
