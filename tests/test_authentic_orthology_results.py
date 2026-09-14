from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from unittest import TestCase


ROOT = Path(__file__).parents[1]
RESULT = ROOT / "results" / "orthology_analysis"


def rows(name: str) -> list[dict[str, str]]:
    with (RESULT / "tables" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class AuthenticOrthologyResultTests(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(
            (RESULT / "run_manifest.json").read_text(encoding="utf-8")
        )
        cls.protein_map = rows("protein_gene_map.csv")
        cls.members = rows("orthogroup_members.csv")
        cls.pairwise = rows("pairwise_orthologs.csv")
        cls.relationships = rows("orthology_relationships.csv")

    def test_manifest_records_complete_six_reference_analysis(self) -> None:
        counts = self.manifest["counts"]
        self.assertEqual(self.manifest["status"], "pass")
        self.assertEqual(counts["references"], 6)
        self.assertEqual(counts["prepared_proteins"], 26_768)
        self.assertEqual(counts["orthogroup_members"], 26_768)
        self.assertEqual(counts["orthogroups_and_unassigned"], 6_170)
        self.assertEqual(counts["pairwise_protein_orthologs"], 56_346)
        self.assertEqual(counts["baseline_alternative_gene_comparisons"], 23_255)
        self.assertTrue(self.manifest["protein_accounting_complete"])

    def test_every_prepared_protein_occurs_in_exactly_one_group(self) -> None:
        prepared = {
            row["prepared_protein_id"]
            for row in self.protein_map
            if row["prepared_protein_id"]
        }
        membership_counts = Counter(
            row["prepared_protein_id"] for row in self.members
        )
        self.assertEqual(len(prepared), 26_768)
        self.assertEqual(set(membership_counts), prepared)
        self.assertTrue(all(count == 1 for count in membership_counts.values()))

    def test_every_pairwise_edge_links_members_of_the_same_group(self) -> None:
        group_by_protein = {
            row["prepared_protein_id"]: row["orthogroup"] for row in self.members
        }
        for edge in self.pairwise:
            self.assertEqual(
                group_by_protein[edge["protein_a"]], edge["orthogroup"]
            )
            self.assertEqual(
                group_by_protein[edge["protein_b"]], edge["orthogroup"]
            )

    def test_each_alternative_has_one_row_per_baseline_gene(self) -> None:
        per_reference = Counter(
            row["alternative_accession"] for row in self.relationships
        )
        self.assertEqual(len(per_reference), 5)
        self.assertEqual(set(per_reference.values()), {4_651})
        keys = [
            (row["baseline_gene_id"], row["alternative_accession"])
            for row in self.relationships
        ]
        self.assertEqual(len(keys), len(set(keys)))

    def test_one_to_one_rows_are_unambiguous_gene_pairs(self) -> None:
        one_to_one = [
            row for row in self.relationships if row["relationship"] == "one_to_one"
        ]
        self.assertTrue(one_to_one)
        for row in one_to_one:
            self.assertEqual(row["baseline_family_gene_count"], "1")
            self.assertEqual(row["alternative_family_gene_count"], "1")
            self.assertNotEqual(row["alternative_gene_ids"], "")
            self.assertNotIn(";", row["alternative_gene_ids"])

    def test_relationship_summary_exactly_matches_comparison_rows(self) -> None:
        observed = Counter(
            (row["alternative_accession"], row["relationship"])
            for row in self.relationships
        )
        with (
            RESULT / "reports" / "orthology_relationship_summary.csv"
        ).open(encoding="utf-8", newline="") as handle:
            reported = {
                (row["alternative_accession"], row["relationship"]): int(
                    row["gene_count"]
                )
                for row in csv.DictReader(handle)
            }
        self.assertEqual(dict(observed), reported)
