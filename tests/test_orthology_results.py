from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from unittest import TestCase


ROOT = Path(__file__).parents[1]
RESULT = ROOT / "results" / "controlled_orthology"


def rows(name: str) -> list[dict[str, str]]:
    with (RESULT / "tables" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class ControlledOrthologyResultTests(TestCase):
    def test_expected_relationship_classes_are_recovered(self) -> None:
        comparisons = {
            row["baseline_gene_id"]: row
            for row in rows("orthology_relationships.csv")
            if row["alternative_accession"] == "CTRL_ALT.1"
        }
        expected = {
            "base_one": "one_to_one",
            "base_old_name": "one_to_one",
            "base_one_many": "one_to_many",
            "base_many_a": "many_to_one",
            "base_many_b": "many_to_one",
            "base_mm_a": "many_to_many",
            "base_mm_b": "many_to_many",
            "base_missing_a": "no_ortholog",
            "base_missing_b": "no_ortholog",
            "base_unassigned": "no_ortholog",
            "base_no_protein": "no_protein",
        }
        self.assertEqual(
            {gene_id: row["relationship"] for gene_id, row in comparisons.items()},
            expected,
        )

    def test_renamed_one_to_one_is_distinguished_from_stable_symbol(self) -> None:
        comparisons = {
            row["baseline_gene_id"]: row
            for row in rows("orthology_relationships.csv")
            if row["alternative_accession"] == "CTRL_ALT.1"
        }
        self.assertEqual(comparisons["base_one"]["renamed_one_to_one"], "False")
        self.assertEqual(comparisons["base_old_name"]["alternative_gene_names"], "new_symbol")
        self.assertEqual(comparisons["base_old_name"]["renamed_one_to_one"], "True")

    def test_every_prepared_protein_is_assigned_once(self) -> None:
        protein_map = [
            row for row in rows("protein_gene_map.csv") if row["prepared_protein_id"]
        ]
        members = rows("orthogroup_members.csv")
        counts = Counter(row["prepared_protein_id"] for row in members)
        self.assertEqual(set(counts), {row["prepared_protein_id"] for row in protein_map})
        self.assertTrue(all(count == 1 for count in counts.values()))
        unassigned = {
            row["gene_id"] for row in members if row["assignment_status"] == "unassigned"
        }
        self.assertIn("base_unassigned", unassigned)
        self.assertIn("alt_unassigned", unassigned)

    def test_no_protein_gene_is_explicit_in_the_mapping(self) -> None:
        matches = [
            row
            for row in rows("protein_gene_map.csv")
            if row["gene_id"] == "base_no_protein"
        ]
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["mapping_status"], "no_protein")
        self.assertEqual(matches[0]["prepared_protein_id"], "")

    def test_manifest_records_deterministic_real_tool_settings(self) -> None:
        manifest = json.loads((RESULT / "run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["parameters"]["method"], "msa")
        self.assertEqual(manifest["parameters"]["sequence_search"], "diamond")
        self.assertEqual(manifest["parameters"]["msa_program"], "famsa")
        self.assertEqual(manifest["parameters"]["tree_program"], "fasttree")
        self.assertEqual(manifest["parameters"]["threads"], 1)
        self.assertEqual(manifest["counts"]["prepared_proteins"], 24)
        self.assertEqual(manifest["counts"]["orthogroup_members"], 24)
        self.assertTrue(manifest["protein_accounting_complete"])

