from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest import TestCase


ROOT = Path(__file__).parents[1]
RESULT = ROOT / "results" / "controlled_core"


def read_csv(name: str) -> list[dict[str, str]]:
    with (RESULT / "tables" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class ControlledCoreResultTests(TestCase):
    def test_plus_and_minus_strand_annotations_are_correct(self) -> None:
        rows = read_csv("variant_annotations.csv")
        self.assertEqual(len(rows), 4)
        by_accession_position = {
            (row["accession"], int(row["position"])): row for row in rows
        }
        for accession in ("CTRL_PLUS.1", "CTRL_MINUS.1"):
            plus = by_accession_position[(accession, 4)]
            minus = by_accession_position[(accession, 387)]
            self.assertEqual(plus["primary_gene_id"], "ctrlA")
            self.assertEqual(plus["primary_consequence"], "missense_variant")
            self.assertEqual(plus["primary_hgvs_p"], "p.Ala2Thr")
            self.assertEqual(minus["primary_gene_id"], "ctrlB")
            self.assertEqual(minus["primary_consequence"], "stop_gained")
            self.assertEqual(minus["primary_hgvs_p"], "p.Glu2*")

    def test_low_callability_is_reference_and_gene_specific(self) -> None:
        rows = read_csv("gene_callability.csv")
        observed = {
            (row["accession"], row["gene_id"]): (
                float(row["callable_fraction"]), row["status"]
            )
            for row in rows
        }
        self.assertEqual(observed[("CTRL_PLUS.1", "ctrlA")], (1.0, "callable"))
        self.assertEqual(observed[("CTRL_PLUS.1", "ctrlB")], (0.5, "low_callability"))
        self.assertEqual(observed[("CTRL_MINUS.1", "ctrlA")], (0.5, "low_callability"))
        self.assertEqual(observed[("CTRL_MINUS.1", "ctrlB")], (1.0, "callable"))

    def test_accounting_and_database_gates_agree(self) -> None:
        manifest = json.loads((RESULT / "run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["counts"]["references"], 2)
        self.assertEqual(manifest["counts"]["genes"], 4)
        self.assertEqual(manifest["counts"]["callable_coding_genes"], 2)
        self.assertEqual(manifest["counts"]["accepted_variants"], 4)
        self.assertEqual(manifest["counts"]["annotation_effects"], 4)
        self.assertTrue(manifest["source_inputs_unchanged"])
        with (RESULT / "reports" / "database_validation.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            checks = list(csv.DictReader(handle))
        self.assertEqual(len(checks), 4)
        self.assertTrue(all(row["passed"] == "True" for row in checks))

    def test_each_reference_has_readable_annotation_artifacts(self) -> None:
        for accession in ("CTRL_PLUS.1", "CTRL_MINUS.1"):
            root = RESULT / "per_reference" / accession
            self.assertTrue((root / "vcf" / "accepted.annotated.vcf").is_file())
            self.assertTrue((root / "vcf" / "accepted.annotated.vcf.gz").is_file())
            self.assertTrue((root / "vcf" / "accepted.annotated.vcf.gz.csi").is_file())
            self.assertTrue((root / "reports" / "snpeff_summary.html").is_file())

    def test_published_artifacts_have_no_machine_or_authorship_markers(self) -> None:
        forbidden = (
            "/" + "home/",
            "/mnt" + "/c/",
            "OneDrive",
            "private-author-marker",
        )
        for path in RESULT.rglob("*"):
            if not path.is_file() or path.suffix in {".gz", ".csi"}:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            self.assertFalse(
                any(marker.lower() in text.lower() for marker in forbidden),
                f"private marker found in {path.relative_to(RESULT)}",
            )
