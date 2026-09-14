from __future__ import annotations

import csv
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.validation import (
    MANIFEST_FIELDS,
    ValidationError,
    parse_gff3,
    validate_inputs,
)


class CoreValidationTests(TestCase):
    def test_rejects_duplicate_manifest_accessions_before_opening_sources(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.csv"
            rows = [
                {
                    **{field: "missing" for field in MANIFEST_FIELDS},
                    "accession": "DUPLICATE.1",
                    "label": label,
                    "organism": "controlled",
                    "role": role,
                }
                for label, role in (("baseline", "baseline"), ("alternative", "alternative"))
            ]
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

            with self.assertRaisesRegex(ValidationError, "duplicate accessions"):
                validate_inputs(
                    WorkflowConfig(
                        reference_manifest=manifest,
                        output_dir=root / "output",
                        baseline_accession="DUPLICATE.1",
                    )
                )

    def test_rejects_invalid_callable_fraction(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.csv"
            manifest.write_text(",".join(MANIFEST_FIELDS) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValidationError, "callable fraction"):
                validate_inputs(
                    WorkflowConfig(
                        reference_manifest=manifest,
                        output_dir=root / "output",
                        minimum_gene_callable_fraction=1.1,
                    )
                )

    def test_parses_multicontig_strands_and_circular_origin_segments(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "annotation.gff3"
            path.write_text(
                "##gff-version 3\n"
                "a\tRefSeq\tregion\t1\t100\t.\t+\t.\tID=a;Is_circular=true\n"
                "a\tRefSeq\tgene\t90\t120\t.\t+\t.\tID=gene-a;locus_tag=a_gene\n"
                "a\tRefSeq\tCDS\t90\t120\t.\t+\t0\tID=cds-a;Parent=gene-a;locus_tag=a_gene;protein_id=A.1\n"
                "b\tRefSeq\tregion\t1\t50\t.\t+\t.\tID=b\n"
                "b\tRefSeq\tgene\t5\t34\t.\t-\t.\tID=gene-b;locus_tag=b_gene\n"
                "b\tRefSeq\tCDS\t5\t34\t.\t-\t0\tID=cds-b;Parent=gene-b;locus_tag=b_gene;protein_id=B.1\n",
                encoding="utf-8",
            )

            genes, circular = parse_gff3(
                path, "CONTROL.1", {"a": "A" * 100, "b": "C" * 50}
            )

            by_id = {gene.gene_id: gene for gene in genes}
            self.assertEqual(circular, {"a"})
            self.assertEqual(by_id["a_gene"].cds_intervals, ((0, 20), (89, 100)))
            self.assertEqual(by_id["a_gene"].strand, "+")
            self.assertEqual(by_id["b_gene"].cds_intervals, ((4, 34),))
            self.assertEqual(by_id["b_gene"].strand, "-")

