from __future__ import annotations

import csv
import gzip
import shutil
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

from comparative_gene_explorer.comparison import (
    _build_gene_burden,
    _select_variant_gene_assignments,
)
from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.external import (
    CommandResult,
    ExternalToolError,
    run_command,
    verify_tool_versions,
)
from comparative_gene_explorer.orthology import (
    OrthologyError,
    OrthologyReference,
    prepare_orthology,
    run_orthofinder,
)
from comparative_gene_explorer.validation import (
    MANIFEST_FIELDS,
    GeneRecord,
    ValidationError,
    inspect_vcf,
    parse_gff3,
    read_callable_bed,
    read_fasta,
    sha256,
    validate_inputs,
)
from comparative_gene_explorer.workflow import (
    WorkflowError,
    _enforce_database_checks,
    _parse_annotated_vcf,
    run_workflow,
)


ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "controlled"


def _write_vcf(path: Path, records: list[str], sample: str = "sample") -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr,length=4>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + "".join(record + "\n" for record in records),
        encoding="utf-8",
        newline="\n",
    )


def _make_manifest(root: Path) -> Path:
    rows: list[dict[str, str]] = []
    for accession, label, role in (
        ("BASE.1", "baseline", "baseline"),
        ("ALT.1", "alternative", "alternative"),
    ):
        directory = root / accession
        directory.mkdir()
        for source, destination in (
            ("reference.fasta", "reference.fasta"),
            ("genomic.gbff", "genomic.gbff"),
            ("annotation.gff3", "annotation.gff3"),
            ("cds_from_genomic.fna", "cds.fna"),
            ("protein.faa", "protein.faa"),
        ):
            shutil.copyfile(FIXTURE / source, directory / destination)
        with (FIXTURE / "accepted.vcf").open("rb") as source, gzip.open(
            directory / "accepted.vcf.gz", "wb"
        ) as destination:
            shutil.copyfileobj(source, destination)
        (directory / "accepted.vcf.gz.csi").write_bytes(b"controlled-index")
        (directory / "callable.bed").write_text(
            "NC_CONTROL.1\t0\t600\n", encoding="utf-8", newline="\n"
        )
        rows.append(
            {
                "accession": accession,
                "label": label,
                "organism": "controlled bacterium",
                "role": role,
                "reference_fasta": f"{accession}/reference.fasta",
                "genbank_annotation": f"{accession}/genomic.gbff",
                "gff3_annotation": f"{accession}/annotation.gff3",
                "cds_fasta": f"{accession}/cds.fna",
                "protein_fasta": f"{accession}/protein.faa",
                "accepted_vcf": f"{accession}/accepted.vcf.gz",
                "callable_bed": f"{accession}/callable.bed",
                "reference_sha256": sha256(directory / "reference.fasta"),
                "accepted_vcf_sha256": sha256(directory / "accepted.vcf.gz"),
            }
        )
    manifest = root / "manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def _gene(
    accession: str,
    gene_id: str,
    protein_ids: tuple[str, ...] = (),
    *,
    feature_type: str = "gene",
) -> GeneRecord:
    return GeneRecord(
        accession=accession,
        gene_id=gene_id,
        gene_name=gene_id,
        feature_id=f"gene-{gene_id}",
        feature_type=feature_type,
        gene_biotype="protein_coding" if protein_ids else "pseudogene",
        contig="chr",
        strand="+",
        start_1based=1,
        end_1based_raw=9,
        feature_intervals=((0, 9),),
        cds_intervals=((0, 9),) if protein_ids else (),
        protein_ids=protein_ids,
    )


class ManifestAndSourceQualityTests(TestCase):
    def test_valid_controlled_manifest_is_accepted(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = validate_inputs(
                WorkflowConfig(
                    reference_manifest=_make_manifest(root),
                    output_dir=root / "output",
                    baseline_accession="BASE.1",
                    sample_id="controlled",
                )
            )
            self.assertEqual(len(inputs.references), 2)

    def test_malformed_manifest_header_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.csv"
            manifest.write_text("accession,label\nBASE.1,baseline\n", encoding="utf-8")
            with self.assertRaisesRegex(ValidationError, "columns"):
                validate_inputs(
                    WorkflowConfig(manifest, root / "output", baseline_accession="BASE.1")
                )

    def test_duplicate_labels_are_rejected_before_sources_are_opened(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.csv"
            rows = [
                {
                    **{field: "missing" for field in MANIFEST_FIELDS},
                    "accession": accession,
                    "label": "duplicate",
                    "organism": "controlled",
                    "role": role,
                }
                for accession, role in (("BASE.1", "baseline"), ("ALT.1", "alternative"))
            ]
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValidationError, "duplicate labels"):
                validate_inputs(
                    WorkflowConfig(manifest, root / "output", baseline_accession="BASE.1")
                )

    def test_reference_and_genbank_sequence_mismatch_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = _make_manifest(root)
            genbank = root / "ALT.1" / "genomic.gbff"
            text = genbank.read_text(encoding="utf-8")
            self.assertIn("atggctgctg", text)
            genbank.write_text(
                text.replace("atggctgctg", "ctggctgctg", 1), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValidationError, "GenBank sequences differ"):
                validate_inputs(
                    WorkflowConfig(
                        manifest,
                        root / "output",
                        baseline_accession="BASE.1",
                        sample_id="controlled",
                    )
                )

    def test_existing_output_directory_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = _make_manifest(root)
            output = root / "output"
            output.mkdir()
            with self.assertRaisesRegex(ValidationError, "already exists"):
                validate_inputs(
                    WorkflowConfig(
                        manifest,
                        output,
                        baseline_accession="BASE.1",
                        sample_id="controlled",
                    )
                )


class SequenceAnnotationAndBedQualityTests(TestCase):
    def test_duplicate_fasta_identifier_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.fasta"
            path.write_text(">same\nACGT\n>same\nACGT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValidationError, "duplicate FASTA"):
                read_fasta(path)

    def test_duplicate_unrelated_gff_feature_identifier_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.gff3"
            path.write_text(
                "chr\tRefSeq\tgene\t1\t9\t.\t+\t.\tID=gene-a;locus_tag=a\n"
                "chr\tRefSeq\tgene\t20\t28\t.\t+\t.\tID=gene-a;locus_tag=b\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValidationError, "duplicate GFF3 feature ID"):
                parse_gff3(path, "CTRL.1", {"chr": "A" * 40})

    def test_discontinuous_cds_may_reuse_its_gff_identifier(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "split.gff3"
            path.write_text(
                "chr\tRefSeq\tgene\t1\t30\t.\t+\t.\tID=gene-a;locus_tag=a\n"
                "chr\tRefSeq\tCDS\t1\t9\t.\t+\t0\tID=cds-a;Parent=gene-a;protein_id=P.1\n"
                "chr\tRefSeq\tCDS\t22\t30\t.\t+\t0\tID=cds-a;Parent=gene-a;protein_id=P.1\n",
                encoding="utf-8",
            )
            genes, _circular = parse_gff3(path, "CTRL.1", {"chr": "A" * 40})
            self.assertEqual(genes[0].cds_intervals, ((0, 9), (21, 30)))

    def test_overlapping_gene_and_pseudogene_are_preserved(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "overlap.gff3"
            path.write_text(
                "chr\tRefSeq\tgene\t1\t30\t.\t+\t.\tID=gene-a;locus_tag=a\n"
                "chr\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=cds-a;Parent=gene-a;protein_id=P.1\n"
                "chr\tRefSeq\tpseudogene\t20\t40\t.\t-\t.\tID=gene-p;locus_tag=p;gene_biotype=pseudogene\n",
                encoding="utf-8",
            )
            genes, _circular = parse_gff3(path, "CTRL.1", {"chr": "A" * 50})
            by_id = {gene.gene_id: gene for gene in genes}
            self.assertEqual(set(by_id), {"a", "p"})
            self.assertEqual(by_id["p"].feature_type, "pseudogene")
            self.assertEqual(by_id["p"].cds_intervals, ())
            self.assertGreater(
                by_id["a"].feature_intervals[0][1],
                by_id["p"].feature_intervals[0][0],
            )

    def test_invalid_bed_coordinates_are_rejected(self) -> None:
        cases = {
            "unknown": "other\t0\t5\n",
            "negative": "chr\t-1\t5\n",
            "empty": "chr\t5\t5\n",
            "past_end": "chr\t0\t11\n",
            "overlap": "chr\t0\t6\nchr\t5\t9\n",
            "unsorted": "chr\t5\t9\nchr\t0\t4\n",
        }
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in cases.items():
                with self.subTest(name=name):
                    path = root / f"{name}.bed"
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaises(ValidationError):
                        read_callable_bed(path, {"chr": "A" * 10})


class VcfQualityTests(TestCase):
    def test_empty_callset_with_valid_sample_header_is_supported(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "empty.vcf"
            _write_vcf(path, [])
            summary = inspect_vcf(path, {"chr": "ACGT"}, "sample")
            self.assertEqual(summary.records, 0)
            self.assertEqual(summary.keys, ())

    def test_malformed_variant_records_are_rejected(self) -> None:
        cases = {
            "duplicate": [
                "chr\t1\t.\tA\tC\t60\tPASS\t.\tGT\t1",
                "chr\t1\t.\tA\tC\t60\tPASS\t.\tGT\t1",
            ],
            "multiallelic": ["chr\t1\t.\tA\tC,G\t60\tPASS\t.\tGT\t1"],
            "out_of_range": ["chr\t5\t.\tA\tC\t60\tPASS\t.\tGT\t1"],
            "reference_mismatch": ["chr\t1\t.\tC\tA\t60\tPASS\t.\tGT\t1"],
            "not_pass": ["chr\t1\t.\tA\tC\t60\tLowQual\t.\tGT\t1"],
        }
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for name, records in cases.items():
                with self.subTest(name=name):
                    path = root / f"{name}.vcf"
                    _write_vcf(path, records)
                    with self.assertRaises(ValidationError):
                        inspect_vcf(path, {"chr": "ACGT"}, "sample")

    def test_wrong_sample_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "sample.vcf"
            _write_vcf(path, [], sample="other")
            with self.assertRaisesRegex(ValidationError, "sample columns"):
                inspect_vcf(path, {"chr": "ACGT"}, "sample")

    def test_empty_annotated_vcf_is_parsed_without_fabricated_effects(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "annotated.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Effects\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                encoding="utf-8",
            )
            variants, effects = _parse_annotated_vcf(
                path, SimpleNamespace(accession="CTRL.1")
            )
            self.assertEqual(variants, [])
            self.assertEqual(effects, [])


class OrthologyAndEffectAccountingQualityTests(TestCase):
    def test_ambiguous_protein_to_gene_mapping_is_explicit(self) -> None:
        reference = OrthologyReference(
            accession="CTRL.1",
            label="controlled",
            role="baseline",
            genes=(
                _gene("CTRL.1", "a", ("P.1",)),
                _gene("CTRL.1", "b", ("P.1",)),
            ),
            proteins={"P.1": "MAA"},
            source_sha256="controlled",
        )
        alternative = OrthologyReference(
            accession="ALT.1",
            label="alternative",
            role="alternative",
            genes=(_gene("ALT.1", "c", ("Q.1",)),),
            proteins={"Q.1": "MAA"},
            source_sha256="controlled-alt",
        )
        prepared = prepare_orthology((reference, alternative))
        row = next(
            row for row in prepared.protein_rows if row["source_protein_id"] == "P.1"
        )
        self.assertEqual(row["mapping_status"], "ambiguous_gene")
        self.assertEqual(row["gene_id"], "")

    def test_sanitized_protein_identifier_collision_is_rejected(self) -> None:
        reference = OrthologyReference(
            accession="CTRL.1",
            label="controlled",
            role="baseline",
            genes=(),
            proteins={"A/B": "MAA", "A?B": "MCC"},
            source_sha256="controlled",
        )
        alternative = OrthologyReference(
            accession="ALT.1",
            label="alternative",
            role="alternative",
            genes=(),
            proteins={"Q.1": "MAA"},
            source_sha256="controlled-alt",
        )
        with self.assertRaisesRegex(OrthologyError, "collision"):
            prepare_orthology((reference, alternative))

    def test_overlapping_gene_effects_are_counted_once_per_variant_gene(self) -> None:
        effects = [
            {
                "accession": "CTRL.1",
                "variant_key": "chr:2:A:C",
                "gene_id": "a",
                "impact": "LOW",
                "effect_order": "1",
                "feature_id": "a.1",
                "consequence": "synonymous_variant",
                "warnings": "",
            },
            {
                "accession": "CTRL.1",
                "variant_key": "chr:2:A:C",
                "gene_id": "a",
                "impact": "HIGH",
                "effect_order": "2",
                "feature_id": "a.2",
                "consequence": "stop_gained",
                "warnings": "",
            },
            {
                "accession": "CTRL.1",
                "variant_key": "chr:2:A:C",
                "gene_id": "b",
                "impact": "MODERATE",
                "effect_order": "3",
                "feature_id": "b.1",
                "consequence": "missense_variant",
                "warnings": "",
            },
            {
                "accession": "CTRL.1",
                "variant_key": "chr:2:A:C",
                "gene_id": "intergenic-a-b",
                "impact": "MODIFIER",
                "effect_order": "4",
                "feature_id": "intergenic-a-b",
                "consequence": "intergenic_region",
                "warnings": "",
            },
        ]
        selected, unlinked = _select_variant_gene_assignments(
            effects, {("CTRL.1", "a"), ("CTRL.1", "b")}
        )
        self.assertEqual(len(selected), 2)
        self.assertEqual(unlinked, 1)
        self.assertEqual(selected[("CTRL.1", "chr:2:A:C", "a")]["impact"], "HIGH")
        self.assertEqual(selected[("CTRL.1", "chr:2:A:C", "b")]["impact"], "MODERATE")

    def test_empty_effect_set_produces_zero_gene_burden(self) -> None:
        catalog = [
            {
                "accession": "CTRL.1",
                "gene_id": "a",
                "gene_name": "a",
                "feature_type": "gene",
                "gene_biotype": "protein_coding",
                "contig": "chr",
                "strand": "+",
            }
        ]
        callability = [
            {
                "accession": "CTRL.1",
                "gene_id": "a",
                "cds_bases": "9",
                "callable_cds_bases": "9",
                "callable_fraction": "1.0",
                "meets_threshold": "True",
            }
        ]
        result = _build_gene_burden(
            catalog,
            callability,
            [],
            {"CTRL.1": {"label": "controlled", "role": "baseline", "project13_eligible": True}},
        )
        burdens, _index, _consequences, consequence_summary, impact_summary, unlinked = result
        self.assertEqual(burdens[0]["distinct_variant_count"], 0)
        self.assertEqual(burdens[0]["variant_density_per_1000_callable_cds"], 0.0)
        self.assertEqual(consequence_summary, [])
        self.assertEqual(impact_summary, [])
        self.assertEqual(unlinked, 0)


class ToolAndAtomicFailureQualityTests(TestCase):
    def test_database_validation_gate_is_inclusive_and_fails_above_limit(self) -> None:
        passing = _enforce_database_checks(
            [{"check": "cds", "error_percentage": 2.0}], "CTRL.1", 0.02
        )
        self.assertTrue(passing[0]["passed"])
        with self.assertRaisesRegex(ValidationError, "exceeds the gate"):
            _enforce_database_checks(
                [{"check": "protein", "error_percentage": 2.001}],
                "CTRL.1",
                0.02,
            )

    def test_success_exit_with_native_failure_marker_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ExternalToolError, "RuntimeException"):
                run_command(
                    [sys.executable, "-c", "print('RuntimeException: controlled')"],
                    cwd=Path(directory),
                )

    @patch("comparative_gene_explorer.external._conda_package_version")
    @patch("comparative_gene_explorer.external.run_command")
    def test_wrong_pinned_tool_version_is_rejected(
        self, mocked_run, mocked_package_version
    ) -> None:
        mocked_run.return_value = CommandResult(
            ("orthofinder", "--version"), "OrthoFinder:v3.1.5", ""
        )
        mocked_package_version.return_value = "0.0"
        with self.assertRaisesRegex(ExternalToolError, "does not satisfy"):
            verify_tool_versions(Path("."))

    @patch("comparative_gene_explorer.workflow._run_reference_annotation")
    @patch("comparative_gene_explorer.workflow.build_gene_tables")
    @patch("comparative_gene_explorer.workflow.verify_tool_versions")
    @patch("comparative_gene_explorer.workflow.validate_inputs")
    def test_subprocess_failure_leaves_no_partial_workflow_output(
        self, mocked_validate, mocked_versions, mocked_gene_tables, mocked_annotation
    ) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "output"
            config = WorkflowConfig(root / "manifest.csv", output)
            mocked_validate.return_value = SimpleNamespace(
                references=(SimpleNamespace(accession="CTRL.1"),)
            )
            mocked_versions.return_value = []
            mocked_gene_tables.return_value = ([{"gene": "a"}], [{"gene": "a"}])
            mocked_annotation.side_effect = ExternalToolError("controlled tool failure")
            with self.assertRaisesRegex(WorkflowError, "controlled tool failure"):
                run_workflow(config)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".output.staging-*")), [])

    @patch("comparative_gene_explorer.orthology.run_command")
    def test_orthofinder_failure_removes_staged_proteomes(self, mocked_run) -> None:
        reference = OrthologyReference(
            accession="CTRL.1",
            label="controlled",
            role="baseline",
            genes=(_gene("CTRL.1", "a", ("P.1",)),),
            proteins={"P.1": "MAA"},
            source_sha256="controlled",
        )
        alternative = OrthologyReference(
            accession="ALT.1",
            label="alternative",
            role="alternative",
            genes=(_gene("ALT.1", "b", ("Q.1",)),),
            proteins={"Q.1": "MAA"},
            source_sha256="controlled-alt",
        )
        prepared = prepare_orthology((reference, alternative))
        mocked_run.side_effect = ExternalToolError("controlled OrthoFinder failure")
        with TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / "orthofinder"
            with self.assertRaisesRegex(ExternalToolError, "controlled OrthoFinder failure"):
                run_orthofinder(prepared, destination, 1)
            self.assertFalse(destination.exists())
            self.assertEqual(list(root.glob(".orthofinder.staging-*")), [])
