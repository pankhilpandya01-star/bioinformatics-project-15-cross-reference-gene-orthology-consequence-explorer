"""Core multi-reference annotation and callable-CDS workflow."""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
from pathlib import Path

from comparative_gene_explorer import __version__
from comparative_gene_explorer.atomic import atomic_output_directory
from comparative_gene_explorer.callability import build_gene_tables
from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.external import (
    ExternalToolError,
    parse_database_check,
    run_command,
    sanitize,
    verify_tool_versions,
)
from comparative_gene_explorer.validation import (
    ReferenceInput,
    ValidationError,
    ValidatedInputs,
    sha256,
    validate_inputs,
)


BACTERIAL_CODON_TABLE = (
    "TTT/F , TTC/F , TTA/L , TTG/L+ , TCT/S , TCC/S , TCA/S , TCG/S , "
    "TAT/Y , TAC/Y , TAA/* , TAG/* , TGT/C , TGC/C , TGA/* , TGG/W , "
    "CTT/L , CTC/L , CTA/L , CTG/L+ , CCT/P , CCC/P , CCA/P , CCG/P , "
    "CAT/H , CAC/H , CAA/Q , CAG/Q , CGT/R , CGC/R , CGA/R , CGG/R , "
    "ATT/I+ , ATC/I+ , ATA/I+ , ATG/M+ , ACT/T , ACC/T , ACA/T , ACG/T , "
    "AAT/N , AAC/N , AAA/K , AAG/K , AGT/S , AGC/S , AGA/R , AGG/R , "
    "GTT/V , GTC/V , GTA/V , GTG/V+ , GCT/A , GCC/A , GCA/A , GCG/A , "
    "GAT/D , GAC/D , GAA/E , GAG/E , GGT/G , GGC/G , GGA/G , GGG/G"
)
LOCUS_TAG = re.compile(r"\[locus_tag=([^\]]+)\]")
MACHINE_COMMAND_HEADER = re.compile(
    r"^##(?:bcftools|samtools|SnpEff)[A-Za-z0-9_]*(?:Command|Cmd)="
)
IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
EFFECT_FIELDS = (
    "allele",
    "consequence",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank",
    "hgvs_c",
    "hgvs_p",
    "cdna_position",
    "cds_position",
    "protein_position",
    "distance",
    "warnings",
)


class WorkflowError(RuntimeError):
    """Raised when workflow validation or orchestration fails."""


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _copy_vcf_without_machine_commands(source: Path, destination: Path) -> None:
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt", encoding="utf-8") as input_handle, destination.open(
        "x", encoding="utf-8", newline="\n"
    ) as output_handle:
        for line in input_handle:
            if not MACHINE_COMMAND_HEADER.match(line):
                output_handle.write(line)


def _reheader_cds(source: Path, destination: Path) -> int:
    occurrences: dict[str, int] = {}
    records = 0
    with source.open("r", encoding="utf-8") as input_handle, destination.open(
        "x", encoding="utf-8", newline="\n"
    ) as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if not line.startswith(">"):
                output_handle.write(line)
                continue
            match = LOCUS_TAG.search(line)
            if match is None:
                raise ValidationError(f"CDS header on line {line_number} has no locus_tag")
            locus_tag = match.group(1)
            occurrences[locus_tag] = occurrences.get(locus_tag, 0) + 1
            suffix = "" if occurrences[locus_tag] == 1 else f".{occurrences[locus_tag]}"
            output_handle.write(f">{locus_tag}{suffix}\n")
            records += 1
    if records == 0:
        raise ValidationError("CDS FASTA contains no records")
    return records


def _write_snpeff_config(path: Path, reference: ReferenceInput, genome_id: str) -> None:
    contigs = list(reference.sequences)
    lines = [
        "# Project-local database; annotation runs without remote downloads.",
        "data.dir = ./database",
        "",
        f"codon.Bacterial_and_Plant_Plastid : {BACTERIAL_CODON_TABLE}",
        "",
        f"{genome_id}.genome : {genome_id}",
        f"{genome_id}.chromosomes : {', '.join(contigs)}",
    ]
    lines.extend(
        f"{genome_id}.{contig}.codonTable : Bacterial_and_Plant_Plastid"
        for contig in contigs
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _vcf_keys(path: Path) -> list[tuple[str, int, str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    keys: list[tuple[str, int, str, str]] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys.append((fields[0], int(fields[1]), fields[3], fields[4]))
    return keys


def _parse_annotated_vcf(
    path: Path, reference: ReferenceInput
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    variants: list[dict[str, object]] = []
    effects: list[dict[str, object]] = []
    ann_header = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=ANN,"):
                ann_header = True
            elif line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                chrom, position, _identifier, ref, alt = fields[:5]
                variant_key = f"{chrom}:{position}:{ref}:{alt}"
                info = fields[7].split(";")
                ann_values = next((value[4:] for value in info if value.startswith("ANN=")), "")
                if not ann_values:
                    raise WorkflowError(f"SnpEff emitted no ANN effect for {reference.accession}:{variant_key}")
                variant_effects: list[dict[str, object]] = []
                for effect_order, value in enumerate(ann_values.split(","), start=1):
                    pieces = value.split("|")
                    pieces.extend([""] * (len(EFFECT_FIELDS) - len(pieces)))
                    effect = {
                        "accession": reference.accession,
                        "variant_key": variant_key,
                        "contig": chrom,
                        "position": int(position),
                        "reference": ref,
                        "alternate": alt,
                        "effect_order": effect_order,
                        **dict(zip(EFFECT_FIELDS, pieces[: len(EFFECT_FIELDS)], strict=True)),
                    }
                    variant_effects.append(effect)
                    effects.append(effect)
                primary = min(
                    variant_effects,
                    key=lambda effect: (
                        IMPACT_ORDER.get(str(effect["impact"]), 99),
                        int(effect["effect_order"]),
                        str(effect["gene_id"]),
                        str(effect["feature_id"]),
                    ),
                )
                variants.append(
                    {
                        "accession": reference.accession,
                        "variant_key": variant_key,
                        "contig": chrom,
                        "position": int(position),
                        "reference": ref,
                        "alternate": alt,
                        "effect_count": len(variant_effects),
                        "primary_consequence": primary["consequence"],
                        "primary_impact": primary["impact"],
                        "primary_gene_id": primary["gene_id"],
                        "primary_gene_name": primary["gene_name"],
                        "primary_feature_id": primary["feature_id"],
                        "primary_hgvs_c": primary["hgvs_c"],
                        "primary_hgvs_p": primary["hgvs_p"],
                        "warnings": primary["warnings"],
                    }
                )
    if not ann_header:
        raise WorkflowError(f"SnpEff output for {reference.accession} lacks the ANN header")
    return variants, effects


def _portable(command: tuple[str, ...]) -> str:
    return " ".join(command)


def _enforce_database_checks(
    checks: list[dict[str, object]],
    accession: str,
    maximum_error_rate: float,
) -> list[dict[str, object]]:
    maximum_percentage = maximum_error_rate * 100
    for check in checks:
        check["accession"] = accession
        check["maximum_error_percentage"] = maximum_percentage
        check["passed"] = float(check["error_percentage"]) <= maximum_percentage
        if not check["passed"]:
            raise ValidationError(
                f"{accession} {check['check']} database validation exceeds the gate"
            )
    return checks


def _run_reference_annotation(
    reference: ReferenceInput,
    staging: Path,
    config: WorkflowConfig,
    commands: list[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    accession = reference.accession
    genome_id = f"P15_{accession.replace('.', '_')}"
    work = staging / "work" / accession
    output = staging / "per_reference" / accession
    logs = output / "logs"
    reports = output / "reports"
    vcf_dir = output / "vcf"
    database_dir = work / "database" / genome_id
    for directory in (work, logs, reports, vcf_dir, database_dir):
        directory.mkdir(parents=True, exist_ok=True)

    replacements = (
        (str(staging), "."),
        (str(reference.reference_fasta), "<REFERENCE_FASTA>"),
        (str(reference.accepted_vcf), "<ACCEPTED_VCF>"),
        (str(reference.genbank_annotation), "<GENBANK_ANNOTATION>"),
        (str(reference.cds_fasta), "<CDS_FASTA>"),
    )
    shutil.copyfile(reference.reference_fasta, work / "reference.fasta")
    shutil.copyfile(reference.genbank_annotation, database_dir / "genes.gbk")
    _copy_vcf_without_machine_commands(reference.accepted_vcf, work / "accepted.vcf")
    cds_records = _reheader_cds(reference.cds_fasta, database_dir / "cds.by_transcript_id.fa")
    _write_snpeff_config(work / "snpEff.config", reference, genome_id)

    norm_command = (
        "bcftools",
        "norm",
        "-f",
        "reference.fasta",
        "-c",
        "e",
        "-m",
        "-any",
        "-d",
        "exact",
        "-Oz",
        "-o",
        "accepted.normalized.vcf.gz",
        "accepted.vcf",
    )
    norm = run_command(
        norm_command,
        cwd=work,
        stderr_path=logs / "bcftools_norm.stderr.txt",
        replacements=replacements,
        reject_markers=(),
    )
    commands.append(f"{accession}: {_portable(norm.command)}")
    index_norm = run_command(
        ("bcftools", "index", "--csi", "accepted.normalized.vcf.gz"),
        cwd=work,
        stderr_path=logs / "bcftools_norm_index.stderr.txt",
        replacements=replacements,
        reject_markers=(),
    )
    commands.append(f"{accession}: {_portable(index_norm.command)}")
    if _vcf_keys(work / "accepted.vcf") != _vcf_keys(work / "accepted.normalized.vcf.gz"):
        raise ValidationError(f"{accession} accepted VCF is not split, unique, and normalized")

    build_command = ("snpEff", "build", "-c", "snpEff.config", "-genbank", "-v", genome_id)
    build = run_command(
        build_command,
        cwd=work,
        stdout_path=logs / "snpeff_build.stdout.txt",
        stderr_path=logs / "snpeff_build.stderr.txt",
        replacements=replacements,
    )
    commands.append(f"{accession}: {_portable(build.command)}")
    protein_check = parse_database_check(build.stdout + "\n" + build.stderr, "protein")

    cds_command = (
        "snpEff",
        "cds",
        "-c",
        "snpEff.config",
        "-v",
        genome_id,
        f"database/{genome_id}/cds.by_transcript_id.fa",
    )
    cds = run_command(
        cds_command,
        cwd=work,
        stdout_path=logs / "snpeff_cds.stdout.txt",
        stderr_path=logs / "snpeff_cds.stderr.txt",
        replacements=replacements,
    )
    commands.append(f"{accession}: {_portable(cds.command)}")
    cds_check = parse_database_check(cds.stdout + "\n" + cds.stderr, "cds")
    checks = _enforce_database_checks(
        sorted([cds_check, protein_check], key=lambda row: str(row["check"])),
        accession,
        config.maximum_database_error_rate,
    )

    annotated = vcf_dir / "accepted.annotated.vcf"
    annotate_command = (
        "snpEff",
        "-c",
        "snpEff.config",
        "-nodownload",
        "-noLog",
        "-stats",
        f"../../per_reference/{accession}/reports/snpeff_summary.html",
        "-csvStats",
        f"../../per_reference/{accession}/reports/snpeff_summary.csv",
        "-ud",
        str(config.upstream_downstream_length),
        genome_id,
        "accepted.normalized.vcf.gz",
    )
    annotation = run_command(
        annotate_command,
        cwd=work,
        stdout_path=annotated,
        stderr_path=logs / "snpeff_annotate.stderr.txt",
        raw_stdout=True,
        replacements=replacements,
    )
    commands.append(f"{accession}: {_portable(annotation.command)}")
    generated_gene_table = work / "snpEff_genes.txt"
    if generated_gene_table.is_file():
        shutil.copyfile(generated_gene_table, reports / "snpeff_genes.txt")

    variants, effects = _parse_annotated_vcf(annotated, reference)
    if len(variants) != reference.vcf.records:
        raise WorkflowError(f"{accession} annotation changed the VCF record count")

    compress_command = (
        "bcftools",
        "view",
        "--threads",
        str(config.threads),
        "-Oz",
        "-o",
        "accepted.annotated.vcf.gz",
        "accepted.annotated.vcf",
    )
    compress = run_command(
        compress_command,
        cwd=vcf_dir,
        stderr_path=logs / "bcftools_compress.stderr.txt",
        replacements=replacements,
        reject_markers=(),
    )
    commands.append(f"{accession}: {_portable(compress.command)}")
    index = run_command(
        ("bcftools", "index", "--csi", "accepted.annotated.vcf.gz"),
        cwd=vcf_dir,
        stderr_path=logs / "bcftools_index.stderr.txt",
        replacements=replacements,
        reject_markers=(),
    )
    commands.append(f"{accession}: {_portable(index.command)}")
    stats = run_command(
        ("bcftools", "stats", "accepted.annotated.vcf.gz"),
        cwd=vcf_dir,
        stdout_path=reports / "annotated.bcftools.stats.txt",
        stderr_path=logs / "bcftools_stats.stderr.txt",
        replacements=replacements,
        reject_markers=(),
    )
    commands.append(f"{accession}: {_portable(stats.command)}")
    if int(run_command(("bcftools", "index", "-n", "accepted.annotated.vcf.gz"), cwd=vcf_dir, reject_markers=()).stdout.strip()) != reference.vcf.records:
        raise WorkflowError(f"{accession} annotated VCF index count disagrees with input")

    validation = [
        {
            "accession": accession,
            "contigs": len(reference.sequences),
            "reference_bases": sum(map(len, reference.sequences.values())),
            "genes": len(reference.genes),
            "cds_records": cds_records,
            "input_variants": reference.vcf.records,
            "annotated_variants": len(variants),
            "annotation_effects": len(effects),
            "reference_genbank_identical": True,
            "normalized_keys_unchanged": True,
            "passed": True,
        }
    ]
    return variants, effects, checks + validation


def _current_hashes(reference: ReferenceInput) -> dict[str, str]:
    paths = {
        "reference_fasta": reference.reference_fasta,
        "genbank_annotation": reference.genbank_annotation,
        "gff3_annotation": reference.gff3_annotation,
        "cds_fasta": reference.cds_fasta,
        "protein_fasta": reference.protein_fasta,
        "accepted_vcf": reference.accepted_vcf,
        "accepted_vcf_index": reference.accepted_vcf_index,
        "callable_bed": reference.callable_bed,
    }
    return {name: sha256(path) for name, path in paths.items()}


def run_workflow(config: WorkflowConfig) -> Path:
    try:
        inputs: ValidatedInputs = validate_inputs(config)
        versions = verify_tool_versions(config.reference_manifest.resolve().parent)
        catalog_rows, callability_rows = build_gene_tables(
            inputs.references, config.minimum_gene_callable_fraction
        )
        commands: list[str] = []
        all_variants: list[dict[str, object]] = []
        all_effects: list[dict[str, object]] = []
        database_rows: list[dict[str, object]] = []
        validation_rows: list[dict[str, object]] = []
        per_reference_counts: list[dict[str, object]] = []

        with atomic_output_directory(config.output_dir.resolve()) as staging:
            tables = staging / "tables"
            reports = staging / "reports"
            tables.mkdir()
            reports.mkdir()
            for position, reference in enumerate(inputs.references, start=1):
                print(
                    f"[{position}/{len(inputs.references)}] {reference.accession}: annotating and validating",
                    flush=True,
                )
                variants, effects, checks_and_validation = _run_reference_annotation(
                    reference, staging, config, commands
                )
                database_rows.extend(checks_and_validation[:-1])
                validation_rows.extend(checks_and_validation[-1:])
                all_variants.extend(variants)
                all_effects.extend(effects)
                matching_callability = [
                    row for row in callability_rows if row["accession"] == reference.accession
                ]
                per_reference_counts.append(
                    {
                        "accession": reference.accession,
                        "role": reference.role,
                        "genes": len(reference.genes),
                        "coding_genes": sum(bool(gene.cds_intervals) for gene in reference.genes),
                        "callable_coding_genes": sum(
                            bool(row["meets_threshold"]) for row in matching_callability
                        ),
                        "accepted_variants": reference.vcf.records,
                        "annotation_effects": len(effects),
                    }
                )

            variant_fields = list(all_variants[0]) if all_variants else [
                "accession", "variant_key", "contig", "position", "reference", "alternate",
                "effect_count", "primary_consequence", "primary_impact", "primary_gene_id",
                "primary_gene_name", "primary_feature_id", "primary_hgvs_c", "primary_hgvs_p", "warnings",
            ]
            effect_fields = list(all_effects[0]) if all_effects else [
                "accession", "variant_key", "contig", "position", "reference", "alternate",
                "effect_order", *EFFECT_FIELDS,
            ]
            write_csv(tables / "gene_catalog.csv", list(catalog_rows[0]), catalog_rows)
            write_csv(tables / "gene_callability.csv", list(callability_rows[0]), callability_rows)
            write_csv(tables / "variant_annotations.csv", variant_fields, all_variants)
            write_csv(tables / "annotation_effects.csv", effect_fields, all_effects)
            write_csv(reports / "database_validation.csv", list(database_rows[0]), database_rows)
            write_csv(reports / "input_validation.csv", list(validation_rows[0]), validation_rows)
            write_csv(reports / "tool_versions.csv", list(versions[0]), versions)
            write_csv(reports / "per_reference_counts.csv", list(per_reference_counts[0]), per_reference_counts)

            for reference in inputs.references:
                if _current_hashes(reference) != reference.hashes:
                    raise ValidationError(f"source input changed during run: {reference.accession}")

            shutil.rmtree(staging / "work")
            manifest = {
                "workflow": "comparative-gene-explorer",
                "workflow_version": __version__,
                "milestone": "core_workflow",
                "status": "pass",
                "parameters": {
                    "baseline_accession": config.baseline_accession,
                    "sample_id": config.sample_id,
                    "minimum_gene_callable_fraction": config.minimum_gene_callable_fraction,
                    "maximum_database_error_rate": config.maximum_database_error_rate,
                    "upstream_downstream_length": config.upstream_downstream_length,
                    "threads": config.threads,
                    "codon_table": "Bacterial_and_Plant_Plastid",
                },
                "references": per_reference_counts,
                "source_checksums": {
                    reference.accession: reference.hashes for reference in inputs.references
                },
                "tools": versions,
                "commands": [sanitize(command, ((str(staging), "."),)) for command in commands],
                "counts": {
                    "references": len(inputs.references),
                    "genes": len(catalog_rows),
                    "coding_genes": sum(int(row["cds_bases"]) > 0 for row in callability_rows),
                    "callable_coding_genes": sum(bool(row["meets_threshold"]) for row in callability_rows),
                    "accepted_variants": len(all_variants),
                    "annotation_effects": len(all_effects),
                    "database_checks": len(database_rows),
                },
                "outputs": {
                    "gene_catalog_sha256": sha256(tables / "gene_catalog.csv"),
                    "gene_callability_sha256": sha256(tables / "gene_callability.csv"),
                    "variant_annotations_sha256": sha256(tables / "variant_annotations.csv"),
                    "annotation_effects_sha256": sha256(tables / "annotation_effects.csv"),
                },
                "source_inputs_unchanged": True,
                "interpretation_boundary": (
                    "Alternative references are sensitivity-analysis cohorts; annotation does not imply "
                    "coordinate-equivalent variants or proven biological effects."
                ),
            }
            (staging / "run_manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
            )
        return config.output_dir.resolve()
    except (ValidationError, ExternalToolError, FileExistsError, OSError, ValueError) as error:
        raise WorkflowError(str(error)) from error
