"""Manifest, sequence, VCF, annotation, and BED validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from comparative_gene_explorer.config import WorkflowConfig


MANIFEST_FIELDS = (
    "accession",
    "label",
    "organism",
    "role",
    "reference_fasta",
    "genbank_annotation",
    "gff3_annotation",
    "cds_fasta",
    "protein_fasta",
    "accepted_vcf",
    "callable_bed",
    "reference_sha256",
    "accepted_vcf_sha256",
)
SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9_.-]+$")


class ValidationError(ValueError):
    """Raised when an input cannot satisfy the workflow contract."""


@dataclass(frozen=True, slots=True)
class GeneRecord:
    accession: str
    gene_id: str
    gene_name: str
    feature_id: str
    feature_type: str
    gene_biotype: str
    contig: str
    strand: str
    start_1based: int
    end_1based_raw: int
    feature_intervals: tuple[tuple[int, int], ...]
    cds_intervals: tuple[tuple[int, int], ...]
    protein_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class VcfSummary:
    sample_id: str
    records: int
    keys: tuple[tuple[str, int, str, str], ...]


@dataclass(frozen=True, slots=True)
class ReferenceInput:
    accession: str
    label: str
    organism: str
    role: str
    reference_fasta: Path
    genbank_annotation: Path
    gff3_annotation: Path
    cds_fasta: Path
    protein_fasta: Path
    accepted_vcf: Path
    accepted_vcf_index: Path
    callable_bed: Path
    sequences: dict[str, str]
    circular_contigs: frozenset[str]
    genes: tuple[GeneRecord, ...]
    callable_intervals: dict[str, tuple[tuple[int, int], ...]]
    vcf: VcfSummary
    hashes: dict[str, str]


@dataclass(frozen=True, slots=True)
class ValidatedInputs:
    config: WorkflowConfig
    references: tuple[ReferenceInput, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def read_fasta(path: Path, *, protein: bool = False) -> dict[str, str]:
    sequences: dict[str, str] = {}
    identifier = ""
    parts: list[str] = []
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if identifier:
                    sequences[identifier] = "".join(parts).upper()
                identifier = line[1:].split()[0]
                if not identifier or identifier in sequences:
                    raise ValidationError(f"invalid or duplicate FASTA identifier on line {line_number}")
                parts = []
            else:
                if not identifier:
                    raise ValidationError(f"sequence precedes FASTA header on line {line_number}")
                sequence = line.strip().upper()
                alphabet = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-") if protein else set("ACGTNRYKMSWBDHVU*-")
                if sequence and not set(sequence) <= alphabet:
                    raise ValidationError(f"invalid FASTA character on line {line_number}")
                parts.append(sequence)
    if identifier:
        sequences[identifier] = "".join(parts).upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValidationError(f"FASTA contains no complete records: {path.name}")
    return sequences


def read_genbank_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    version = ""
    parts: list[str] = []
    in_origin = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("VERSION"):
                version = line.split()[1]
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                if not version or version in sequences:
                    raise ValidationError("GenBank records require unique VERSION identifiers")
                sequences[version] = "".join(parts).upper()
                version, parts, in_origin = "", [], False
            elif in_origin:
                parts.append("".join(character for character in line if character.isalpha()))
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValidationError(f"GenBank file contains no complete records: {path.name}")
    return sequences


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValidationError(f"malformed GFF3 attribute: {item}")
        key, value = item.split("=", 1)
        attributes[unquote(key)] = unquote(value)
    return attributes


def _normalized_intervals(
    start: int, end: int, length: int, circular: bool
) -> tuple[tuple[int, int], ...]:
    if start < 1 or end < start or start > length:
        raise ValidationError("GFF3 feature has invalid coordinates")
    if end <= length:
        return ((start - 1, end),)
    # NCBI linearizes an origin-spanning feature as L-n..L+m.
    if not circular or end - length > length:
        raise ValidationError("GFF3 feature exceeds a non-circular sequence")
    return ((start - 1, length), (0, end - length))


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    if not intervals:
        return ()
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return tuple((start, end) for start, end in merged)


def parse_gff3(path: Path, accession: str, sequences: dict[str, str]) -> tuple[tuple[GeneRecord, ...], frozenset[str]]:
    rows: list[tuple[str, str, int, int, str, dict[str, str]]] = []
    circular_contigs: set[str] = set()
    feature_ids: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValidationError(f"GFF3 line {line_number} does not contain nine columns")
            contig, _source, feature_type, start_text, end_text, _score, strand, _phase, attribute_text = fields
            if contig not in sequences:
                raise ValidationError(f"GFF3 line {line_number} uses unknown contig {contig}")
            if strand not in {"+", "-", ".", "?"}:
                raise ValidationError(f"GFF3 line {line_number} has invalid strand {strand}")
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise ValidationError(f"GFF3 line {line_number} has non-integer coordinates") from error
            attributes = parse_attributes(attribute_text)
            feature_id = attributes.get("ID", "")
            if feature_id:
                identity = (feature_type, attributes.get("Parent", ""))
                previous_identity = feature_ids.get(feature_id)
                # GFF3 permits one discontinuous CDS to use the same ID on
                # multiple rows, but unrelated features must remain unique.
                if previous_identity is not None and not (
                    feature_type == "CDS" and previous_identity == identity
                ):
                    raise ValidationError(f"duplicate GFF3 feature ID: {feature_id}")
                feature_ids[feature_id] = identity
            if feature_type == "region" and attributes.get("Is_circular", "").lower() == "true":
                circular_contigs.add(contig)
            rows.append((contig, feature_type, start, end, strand, attributes))

    genes: dict[str, dict[str, object]] = {}
    gene_ids: set[str] = set()
    locus_to_feature: dict[str, str] = {}
    feature_parents: dict[str, tuple[str, ...]] = {}
    for contig, feature_type, start, end, strand, attributes in rows:
        feature_id = attributes.get("ID", "")
        if feature_id:
            feature_parents[feature_id] = tuple(filter(None, attributes.get("Parent", "").split(",")))
        if feature_type not in {"gene", "pseudogene"}:
            continue
        if not feature_id or feature_id in genes:
            raise ValidationError("GFF3 gene and pseudogene features require unique ID values")
        locus_tag = attributes.get("locus_tag", "")
        gene_id = locus_tag or attributes.get("gene", "") or feature_id.removeprefix("gene-")
        if not gene_id or not SAFE_IDENTIFIER.fullmatch(gene_id):
            raise ValidationError(f"invalid gene identifier: {gene_id!r}")
        if gene_id in gene_ids:
            raise ValidationError(f"duplicate gene identifier in GFF3: {gene_id}")
        gene_ids.add(gene_id)
        if locus_tag:
            if locus_tag in locus_to_feature:
                raise ValidationError(f"duplicate locus_tag in GFF3: {locus_tag}")
            locus_to_feature[locus_tag] = feature_id
        feature_intervals = _normalized_intervals(
            start, end, len(sequences[contig]), contig in circular_contigs
        )
        genes[feature_id] = {
            "gene_id": gene_id,
            "gene_name": attributes.get("gene") or attributes.get("Name") or gene_id,
            "feature_id": feature_id,
            "feature_type": feature_type,
            "gene_biotype": attributes.get("gene_biotype", ""),
            "contig": contig,
            "strand": strand,
            "start": start,
            "end": end,
            "feature_intervals": feature_intervals,
            "cds_intervals": [],
            "protein_ids": set(),
        }

    def parent_gene(feature_id: str) -> str | None:
        seen: set[str] = set()
        pending = [feature_id]
        while pending:
            current = pending.pop()
            if current in seen:
                continue
            seen.add(current)
            if current in genes:
                return current
            pending.extend(feature_parents.get(current, ()))
        return None

    for contig, feature_type, start, end, _strand, attributes in rows:
        if feature_type != "CDS":
            continue
        parent_candidates = tuple(filter(None, attributes.get("Parent", "").split(",")))
        parent = None
        for candidate in parent_candidates:
            resolved = parent_gene(candidate)
            if resolved is not None:
                parent = resolved
                break
        if parent is None and attributes.get("locus_tag"):
            parent = locus_to_feature.get(attributes["locus_tag"])
        if parent is None:
            raise ValidationError("CDS feature cannot be mapped unambiguously to a gene")
        gene = genes[parent]
        if gene["contig"] != contig:
            raise ValidationError("CDS and parent gene use different contigs")
        intervals = _normalized_intervals(
            start, end, len(sequences[contig]), contig in circular_contigs
        )
        gene["cds_intervals"].extend(intervals)  # type: ignore[union-attr]
        protein_id = attributes.get("protein_id", "")
        if protein_id:
            gene["protein_ids"].add(protein_id)  # type: ignore[union-attr]

    if not genes:
        raise ValidationError("GFF3 contains no gene or pseudogene features")
    records = tuple(
        GeneRecord(
            accession=accession,
            gene_id=str(gene["gene_id"]),
            gene_name=str(gene["gene_name"]),
            feature_id=str(gene["feature_id"]),
            feature_type=str(gene["feature_type"]),
            gene_biotype=str(gene["gene_biotype"]),
            contig=str(gene["contig"]),
            strand=str(gene["strand"]),
            start_1based=int(gene["start"]),
            end_1based_raw=int(gene["end"]),
            feature_intervals=tuple(gene["feature_intervals"]),  # type: ignore[arg-type]
            cds_intervals=merge_intervals(list(gene["cds_intervals"])),  # type: ignore[arg-type]
            protein_ids=tuple(sorted(gene["protein_ids"])),  # type: ignore[arg-type]
        )
        for gene in genes.values()
    )
    return records, frozenset(circular_contigs)


def read_callable_bed(path: Path, sequences: dict[str, str]) -> dict[str, tuple[tuple[int, int], ...]]:
    intervals: dict[str, list[tuple[int, int]]] = {contig: [] for contig in sequences}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValidationError(f"BED line {line_number} has fewer than three columns")
            contig, start_text, end_text = fields[:3]
            if contig not in sequences:
                raise ValidationError(f"BED line {line_number} uses unknown contig {contig}")
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise ValidationError(f"BED line {line_number} has non-integer coordinates") from error
            if start < 0 or end <= start or end > len(sequences[contig]):
                raise ValidationError(f"BED line {line_number} has invalid half-open coordinates")
            if intervals[contig] and start < intervals[contig][-1][1]:
                raise ValidationError(f"BED intervals overlap or are unsorted on {contig}")
            intervals[contig].append((start, end))
    return {contig: tuple(values) for contig, values in intervals.items()}


def inspect_vcf(path: Path, sequences: dict[str, str], sample_id: str) -> VcfSummary:
    samples: list[str] = []
    keys: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    contig_order = {contig: index for index, contig in enumerate(sequences)}
    previous: tuple[int, int, str, str] | None = None
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValidationError(f"VCF line {line_number} has fewer than eight columns")
                chrom, position_text, _identifier, ref, alt, _quality, filter_value = fields[:7]
                if chrom not in sequences:
                    raise ValidationError(f"VCF line {line_number} uses unknown contig {chrom}")
                if "," in alt:
                    raise ValidationError("accepted VCF contains an unsplit multiallelic record")
                try:
                    position = int(position_text)
                except ValueError as error:
                    raise ValidationError(
                        f"VCF line {line_number} has a non-integer position"
                    ) from error
                if position < 1 or position + len(ref) - 1 > len(sequences[chrom]):
                    raise ValidationError(
                        f"VCF position is outside the reference at {chrom}:{position}"
                    )
                if not ref or not alt or ref == "." or alt in {".", "*"}:
                    raise ValidationError(f"VCF line {line_number} has an empty allele")
                key = (chrom, position, ref, alt)
                if key in seen:
                    raise ValidationError(f"duplicate VCF variant key: {key}")
                seen.add(key)
                sort_key = (contig_order[chrom], position, ref, alt)
                if previous is not None and sort_key < previous:
                    raise ValidationError("accepted VCF is not sorted in reference order")
                previous = sort_key
                if filter_value != "PASS":
                    raise ValidationError("accepted VCF contains a non-PASS record")
                observed_ref = sequences[chrom][position - 1 : position - 1 + len(ref)]
                if observed_ref != ref.upper():
                    raise ValidationError(f"VCF reference mismatch at {chrom}:{position}")
                keys.append(key)
    if samples != [sample_id]:
        raise ValidationError(f"VCF sample columns are {samples}; expected [{sample_id}]")
    return VcfSummary(sample_id, len(keys), tuple(keys))


def _input_index(path: Path) -> Path:
    indexes = [candidate for suffix in (".csi", ".tbi") if (candidate := Path(f"{path}{suffix}")).is_file()]
    if len(indexes) != 1:
        raise ValidationError(f"accepted VCF requires exactly one CSI or TBI index: {path.name}")
    return indexes[0]


def _resolve(manifest: Path, value: str) -> Path:
    return (manifest.parent / value).resolve()


def validate_inputs(config: WorkflowConfig) -> ValidatedInputs:
    manifest = config.reference_manifest.resolve()
    output = config.output_dir.resolve()
    if output.exists():
        raise ValidationError(f"output directory already exists: {output}")
    if not manifest.is_file():
        raise ValidationError(f"reference manifest does not exist: {manifest}")
    if not SAFE_IDENTIFIER.fullmatch(config.baseline_accession):
        raise ValidationError("baseline accession contains unsafe characters")
    if not SAFE_IDENTIFIER.fullmatch(config.sample_id):
        raise ValidationError("sample ID contains unsafe characters")
    if not 0 < config.minimum_gene_callable_fraction <= 1:
        raise ValidationError("minimum gene callable fraction must be greater than zero and at most one")
    if not 0 <= config.maximum_database_error_rate <= 1:
        raise ValidationError("maximum database error rate must be between zero and one")
    if config.upstream_downstream_length < 0:
        raise ValidationError("upstream/downstream length cannot be negative")
    if config.threads < 1:
        raise ValidationError("threads must be at least one")

    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise ValidationError("reference manifest columns do not match the required interface")
        manifest_rows = list(reader)
    if len(manifest_rows) < 2:
        raise ValidationError("reference manifest requires a baseline and at least one alternative")
    accessions = [row["accession"] for row in manifest_rows]
    labels = [row["label"] for row in manifest_rows]
    if len(accessions) != len(set(accessions)):
        raise ValidationError("reference manifest contains duplicate accessions")
    if len(labels) != len(set(labels)):
        raise ValidationError("reference manifest contains duplicate labels")
    baseline_rows = [row for row in manifest_rows if row["role"] == "baseline"]
    if len(baseline_rows) != 1 or baseline_rows[0]["accession"] != config.baseline_accession:
        raise ValidationError("manifest must contain the selected accession as its only baseline")
    if any(row["role"] not in {"baseline", "alternative"} for row in manifest_rows):
        raise ValidationError("manifest role must be baseline or alternative")

    references: list[ReferenceInput] = []
    for row in manifest_rows:
        accession = row["accession"]
        if not SAFE_IDENTIFIER.fullmatch(accession):
            raise ValidationError(f"invalid accession: {accession}")
        paths = {field: _resolve(manifest, row[field]) for field in MANIFEST_FIELDS[4:11]}
        for field, path in paths.items():
            if not path.is_file():
                raise ValidationError(f"{accession} {field} is missing: {path}")
        reference_hash = sha256(paths["reference_fasta"])
        accepted_hash = sha256(paths["accepted_vcf"])
        if reference_hash != row["reference_sha256"]:
            raise ValidationError(f"{accession} reference checksum does not match the manifest")
        if accepted_hash != row["accepted_vcf_sha256"]:
            raise ValidationError(f"{accession} accepted VCF checksum does not match the manifest")
        sequences = read_fasta(paths["reference_fasta"])
        if read_genbank_sequences(paths["genbank_annotation"]) != sequences:
            raise ValidationError(f"{accession} GenBank sequences differ from the reference FASTA")
        genes, circular_contigs = parse_gff3(paths["gff3_annotation"], accession, sequences)
        read_fasta(paths["cds_fasta"])
        read_fasta(paths["protein_fasta"], protein=True)
        callable_intervals = read_callable_bed(paths["callable_bed"], sequences)
        vcf = inspect_vcf(paths["accepted_vcf"], sequences, config.sample_id)
        index = _input_index(paths["accepted_vcf"])
        hashes = {
            "reference_fasta": reference_hash,
            "genbank_annotation": sha256(paths["genbank_annotation"]),
            "gff3_annotation": sha256(paths["gff3_annotation"]),
            "cds_fasta": sha256(paths["cds_fasta"]),
            "protein_fasta": sha256(paths["protein_fasta"]),
            "accepted_vcf": accepted_hash,
            "accepted_vcf_index": sha256(index),
            "callable_bed": sha256(paths["callable_bed"]),
        }
        references.append(
            ReferenceInput(
                accession=accession,
                label=row["label"],
                organism=row["organism"],
                role=row["role"],
                reference_fasta=paths["reference_fasta"],
                genbank_annotation=paths["genbank_annotation"],
                gff3_annotation=paths["gff3_annotation"],
                cds_fasta=paths["cds_fasta"],
                protein_fasta=paths["protein_fasta"],
                accepted_vcf=paths["accepted_vcf"],
                accepted_vcf_index=index,
                callable_bed=paths["callable_bed"],
                sequences=sequences,
                circular_contigs=circular_contigs,
                genes=genes,
                callable_intervals=callable_intervals,
                vcf=vcf,
                hashes=hashes,
            )
        )
    return ValidatedInputs(config, tuple(references))
