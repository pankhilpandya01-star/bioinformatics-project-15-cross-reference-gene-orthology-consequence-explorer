"""Acquire and validate Project 15 annotation sources and local SnpEff databases."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_ACCESSIONS = {
    "GCF_000005845.2": {"variants": 30_973, "callable_bases": 4_199_245},
    "GCF_003697165.2": {"variants": 92_665, "callable_bases": 3_888_294},
    "GCF_036503815.1": {"variants": 125_327, "callable_bases": 3_910_557},
    "GCF_005843885.1": {"variants": 49_201, "callable_bases": 777_722},
    "GCF_000026225.1": {"variants": 15_264, "callable_bases": 239_191},
    "GCF_002900365.1": {"variants": 44_151, "callable_bases": 776_835},
}
EXPECTED_FIELDS = (
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
PACKAGE_FILES = {
    "reference_fasta": "reference.fna",
    "genbank_annotation": "annotation.gbff",
    "gff3_annotation": "annotation.gff3",
    "cds_fasta": "cds.fna",
    "protein_fasta": "protein.faa",
}
DATABASE_CHECK = re.compile(
    r"(?P<name>CDS|Protein) check:\s+\S+\s+OK:\s+(?P<ok>\d+)"
    r"(?:\s+Warnings:\s+(?P<warnings>\d+))?"
    r"\s+Not found:\s+(?P<not_found>\d+)\s+Errors:\s+(?P<errors>\d+)"
    r"\s+Error percentage:\s+(?P<percentage>[\d.]+)%"
)
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")
LOCUS_TAG = re.compile(r"\[locus_tag=([^\]]+)\]")
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


class ReviewError(RuntimeError):
    """Raised when a source, tool, or validation gate is invalid."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ReviewError(f"refusing to write empty table: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(
    command: Sequence[str],
    *,
    cwd: Path,
    log_prefix: Path | None = None,
    allow_nonzero: bool = False,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command), cwd=cwd, capture_output=True, text=True, check=False
    )
    if log_prefix is not None:
        log_prefix.parent.mkdir(parents=True, exist_ok=True)
        log_prefix.with_suffix(".stdout.txt").write_text(
            completed.stdout, encoding="utf-8", newline="\n"
        )
        log_prefix.with_suffix(".stderr.txt").write_text(
            completed.stderr, encoding="utf-8", newline="\n"
        )
    if completed.returncode and not allow_nonzero:
        tail = (completed.stdout + "\n" + completed.stderr)[-3000:].strip()
        raise ReviewError(f"{command[0]} failed with status {completed.returncode}: {tail}")
    return completed


def conda_package_version(package: str) -> str:
    metadata_dir = Path(sys.prefix) / "conda-meta"
    versions: list[str] = []
    for path in metadata_dir.glob(f"{package}-*.json"):
        metadata = json.loads(path.read_text(encoding="utf-8"))
        if metadata.get("name") == package:
            versions.append(str(metadata["version"]))
    if len(versions) != 1:
        raise ReviewError(f"expected one installed {package} package; found {versions}")
    return versions[0]


def verify_tool_versions(cwd: Path) -> list[dict[str, object]]:
    definitions = (
        ("OrthoFinder", "orthofinder", "3.1.5", ["orthofinder", "--version"], r"OrthoFinder:v([\w.]+)"),
        ("SnpEff", "snpeff", "5.4.0c", ["snpEff", "-version"], r"SnpEff\s+([\w.]+)"),
        ("BCFtools", "bcftools", "1.24", ["bcftools", "--version"], r"bcftools\s+([\w.]+)"),
        ("NCBI Datasets CLI", "ncbi-datasets-cli", "18.36.0", ["datasets", "version"], r"datasets version:\s+([\w.]+)"),
        ("OpenJDK", "openjdk", "21", ["java", "-version"], r'openjdk version "([^"]+)"'),
    )
    rows: list[dict[str, object]] = []
    for label, package, required, command, pattern in definitions:
        completed = run(command, cwd=cwd)
        version_text = ANSI_ESCAPE.sub("", completed.stdout + completed.stderr)
        match = re.search(pattern, version_text)
        if match is None:
            raise ReviewError(f"could not parse {label} CLI version")
        installed = conda_package_version(package)
        passed = installed == required if label != "OpenJDK" else installed.startswith("21.")
        if not passed:
            raise ReviewError(f"{label} package {installed} does not satisfy {required}")
        rows.append(
            {
                "tool": label,
                "required": required,
                "installed_package": installed,
                "cli_reported": match.group(1),
                "passed": True,
            }
        )

    python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    matplotlib_version = __import__("matplotlib").__version__
    for label, package, required, reported, passed in (
        ("Python", "python", "3.12", python_version, python_version.startswith("3.12.")),
        ("Matplotlib", "matplotlib", "3.11.1", matplotlib_version, matplotlib_version == "3.11.1"),
    ):
        installed = conda_package_version(package)
        if not passed:
            raise ReviewError(f"{label} version {reported} does not satisfy {required}")
        rows.append(
            {
                "tool": label,
                "required": required,
                "installed_package": installed,
                "cli_reported": reported,
                "passed": True,
            }
        )
    return rows


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != EXPECTED_FIELDS:
            raise ReviewError("reference manifest columns do not match the required interface")
        rows = list(reader)
    accessions = [row["accession"] for row in rows]
    if len(rows) != 6 or set(accessions) != set(EXPECTED_ACCESSIONS):
        raise ReviewError("manifest must contain each of the six approved references exactly once")
    if len(accessions) != len(set(accessions)):
        raise ReviewError("manifest contains duplicate accessions")
    if [row["accession"] for row in rows if row["role"] == "baseline"] != ["GCF_000005845.2"]:
        raise ReviewError("MG1655 must be the only baseline reference")
    if any(row["role"] not in {"baseline", "alternative"} for row in rows):
        raise ReviewError("manifest role must be baseline or alternative")
    return rows


def resolve_manifest_path(manifest: Path, value: str) -> Path:
    return (manifest.parent / value).resolve()


def one_match(root: Path, pattern: str) -> Path:
    matches = [path for path in root.rglob(pattern) if path.is_file()]
    if len(matches) != 1:
        raise ReviewError(f"expected one {pattern} in NCBI package; found {len(matches)}")
    return matches[0]


def one_reference_fasta(root: Path) -> Path:
    matches = [
        path
        for path in root.rglob("*_genomic.fna")
        if "cds_from_genomic" not in path.name
        and "rna_from_genomic" not in path.name
        and "_rna.fna" not in path.name
    ]
    if len(matches) != 1:
        observed = ", ".join(str(path.relative_to(root)) for path in matches)
        raise ReviewError(
            f"expected one genomic FASTA in NCBI package; found {len(matches)}: {observed}"
        )
    return matches[0]


def acquire_package(accession: str, source_root: Path, commands: list[str]) -> Path:
    commands.append(
        "datasets download genome accession "
        f"{accession} --include genome,gbff,gff3,cds,protein,seq-report "
        "--filename <PACKAGE_ZIP>"
    )
    destination = source_root / accession
    required = [destination / name for name in (*PACKAGE_FILES.values(), "sequence_report.jsonl", "assembly_data_report.jsonl", "dataset_catalog.json", "ncbi_dataset.zip")]
    if destination.is_dir():
        missing = [path.name for path in required if not path.is_file()]
        if missing:
            raise ReviewError(f"existing local package {accession} is incomplete: {missing}")
        return destination
    if destination.exists():
        raise ReviewError(f"source package target is not a directory: {destination}")

    source_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{accession}.staging-", dir=source_root))
    try:
        archive = staging / "ncbi_dataset.zip"
        command = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome,gbff,gff3,cds,protein,seq-report",
            "--filename",
            str(archive),
        ]
        run(command, cwd=staging, log_prefix=staging / "download")
        extract_root = staging / "expanded"
        with zipfile.ZipFile(archive) as package:
            package.extractall(extract_root)
        accession_root = extract_root / "ncbi_dataset" / "data" / accession
        if not accession_root.is_dir():
            raise ReviewError(f"downloaded package does not contain exact accession {accession}")
        selections = {
            "reference.fna": one_reference_fasta(accession_root),
            "annotation.gbff": one_match(accession_root, "genomic.gbff"),
            "annotation.gff3": one_match(accession_root, "genomic.gff"),
            "cds.fna": one_match(accession_root, "cds_from_genomic.fna"),
            "protein.faa": one_match(accession_root, "protein.faa"),
            "sequence_report.jsonl": one_match(accession_root, "sequence_report.jsonl"),
            "assembly_data_report.jsonl": one_match(extract_root / "ncbi_dataset" / "data", "assembly_data_report.jsonl"),
            "dataset_catalog.json": one_match(extract_root / "ncbi_dataset" / "data", "dataset_catalog.json"),
        }
        for output_name, source in selections.items():
            shutil.copyfile(source, staging / output_name)
        shutil.rmtree(extract_root)
        for temporary_log in staging.glob("download.*.txt"):
            temporary_log.unlink()
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return destination


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    identifier = ""
    parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if identifier:
                    sequences[identifier] = "".join(parts).upper()
                identifier = line[1:].split()[0]
                if not identifier or identifier in sequences:
                    raise ReviewError(f"invalid or duplicate FASTA identifier on line {line_number}")
                parts = []
            else:
                if not identifier:
                    raise ReviewError(f"sequence precedes FASTA header on line {line_number}")
                parts.append(line.strip())
    if identifier:
        sequences[identifier] = "".join(parts).upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ReviewError(f"empty FASTA content: {path}")
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
                    raise ReviewError("GenBank records require unique versioned identifiers")
                sequences[version] = "".join(parts).upper()
                version, parts, in_origin = "", [], False
            elif in_origin:
                parts.append("".join(character for character in line if character.isalpha()))
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ReviewError(f"empty GenBank sequence content: {path}")
    return sequences


def validate_gff(
    path: Path, sequences: dict[str, str]
) -> tuple[Counter[str], set[str], int]:
    counts: Counter[str] = Counter()
    observed_seqids: set[str] = set()
    circular_contigs: set[str] = set()
    origin_spanning_features = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ReviewError(f"GFF3 line {line_number} does not contain nine columns")
            seqid, _source, feature_type, start_text, end_text = fields[:5]
            if seqid not in sequences:
                raise ReviewError(f"GFF3 line {line_number} uses unknown sequence {seqid}")
            start, end = int(start_text), int(end_text)
            if feature_type == "region" and "Is_circular=true" in fields[8]:
                circular_contigs.add(seqid)
            if start < 1 or end < start or start > len(sequences[seqid]):
                raise ReviewError(f"GFF3 line {line_number} has invalid coordinates")
            if end > len(sequences[seqid]):
                # NCBI linearizes circular-origin joins by extending the end
                # coordinate beyond the sequence length (for example L-10..L+20).
                if seqid not in circular_contigs or end - len(sequences[seqid]) > len(sequences[seqid]):
                    raise ReviewError(f"GFF3 line {line_number} exceeds a non-circular sequence")
                origin_spanning_features += 1
            observed_seqids.add(seqid)
            counts[feature_type] += 1
    if not counts or not observed_seqids:
        raise ReviewError("GFF3 contains no feature rows")
    return counts, circular_contigs, origin_spanning_features


def inspect_vcf(path: Path, sequences: dict[str, str], sample_id: str) -> dict[str, object]:
    samples: list[str] = []
    records = 0
    keys: set[tuple[str, int, str, str]] = set()
    contig_order = {contig: order for order, contig in enumerate(sequences)}
    previous: tuple[int, int] | None = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ReviewError(f"VCF line {line_number} has fewer than eight columns")
                chrom, position_text, _identifier, ref, alt, _qual, filter_value = fields[:7]
                if chrom not in sequences:
                    raise ReviewError(f"VCF line {line_number} uses unknown sequence {chrom}")
                if "," in alt:
                    raise ReviewError("accepted VCF contains an unsplit multiallelic record")
                position = int(position_text)
                key = (chrom, position, ref, alt)
                if key in keys:
                    raise ReviewError(f"duplicate VCF variant key: {key}")
                keys.add(key)
                sort_key = (contig_order[chrom], position)
                if previous is not None and sort_key < previous:
                    raise ReviewError("accepted VCF is not sorted in reference order")
                previous = sort_key
                if filter_value != "PASS":
                    raise ReviewError("accepted VCF contains a non-PASS record")
                observed_ref = sequences[chrom][position - 1 : position - 1 + len(ref)]
                if observed_ref != ref.upper():
                    raise ReviewError(f"VCF reference allele mismatch at {chrom}:{position}")
                records += 1
    if samples != [sample_id]:
        raise ReviewError(f"VCF sample columns are {samples}, expected [{sample_id}]")
    return {"records": records, "sample": sample_id}


def inspect_bed(path: Path, sequences: dict[str, str]) -> dict[str, int]:
    total = 0
    intervals = 0
    previous: dict[str, tuple[int, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ReviewError(f"BED line {line_number} has fewer than three columns")
            contig, start_text, end_text = fields[:3]
            if contig not in sequences:
                raise ReviewError(f"BED line {line_number} uses unknown sequence {contig}")
            start, end = int(start_text), int(end_text)
            if start < 0 or end <= start or end > len(sequences[contig]):
                raise ReviewError(f"BED line {line_number} has invalid half-open coordinates")
            if contig in previous and start < previous[contig][1]:
                raise ReviewError(f"BED intervals overlap or are unsorted on {contig}")
            previous[contig] = (start, end)
            total += end - start
            intervals += 1
    return {"callable_intervals": intervals, "callable_bases": total}


def reheader_cds(source: Path, destination: Path) -> int:
    occurrences: Counter[str] = Counter()
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
                raise ReviewError(f"CDS header on line {line_number} has no locus_tag")
            locus_tag = match.group(1)
            occurrences[locus_tag] += 1
            suffix = "" if occurrences[locus_tag] == 1 else f".{occurrences[locus_tag]}"
            output_handle.write(f">{locus_tag}{suffix}\n")
            records += 1
    if not records:
        raise ReviewError("CDS FASTA has no records")
    return records


def write_snpeff_config(path: Path, genome_id: str, contigs: Iterable[str]) -> None:
    contig_list = list(contigs)
    lines = [
        "# Project-local database; remote lookup is disabled during annotation.",
        "data.dir = ./database",
        "",
        f"codon.Bacterial_and_Plant_Plastid : {BACTERIAL_CODON_TABLE}",
        "",
        f"{genome_id}.genome : {genome_id}",
        f"{genome_id}.chromosomes : {', '.join(contig_list)}",
    ]
    lines.extend(
        f"{genome_id}.{contig}.codonTable : Bacterial_and_Plant_Plastid"
        for contig in contig_list
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def parse_database_check(text: str, name: str) -> dict[str, object]:
    matches = [match for match in DATABASE_CHECK.finditer(text) if match.group("name").lower() == name]
    if len(matches) != 1:
        raise ReviewError(f"expected exactly one {name} database validation summary")
    match = matches[0]
    return {
        "check": name,
        "ok": int(match.group("ok")),
        "warnings": int(match.group("warnings") or 0),
        "not_found": int(match.group("not_found")),
        "errors": int(match.group("errors")),
        "error_percentage": float(match.group("percentage")),
    }


def build_database(
    row: dict[str, str],
    source: Path,
    database_root: Path,
    sequences: dict[str, str],
    maximum_error_rate: float,
    commands: list[str],
) -> tuple[list[dict[str, object]], int]:
    accession = row["accession"]
    genome_id = f"P15_{accession.replace('.', '_')}"
    commands.extend(
        [
            f"snpEff build -c <SNPEFF_CONFIG> -genbank -v {genome_id}",
            f"snpEff cds -c <SNPEFF_CONFIG> -v {genome_id} <REHEADERED_CDS>",
        ]
    )
    destination = database_root / accession
    validation_path = destination / "validation.json"
    if destination.is_dir():
        if not validation_path.is_file():
            raise ReviewError(f"existing SnpEff database is incomplete: {accession}")
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if validation.get("genbank_sha256") != sha256(source / "annotation.gbff"):
            raise ReviewError(f"existing SnpEff database has stale GenBank input: {accession}")
        if validation.get("cds_sha256") != sha256(source / "cds.fna"):
            raise ReviewError(f"existing SnpEff database has stale CDS input: {accession}")
        checks = validation.get("checks", [])
        if len(checks) != 2 or not all(bool(check.get("passed")) for check in checks):
            raise ReviewError(f"existing SnpEff database lacks both checks: {accession}")
        if any(float(check["error_percentage"]) > maximum_error_rate * 100 for check in checks):
            raise ReviewError(f"existing SnpEff database no longer satisfies the configured gate: {accession}")
        return checks, int(validation["cds_records"])

    database_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{accession}.staging-", dir=database_root))
    try:
        database_dir = staging / "database" / genome_id
        database_dir.mkdir(parents=True)
        shutil.copyfile(source / "annotation.gbff", database_dir / "genes.gbk")
        cds_records = reheader_cds(source / "cds.fna", database_dir / "cds.by_transcript_id.fa")
        config = staging / "snpEff.config"
        write_snpeff_config(config, genome_id, sequences)
        build_command = ["snpEff", "build", "-c", "snpEff.config", "-genbank", "-v", genome_id]
        build = run(build_command, cwd=staging, log_prefix=staging / "logs" / "snpeff_build")
        cds_command = [
            "snpEff",
            "cds",
            "-c",
            "snpEff.config",
            "-v",
            genome_id,
            f"database/{genome_id}/cds.by_transcript_id.fa",
        ]
        cds = run(cds_command, cwd=staging, log_prefix=staging / "logs" / "snpeff_cds")
        checks = [
            parse_database_check(build.stdout + "\n" + build.stderr, "protein"),
            parse_database_check(cds.stdout + "\n" + cds.stderr, "cds"),
        ]
        maximum_percentage = maximum_error_rate * 100
        for check in checks:
            check["maximum_error_percentage"] = maximum_percentage
            check["passed"] = float(check["error_percentage"]) <= maximum_percentage
            if not check["passed"]:
                raise ReviewError(f"{accession} {check['check']} error rate exceeds the 2% gate")
        checks.sort(key=lambda item: str(item["check"]))
        validation = {
            "accession": accession,
            "genome_id": genome_id,
            "genbank_sha256": sha256(source / "annotation.gbff"),
            "cds_sha256": sha256(source / "cds.fna"),
            "cds_records": cds_records,
            "checks": checks,
        }
        validation_path = staging / "validation.json"
        validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8", newline="\n")
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return checks, cds_records


def source_manifest_rows(accession: str, source: Path, retrieved: str) -> list[dict[str, object]]:
    roles = {
        "ncbi_package_zip": "ncbi_dataset.zip",
        **PACKAGE_FILES,
        "sequence_report": "sequence_report.jsonl",
        "assembly_data_report": "assembly_data_report.jsonl",
        "dataset_catalog": "dataset_catalog.json",
    }
    return [
        {
            "accession": accession,
            "role": role,
            "file": filename,
            "bytes": (source / filename).stat().st_size,
            "sha256": sha256(source / filename),
            "source_url": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/",
            "retrieved": retrieved,
        }
        for role, filename in roles.items()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--maximum-database-error-rate", type=float, default=0.02)
    parser.add_argument("--retrieved", required=True)
    args = parser.parse_args()

    manifest = args.reference_manifest.resolve()
    source_root = args.source_dir.resolve()
    database_root = args.database_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise ReviewError(f"output directory already exists: {output}")
    if not 0 <= args.maximum_database_error_rate <= 1:
        raise ReviewError("maximum database error rate must be between zero and one")

    rows = load_manifest(manifest)
    tool_rows = verify_tool_versions(manifest.parent)
    commands: list[str] = []
    acquisition_rows: list[dict[str, object]] = []
    validation_rows: list[dict[str, object]] = []
    database_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    inherited_hashes: dict[Path, str] = {}

    for position, row in enumerate(rows, start=1):
        accession = row["accession"]
        print(f"[{position}/6] {accession}: acquiring exact NCBI package", flush=True)
        source = acquire_package(accession, source_root, commands)
        acquisition_rows.extend(source_manifest_rows(accession, source, args.retrieved))
        for field, filename in PACKAGE_FILES.items():
            expected = resolve_manifest_path(manifest, row[field])
            observed = source / filename
            if expected != observed:
                raise ReviewError(f"manifest {field} path does not resolve to the acquired {accession} file")

        reference = source / "reference.fna"
        if sha256(reference) != row["reference_sha256"]:
            raise ReviewError(f"{accession} reference checksum differs from Project 13")
        sequences = read_fasta(reference)
        genbank_sequences = read_genbank_sequences(source / "annotation.gbff")
        if genbank_sequences != sequences:
            raise ReviewError(f"{accession} GenBank sequences differ from the reference FASTA")
        counts, circular_contigs, origin_spanning_features = validate_gff(
            source / "annotation.gff3", sequences
        )

        accepted_vcf = resolve_manifest_path(manifest, row["accepted_vcf"])
        callable_bed = resolve_manifest_path(manifest, row["callable_bed"])
        for inherited in (accepted_vcf, callable_bed):
            if not inherited.is_file():
                raise ReviewError(f"missing inherited Project 13 input: {inherited}")
            inherited_hashes[inherited] = sha256(inherited)
        if inherited_hashes[accepted_vcf] != row["accepted_vcf_sha256"]:
            raise ReviewError(f"{accession} accepted VCF checksum differs from Project 13")

        vcf = inspect_vcf(accepted_vcf, sequences, args.sample_id)
        bed = inspect_bed(callable_bed, sequences)
        expected = EXPECTED_ACCESSIONS[accession]
        if vcf["records"] != expected["variants"]:
            raise ReviewError(f"{accession} accepted VCF count does not reproduce Project 13")
        if bed["callable_bases"] != expected["callable_bases"]:
            raise ReviewError(f"{accession} callable-base count does not reproduce Project 13")
        print(f"[{position}/6] {accession}: source, VCF, and callable BED validated", flush=True)
        index_count = run(["bcftools", "index", "-n", str(accepted_vcf)], cwd=manifest.parent)
        if int(index_count.stdout.strip()) != vcf["records"]:
            raise ReviewError(f"{accession} VCF index count disagrees with the VCF")

        checks, cds_records = build_database(
            row,
            source,
            database_root,
            sequences,
            args.maximum_database_error_rate,
            commands,
        )
        print(f"[{position}/6] {accession}: SnpEff CDS/protein gates passed", flush=True)
        for check in checks:
            database_rows.append({"accession": accession, **check})
        for feature_type, count in sorted(counts.items()):
            feature_rows.append({"accession": accession, "feature_type": feature_type, "count": count})
        validation_rows.append(
            {
                "accession": accession,
                "label": row["label"],
                "role": row["role"],
                "contigs": len(sequences),
                "reference_bases": sum(map(len, sequences.values())),
                "reference_sha256": sha256(reference),
                "genbank_sequence_identical": True,
                "gff3_coordinates_valid": True,
                "circular_contigs": len(circular_contigs),
                "origin_spanning_features": origin_spanning_features,
                "cds_records": cds_records,
                "protein_records": len(read_fasta(source / "protein.faa")),
                "accepted_variants": vcf["records"],
                "accepted_vcf_sha256": inherited_hashes[accepted_vcf],
                "callable_intervals": bed["callable_intervals"],
                "callable_bases": bed["callable_bases"],
                "callable_bed_sha256": inherited_hashes[callable_bed],
                "passed": True,
            }
        )

    for path, before in inherited_hashes.items():
        if sha256(path) != before:
            raise ReviewError(f"inherited input changed during review: {path.name}")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging-", dir=output.parent))
    try:
        write_csv(staging / "tool_versions.csv", tool_rows)
        write_csv(staging / "acquisition_manifest.csv", acquisition_rows)
        write_csv(staging / "source_validation.csv", validation_rows)
        write_csv(staging / "database_validation.csv", database_rows)
        write_csv(staging / "feature_counts.csv", feature_rows)
        (staging / "portable_commands.txt").write_text(
            "\n".join(commands) + "\n", encoding="utf-8", newline="\n"
        )
        summary = {
            "workflow": "project-15-source-review",
            "status": "pass",
            "sample_id": args.sample_id,
            "baseline_accession": "GCF_000005845.2",
            "references": len(validation_rows),
            "source_checks_passed": len(validation_rows),
            "database_checks_passed": sum(bool(row["passed"]) for row in database_rows),
            "maximum_database_error_rate": args.maximum_database_error_rate,
            "reference_sequence_identity": "FASTA and GenBank sequences match for all references",
            "inherited_inputs_unchanged": True,
            "comparison_boundary": (
                "MG1655 is primary; five alternatives are sensitivity-analysis cohorts. "
                "Gene-level orthology will not imply coordinate-equivalent variants."
            ),
        }
        (staging / "review_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )
        staging.rename(output)
        print(f"Published source-review evidence to {output}", flush=True)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
