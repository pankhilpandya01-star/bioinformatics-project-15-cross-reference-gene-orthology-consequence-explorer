"""Create and run a two-reference controlled Milestone 2 demonstration."""

from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path

from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.workflow import run_workflow


ROOT = Path(__file__).parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "controlled"
INPUTS = ROOT / "local" / "controlled_core_inputs"
OUTPUT = ROOT / "results" / "controlled_core"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_identifier(source: Path, destination: Path, sequence_id: str) -> None:
    text = source.read_text(encoding="utf-8")
    text = text.replace("NC_CONTROL.1", sequence_id)
    text = text.replace("NC_CONTROL", sequence_id.rsplit(".", 1)[0])
    destination.write_text(text, encoding="utf-8", newline="\n")


def prepare_inputs() -> Path:
    if INPUTS.exists():
        return INPUTS / "reference_manifest.csv"
    INPUTS.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".controlled_core_inputs.staging-", dir=INPUTS.parent))
    specifications = (
        ("CTRL_PLUS.1", "plus-primary", "baseline", "NC_CTRL_PLUS.1", ((0, 90), (300, 345))),
        ("CTRL_MINUS.1", "minus-primary", "alternative", "NC_CTRL_MINUS.1", ((0, 45), (300, 390))),
    )
    rows: list[dict[str, str]] = []
    try:
        for accession, label, role, sequence_id, callable_intervals in specifications:
            directory = staging / accession
            directory.mkdir()
            for source_name, destination_name in (
                ("reference.fasta", "reference.fasta"),
                ("genomic.gbff", "genomic.gbff"),
                ("annotation.gff3", "annotation.gff3"),
                ("cds_from_genomic.fna", "cds_from_genomic.fna"),
                ("accepted.vcf", "accepted.vcf"),
            ):
                replace_identifier(FIXTURE / source_name, directory / destination_name, sequence_id)
            shutil.copyfile(FIXTURE / "protein.faa", directory / "protein.faa")
            with (directory / "callable_regions.bed").open("w", encoding="utf-8", newline="\n") as handle:
                for start, end in callable_intervals:
                    handle.write(f"{sequence_id}\t{start}\t{end}\n")
            subprocess.run(
                ["bcftools", "view", "-Oz", "-o", "accepted.vcf.gz", "accepted.vcf"],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                ["bcftools", "index", "--csi", "accepted.vcf.gz"],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
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
                    "cds_fasta": f"{accession}/cds_from_genomic.fna",
                    "protein_fasta": f"{accession}/protein.faa",
                    "accepted_vcf": f"{accession}/accepted.vcf.gz",
                    "callable_bed": f"{accession}/callable_regions.bed",
                    "reference_sha256": sha256(directory / "reference.fasta"),
                    "accepted_vcf_sha256": sha256(directory / "accepted.vcf.gz"),
                }
            )
        manifest = staging / "reference_manifest.csv"
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        staging.rename(INPUTS)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return INPUTS / "reference_manifest.csv"


def main() -> None:
    manifest = prepare_inputs()
    output = run_workflow(
        WorkflowConfig(
            reference_manifest=manifest,
            output_dir=OUTPUT,
            baseline_accession="CTRL_PLUS.1",
            sample_id="controlled",
        )
    )
    print(f"Controlled Milestone 2 demonstration: {output}")


if __name__ == "__main__":
    main()

