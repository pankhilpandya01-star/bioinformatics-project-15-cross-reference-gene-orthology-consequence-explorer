"""Run controlled orthology cases through the real OrthoFinder workflow."""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

from comparative_gene_explorer.orthology import OrthologyReference, run_orthology_workflow
from comparative_gene_explorer.validation import GeneRecord


ROOT = Path(__file__).parents[1]
OUTPUT = ROOT / "results" / "controlled_orthology"
WORK = ROOT / "local" / "controlled_orthology_work"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def family_sequence(seed: int, length: int = 180) -> str:
    generator = random.Random(seed)
    return "M" + "".join(generator.choice(AMINO_ACIDS) for _ in range(length - 1))


def variant(sequence: str, seed: int, changes: int = 5) -> str:
    generator = random.Random(seed)
    values = list(sequence)
    for position in generator.sample(range(1, len(values)), changes):
        alternatives = AMINO_ACIDS.replace(values[position], "")
        values[position] = generator.choice(alternatives)
    return "".join(values)


def gene(
    accession: str,
    gene_id: str,
    protein_id: str | None,
    gene_name: str | None = None,
) -> GeneRecord:
    return GeneRecord(
        accession=accession,
        gene_id=gene_id,
        gene_name=gene_name or gene_id,
        feature_id=f"gene-{gene_id}",
        feature_type="gene",
        gene_biotype="protein_coding" if protein_id else "pseudogene",
        contig=f"{accession}_contig",
        strand="+",
        start_1based=1,
        end_1based_raw=90,
        feature_intervals=((0, 90),),
        cds_intervals=((0, 90),) if protein_id else (),
        protein_ids=(protein_id,) if protein_id else (),
    )


def digest(proteins: dict[str, str]) -> str:
    payload = "".join(f">{key}\n{value}\n" for key, value in proteins.items())
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def build_references() -> tuple[OrthologyReference, OrthologyReference, OrthologyReference]:
    families = {name: family_sequence(seed) for seed, name in enumerate(
        ("one", "renamed", "one_many", "many_one", "many_many", "missing", "unassigned"),
        start=101,
    )}
    baseline_proteins = {
        "BASE_ONE.1": variant(families["one"], 1),
        "BASE_RENAMED.1": variant(families["renamed"], 2),
        "BASE_ONE_MANY.1": variant(families["one_many"], 3),
        "BASE_MANY_A.1": variant(families["many_one"], 4),
        "BASE_MANY_B.1": variant(families["many_one"], 5),
        "BASE_MM_A.1": variant(families["many_many"], 6),
        "BASE_MM_B.1": variant(families["many_many"], 7),
        "BASE_MISSING_A.1": variant(families["missing"], 8),
        "BASE_MISSING_B.1": variant(families["missing"], 9),
        "BASE_UNASSIGNED.1": families["unassigned"],
    }
    alternative_proteins = {
        "ALT_ONE.1": variant(families["one"], 11),
        "ALT_NEW_NAME.1": variant(families["renamed"], 12),
        "ALT_ONE_MANY_A.1": variant(families["one_many"], 13),
        "ALT_ONE_MANY_B.1": variant(families["one_many"], 14),
        "ALT_MANY.1": variant(families["many_one"], 15),
        "ALT_MM_A.1": variant(families["many_many"], 16),
        "ALT_MM_B.1": variant(families["many_many"], 17),
        "ALT_UNASSIGNED.1": family_sequence(999),
    }
    outgroup_proteins = {
        "OUT_ONE.1": variant(families["one"], 21, 10),
        "OUT_RENAMED.1": variant(families["renamed"], 22, 10),
        "OUT_ONE_MANY.1": variant(families["one_many"], 23, 10),
        "OUT_MANY.1": variant(families["many_one"], 24, 10),
        "OUT_MM.1": variant(families["many_many"], 25, 10),
        "OUT_MISSING.1": variant(families["missing"], 26, 10),
    }
    baseline_genes = tuple(
        gene("CTRL_BASE.1", gene_id, protein_id, gene_name)
        for gene_id, protein_id, gene_name in (
            ("base_one", "BASE_ONE.1", "stable_symbol"),
            ("base_old_name", "BASE_RENAMED.1", "old_symbol"),
            ("base_one_many", "BASE_ONE_MANY.1", None),
            ("base_many_a", "BASE_MANY_A.1", None),
            ("base_many_b", "BASE_MANY_B.1", None),
            ("base_mm_a", "BASE_MM_A.1", None),
            ("base_mm_b", "BASE_MM_B.1", None),
            ("base_missing_a", "BASE_MISSING_A.1", None),
            ("base_missing_b", "BASE_MISSING_B.1", None),
            ("base_unassigned", "BASE_UNASSIGNED.1", None),
            ("base_no_protein", None, None),
        )
    )
    alternative_genes = tuple(
        gene("CTRL_ALT.1", gene_id, protein_id, gene_name)
        for gene_id, protein_id, gene_name in (
            ("alt_one", "ALT_ONE.1", "stable_symbol"),
            ("alt_new_name", "ALT_NEW_NAME.1", "new_symbol"),
            ("alt_one_many_a", "ALT_ONE_MANY_A.1", None),
            ("alt_one_many_b", "ALT_ONE_MANY_B.1", None),
            ("alt_many", "ALT_MANY.1", None),
            ("alt_mm_a", "ALT_MM_A.1", None),
            ("alt_mm_b", "ALT_MM_B.1", None),
            ("alt_unassigned", "ALT_UNASSIGNED.1", None),
        )
    )
    outgroup_genes = tuple(
        gene("CTRL_OUT.1", gene_id, protein_id, gene_name)
        for gene_id, protein_id, gene_name in (
            ("out_one", "OUT_ONE.1", "stable_symbol"),
            ("out_renamed", "OUT_RENAMED.1", "out_symbol"),
            ("out_one_many", "OUT_ONE_MANY.1", None),
            ("out_many", "OUT_MANY.1", None),
            ("out_mm", "OUT_MM.1", None),
            ("out_missing", "OUT_MISSING.1", None),
        )
    )
    return (
        OrthologyReference(
            "CTRL_BASE.1", "controlled baseline", "baseline", baseline_genes,
            baseline_proteins, digest(baseline_proteins),
        ),
        OrthologyReference(
            "CTRL_ALT.1", "controlled alternative", "alternative", alternative_genes,
            alternative_proteins, digest(alternative_proteins),
        ),
        OrthologyReference(
            "CTRL_OUT.1", "controlled outgroup", "alternative", outgroup_genes,
            outgroup_proteins, digest(outgroup_proteins),
        ),
    )


def main() -> None:
    output = run_orthology_workflow(
        build_references(),
        baseline_accession="CTRL_BASE.1",
        threads=1,
        work_directory=WORK,
        output_directory=OUTPUT,
    )
    print(f"Controlled orthology demonstration: {output}")


if __name__ == "__main__":
    main()
