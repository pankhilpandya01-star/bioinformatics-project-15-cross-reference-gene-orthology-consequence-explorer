"""Protein preparation, OrthoFinder execution, parsing, and relationship classification."""

from __future__ import annotations

import csv
import json
import re
import shutil
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from comparative_gene_explorer.atomic import atomic_output_directory
from comparative_gene_explorer.external import run_command, verify_tool_versions
from comparative_gene_explorer.validation import (
    GeneRecord,
    ReferenceInput,
    ValidationError,
    read_fasta,
    sha256,
)


SAFE_TOKEN = re.compile(r"[^A-Za-z0-9_.-]+")


class OrthologyError(RuntimeError):
    """Raised when orthology inputs or OrthoFinder outputs are inconsistent."""


@dataclass(frozen=True, slots=True)
class OrthologyReference:
    accession: str
    label: str
    role: str
    genes: tuple[GeneRecord, ...]
    proteins: dict[str, str]
    source_sha256: str


@dataclass(frozen=True, slots=True)
class PreparedOrthology:
    references: tuple[OrthologyReference, ...]
    protein_rows: tuple[dict[str, object], ...]
    prepared_sequences: dict[str, dict[str, str]]
    species_to_accession: dict[str, str]


def species_key(accession: str) -> str:
    key = SAFE_TOKEN.sub("_", accession).strip("_")
    if not key:
        raise OrthologyError(f"cannot create a safe species key for {accession!r}")
    return key


def prepared_protein_id(accession: str, source_protein_id: str) -> str:
    protein = SAFE_TOKEN.sub("_", source_protein_id).strip("_")
    if not protein:
        raise OrthologyError(f"cannot create a safe protein ID for {source_protein_id!r}")
    return f"{species_key(accession)}__{protein}"


def references_from_validated(
    references: tuple[ReferenceInput, ...]
) -> tuple[OrthologyReference, ...]:
    return tuple(
        OrthologyReference(
            accession=reference.accession,
            label=reference.label,
            role=reference.role,
            genes=reference.genes,
            proteins=read_fasta(reference.protein_fasta, protein=True),
            source_sha256=reference.hashes["protein_fasta"],
        )
        for reference in references
    )


def prepare_orthology(references: tuple[OrthologyReference, ...]) -> PreparedOrthology:
    if len(references) < 2:
        raise OrthologyError("OrthoFinder requires at least two reference proteomes")
    accessions = [reference.accession for reference in references]
    if len(accessions) != len(set(accessions)):
        raise OrthologyError("orthology references contain duplicate accessions")

    protein_rows: list[dict[str, object]] = []
    prepared_sequences: dict[str, dict[str, str]] = {}
    species_to_accession: dict[str, str] = {}
    seen_prepared: set[str] = set()
    for reference in references:
        key = species_key(reference.accession)
        if key in species_to_accession:
            raise OrthologyError(f"species-key collision: {key}")
        species_to_accession[key] = reference.accession
        protein_to_genes: dict[str, list[GeneRecord]] = {}
        for gene in reference.genes:
            for protein_id in gene.protein_ids:
                protein_to_genes.setdefault(protein_id, []).append(gene)

        prepared_sequences[key] = {}
        for source_id, sequence in reference.proteins.items():
            prepared_id = prepared_protein_id(reference.accession, source_id)
            if prepared_id in seen_prepared:
                raise OrthologyError(f"prepared protein identifier collision: {prepared_id}")
            seen_prepared.add(prepared_id)
            prepared_sequences[key][prepared_id] = sequence
            genes = protein_to_genes.get(source_id, [])
            if len(genes) == 1:
                gene = genes[0]
                mapping_status = "mapped"
                gene_id, gene_name, feature_id = gene.gene_id, gene.gene_name, gene.feature_id
            elif not genes:
                mapping_status = "no_gene"
                gene_id = gene_name = feature_id = ""
            else:
                mapping_status = "ambiguous_gene"
                gene_id = gene_name = feature_id = ""
            protein_rows.append(
                {
                    "accession": reference.accession,
                    "species_key": key,
                    "source_protein_id": source_id,
                    "prepared_protein_id": prepared_id,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "feature_id": feature_id,
                    "mapping_status": mapping_status,
                }
            )

        for source_id, genes in sorted(protein_to_genes.items()):
            if source_id not in reference.proteins:
                gene = genes[0] if len(genes) == 1 else None
                protein_rows.append(
                    {
                        "accession": reference.accession,
                        "species_key": key,
                        "source_protein_id": source_id,
                        "prepared_protein_id": "",
                        "gene_id": gene.gene_id if gene else "",
                        "gene_name": gene.gene_name if gene else "",
                        "feature_id": gene.feature_id if gene else "",
                        "mapping_status": "missing_sequence" if gene else "ambiguous_gene",
                    }
                )
        for gene in reference.genes:
            if not gene.protein_ids:
                protein_rows.append(
                    {
                        "accession": reference.accession,
                        "species_key": key,
                        "source_protein_id": "",
                        "prepared_protein_id": "",
                        "gene_id": gene.gene_id,
                        "gene_name": gene.gene_name,
                        "feature_id": gene.feature_id,
                        "mapping_status": "no_protein",
                    }
                )

    return PreparedOrthology(
        references,
        tuple(protein_rows),
        prepared_sequences,
        species_to_accession,
    )


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in sequences.items():
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def _source_fingerprint(prepared: PreparedOrthology) -> dict[str, str]:
    return {
        reference.accession: reference.source_sha256 for reference in prepared.references
    }


def run_orthofinder(
    prepared: PreparedOrthology,
    work_directory: Path,
    threads: int,
) -> tuple[Path, str]:
    destination = work_directory.resolve()
    marker = destination / "completed.json"
    fingerprint = _source_fingerprint(prepared)
    command_text = (
        "orthofinder -f proteomes -M msa -S diamond -A famsa -T fasttree "
        f"-t {threads} -a {threads} -X -o results"
    )
    if destination.is_dir():
        if not marker.is_file():
            raise OrthologyError(f"existing OrthoFinder work directory is incomplete: {destination}")
        metadata = json.loads(marker.read_text(encoding="utf-8"))
        if metadata.get("protein_source_sha256") != fingerprint:
            raise OrthologyError("existing OrthoFinder work directory has stale protein inputs")
        if metadata.get("command") != command_text:
            raise OrthologyError("existing OrthoFinder work directory used different parameters")
        return _find_results_root(destination), command_text
    if destination.exists():
        raise OrthologyError(f"OrthoFinder work target is not a directory: {destination}")
    if threads < 1:
        raise OrthologyError("OrthoFinder threads must be at least one")

    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent)
    )
    try:
        proteomes = staging / "proteomes"
        proteomes.mkdir()
        for key, sequences in prepared.prepared_sequences.items():
            _write_fasta(proteomes / f"{key}.faa", sequences)
        command = (
            "orthofinder",
            "-f",
            "proteomes",
            "-M",
            "msa",
            "-S",
            "diamond",
            "-A",
            "famsa",
            "-T",
            "fasttree",
            "-t",
            str(threads),
            "-a",
            str(threads),
            "-X",
            "-o",
            "results",
        )
        run_command(
            command,
            cwd=staging,
            stdout_path=staging / "orthofinder.stdout.txt",
            stderr_path=staging / "orthofinder.stderr.txt",
            replacements=((str(staging), "."),),
            reject_markers=(),
        )
        _find_results_root(staging)
        marker_data = {
            "status": "complete",
            "command": command_text,
            "protein_source_sha256": fingerprint,
            "prepared_proteins": sum(map(len, prepared.prepared_sequences.values())),
        }
        (staging / "completed.json").write_text(
            json.dumps(marker_data, indent=2) + "\n", encoding="utf-8", newline="\n"
        )
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return _find_results_root(destination), command_text


def _find_results_root(work_directory: Path) -> Path:
    direct = work_directory / "results"
    if (direct / "Orthogroups").is_dir():
        return direct
    candidates = [
        path.parent
        for path in direct.rglob("Orthogroups")
        if path.is_dir()
    ]
    if len(candidates) != 1:
        raise OrthologyError(f"expected one complete OrthoFinder result; found {len(candidates)}")
    return candidates[0]


def _split_members(value: str) -> list[str]:
    return [member.strip() for member in value.split(",") if member.strip()]


def _accession_for_column(column: str, species_map: dict[str, str]) -> str:
    normalized = column.removesuffix(".faa")
    if normalized not in species_map:
        raise OrthologyError(f"unknown OrthoFinder species column: {column}")
    return species_map[normalized]


def parse_orthogroups(
    results: Path,
    prepared: PreparedOrthology,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, str]]:
    mapped = {
        str(row["prepared_protein_id"]): row
        for row in prepared.protein_rows
        if row["prepared_protein_id"]
    }
    catalog: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    protein_to_group: dict[str, str] = {}

    def consume(path: Path, assignment_status: str) -> None:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
                raise OrthologyError(f"unexpected OrthoFinder table header: {path.name}")
            for row in reader:
                orthogroup = row["Orthogroup"]
                group_members: list[tuple[str, str]] = []
                for column in reader.fieldnames[1:]:
                    accession = _accession_for_column(column, prepared.species_to_accession)
                    group_members.extend((accession, member) for member in _split_members(row[column]))
                if not group_members:
                    raise OrthologyError(f"empty OrthoFinder group: {orthogroup}")
                catalog.append(
                    {
                        "orthogroup": orthogroup,
                        "assignment_status": assignment_status,
                        "protein_count": len(group_members),
                        "species_count": len({accession for accession, _member in group_members}),
                    }
                )
                for accession, prepared_id in group_members:
                    if prepared_id not in mapped:
                        raise OrthologyError(f"unknown prepared protein in OrthoFinder output: {prepared_id}")
                    if prepared_id in protein_to_group:
                        raise OrthologyError(f"protein occurs in multiple OrthoFinder groups: {prepared_id}")
                    protein_to_group[prepared_id] = orthogroup
                    protein = mapped[prepared_id]
                    members.append(
                        {
                            "orthogroup": orthogroup,
                            "assignment_status": assignment_status,
                            "accession": accession,
                            "prepared_protein_id": prepared_id,
                            "source_protein_id": protein["source_protein_id"],
                            "gene_id": protein["gene_id"],
                            "gene_name": protein["gene_name"],
                            "mapping_status": protein["mapping_status"],
                        }
                    )

    consume(results / "Orthogroups" / "Orthogroups.tsv", "orthogroup")
    unassigned = results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
    if unassigned.is_file() and unassigned.stat().st_size:
        consume(unassigned, "unassigned")
    expected = set(mapped)
    observed = set(protein_to_group)
    if observed != expected:
        missing = len(expected - observed)
        unexpected = len(observed - expected)
        raise OrthologyError(
            f"OrthoFinder protein accounting differs from prepared proteomes: missing={missing}, unexpected={unexpected}"
        )
    return catalog, members, protein_to_group


def parse_pairwise_orthologs(
    results: Path,
    prepared: PreparedOrthology,
    protein_to_group: dict[str, str],
) -> list[dict[str, object]]:
    mapped = {
        str(row["prepared_protein_id"]): row
        for row in prepared.protein_rows
        if row["prepared_protein_id"]
    }
    order = {reference.accession: index for index, reference in enumerate(prepared.references)}
    edges: dict[tuple[str, str], dict[str, object]] = {}
    orthologues_root = results / "Orthologues"
    if not orthologues_root.is_dir():
        raise OrthologyError("OrthoFinder result is missing pairwise orthologue tables")
    pairwise_paths = sorted(orthologues_root.glob("Orthologues_*/*__v__*.tsv"))
    if not pairwise_paths:
        raise OrthologyError("OrthoFinder result contains no pairwise orthologue tables")
    for path in pairwise_paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or len(reader.fieldnames) != 3:
                raise OrthologyError(f"unexpected pairwise orthologue table: {path.name}")
            left_column, right_column = reader.fieldnames[1:]
            left_accession = _accession_for_column(left_column, prepared.species_to_accession)
            right_accession = _accession_for_column(right_column, prepared.species_to_accession)
            for row in reader:
                for left in _split_members(row[left_column]):
                    for right in _split_members(row[right_column]):
                        if left not in mapped or right not in mapped:
                            raise OrthologyError("pairwise table contains an unknown prepared protein")
                        if order[left_accession] <= order[right_accession]:
                            protein_a, protein_b = left, right
                            accession_a, accession_b = left_accession, right_accession
                        else:
                            protein_a, protein_b = right, left
                            accession_a, accession_b = right_accession, left_accession
                        key = (protein_a, protein_b)
                        group_a = protein_to_group.get(protein_a, "")
                        group_b = protein_to_group.get(protein_b, "")
                        if group_a != group_b or (row[reader.fieldnames[0]] and row[reader.fieldnames[0]] != group_a):
                            raise OrthologyError("pairwise orthologue and orthogroup identifiers disagree")
                        map_a, map_b = mapped[protein_a], mapped[protein_b]
                        edges[key] = {
                            "orthogroup": group_a,
                            "accession_a": accession_a,
                            "protein_a": protein_a,
                            "gene_a": map_a["gene_id"],
                            "accession_b": accession_b,
                            "protein_b": protein_b,
                            "gene_b": map_b["gene_id"],
                        }
    return sorted(
        edges.values(),
        key=lambda row: (
            order[str(row["accession_a"])],
            order[str(row["accession_b"])],
            str(row["orthogroup"]),
            str(row["protein_a"]),
            str(row["protein_b"]),
        ),
    )


def classify_baseline_relationships(
    prepared: PreparedOrthology,
    members: list[dict[str, object]],
    pairwise: list[dict[str, object]],
    baseline_accession: str,
) -> list[dict[str, object]]:
    reference_by_accession = {reference.accession: reference for reference in prepared.references}
    if baseline_accession not in reference_by_accession:
        raise OrthologyError(f"baseline accession is absent: {baseline_accession}")
    alternatives = [
        reference for reference in prepared.references if reference.accession != baseline_accession
    ]
    mapped_rows = [row for row in prepared.protein_rows if row["mapping_status"] == "mapped"]
    proteins_by_gene: dict[tuple[str, str], set[str]] = {}
    for row in mapped_rows:
        proteins_by_gene.setdefault(
            (str(row["accession"]), str(row["gene_id"])), set()
        ).add(str(row["prepared_protein_id"]))
    group_by_protein = {
        str(row["prepared_protein_id"]): str(row["orthogroup"]) for row in members
    }
    mapping_by_protein = {
        str(row["prepared_protein_id"]): row
        for row in mapped_rows
    }
    # Index the undirected pairwise graph once; authentic comparisons contain
    # thousands of genes, so rescanning every edge for every gene is avoidable.
    neighbours: dict[str, set[str]] = defaultdict(set)
    for row in pairwise:
        protein_a = str(row["protein_a"])
        protein_b = str(row["protein_b"])
        neighbours[protein_a].add(protein_b)
        neighbours[protein_b].add(protein_a)
    proteins_by_accession_group: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in members:
        if row["mapping_status"] == "mapped":
            proteins_by_accession_group[
                (str(row["accession"]), str(row["orthogroup"]))
            ].add(str(row["prepared_protein_id"]))

    baseline = reference_by_accession[baseline_accession]
    rows: list[dict[str, object]] = []
    for alternative in alternatives:
        for gene in baseline.genes:
            baseline_proteins = proteins_by_gene.get((baseline_accession, gene.gene_id), set())
            linked_target_proteins = {
                target
                for source in baseline_proteins
                for target in neighbours.get(source, set())
                if mapping_by_protein.get(target, {}).get("accession") == alternative.accession
            }
            if not baseline_proteins:
                relationship = "no_protein"
                target_genes: set[str] = set()
                source_genes = {gene.gene_id}
                groups: set[str] = set()
            elif not linked_target_proteins:
                relationship = "no_ortholog"
                target_genes = set()
                source_genes = {gene.gene_id}
                groups = {group_by_protein[p] for p in baseline_proteins if p in group_by_protein}
            else:
                target_genes = {
                    str(mapping_by_protein[protein]["gene_id"])
                    for protein in linked_target_proteins
                }
                groups = {
                    group_by_protein[protein]
                    for protein in baseline_proteins | linked_target_proteins
                    if protein in group_by_protein
                }
                family_target_proteins = set().union(
                    *(
                        proteins_by_accession_group.get((alternative.accession, group), set())
                        for group in groups
                    )
                )
                source_proteins = {
                    source
                    for target in family_target_proteins
                    for source in neighbours.get(target, set())
                    if mapping_by_protein.get(source, {}).get("accession") == baseline_accession
                }
                source_genes = {
                    str(mapping_by_protein[protein]["gene_id"]) for protein in source_proteins
                }
                source_count, target_count = len(source_genes), len(target_genes)
                if source_count == 1 and target_count == 1:
                    relationship = "one_to_one"
                elif source_count == 1 and target_count > 1:
                    relationship = "one_to_many"
                elif source_count > 1 and target_count == 1:
                    relationship = "many_to_one"
                else:
                    relationship = "many_to_many"

            target_gene_lookup = {item.gene_id: item for item in alternative.genes}
            target_names = {
                target_gene_lookup[target].gene_name
                for target in target_genes
                if target in target_gene_lookup
            }
            renamed = (
                relationship == "one_to_one"
                and len(target_genes) == 1
                and target_names != {gene.gene_name}
            )
            rows.append(
                {
                    "baseline_accession": baseline_accession,
                    "baseline_gene_id": gene.gene_id,
                    "baseline_gene_name": gene.gene_name,
                    "alternative_accession": alternative.accession,
                    "relationship": relationship,
                    "baseline_family_gene_count": len(source_genes),
                    "alternative_family_gene_count": len(target_genes),
                    "orthogroups": ";".join(sorted(groups)),
                    "alternative_gene_ids": ";".join(sorted(target_genes)),
                    "alternative_gene_names": ";".join(sorted(target_names)),
                    "renamed_one_to_one": renamed,
                }
            )
    return rows


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _file_digest(path: Path) -> dict[str, object]:
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def publish_orthology_results(
    prepared: PreparedOrthology,
    results_root: Path,
    output_directory: Path,
    baseline_accession: str,
    threads: int,
    command_text: str,
) -> Path:
    catalog, members, protein_to_group = parse_orthogroups(results_root, prepared)
    pairwise = parse_pairwise_orthologs(results_root, prepared, protein_to_group)
    relationships = classify_baseline_relationships(
        prepared, members, pairwise, baseline_accession
    )
    relationship_counts: dict[tuple[str, str], int] = {}
    for row in relationships:
        key = (str(row["alternative_accession"]), str(row["relationship"]))
        relationship_counts[key] = relationship_counts.get(key, 0) + 1
    summary_rows = [
        {"alternative_accession": accession, "relationship": relationship, "gene_count": count}
        for (accession, relationship), count in sorted(relationship_counts.items())
    ]
    source_files = [
        results_root / "Orthogroups" / "Orthogroups.tsv",
        results_root / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv",
    ]
    source_manifest = [_file_digest(path) for path in source_files if path.is_file()]

    with atomic_output_directory(output_directory.resolve()) as staging:
        tables = staging / "tables"
        reports = staging / "reports"
        tables.mkdir()
        reports.mkdir()
        protein_fields = [
            "accession", "species_key", "source_protein_id", "prepared_protein_id",
            "gene_id", "gene_name", "feature_id", "mapping_status",
        ]
        _write_csv(tables / "protein_gene_map.csv", list(prepared.protein_rows), protein_fields)
        _write_csv(
            tables / "orthogroup_catalog.csv",
            catalog,
            ["orthogroup", "assignment_status", "protein_count", "species_count"],
        )
        _write_csv(
            tables / "orthogroup_members.csv",
            members,
            [
                "orthogroup", "assignment_status", "accession", "prepared_protein_id",
                "source_protein_id", "gene_id", "gene_name", "mapping_status",
            ],
        )
        _write_csv(
            tables / "pairwise_orthologs.csv",
            pairwise,
            [
                "orthogroup", "accession_a", "protein_a", "gene_a",
                "accession_b", "protein_b", "gene_b",
            ],
        )
        _write_csv(
            tables / "orthology_relationships.csv",
            relationships,
            [
                "baseline_accession", "baseline_gene_id", "baseline_gene_name",
                "alternative_accession", "relationship", "baseline_family_gene_count",
                "alternative_family_gene_count", "orthogroups", "alternative_gene_ids",
                "alternative_gene_names", "renamed_one_to_one",
            ],
        )
        _write_csv(
            reports / "orthology_relationship_summary.csv",
            summary_rows,
            ["alternative_accession", "relationship", "gene_count"],
        )
        _write_csv(
            reports / "orthofinder_source_files.csv",
            source_manifest,
            ["file", "bytes", "sha256"],
        )
        versions = verify_tool_versions(staging)
        _write_csv(
            reports / "tool_versions.csv",
            versions,
            ["tool", "required", "installed_package", "cli_reported", "passed"],
        )
        total_prepared = sum(map(len, prepared.prepared_sequences.values()))
        if len(members) != total_prepared or len(protein_to_group) != total_prepared:
            raise OrthologyError("published orthogroup membership does not account for every protein")
        manifest = {
            "workflow": "comparative-gene-explorer",
            "milestone": "orthology_workflow",
            "status": "pass",
            "baseline_accession": baseline_accession,
            "parameters": {
                "method": "msa",
                "sequence_search": "diamond",
                "msa_program": "famsa",
                "tree_program": "fasttree",
                "threads": threads,
                "preserve_prepared_ids": True,
            },
            "command": command_text,
            "protein_source_sha256": _source_fingerprint(prepared),
            "counts": {
                "references": len(prepared.references),
                "prepared_proteins": total_prepared,
                "protein_map_rows": len(prepared.protein_rows),
                "orthogroups_and_unassigned": len(catalog),
                "orthogroup_members": len(members),
                "pairwise_protein_orthologs": len(pairwise),
                "baseline_alternative_gene_comparisons": len(relationships),
            },
            "protein_accounting_complete": True,
            "interpretation_boundary": (
                "Orthogroup membership suggests evolutionary relationship; it does not prove "
                "identical function or coordinate-equivalent variants."
            ),
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
        )
    return output_directory.resolve()


def run_orthology_workflow(
    references: tuple[OrthologyReference, ...],
    *,
    baseline_accession: str,
    threads: int,
    work_directory: Path,
    output_directory: Path,
) -> Path:
    prepared = prepare_orthology(references)
    results, command = run_orthofinder(prepared, work_directory, threads)
    return publish_orthology_results(
        prepared, results, output_directory, baseline_accession, threads, command
    )
