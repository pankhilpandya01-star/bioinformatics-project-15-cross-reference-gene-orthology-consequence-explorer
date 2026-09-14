"""Gene-level consequence burden and cross-reference comparison analysis."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from comparative_gene_explorer import __version__
from comparative_gene_explorer.atomic import atomic_output_directory


IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
RELATIONSHIP_ORDER = (
    "one_to_one",
    "one_to_many",
    "many_to_one",
    "many_to_many",
    "no_ortholog",
    "no_protein",
)


class ComparisonError(RuntimeError):
    """Raised when compact comparison inputs fail an accounting rule."""


@dataclass(frozen=True, slots=True)
class ComparisonConfig:
    annotation_directory: Path
    orthology_directory: Path
    project14_review_candidates: Path
    project14_annotation_effects: Path
    project13_reference_concordance: Path
    output_directory: Path
    baseline_accession: str = "GCF_000005845.2"
    minimum_gene_callable_fraction: float = 0.90


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ComparisonError(f"required comparison input is missing: {path.name}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ComparisonError(f"CSV has no header: {path.name}")
        return list(reader)


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bool(value: str) -> bool:
    if value == "True" or value.lower() == "true":
        return True
    if value == "False" or value.lower() == "false":
        return False
    raise ComparisonError(f"expected a Boolean CSV value; found {value!r}")


def _float_or_blank(value: object) -> float | str:
    return "" if value == "" else float(value)


def _load_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise ComparisonError(f"required manifest is missing: {path.name}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("status") != "pass":
        raise ComparisonError(f"input manifest is not a passing workflow: {path.name}")
    return data


def _prepare_reference_metadata(
    rows: list[dict[str, str]], baseline_accession: str
) -> tuple[list[str], dict[str, dict[str, object]]]:
    accessions = [row["accession"] for row in rows]
    if len(accessions) != 6 or len(set(accessions)) != 6:
        raise ComparisonError("Project 13 concordance must contain six unique references")
    if baseline_accession not in accessions:
        raise ComparisonError("baseline accession is absent from Project 13 concordance")
    metadata: dict[str, dict[str, object]] = {}
    for row in rows:
        eligible = row["eligible"].lower() == "true"
        metadata[row["accession"]] = {
            "label": row["label"],
            "role": row["role"],
            "project13_eligible": eligible,
            "project13_mapping_eligible": row["mapping_eligible"].lower() == "true",
            "project13_callability_eligible": row["callability_eligible"].lower() == "true",
        }
    eligible = [accession for accession in accessions if metadata[accession]["project13_eligible"]]
    if eligible != [baseline_accession]:
        raise ComparisonError("Project 13 eligibility does not match the inherited decision")
    return accessions, metadata


def _select_variant_gene_assignments(
    effects: list[dict[str, str]], valid_genes: set[tuple[str, str]]
) -> tuple[dict[tuple[str, str, str], dict[str, str]], int]:
    selected: dict[tuple[str, str, str], dict[str, str]] = {}
    unlinked_effects = 0
    for effect in effects:
        gene_key = (effect["accession"], effect["gene_id"])
        if gene_key not in valid_genes:
            unlinked_effects += 1
            continue
        key = (effect["accession"], effect["variant_key"], effect["gene_id"])
        existing = selected.get(key)
        rank = (
            IMPACT_ORDER.get(effect["impact"], 99),
            int(effect["effect_order"]),
            effect["feature_id"],
        )
        if existing is None:
            selected[key] = effect
            continue
        existing_rank = (
            IMPACT_ORDER.get(existing["impact"], 99),
            int(existing["effect_order"]),
            existing["feature_id"],
        )
        if rank < existing_rank:
            selected[key] = effect
    return selected, unlinked_effects


def _validate_project14_reproduction(
    project14_effects: list[dict[str, str]],
    current_effects: list[dict[str, str]],
    baseline_accession: str,
) -> list[dict[str, object]]:
    current_baseline = [
        row for row in current_effects if row["accession"] == baseline_accession
    ]
    project14_fields = (
        "variant_key",
        "effect_index",
        "allele",
        "consequence_terms",
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
        "messages",
    )
    current_fields = (
        "variant_key",
        "effect_order",
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
    expected = Counter(
        tuple(row[field] for field in project14_fields) for row in project14_effects
    )
    observed = Counter(
        tuple(row[field] for field in current_fields) for row in current_baseline
    )
    if len(project14_effects) != 31_059 or len(current_baseline) != 31_059:
        raise ComparisonError("MG1655 annotation does not contain 31,059 effects")
    if expected != observed:
        raise ComparisonError("MG1655 effects do not exactly reproduce Project 14")
    return [
        {
            "baseline_accession": baseline_accession,
            "project14_effects": len(project14_effects),
            "project15_effects": len(current_baseline),
            "canonical_effect_multiset_identical": True,
            "passed": True,
        }
    ]


def _build_gene_burden(
    gene_catalog: list[dict[str, str]],
    callability: list[dict[str, str]],
    effects: list[dict[str, str]],
    metadata: dict[str, dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    dict[tuple[str, str], dict[str, object]],
    dict[tuple[str, str], Counter[str]],
    list[dict[str, object]],
    list[dict[str, object]],
    int,
]:
    catalog_by_gene = {
        (row["accession"], row["gene_id"]): row for row in gene_catalog
    }
    if len(catalog_by_gene) != len(gene_catalog):
        raise ComparisonError("gene catalog contains duplicate accession/gene keys")
    callable_by_gene = {
        (row["accession"], row["gene_id"]): row for row in callability
    }
    if set(callable_by_gene) != set(catalog_by_gene):
        raise ComparisonError("gene catalog and callability keys disagree")

    assignments, unlinked_effects = _select_variant_gene_assignments(
        effects, set(catalog_by_gene)
    )
    assignments_by_gene: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for (accession, _variant_key, gene_id), effect in assignments.items():
        assignments_by_gene[(accession, gene_id)].append(effect)

    burden_rows: list[dict[str, object]] = []
    burden_by_gene: dict[tuple[str, str], dict[str, object]] = {}
    consequence_by_gene: dict[tuple[str, str], Counter[str]] = {}
    consequence_summary: Counter[tuple[str, str]] = Counter()
    impact_summary: Counter[tuple[str, str]] = Counter()
    for key, catalog in catalog_by_gene.items():
        accession, gene_id = key
        gene_effects = assignments_by_gene.get(key, [])
        variant_keys = {effect["variant_key"] for effect in gene_effects}
        if len(variant_keys) != len(gene_effects):
            raise ComparisonError(f"duplicate variant/gene assignment: {accession}:{gene_id}")
        consequences = Counter(effect["consequence"] for effect in gene_effects)
        impacts = Counter(effect["impact"] for effect in gene_effects)
        for consequence, count in consequences.items():
            consequence_summary[(accession, consequence)] += count
        for impact, count in impacts.items():
            impact_summary[(accession, impact)] += count
        callable_row = callable_by_gene[key]
        cds_bases = int(callable_row["cds_bases"])
        callable_bases = int(callable_row["callable_cds_bases"])
        density: float | str = (
            len(gene_effects) * 1000 / callable_bases if callable_bases else ""
        )
        row: dict[str, object] = {
            "accession": accession,
            "label": metadata[accession]["label"],
            "role": metadata[accession]["role"],
            "project13_eligible": metadata[accession]["project13_eligible"],
            "gene_id": gene_id,
            "gene_name": catalog["gene_name"],
            "feature_type": catalog["feature_type"],
            "gene_biotype": catalog["gene_biotype"],
            "contig": catalog["contig"],
            "strand": catalog["strand"],
            "cds_bases": cds_bases,
            "callable_cds_bases": callable_bases,
            "callable_fraction": callable_row["callable_fraction"],
            "meets_callable_threshold": _bool(callable_row["meets_threshold"]),
            "distinct_variant_count": len(gene_effects),
            "high_impact_variants": impacts["HIGH"],
            "moderate_impact_variants": impacts["MODERATE"],
            "low_impact_variants": impacts["LOW"],
            "modifier_impact_variants": impacts["MODIFIER"],
            "warning_variant_count": sum(bool(effect["warnings"]) for effect in gene_effects),
            "variant_density_per_1000_callable_cds": density,
            "consequence_profile": ";".join(
                f"{name}:{count}" for name, count in sorted(consequences.items())
            ),
        }
        burden_rows.append(row)
        burden_by_gene[key] = row
        consequence_by_gene[key] = consequences

    consequence_rows = [
        {
            "accession": accession,
            "label": metadata[accession]["label"],
            "consequence": consequence,
            "variant_gene_count": count,
        }
        for (accession, consequence), count in sorted(consequence_summary.items())
    ]
    impact_rows = [
        {
            "accession": accession,
            "label": metadata[accession]["label"],
            "impact": impact,
            "variant_gene_count": count,
        }
        for (accession, impact), count in sorted(impact_summary.items())
    ]
    return (
        burden_rows,
        burden_by_gene,
        consequence_by_gene,
        consequence_rows,
        impact_rows,
        unlinked_effects,
    )


def _comparison_status(
    relationship: str,
    baseline: dict[str, object],
    alternative: dict[str, object] | None,
) -> str:
    if relationship != "one_to_one":
        return relationship
    if alternative is None:
        raise ComparisonError("one-to-one relationship lacks an alternative gene")
    baseline_callable = bool(baseline["meets_callable_threshold"])
    alternative_callable = bool(alternative["meets_callable_threshold"])
    if baseline_callable and alternative_callable:
        return "comparable"
    if not baseline_callable and not alternative_callable:
        return "both_low_callability"
    if not baseline_callable:
        return "baseline_low_callability"
    return "alternative_low_callability"


def _build_baseline_comparison(
    relationships: list[dict[str, str]],
    burden_by_gene: dict[tuple[str, str], dict[str, object]],
    metadata: dict[str, dict[str, object]],
    baseline_accession: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    statuses: Counter[tuple[str, str]] = Counter()
    for relationship_row in relationships:
        if relationship_row["baseline_accession"] != baseline_accession:
            raise ComparisonError("orthology relationship table uses an unexpected baseline")
        baseline_gene = relationship_row["baseline_gene_id"]
        alternative_accession = relationship_row["alternative_accession"]
        baseline = burden_by_gene[(baseline_accession, baseline_gene)]
        alternative_gene_ids = [
            item for item in relationship_row["alternative_gene_ids"].split(";") if item
        ]
        alternative: dict[str, object] | None = None
        if relationship_row["relationship"] == "one_to_one":
            if len(alternative_gene_ids) != 1:
                raise ComparisonError("one-to-one relationship does not contain one target gene")
            alternative = burden_by_gene.get(
                (alternative_accession, alternative_gene_ids[0])
            )
            if alternative is None:
                raise ComparisonError("one-to-one target gene is absent from the burden table")
        status = _comparison_status(
            relationship_row["relationship"], baseline, alternative
        )
        comparable = status == "comparable"
        baseline_density = _float_or_blank(
            baseline["variant_density_per_1000_callable_cds"]
        )
        alternative_density = (
            _float_or_blank(alternative["variant_density_per_1000_callable_cds"])
            if alternative is not None
            else ""
        )
        density_difference: float | str = ""
        density_ratio: float | str = ""
        density_direction = "not_comparable"
        if comparable:
            density_difference = float(alternative_density) - float(baseline_density)
            if float(baseline_density) > 0:
                density_ratio = float(alternative_density) / float(baseline_density)
            if float(alternative_density) > float(baseline_density):
                density_direction = "higher_in_alternative"
            elif float(alternative_density) < float(baseline_density):
                density_direction = "lower_in_alternative"
            else:
                density_direction = "equal"
        row: dict[str, object] = {
            "baseline_accession": baseline_accession,
            "baseline_gene_id": baseline_gene,
            "baseline_gene_name": baseline["gene_name"],
            "alternative_accession": alternative_accession,
            "alternative_label": metadata[alternative_accession]["label"],
            "project13_alternative_eligible": metadata[alternative_accession][
                "project13_eligible"
            ],
            "relationship": relationship_row["relationship"],
            "orthogroups": relationship_row["orthogroups"],
            "alternative_gene_ids": relationship_row["alternative_gene_ids"],
            "alternative_gene_names": relationship_row["alternative_gene_names"],
            "renamed_one_to_one": _bool(relationship_row["renamed_one_to_one"]),
            "comparison_status": status,
            "comparable": comparable,
            "baseline_callable_fraction": baseline["callable_fraction"],
            "alternative_callable_fraction": (
                alternative["callable_fraction"] if alternative is not None else ""
            ),
            "baseline_callable_cds_bases": baseline["callable_cds_bases"],
            "alternative_callable_cds_bases": (
                alternative["callable_cds_bases"] if alternative is not None else ""
            ),
            "baseline_distinct_variant_count": baseline["distinct_variant_count"],
            "alternative_distinct_variant_count": (
                alternative["distinct_variant_count"] if alternative is not None else ""
            ),
            "baseline_high_impact_variants": baseline["high_impact_variants"],
            "alternative_high_impact_variants": (
                alternative["high_impact_variants"] if alternative is not None else ""
            ),
            "baseline_variant_density_per_1000_callable_cds": baseline_density,
            "alternative_variant_density_per_1000_callable_cds": alternative_density,
            "alternative_minus_baseline_density": density_difference,
            "alternative_to_baseline_density_ratio": density_ratio,
            "density_direction": density_direction,
        }
        rows.append(row)
        statuses[(alternative_accession, status)] += 1

    keys = [
        (row["baseline_gene_id"], row["alternative_accession"]) for row in rows
    ]
    if len(keys) != len(set(keys)):
        raise ComparisonError("baseline comparison contains duplicate gene/reference keys")
    status_rows = [
        {
            "alternative_accession": accession,
            "alternative_label": metadata[accession]["label"],
            "comparison_status": status,
            "gene_count": count,
        }
        for (accession, status), count in sorted(statuses.items())
    ]

    density_rows: list[dict[str, object]] = []
    alternatives = [
        accession for accession in metadata if accession != baseline_accession
    ]
    for accession in alternatives:
        comparable_rows = [
            row
            for row in rows
            if row["alternative_accession"] == accession and row["comparable"]
        ]
        baseline_variants = sum(
            int(row["baseline_distinct_variant_count"]) for row in comparable_rows
        )
        alternative_variants = sum(
            int(row["alternative_distinct_variant_count"]) for row in comparable_rows
        )
        baseline_bases = sum(
            int(row["baseline_callable_cds_bases"]) for row in comparable_rows
        )
        alternative_bases = sum(
            int(row["alternative_callable_cds_bases"]) for row in comparable_rows
        )
        density_rows.append(
            {
                "alternative_accession": accession,
                "alternative_label": metadata[accession]["label"],
                "project13_alternative_eligible": metadata[accession][
                    "project13_eligible"
                ],
                "comparable_gene_pairs": len(comparable_rows),
                "baseline_variant_gene_assignments": baseline_variants,
                "alternative_variant_gene_assignments": alternative_variants,
                "baseline_callable_cds_bases": baseline_bases,
                "alternative_callable_cds_bases": alternative_bases,
                "baseline_variants_per_1000_callable_cds": (
                    baseline_variants * 1000 / baseline_bases if baseline_bases else ""
                ),
                "alternative_variants_per_1000_callable_cds": (
                    alternative_variants * 1000 / alternative_bases
                    if alternative_bases
                    else ""
                ),
            }
        )
    return rows, status_rows, density_rows


def _build_review_contexts(
    review_candidates: list[dict[str, str]],
    comparison_rows: list[dict[str, object]],
    burden_by_gene: dict[tuple[str, str], dict[str, object]],
    consequence_by_gene: dict[tuple[str, str], Counter[str]],
    metadata: dict[str, dict[str, object]],
    baseline_accession: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if len(review_candidates) != 80:
        raise ComparisonError("Project 14 review-candidate count is not 80")
    if len({row["variant_key"] for row in review_candidates}) != 80:
        raise ComparisonError("Project 14 review candidates contain duplicate variant keys")
    comparisons = {
        (str(row["baseline_gene_id"]), str(row["alternative_accession"])): row
        for row in comparison_rows
    }
    alternatives = [
        accession for accession in metadata if accession != baseline_accession
    ]
    rows: list[dict[str, object]] = []
    context_counts: Counter[tuple[str, str]] = Counter()
    for candidate in review_candidates:
        for accession in alternatives:
            comparison = comparisons.get((candidate["gene_id"], accession))
            if comparison is None:
                raise ComparisonError(
                    f"review-candidate gene is absent from orthology table: {candidate['gene_id']}"
                )
            alternative_gene_ids = [
                item
                for item in str(comparison["alternative_gene_ids"]).split(";")
                if item
            ]
            alternative: dict[str, object] | None = None
            if comparison["relationship"] == "one_to_one":
                alternative = burden_by_gene[(accession, alternative_gene_ids[0])]
            comparable = bool(comparison["comparable"])
            has_high = bool(alternative and int(alternative["high_impact_variants"]) > 0)
            has_warning = bool(alternative and int(alternative["warning_variant_count"]) > 0)
            same_consequence = bool(
                alternative
                and candidate["primary_consequence"]
                in consequence_by_gene[(accession, alternative_gene_ids[0])]
            )
            if not comparable:
                context = "not_comparable"
            elif has_high:
                context = "comparable_with_high_impact_signal"
            else:
                context = "comparable_without_high_impact_signal"
            context_counts[(accession, context)] += 1
            rows.append(
                {
                    "baseline_variant_key": candidate["variant_key"],
                    "baseline_gene_id": candidate["gene_id"],
                    "baseline_gene_name": candidate["gene_name"],
                    "baseline_primary_consequence": candidate["primary_consequence"],
                    "baseline_primary_impact": candidate["primary_impact"],
                    "baseline_review_reasons": candidate["review_reasons"],
                    "alternative_accession": accession,
                    "alternative_label": metadata[accession]["label"],
                    "relationship": comparison["relationship"],
                    "comparison_status": comparison["comparison_status"],
                    "comparable": comparable,
                    "alternative_gene_ids": comparison["alternative_gene_ids"],
                    "alternative_gene_names": comparison["alternative_gene_names"],
                    "alternative_distinct_variant_count": (
                        alternative["distinct_variant_count"] if alternative else ""
                    ),
                    "alternative_high_impact_variants": (
                        alternative["high_impact_variants"] if alternative else ""
                    ),
                    "alternative_has_high_impact_signal": has_high,
                    "alternative_has_annotation_warning": has_warning,
                    "alternative_has_same_consequence_term": same_consequence,
                    "comparison_scope": "gene_level_only_no_variant_equivalence",
                    "context_category": context,
                }
            )
    summary = [
        {
            "alternative_accession": accession,
            "alternative_label": metadata[accession]["label"],
            "context_category": context,
            "review_candidate_count": count,
        }
        for (accession, context), count in sorted(context_counts.items())
    ]
    if len(rows) != 400:
        raise ComparisonError("review-candidate orthology table must contain 400 rows")
    return rows, summary


def _build_callability_summary(
    burdens: list[dict[str, object]],
    accessions: list[str],
    metadata: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for accession in accessions:
        genes = [row for row in burdens if row["accession"] == accession]
        coding = [row for row in genes if int(row["cds_bases"]) > 0]
        callable_genes = [row for row in coding if row["meets_callable_threshold"]]
        rows.append(
            {
                "accession": accession,
                "label": metadata[accession]["label"],
                "role": metadata[accession]["role"],
                "project13_eligible": metadata[accession]["project13_eligible"],
                "genes": len(genes),
                "coding_genes": len(coding),
                "callable_coding_genes": len(callable_genes),
                "callable_coding_gene_percent": (
                    len(callable_genes) * 100 / len(coding) if coding else 0
                ),
                "total_cds_bases": sum(int(row["cds_bases"]) for row in coding),
                "callable_cds_bases": sum(
                    int(row["callable_cds_bases"]) for row in coding
                ),
            }
        )
    return rows


def _draw_dashboard(
    path: Path,
    comparison_rows: list[dict[str, object]],
    density_rows: list[dict[str, object]],
    review_summary: list[dict[str, object]],
) -> None:
    alternatives = [row["alternative_accession"] for row in density_rows]
    labels = [str(row["alternative_label"]) for row in density_rows]
    relationship_counts = Counter(
        (str(row["alternative_accession"]), str(row["relationship"]))
        for row in comparison_rows
    )
    status_counts = Counter(
        (str(row["alternative_accession"]), str(row["comparison_status"]))
        for row in comparison_rows
    )
    review_counts = {
        (row["alternative_accession"], row["context_category"]): int(
            row["review_candidate_count"]
        )
        for row in review_summary
    }

    colors = {
        "one_to_one": "#2F6690",
        "one_to_many": "#3A7D44",
        "many_to_one": "#70A37F",
        "many_to_many": "#F2B134",
        "no_ortholog": "#C8553D",
        "no_protein": "#8D99AE",
    }
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    fig.suptitle(
        "Cross-reference gene orthology and consequence burden",
        fontsize=18,
        fontweight="bold",
    )

    ax = axes[0, 0]
    bottoms = [0] * len(alternatives)
    for relationship in RELATIONSHIP_ORDER:
        values = [relationship_counts[(accession, relationship)] for accession in alternatives]
        ax.bar(labels, values, bottom=bottoms, label=relationship.replace("_", " "), color=colors[relationship])
        bottoms = [bottom + value for bottom, value in zip(bottoms, values, strict=True)]
    ax.set_title("A. MG1655 gene relationships")
    ax.set_ylabel("MG1655 genes")
    ax.legend(fontsize=8, ncol=2, frameon=False)

    ax = axes[0, 1]
    one_to_one = [relationship_counts[(accession, "one_to_one")] for accession in alternatives]
    comparable = [status_counts[(accession, "comparable")] for accession in alternatives]
    positions = range(len(alternatives))
    ax.bar([position - 0.2 for position in positions], one_to_one, width=0.4, label="one-to-one", color="#7FA8C9")
    ax.bar([position + 0.2 for position in positions], comparable, width=0.4, label="one-to-one + both callable", color="#2F6690")
    for position, value in zip(positions, comparable, strict=True):
        ax.text(
            position + 0.2,
            value + max(one_to_one) * 0.015,
            f"n={value:,}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#1F3F5B",
        )
    ax.set_xticks(list(positions), labels)
    ax.set_title("B. Comparable callable genes")
    ax.set_ylabel("Gene pairs")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    baseline_density = [float(row["baseline_variants_per_1000_callable_cds"]) if row["baseline_variants_per_1000_callable_cds"] != "" else 0 for row in density_rows]
    alternative_density = [float(row["alternative_variants_per_1000_callable_cds"]) if row["alternative_variants_per_1000_callable_cds"] != "" else 0 for row in density_rows]
    ax.bar([position - 0.2 for position in positions], baseline_density, width=0.4, label="MG1655", color="#2F6690")
    ax.bar([position + 0.2 for position in positions], alternative_density, width=0.4, label="alternative", color="#C8553D")
    for position, row, left, right in zip(
        positions, density_rows, baseline_density, alternative_density, strict=True
    ):
        ax.text(
            position,
            max(left, right) + max(alternative_density) * 0.025,
            f"{int(row['comparable_gene_pairs']):,} pairs",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#444444",
        )
    ax.set_xticks(list(positions), labels)
    ax.set_title("C. Burden in comparable genes")
    ax.set_ylabel("Variant–gene assignments per 1,000 callable CDS bases")
    ax.set_ylim(0, max(1.0, max(alternative_density + baseline_density) * 1.13))
    ax.legend(frameon=False)

    ax = axes[1, 1]
    categories = (
        "comparable_with_high_impact_signal",
        "comparable_without_high_impact_signal",
        "not_comparable",
    )
    category_colors = ("#C8553D", "#3A7D44", "#8D99AE")
    bottoms = [0] * len(alternatives)
    for category, color in zip(categories, category_colors, strict=True):
        values = [review_counts.get((accession, category), 0) for accession in alternatives]
        ax.bar(labels, values, bottom=bottoms, label=category.replace("_", " "), color=color)
        bottoms = [bottom + value for bottom, value in zip(bottoms, values, strict=True)]
    ax.set_title("D. Project 14 review-candidate gene contexts")
    ax.set_ylabel("Review candidates (80 per reference)")
    ax.legend(fontsize=8, frameon=False)

    for ax in axes.flat:
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelrotation=24, labelsize=9)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("right")
    fig.text(
        0.5,
        0.004,
        "Alternative references are sensitivity analyses; gene-level similarity does not establish variant equivalence.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def run_comparative_analysis(config: ComparisonConfig) -> Path:
    output = config.output_directory.resolve()
    if output.exists():
        raise ComparisonError(f"output directory already exists: {output}")
    if not 0 < config.minimum_gene_callable_fraction <= 1:
        raise ComparisonError("minimum callable fraction must be greater than zero and at most one")

    annotation = config.annotation_directory.resolve()
    orthology = config.orthology_directory.resolve()
    annotation_manifest_path = annotation / "run_manifest.json"
    orthology_manifest_path = orthology / "run_manifest.json"
    annotation_manifest = _load_json(annotation_manifest_path)
    orthology_manifest = _load_json(orthology_manifest_path)
    if annotation_manifest.get("parameters", {}).get("baseline_accession") != config.baseline_accession:  # type: ignore[union-attr]
        raise ComparisonError("annotation manifest baseline does not match")
    if orthology_manifest.get("baseline_accession") != config.baseline_accession:
        raise ComparisonError("orthology manifest baseline does not match")

    input_paths = {
        "annotation_manifest": annotation_manifest_path,
        "gene_catalog": annotation / "tables" / "gene_catalog.csv",
        "gene_callability": annotation / "tables" / "gene_callability.csv",
        "variant_annotations": annotation / "tables" / "variant_annotations.csv",
        "annotation_effects": annotation / "tables" / "annotation_effects.csv",
        "orthology_manifest": orthology_manifest_path,
        "orthology_relationships": orthology / "tables" / "orthology_relationships.csv",
        "orthogroup_catalog": orthology / "tables" / "orthogroup_catalog.csv",
        "project14_review_candidates": config.project14_review_candidates.resolve(),
        "project14_annotation_effects": config.project14_annotation_effects.resolve(),
        "project13_reference_concordance": config.project13_reference_concordance.resolve(),
    }
    source_hashes = {name: _sha256(path) for name, path in input_paths.items()}
    gene_catalog = _read_csv(input_paths["gene_catalog"])
    callability = _read_csv(input_paths["gene_callability"])
    variant_annotations = _read_csv(input_paths["variant_annotations"])
    effects = _read_csv(input_paths["annotation_effects"])
    relationships = _read_csv(input_paths["orthology_relationships"])
    review_candidates = _read_csv(input_paths["project14_review_candidates"])
    project14_effects = _read_csv(input_paths["project14_annotation_effects"])
    concordance = _read_csv(input_paths["project13_reference_concordance"])
    accessions, metadata = _prepare_reference_metadata(
        concordance, config.baseline_accession
    )
    if {row["accession"] for row in gene_catalog} != set(accessions):
        raise ComparisonError("annotation and Project 13 reference sets disagree")
    if len(variant_annotations) != int(annotation_manifest["counts"]["accepted_variants"]):  # type: ignore[index]
        raise ComparisonError("variant annotation count disagrees with its manifest")
    if len(effects) != int(annotation_manifest["counts"]["annotation_effects"]):  # type: ignore[index]
        raise ComparisonError("annotation effect count disagrees with its manifest")
    variant_keys = [
        (row["accession"], row["variant_key"]) for row in variant_annotations
    ]
    if len(variant_keys) != len(set(variant_keys)):
        raise ComparisonError("variant annotation table contains duplicate keys")
    if len(relationships) != int(
        orthology_manifest["counts"]["baseline_alternative_gene_comparisons"]  # type: ignore[index]
    ):
        raise ComparisonError("orthology relationship count disagrees with its manifest")
    project14_reproduction = _validate_project14_reproduction(
        project14_effects, effects, config.baseline_accession
    )

    (
        burdens,
        burden_by_gene,
        consequence_by_gene,
        consequence_summary,
        impact_summary,
        unlinked_effects,
    ) = _build_gene_burden(gene_catalog, callability, effects, metadata)
    comparison_rows, status_summary, density_summary = _build_baseline_comparison(
        relationships, burden_by_gene, metadata, config.baseline_accession
    )
    review_rows, review_summary = _build_review_contexts(
        review_candidates,
        comparison_rows,
        burden_by_gene,
        consequence_by_gene,
        metadata,
        config.baseline_accession,
    )
    callability_summary = _build_callability_summary(
        burdens, accessions, metadata
    )
    orthology_summary = [
        {
            "alternative_accession": accession,
            "alternative_label": metadata[accession]["label"],
            "relationship": relationship,
            "gene_count": sum(
                row["alternative_accession"] == accession
                and row["relationship"] == relationship
                for row in comparison_rows
            ),
        }
        for accession in accessions
        if accession != config.baseline_accession
        for relationship in RELATIONSHIP_ORDER
    ]

    with atomic_output_directory(output) as staging:
        tables = staging / "tables"
        summaries = staging / "summaries"
        dashboard = staging / "dashboard"
        tables.mkdir()
        summaries.mkdir()
        dashboard.mkdir()
        burden_fields = list(burdens[0])
        comparison_fields = list(comparison_rows[0])
        review_fields = list(review_rows[0])
        _write_csv(tables / "gene_consequence_burden.csv", burdens, burden_fields)
        _write_csv(
            tables / "baseline_ortholog_comparison.csv",
            comparison_rows,
            comparison_fields,
        )
        _write_csv(
            tables / "review_candidate_orthology.csv", review_rows, review_fields
        )
        _write_csv(
            summaries / "orthology_summary.csv",
            orthology_summary,
            list(orthology_summary[0]),
        )
        _write_csv(
            summaries / "callability_summary.csv",
            callability_summary,
            list(callability_summary[0]),
        )
        _write_csv(
            summaries / "consequence_summary.csv",
            consequence_summary,
            ["accession", "label", "consequence", "variant_gene_count"],
        )
        _write_csv(
            summaries / "impact_summary.csv",
            impact_summary,
            ["accession", "label", "impact", "variant_gene_count"],
        )
        _write_csv(
            summaries / "comparison_status_summary.csv",
            status_summary,
            list(status_summary[0]),
        )
        _write_csv(
            summaries / "comparable_density_summary.csv",
            density_summary,
            list(density_summary[0]),
        )
        _write_csv(
            summaries / "review_context_summary.csv",
            review_summary,
            list(review_summary[0]),
        )
        _write_csv(
            summaries / "project14_reproduction.csv",
            project14_reproduction,
            list(project14_reproduction[0]),
        )
        dashboard_path = dashboard / "comparative_gene_dashboard.png"
        _draw_dashboard(
            dashboard_path, comparison_rows, density_summary, review_summary
        )

        if {name: _sha256(path) for name, path in input_paths.items()} != source_hashes:
            raise ComparisonError("a source input changed during comparative analysis")
        comparable_pairs = sum(bool(row["comparable"]) for row in comparison_rows)
        variant_gene_assignments = sum(
            int(row["distinct_variant_count"]) for row in burdens
        )
        manifest = {
            "workflow": "comparative-gene-explorer",
            "workflow_version": __version__,
            "milestone": "authentic_comparison",
            "status": "pass",
            "baseline_accession": config.baseline_accession,
            "parameters": {
                "minimum_gene_callable_fraction": config.minimum_gene_callable_fraction,
                "variant_gene_rule": "one row per distinct variant and gene using the highest-impact effect",
                "density_denominator": "callable CDS bases",
            },
            "counts": {
                "references": len(accessions),
                "accepted_variants": len(variant_annotations),
                "annotation_effects": len(effects),
                "genes": len(burdens),
                "variant_gene_assignments": variant_gene_assignments,
                "effects_without_exact_gene_key": unlinked_effects,
                "effects_collapsed_within_variant_gene": (
                    len(effects) - unlinked_effects - variant_gene_assignments
                ),
                "baseline_alternative_gene_comparisons": len(comparison_rows),
                "comparable_gene_pairs": comparable_pairs,
                "project14_review_candidates": len(review_candidates),
                "review_candidate_reference_contexts": len(review_rows),
                "project14_effects_reproduced": len(project14_effects),
            },
            "project13_eligibility": {
                accession: metadata[accession]["project13_eligible"]
                for accession in accessions
            },
            "per_reference_annotations": annotation_manifest["references"],
            "source_files": {
                name: {"file": path.name, "sha256": source_hashes[name]}
                for name, path in input_paths.items()
            },
            "outputs": {
                "gene_consequence_burden_sha256": _sha256(
                    tables / "gene_consequence_burden.csv"
                ),
                "baseline_ortholog_comparison_sha256": _sha256(
                    tables / "baseline_ortholog_comparison.csv"
                ),
                "review_candidate_orthology_sha256": _sha256(
                    tables / "review_candidate_orthology.csv"
                ),
                "dashboard_sha256": _sha256(dashboard_path),
            },
            "tools": {
                "python": "3.12",
                "matplotlib": matplotlib.__version__,
                "orthofinder": "3.1.5",
                "snpeff": "5.4.0c",
                "bcftools": "1.24",
            },
            "source_inputs_unchanged": True,
            "project14_effects_exactly_reproduced": True,
            "interpretation_boundary": (
                "Alternative references are sensitivity-analysis cohorts. Orthology and gene-level burden do not establish coordinate-equivalent variants or identical function."
            ),
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    return output
