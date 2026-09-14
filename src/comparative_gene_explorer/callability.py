"""Gene catalog and zero-based callable-CDS intersection logic."""

from __future__ import annotations

from bisect import bisect_right

from comparative_gene_explorer.validation import GeneRecord, ReferenceInput


def interval_bases(intervals: tuple[tuple[int, int], ...]) -> int:
    return sum(end - start for start, end in intervals)


def callable_overlap_bases(
    query: tuple[tuple[int, int], ...], targets: tuple[tuple[int, int], ...]
) -> int:
    if not query or not targets:
        return 0
    target_ends = [end for _start, end in targets]
    total = 0
    for query_start, query_end in query:
        index = bisect_right(target_ends, query_start)
        while index < len(targets) and targets[index][0] < query_end:
            target_start, target_end = targets[index]
            total += max(0, min(query_end, target_end) - max(query_start, target_start))
            index += 1
    return total


def gene_catalog_row(reference: ReferenceInput, gene: GeneRecord) -> dict[str, object]:
    return {
        "accession": reference.accession,
        "label": reference.label,
        "role": reference.role,
        "gene_id": gene.gene_id,
        "gene_name": gene.gene_name,
        "feature_id": gene.feature_id,
        "feature_type": gene.feature_type,
        "gene_biotype": gene.gene_biotype,
        "contig": gene.contig,
        "strand": gene.strand,
        "start_1based": gene.start_1based,
        "end_1based_raw": gene.end_1based_raw,
        "origin_spanning": len(gene.feature_intervals) > 1,
        "feature_segments": len(gene.feature_intervals),
        "cds_segments": len(gene.cds_intervals),
        "cds_bases": interval_bases(gene.cds_intervals),
        "protein_ids": ";".join(gene.protein_ids),
    }


def gene_callability_row(
    reference: ReferenceInput,
    gene: GeneRecord,
    minimum_fraction: float,
) -> dict[str, object]:
    cds_bases = interval_bases(gene.cds_intervals)
    callable_bases = callable_overlap_bases(
        gene.cds_intervals, reference.callable_intervals[gene.contig]
    )
    if callable_bases > cds_bases:
        raise ValueError(f"callable CDS exceeds total CDS for {reference.accession}:{gene.gene_id}")
    fraction = callable_bases / cds_bases if cds_bases else None
    if fraction is None:
        status = "no_cds"
        meets_threshold = False
    elif fraction >= minimum_fraction:
        status = "callable"
        meets_threshold = True
    else:
        status = "low_callability"
        meets_threshold = False
    return {
        "accession": reference.accession,
        "gene_id": gene.gene_id,
        "contig": gene.contig,
        "strand": gene.strand,
        "cds_bases": cds_bases,
        "callable_cds_bases": callable_bases,
        "callable_fraction": "" if fraction is None else fraction,
        "minimum_callable_fraction": minimum_fraction,
        "meets_threshold": meets_threshold,
        "status": status,
    }


def build_gene_tables(
    references: tuple[ReferenceInput, ...], minimum_fraction: float
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    catalog: list[dict[str, object]] = []
    callability: list[dict[str, object]] = []
    for reference in references:
        for gene in reference.genes:
            catalog.append(gene_catalog_row(reference, gene))
            callability.append(gene_callability_row(reference, gene, minimum_fraction))
    return catalog, callability

