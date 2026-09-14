# Methods

## Study design

Project 15 compares gene-level consequence burden across the six references
evaluated in Project 13. MG1655 (`GCF_000005845.2`) is the baseline. The other
five references are sensitivity-analysis cohorts because they did not satisfy
Project 13's combined mapping and callability eligibility rules.

All comparisons use the same public sample, `SRR13921545`, and the accepted
reference-specific callsets from Project 13. No reads are remapped and no
variants are recalled in this project.

## Input validation

The workflow requires one baseline and at least one alternative in
`reference_manifest.csv`. For every reference it verifies:

- exact accession version and unique label;
- reference FASTA and accepted VCF SHA-256;
- FASTA, GenBank, and GFF3 sequence/contig agreement;
- valid, unique gene and protein identifiers;
- one expected VCF sample, sorted unique variant keys, `PASS` status, and exact
  reference alleles;
- zero-based, half-open, ordered callable BED intervals within contig bounds;
  and
- exactly one CSI or TBI index for each accepted VCF.

Source hashes are calculated again after each workflow completes.

## Reference-specific annotation

One private SnpEff database is built for each exact NCBI GenBank record. The
project-local configuration uses the bacterial and plant-plastid codon table.
Independent CDS and protein checks must each have an error rate no higher than
2%.

Accepted VCFs are checked and normalized with BCFtools before annotation. SnpEff
runs offline with downloads and telemetry disabled and the upstream/downstream
window set to zero. Every `ANN` effect is preserved. This stage does not alter
the accepted/rejected status inherited from Project 13.

## Proteome preparation and orthology

NCBI protein accessions are mapped to GFF3 genes. Each protein identifier is
prefixed with its exact assembly accession to prevent cross-reference name
collisions. Ambiguous protein-to-gene mappings are reported rather than forced.
Genes without a protein product remain explicit.

OrthoFinder 3.1.5 runs once across all six proteomes with:

```text
orthofinder -f proteomes -M msa -S diamond -A famsa -T fasttree -t 1 -a 1 -X -o results
```

The run uses DIAMOND sequence search, FAMSA multiple-sequence alignment,
FastTree gene/species trees, preserved input identifiers, and one thread.
Proteins must occur exactly once in either an orthogroup or the unassigned set.

Pairwise relationships are classified as `one_to_one`, `one_to_many`,
`many_to_one`, `many_to_many`, `no_ortholog`, or `no_protein`. A renamed
one-to-one pair is recorded when the gene identifiers map unambiguously but the
displayed gene symbols differ.

## Callable CDS calculation

GFF3 uses one-based inclusive coordinates; BED uses zero-based half-open
coordinates. Gene and CDS coordinates are converted before intersection.
Overlapping or adjacent CDS segments within one gene are merged so bases cannot
be counted twice. Circular origin-spanning features are split at the reference
boundary.

For each gene:

```text
callable fraction = callable CDS bases / total combined CDS bases
```

Proteinless or noncoding features have no callable CDS fraction. A coding gene
passes when at least 90% of its combined CDS bases are callable.

## Variant–gene assignments

Each SnpEff effect whose gene identifier exactly matches the reference's gene
catalog is grouped by `(reference, variant, gene)`. One deterministic effect is
retained for burden counting:

1. highest impact: `HIGH`, `MODERATE`, `LOW`, then `MODIFIER`;
2. SnpEff annotation order; and
3. feature identifier.

All original effects remain in the local annotation table. This reduction only
prevents multiple transcripts or overlapping effect labels from inflating the
same variant–gene pair. A variant overlapping two distinct genes contributes
once to each gene.

## Cross-reference comparability

An MG1655/alternative gene pair is comparable when:

1. OrthoFinder reports a one-to-one relationship;
2. both proteins map unambiguously to one gene; and
3. both genes meet the 90% callable-CDS threshold.

The burden for a gene or matched gene set is:

```text
1,000 × distinct variant–gene assignments / callable CDS bases
```

Raw counts are retained alongside the normalized density. Alternative-reference
eligibility from Project 13 is carried into every comparison.

## Review-candidate context

Project 14's 80 high-impact or warning-bearing review candidates are joined to
all five alternative relationships by MG1655 gene. For a one-to-one target, the
table reports whether that alternative gene contains any high-impact or warning
annotation and whether it contains the same consequence term.

These are gene-level contexts only. They do not imply coordinate or allele
equivalence.

## Reporting and reproducibility

The workflow writes detailed burden/comparison tables, compact summaries, a
four-panel Matplotlib dashboard, and JSON manifests. Outputs are published by a
same-parent atomic rename only after validation succeeds. Injected failures are
tested to leave no final directory or staging residue.
