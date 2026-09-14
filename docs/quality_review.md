# Milestone 5 — Quality review

## Decision

Milestone 5 passed. The expanded Python 3.12 suite contains 58 tests, all of
which pass. A separate read-only validation of the authentic sources also
passes with six references, 28,859 genes, and 357,581 accepted variants.

No authentic result totals changed during the review.

## Hardening changes

The quality review identified and corrected three boundary cases:

1. GFF3 parsing now rejects a repeated identifier when it belongs to unrelated
   features. A legitimate discontinuous CDS may still reuse one identifier on
   multiple rows when its feature type and parent agree.
2. VCF validation now rejects non-integer or out-of-range positions and empty
   alleles before attempting a reference comparison.
3. Empty callsets retain stable consequence/impact CSV headers and a valid
   dashboard scale instead of depending on the first result row.

The SnpEff database gate was also factored into one testable rule. An error rate
of exactly 2% passes; any value above 2% fails.

## Test coverage

| Area | Cases reviewed |
|---|---|
| Manifests | malformed columns, duplicate accessions and labels, unsafe parameters, existing outputs |
| Reference sources | FASTA/GenBank disagreement, unknown contigs, invalid coordinates, checksum preservation |
| Identifiers | duplicate FASTA and GFF3 identifiers, accession-prefix collisions, ambiguous protein mappings |
| Annotation structures | plus/minus strands, circular and discontinuous CDS features, pseudogenes, overlapping genes |
| Callable BED | negative, empty, out-of-range, overlapping, unsorted, and unknown-contig intervals |
| VCF | duplicate and multiallelic keys, wrong sample, non-PASS records, reference mismatch, empty callsets |
| Orthology | one-to-one, one-to-many, many-to-one, many-to-many, missing, renamed, unassigned, and proteinless cases |
| Consequences | highest-impact effect selection, overlapping genes, exact `(variant, gene)` accounting, empty effects |
| Native tools | failure markers, wrong pinned versions, sanitized diagnostics, database validation failures |
| Atomic output | generic failure, annotation subprocess failure, and OrthoFinder failure with complete cleanup |
| Authentic artifacts | exact callset totals, Project 14 reproduction, protein membership, comparison and dashboard agreement |

The machine-readable matrix is in `results/quality_review/test_matrix.csv`.

## Tool verification

The live WSL2 environment reports:

- Python 3.12.14
- OrthoFinder 3.1.5
- SnpEff 5.4.0c
- BCFtools 1.24
- NCBI Datasets CLI 18.36.0
- OpenJDK 21.0.10
- Matplotlib 3.11.1

All match the pinned project environment.

## Atomic and accounting checks

Injected failures in both annotation and OrthoFinder execution leave no final
output and no sibling staging directory. Existing output directories are
rejected rather than overwritten.

Authentic-result tests continue to verify:

- all 26,768 proteins occur in exactly one orthogroup or unassigned set;
- all 23,255 MG1655/alternative comparison keys are unique;
- callable CDS never exceeds total CDS;
- all 357,959 effects are accounted for as a selected variant–gene record,
  a non-gene effect, or an effect collapsed within an existing pair;
- all 80 Project 14 candidates have five gene-level contexts; and
- the dashboard has the expected full-resolution dimensions.

## Publication boundary

The complete annotation workspace remains under ignored local storage. No
repository was created or published. Documentation and publication checks are
reserved for Milestone 6.
