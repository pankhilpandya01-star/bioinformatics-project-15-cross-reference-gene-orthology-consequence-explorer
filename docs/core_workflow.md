# Core workflow review

## Decision

Milestone 2 passed. The installable `comparative-gene-explorer` CLI now validates
the reference manifest and every linked source, normalizes and annotates each
reference-specific callset offline, converts gene/CDS coordinates correctly,
intersects CDS intervals with callable BED regions, and publishes results only
after all reference and accounting checks pass.

The authentic six-reference manifest also passed the new core validator without
starting the deferred authentic annotation run:

| Accession | Contigs | Parsed genes/features | Accepted variants |
| --- | ---: | ---: | ---: |
| `GCF_000005845.2` | 1 | 4,651 | 30,973 |
| `GCF_003697165.2` | 2 | 4,938 | 92,665 |
| `GCF_036503815.1` | 5 | 4,992 | 125,327 |
| `GCF_005843885.1` | 1 | 5,002 | 49,201 |
| `GCF_000026225.1` | 2 | 4,477 | 15,264 |
| `GCF_002900365.1` | 3 | 4,799 | 44,151 |

The parsed totals include protein-coding, RNA, pseudogene, and other gene-level
features. They are not the same denominator as CDS or protein records.

## Controlled demonstration

Two 600-bp controlled references each contain a plus-strand gene (`ctrlA`) and
a minus-strand gene (`ctrlB`). The same two accepted variants were annotated in
both references:

| Position | Strand context | Expected result | Observed result |
| ---: | --- | --- | --- |
| 4 | Plus-strand `ctrlA` | Missense, p.Ala2Thr | Pass |
| 387 | Minus-strand `ctrlB` | Stop gained, p.Glu2* | Pass |

The callable BED inputs deliberately reverse the callable gene between the two
references:

| Reference | Gene | Callable CDS | Status at 90% |
| --- | --- | ---: | --- |
| `CTRL_PLUS.1` | `ctrlA` | 90/90 (100%) | Callable |
| `CTRL_PLUS.1` | `ctrlB` | 45/90 (50%) | Low callability |
| `CTRL_MINUS.1` | `ctrlA` | 45/90 (50%) | Low callability |
| `CTRL_MINUS.1` | `ctrlB` | 90/90 (100%) | Callable |

All four controlled SnpEff validation checks reported 0% error. Four input
variants produced four annotated variants and four effects, and each compressed
VCF has a CSI index.

## Coordinate and denominator rules

GFF3's one-based inclusive CDS coordinates are converted to zero-based
half-open intervals before BED intersection. Adjacent or overlapping CDS pieces
for the same gene are merged so bases cannot be counted twice. A verified
circular-origin feature is split into end-of-contig and start-of-contig pieces.
Callable CDS bases are required never to exceed total CDS bases.

Genes without CDS remain in the catalog with `no_cds` status. A coding gene is
`callable` only when its combined callable CDS fraction meets the configured
threshold; otherwise it is `low_callability`.

## Failure and publication behavior

- Manifest accessions and labels must be unique, with exactly one selected
  baseline and at least one alternative.
- Checksums, sample names, reference alleles, VCF sorting, indexes, contigs,
  GFF3 relationships, and BED intervals are validated before annotation.
- BCFtools confirms that normalization does not change the accepted variant-key
  set.
- Both CDS and protein database checks must pass before SnpEff annotation.
- Subprocess failures return sanitized diagnostics.
- A failed run removes its staging directory and leaves the requested output
  absent.
- Source checksums are recomputed before atomic publication.

Seventeen Python tests pass, covering the CLI contract, interval arithmetic,
multi-contig and circular-origin parsing, plus/minus-strand annotation,
low-callability classification, duplicate-manifest rejection, sanitized native
failures, result accounting, and atomic cleanup.

Authentic annotation and orthology comparison have not been run yet. Milestone
3 will prepare uniquely identified proteomes, map proteins to genes, run
OrthoFinder, and classify orthology relationships.

