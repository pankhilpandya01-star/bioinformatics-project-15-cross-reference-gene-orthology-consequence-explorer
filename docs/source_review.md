# Source review

## Decision

All six exact NCBI annotation packages and project-local SnpEff databases are
accepted for Project 15. Every reference FASTA reproduced its Project 13
SHA-256, every GenBank sequence matched its corresponding FASTA, every GFF3
feature used a known reference sequence, and all inherited accepted VCF and
callable-region totals reproduced Project 13.

The 12 predeclared SnpEff validation checks all passed the 2% maximum error
rate. This permits the project to proceed to manifest, coordinate, annotation,
and callability implementation. It does not make alternative-reference variants
coordinate-equivalent to MG1655 variants.

## Source and cohort checks

| Accession | Contigs | Reference bp | Accepted variants | Callable bp | CDS FASTA | Protein FASTA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `GCF_000005845.2` | 1 | 4,641,652 | 30,973 | 4,199,245 | 4,318 | 4,300 |
| `GCF_003697165.2` | 2 | 5,034,834 | 92,665 | 3,888,294 | 4,823 | 4,609 |
| `GCF_036503815.1` | 5 | 5,199,626 | 125,327 | 3,910,557 | 4,868 | 4,705 |
| `GCF_005843885.1` | 1 | 5,120,753 | 49,201 | 777,722 | 4,873 | 4,672 |
| `GCF_000026225.1` | 2 | 4,643,861 | 15,264 | 239,191 | 4,362 | 4,265 |
| `GCF_002900365.1` | 3 | 4,896,291 | 44,151 | 776,835 | 4,688 | 4,217 |

The VCF check required one sample (`SRR13921545`), only split `PASS` records,
unique variant keys, reference-order sorting, matching reference alleles, and an
index count equal to the parsed record count. BED checks required known contigs,
valid zero-based half-open intervals, sorted non-overlapping intervals, and the
exact callable-base totals above.

## Annotation contents

| Accession | GFF3 genes | Pseudogenes | CDS | ncRNA | tRNA | rRNA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `GCF_000005845.2` | 4,506 | 145 | 4,340 | 108 | 86 | 22 |
| `GCF_003697165.2` | 4,736 | 202 | 4,834 | 2 | 87 | 22 |
| `GCF_036503815.1` | 4,852 | 140 | 4,871 | 6 | 92 | 22 |
| `GCF_005843885.1` | 4,838 | 164 | 4,875 | 8 | 95 | 22 |
| `GCF_000026225.1` | 4,395 | 82 | 4,367 | 3 | 86 | 22 |
| `GCF_002900365.1` | 4,473 | 326 | 4,707 | 2 | 84 | 21 |

NCBI represents some genes crossing a circular replicon origin by extending the
GFF3 end coordinate past the linearized sequence length. The validator accepted
two such feature rows in `GCF_036503815.1` and four in `GCF_002900365.1` only
after confirming the corresponding `Is_circular=true` region declaration. The
same coordinate pattern remains invalid on a non-circular contig.

## SnpEff database validation

| Accession | CDS error | Protein error | Result |
| --- | ---: | ---: | --- |
| `GCF_000005845.2` | 0.232% | 0.000% | Pass |
| `GCF_003697165.2` | 0.000% | 0.000% | Pass |
| `GCF_036503815.1` | 0.021% | 0.021% | Pass |
| `GCF_005843885.1` | 0.000% | 0.000% | Pass |
| `GCF_000026225.1` | 0.000% | 0.000% | Pass |
| `GCF_002900365.1` | 0.043% | 0.046% | Pass |

Each database was built from the exact GenBank record with the bacterial and
plant-plastid codon table. NCBI CDS identifiers were rewritten in a
validation-only copy to match SnpEff's deterministic GenBank transcript IDs;
the downloaded CDS files were not changed.

## Tool environment

The isolated WSL2 environment passed these version gates: Python 3.12.14,
OrthoFinder 3.1.5, SnpEff package 5.4.0c (CLI 5.4c), BCFtools 1.24, NCBI
Datasets CLI 18.36.0, Matplotlib 3.11.1, and OpenJDK 21.0.10. Java 21 is used
because the pinned Bioconda SnpEff package requires it.

Machine-readable checksums, feature counts, tool versions, database results,
and portable commands are in `results/source_review/`. Large source packages,
raw SnpEff logs, and generated databases remain under ignored `local/` paths.

