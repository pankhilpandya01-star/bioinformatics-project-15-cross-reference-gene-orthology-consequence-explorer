# Data provenance

## Public sample

The read cohort originates from NCBI SRA experiment
[`SRX10301019`](https://www.ncbi.nlm.nih.gov/sra/SRX10301019), run
`SRR13921545`. The source record describes paired NovaSeq whole-genome
sequencing associated with *Escherichia coli* K-12 MG1655.

Project 15 does not independently assert the sample's strain identity. It reuses
Project 13's complete trimmed cohort and its six reference-specific accepted
callsets to test reference sensitivity.

## Inherited lineage

```text
SRR13921545 raw paired reads
        |
        v
Project 11 paired preprocessing design
        |
        v
Project 12 accepted MG1655-relative small variants
        |
        v
Project 13 six-reference accepted VCFs + callable BED files
        |
        v
Project 14 MG1655 annotations + 80 review candidates
        |
        v
Project 15 gene orthology and consequence-burden comparison
```

The MG1655 accepted VCF SHA-256 remains
`5a15f3a7df057a686c8f9eebc805be804d02d217c547c82c0837711572063c0f`.
Project 15 reproduces Project 14's 30,973 variants and 31,059 effects exactly.

## Reference assemblies

| Role | Accession | Label |
|---|---|---|
| Baseline | `GCF_000005845.2` | MG1655 |
| Alternative | `GCF_003697165.2` | DSM 30083 / ATCC 11775 |
| Alternative | `GCF_036503815.1` | 2017-02-2CC |
| Alternative | `GCF_005843885.1` | E4742 |
| Alternative | `GCF_000026225.1` | ATCC 35469 |
| Alternative | `GCF_002900365.1` | HT073016 |

Each NCBI Datasets package contains the exact genome FASTA, GenBank, GFF3, CDS,
protein, and sequence-report files for its accession version. Acquisition was
recorded on 2026-09-13.

The authoritative machine-readable records are:

- [`data/reference_manifest.csv`](../data/reference_manifest.csv) for inherited
  reference and VCF hashes;
- [`results/source_review/acquisition_manifest.csv`](../results/source_review/acquisition_manifest.csv)
  for package files, URLs, retrieval date, sizes, and hashes;
- [`results/source_review/source_validation.csv`](../results/source_review/source_validation.csv)
  for sequence and annotation agreement; and
- [`results/authentic_comparison/run_manifest.json`](../results/authentic_comparison/run_manifest.json)
  for the final source/output hashes and totals.

## Eligibility carried from Project 13

MG1655 was the only reference that met both whole-genome mapping and callability
eligibility rules. The formal Project 13 decision remained inconclusive because
there was no eligible runner-up for the required density comparison. All five
alternatives are therefore labelled sensitivity-only in Project 15.

## Storage boundary

The public tree contains compact checksums, controlled fixtures, summaries,
comparison tables, test evidence, and the dashboard. The following remain under
ignored local storage:

- downloaded NCBI packages;
- complete accepted and annotated VCF working copies;
- complete annotation-effect tables and SnpEff reports;
- generated SnpEff databases; and
- OrthoFinder alignments, trees, and working files.

No source file is modified by the workflow. Hashes are checked before and after
processing.
