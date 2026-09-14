# Milestone 3 — Orthology workflow review

## Decision

Milestone 3 passed. The six exact reference proteomes were prepared with
accession-prefixed identifiers and analysed in one complete OrthoFinder 3.1.5
run. All 26,768 prepared proteins occur exactly once in either an assigned
orthogroup or the unassigned set. The workflow is paused before consequence
burden comparison.

Orthology identifies evolutionary relationships between proteins. It does not
prove identical function, and it does not establish that variants at different
reference coordinates are the same event.

## Reproducible command

The analysis ran in the pinned WSL2 environment with one search and analysis
thread:

```bash
python scripts/run_orthology.py \
  --reference-manifest data/reference_manifest.csv \
  --baseline-accession GCF_000005845.2 \
  --sample-id SRR13921545 \
  --threads 1 \
  --work-dir local/orthofinder_authentic \
  --output-dir results/orthology_analysis
```

The underlying OrthoFinder configuration was:

```text
orthofinder -f proteomes -M msa -S diamond -A famsa -T fasttree -t 1 -a 1 -X -o results
```

`-X` preserves the prepared identifiers. The full OrthoFinder working tree
remains under the ignored `local/` directory; the public result contains only
the parsed tables, checksums, version evidence, and run manifest.

## Identifier and protein-to-gene handling

Each prepared protein identifier begins with its exact assembly accession, so
proteins from different references cannot collide. NCBI protein accessions are
then mapped back to GFF3 gene records.

Of 26,768 proteins, 26,656 map unambiguously to one gene. A further 112 proteins
in the alternative references link to more than one gene record and are marked
`ambiguous_gene`; none is silently forced onto a gene. The MG1655 baseline has
no ambiguous protein mapping. Another 1,869 annotated genes have no protein and
remain explicit as `no_protein` records.

Only unambiguous protein-to-gene mappings can contribute to a reported
one-to-one gene relationship.

## Authentic results

OrthoFinder produced 4,935 assigned orthogroups containing 25,533 proteins and
1,235 unassigned single-protein groups. Parsing yielded 56,346 pairwise protein
orthology edges and 23,255 MG1655-to-alternative gene comparisons—one row for
each of 4,651 baseline genes against each of five alternatives.

| Alternative reference | One-to-one | One-to-many | Many-to-one | Many-to-many | No ortholog | No protein | Renamed one-to-one |
|---|---:|---:|---:|---:|---:|---:|---:|
| GCF_003697165.2 | 3,443 | 80 | 195 | 14 | 558 | 361 | 219 |
| GCF_036503815.1 | 3,492 | 86 | 167 | 20 | 525 | 361 | 233 |
| GCF_005843885.1 | 3,346 | 64 | 136 | 20 | 724 | 361 | 930 |
| GCF_000026225.1 | 3,163 | 91 | 207 | 11 | 818 | 361 | 971 |
| GCF_002900365.1 | 3,223 | 31 | 162 | 15 | 859 | 361 | 882 |

The one-to-one share ranges from 68.0% to 75.1% of MG1655 genes. In total,
2,716 MG1655 genes have a one-to-one relationship in all five alternatives.
Different gene symbols are recorded separately as renamed one-to-one cases;
a name difference does not imply a functional difference.

The remaining relationships are not treated as failures. Duplicated families,
missing orthologues, ambiguous mappings, and proteinless features are retained
so Milestone 4 can apply its comparability rules without hiding excluded genes.

## Controlled demonstration

The controlled three-proteome example runs the real pinned OrthoFinder binary
and recovers all required cases:

- a stable one-to-one gene pair;
- a one-to-one pair whose gene symbol changed;
- one-to-many, many-to-one, and many-to-many duplicated families;
- genes with no ortholog;
- proteins in the unassigned set; and
- a gene with no protein product.

Its 24 prepared proteins are all assigned exactly once. Compact evidence is in
`results/controlled_orthology/`.

## Verification

The complete Python 3.12 suite contains 28 passing tests. Authentic-result tests
check exact protein membership, pairwise orthogroup consistency, one row per
baseline-gene/reference combination, unambiguous one-to-one rows, and agreement
between the detailed and summary tables.

Source protein checksums remain recorded in the run manifest. No repository was
created or published during this milestone.
