# Milestone 4 — Authentic comparison review

## Decision

Milestone 4 passed. All six Project 13 accepted callsets were annotated against
their exact reference annotations, joined to the Milestone 3 orthology results,
and compared using callable CDS bases as the denominator.

The result supports keeping MG1655 as the primary interpretation reference. It
does not establish coordinate-equivalent variants between references, and the
five alternatives remain sensitivity analyses because none passed Project 13's
whole-genome eligibility rule.

## Annotation accounting

| Reference | Accepted variants | Predicted effects | Coding genes meeting 90% callability |
|---|---:|---:|---:|
| MG1655, GCF_000005845.2 | 30,973 | 31,059 | 3,857 of 4,308 |
| DSM 30083 / ATCC 11775, GCF_003697165.2 | 92,665 | 92,706 | 3,624 of 4,823 |
| 2017-02-2CC, GCF_036503815.1 | 125,327 | 125,378 | 3,501 of 4,868 |
| E4742, GCF_005843885.1 | 49,201 | 49,230 | 13 of 4,873 |
| ATCC 35469, GCF_000026225.1 | 15,264 | 15,266 | 1 of 4,362 |
| HT073016, GCF_002900365.1 | 44,151 | 44,320 | 3 of 4,688 |

The combined callsets contain 357,581 accepted variants and 357,959 predicted
effects. All twelve CDS/protein database checks remain below the 2% error gate.
MG1655 reproduces Project 14 exactly: all 31,059 canonical effect records are
identical, not merely equal in count.

The analysis creates one record per distinct `(variant, gene)` and retains that
pair's highest-impact effect. This yields 327,542 variant–gene assignments.
There are 30,382 effects without an exact gene-catalog key, primarily
intergenic annotations, and 35 additional effects collapse into an already
counted variant–gene pair. These three categories account for all 357,959
effects exactly.

## Comparable gene burden

A gene pair is comparable only when it is one-to-one and at least 90% of the
combined CDS bases are callable in both references.

| Alternative | Comparable gene pairs | MG1655 burden | Alternative burden | Interpretation |
|---|---:|---:|---:|---|
| DSM 30083 / ATCC 11775 | 3,206 | 7.42 | 23.90 | Alternative is 3.2× higher |
| 2017-02-2CC | 3,137 | 7.31 | 32.00 | Alternative is 4.4× higher |
| E4742 | 9 | 2.99 | 46.46 | Too few comparable genes |
| ATCC 35469 | 1 | 0.00 | 0.00 | Too few comparable genes |
| HT073016 | 3 | 1.71 | 37.10 | Too few comparable genes |

Burden is the number of distinct variant–gene assignments per 1,000 callable
CDS bases. The first two alternatives provide thousands of comparable gene
pairs and show substantially more reference-relative differences than MG1655
over matched gene sets. This is compatible with MG1655 being the least-divergent
tested reference, but it is not proof of sample identity.

The final three alternatives cannot support a broad comparison. Although they
contain thousands of one-to-one orthologs, Project 13's weak mapping and
callability leave only 9, 1, and 3 gene pairs that pass the fixed comparability
gate. Their plotted burdens are retained as audit evidence and must not be
generalized to the genome.

## Project 14 review candidates

Each of Project 14's 80 review candidates is connected to the corresponding
gene context in all five alternatives, producing 400 audit rows.

- DSM 30083 / ATCC 11775: 28 candidate contexts are comparable; 15 corresponding
  genes contain at least one high-impact alternative-reference annotation.
- 2017-02-2CC: 34 candidate contexts are comparable; 23 corresponding genes
  contain at least one high-impact alternative-reference annotation.
- E4742, ATCC 35469, and HT073016: none of the 80 candidate contexts passes the
  complete orthology-plus-callability gate.

These are gene-level signals. A high-impact annotation in a corresponding gene
does not mean the original MG1655-relative variant is present at an equivalent
coordinate.

## Outputs

Compact public-facing results are in `results/authentic_comparison/`:

- `gene_consequence_burden.csv`, one row per gene and reference;
- `baseline_ortholog_comparison.csv`, one row per MG1655 gene and alternative;
- `review_candidate_orthology.csv`, five contexts for each Project 14 candidate;
- eight summary tables, exact source/output checksums, and the run manifest; and
- the four-panel comparison dashboard.

The 266 MB annotation workspace, annotated VCFs, complete effect table, and
SnpEff reports remain under ignored `local/` storage as required.

## Verification

The Python 3.12 suite contains 35 passing tests. Authentic checks cover Project
13 callset totals, exact Project 14 reproduction, variant–gene accounting,
callable-base bounds, unique comparison keys, strict comparability, density
summary agreement, five contexts per review candidate, and dashboard dimensions.

All source hashes were rechecked after analysis and remain unchanged. No
repository was created or published during this milestone.
