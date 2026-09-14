# Limitations

## Gene-level comparison, not variant liftover

Project 15 compares burden in corresponding genes. It does not align the six
genomes or transform exact variant coordinates and alleles between references.
Two variants associated with orthologous genes must not be described as the
same mutation without a separate whole-genome alignment and allele check.

## Reference dependence

Mapping, callability, variant calling, annotation, and gene burden all depend on
the selected reference. The alternative callsets contain different coordinate
systems and different callable regions. Normalizing by callable CDS reduces one
bias but cannot make the analyses reference-independent.

Only MG1655 passed Project 13's whole-genome mapping and callability eligibility
rules. The five alternatives are sensitivity analyses, not equally supported
primary references.

## Orthology uncertainty

OrthoFinder infers evolutionary relationships from the supplied protein sets
and gene trees. A one-to-one ortholog is not guaranteed to retain identical
function, regulation, expression, or phenotype. One-to-many and many-to-many
families are reported but excluded from direct burden comparison.

RNA genes, pseudogenes, and other proteinless features cannot receive
protein-based orthology. Ambiguous protein-to-gene mappings are retained in the
audit tables and cannot qualify as one-to-one comparisons.

## Annotation uncertainty

SnpEff consequences are sequence-based predictions from one annotation release.
`HIGH`, `MODERATE`, `LOW`, and `MODIFIER` are prioritization categories, not
measurements of gene function or organism phenotype. Annotation warnings remain
visible and are not silently treated as validated biology.

This project performs no resistance prediction, clinical interpretation,
functional validation, pathway enrichment, or gene-expression analysis.

## Callability threshold

The 90% callable-CDS rule is deliberately strict but still arbitrary. A gene
that passes can contain uncallable bases, and a gene that fails can contain
useful evidence. Results near the threshold may change under a different depth,
mapping-quality, or base-quality definition.

E4742, ATCC 35469, and HT073016 have only 9, 1, and 3 comparable gene pairs.
Their normalized burdens are displayed for completeness but cannot support a
broad genomic conclusion.

## Counting rules

One variant can be counted for two genuinely overlapping genes. Within one
gene, multiple SnpEff effects collapse to the highest-impact effect for the
burden table. Consequently, total variant–gene assignments are not expected to
equal either VCF rows or raw annotation-effect rows.

Density per 1,000 callable CDS bases describes reference-relative sequence
difference burden. It does not estimate mutation rate, selection, or biological
damage.

## Scope exclusions

The project does not include remapping, variant recalling, exact coordinate
liftover, pangenome graphs, structural variants, copy-number analysis, pathway
analysis, resistance prediction, phenotype inference, or phylogenetic claims.

Project 16 can add whole-genome alignment and exact variant-coordinate liftover
with explicit handling of non-alignable and structurally rearranged regions.
