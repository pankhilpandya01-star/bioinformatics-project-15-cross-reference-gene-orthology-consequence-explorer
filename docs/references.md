# References

Accessed sources and software documentation used to design and reproduce this
project are listed below. Software versions are pinned in `environment.yml`.

## Biological data

- NCBI Sequence Read Archive. SRX10301019, paired-end whole-genome sequencing
  experiment for SRR13921545:
  https://www.ncbi.nlm.nih.gov/sra/SRX10301019
- NCBI RefSeq. *Escherichia coli* K-12 substrain MG1655 assembly
  GCF_000005845.2:
  https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/
- NCBI Datasets. Genome download command reference:
  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/
- NCBI Datasets. Genome data-package contents:
  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/

Exact accessions, retrieval dates, file checksums, and inherited Project 13
eligibility labels are recorded in `data/reference_manifest.csv` and
`results/authentic_comparison/run_manifest.json`.

## Orthology

- Emms, D. M., and Kelly, S. (2019). OrthoFinder: phylogenetic orthology
  inference for comparative genomics. *Genome Biology*, 20, 238.
  https://doi.org/10.1186/s13059-019-1832-y
- OrthoFinder installation guide:
  https://orthofinder.github.io/OrthoFinder/download_and_install/
- OrthoFinder results guide:
  https://orthofinder.github.io/OrthoFinder/tutorials/guide-to-results/

## Variant annotation and validation

- Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L.,
  Land, S. J., Lu, X., and Ruden, D. M. (2012). A program for annotating and
  predicting the effects of single nucleotide polymorphisms, SnpEff.
  *Fly*, 6(2), 80–92. https://doi.org/10.4161/fly.19695
- SnpEff database-building guide:
  https://pcingola.github.io/SnpEff/snpeff/build_db/
- SnpEff command-line options:
  https://pcingola.github.io/SnpEff/snpeff/commandline/
- SnpEff ANN field specification:
  https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
- Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools.
  *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
- BCFtools 1.24 manual:
  https://www.htslib.org/doc/1.24/bcftools.html

## Runtime and plotting

- Python 3.12 documentation: https://docs.python.org/3.12/
- Matplotlib documentation: https://matplotlib.org/stable/
- GitHub Actions checkout action:
  https://github.com/actions/checkout
- setup-micromamba action:
  https://github.com/mamba-org/setup-micromamba
