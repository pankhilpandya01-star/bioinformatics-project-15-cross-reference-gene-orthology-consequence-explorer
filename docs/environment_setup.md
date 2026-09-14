# Linux and WSL2 setup

## Supported execution model

The workflow uses Linux-native OrthoFinder, SnpEff, BCFtools, and NCBI Datasets
commands. It was run in Ubuntu through WSL2 on Windows and is also suitable for
a regular Linux host. GitHub Actions uses Ubuntu.

The Python source can be edited from Windows, but the scientific commands and
tests should run inside the pinned Linux environment.

## One-time WSL2 preparation

If WSL2 is not installed, open an administrator PowerShell session and run:

```powershell
wsl --install -d Ubuntu
```

Restart if requested, open Ubuntu once, and complete the Linux account setup.
Install micromamba using its official installation instructions.

## Create the environment

From the project root in a Linux shell:

```bash
micromamba create -f environment.yml
micromamba activate comparative-gene-explorer
python -m pip install --no-deps -e .
```

Verify every pinned tool:

```bash
python --version
orthofinder --version
snpEff -version
bcftools --version
datasets version
java -version
python -c "import matplotlib; print(matplotlib.__version__)"
```

Expected versions are Python 3.12, OrthoFinder 3.1.5, SnpEff 5.4.0c (CLI may
report 5.4c), BCFtools 1.24, NCBI Datasets CLI 18.36.0, OpenJDK 21, and
Matplotlib 3.11.1.

OpenJDK 21 is intentional because the pinned Bioconda SnpEff package uses the
Java 21 runtime.

## Portfolio directory layout

The included reference manifest uses relative paths to Project 13 outputs. The
authentic comparison also reads the compact Project 14 effect and review tables.
The simplest layout is:

```text
portfolio/
  bioinformatics-project-13-sample-identity-reference-concordance-explorer/
  bioinformatics-project-14-variant-annotation-gene-consequence-explorer/
  bioinformatics-project-15-cross-reference-gene-orthology-consequence-explorer/
```

If the projects are stored elsewhere, update only the path columns in a copy of
`data/reference_manifest.csv` and pass the correct upstream paths to
`scripts/run_comparative_analysis.py`. Do not change the recorded hashes.

## Acquire exact annotation packages

The source-review script downloads the exact NCBI accession versions, extracts
the required files, verifies them, and builds the six local databases:

```bash
python scripts/run_source_review.py \
  --reference-manifest data/reference_manifest.csv \
  --source-dir local/source_packages \
  --database-dir local/snpeff_databases \
  --output-dir results/source_review \
  --retrieved YYYY-MM-DD
```

Compare newly retrieved files with
`results/source_review/acquisition_manifest.csv`. NCBI annotation content may be
updated while an assembly accession remains unchanged, so a hash difference
requires a new source review rather than silent substitution.

## Run the workflow stages

Annotation and callable-CDS analysis:

```bash
comparative-gene-explorer \
  --reference-manifest data/reference_manifest.csv \
  --baseline-accession GCF_000005845.2 \
  --sample-id SRR13921545 \
  --minimum-gene-callable-fraction 0.90 \
  --maximum-database-error-rate 0.02 \
  --upstream-downstream-length 0 \
  --threads 1 \
  --output-dir local/authentic_annotation
```

Orthology:

```bash
python scripts/run_orthology.py \
  --reference-manifest data/reference_manifest.csv \
  --baseline-accession GCF_000005845.2 \
  --sample-id SRR13921545 \
  --threads 1 \
  --work-dir local/orthofinder_authentic \
  --output-dir results/orthology_analysis
```

Final compact comparison:

```bash
python scripts/run_comparative_analysis.py \
  --annotation-dir local/authentic_annotation \
  --orthology-dir results/orthology_analysis \
  --project14-review-candidates ../bioinformatics-project-14-variant-annotation-gene-consequence-explorer/results/authentic_annotation/analysis/review_candidates.csv \
  --project14-annotation-effects ../bioinformatics-project-14-variant-annotation-gene-consequence-explorer/results/authentic_annotation/analysis/annotation_effects.csv \
  --project13-reference-concordance ../bioinformatics-project-13-sample-identity-reference-concordance-explorer/results/concordance_analysis/reference_concordance.csv \
  --output-dir results/authentic_comparison
```

Every output directory must be new. The workflow will not overwrite a previous
run.

## Run tests

The committed suite is offline and does not download public data:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp/matplotlib
export PYTHONHASHSEED=0
python -m unittest discover -s tests -v
```

## Troubleshooting

- If a checksum fails, restore the exact inherited file or perform a documented
  source review; do not replace the recorded hash silently.
- If a VCF index is missing, create exactly one CSI or TBI index with BCFtools.
- If FASTA, GenBank, GFF3, VCF, or BED contigs disagree, verify that every file
  belongs to the same accession version.
- If SnpEff validation exceeds 2%, inspect the GenBank/CDS/protein package and
  do not continue with annotation.
- If OrthoFinder finds an incomplete existing work directory, move that failed
  directory aside and rerun with a new path.
- If the output directory exists, choose a new directory rather than deleting
  evidence from an earlier run.
- If a native version differs, recreate the environment from `environment.yml`.
