# Publication checklist

Local preparation date: 2026-09-14.

The public repository was created and the audited candidate was published after
explicit authorization.

## Completed locally

- [x] Verified all six reference accessions and inherited checksums.
- [x] Confirmed FASTA, GenBank, GFF3, protein, VCF, and callable-BED
      consistency for every reference.
- [x] Passed CDS and protein validation for all six local SnpEff databases.
- [x] Annotated all 357,581 accepted variants and retained all 357,959 effects.
- [x] Reproduced the Project 14 MG1655 totals: 30,973 variants, 31,059
      canonical effects, and 80 review candidates.
- [x] Accounted for all 26,768 proteins exactly once in an orthogroup or the
      unassigned set.
- [x] Required unambiguous protein-to-gene mappings for every reported
      one-to-one comparison.
- [x] Confirmed callable CDS bases never exceed total CDS bases.
- [x] Produced 23,255 unique MG1655-to-alternative gene comparisons.
- [x] Produced five alternative-reference contexts for each of the 80
      Project 14 review candidates.
- [x] Reconciled VCF, CSV, OrthoFinder, dashboard, and manifest totals.
- [x] Passed the complete 58-test offline Python 3.12 suite.
- [x] Tested cleanup after controlled workflow failures.
- [x] Kept large NCBI packages, SnpEff databases, full annotated VCFs, and
      OrthoFinder working files under the ignored `local/` directory.
- [x] Used portable commands and relative public manifest paths.
- [x] Documented that orthologous genes do not imply equivalent variants or
      identical function.
- [x] Added the MIT license and a read-only GitHub Actions test workflow.
- [x] Connected Projects 11–14 and the planned Project 16 in the README.

## Final local audit

- [x] Recorded the exact public candidate file count, total size, and largest
      file in `results/publication_audit/audit_summary.json`.
- [x] Confirmed the candidate tree contains no private paths, credentials,
      personal metadata, or unwanted authorship wording.
- [x] Confirmed every local README and documentation link resolves; repository
      links name the intended future public URL.

## Publication status

- [x] Created the public GitHub repository.
- [x] Reviewed the exact staged file list.
- [x] Committed and pushed the publication candidate.
- [x] Confirmed the hosted GitHub Actions run passed.
- [x] Checked README tables, links, and dashboard rendering on GitHub.
- [ ] Add Project 15 to the portfolio and earlier-project navigation.
