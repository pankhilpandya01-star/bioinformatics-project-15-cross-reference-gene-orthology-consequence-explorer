"""Run Project 15's deterministic OrthoFinder milestone."""

from __future__ import annotations

import argparse
from pathlib import Path

from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.orthology import (
    references_from_validated,
    run_orthology_workflow,
)
from comparative_gene_explorer.validation import validate_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    config = WorkflowConfig(
        reference_manifest=args.reference_manifest,
        output_dir=args.output_dir,
        baseline_accession=args.baseline_accession,
        sample_id=args.sample_id,
        threads=args.threads,
    )
    inputs = validate_inputs(config)
    output = run_orthology_workflow(
        references_from_validated(inputs.references),
        baseline_accession=args.baseline_accession,
        threads=args.threads,
        work_directory=args.work_dir,
        output_directory=args.output_dir,
    )
    print(f"Completed orthology workflow: {output}")


if __name__ == "__main__":
    main()

