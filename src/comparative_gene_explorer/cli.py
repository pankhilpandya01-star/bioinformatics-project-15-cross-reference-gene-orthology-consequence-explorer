"""Command-line interface for cross-reference gene comparison."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from comparative_gene_explorer import __version__
from comparative_gene_explorer.config import WorkflowConfig
from comparative_gene_explorer.workflow import WorkflowError, run_workflow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="comparative-gene-explorer",
        description=(
            "Annotate reference-specific callsets and measure gene-level callable CDS "
            "coverage before orthology comparison."
        ),
    )
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--minimum-gene-callable-fraction", type=float, default=0.90)
    parser.add_argument("--maximum-database-error-rate", type=float, default=0.02)
    parser.add_argument("--upstream-downstream-length", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = run_workflow(
            WorkflowConfig(
                reference_manifest=args.reference_manifest,
                output_dir=args.output_dir,
                baseline_accession=args.baseline_accession,
                sample_id=args.sample_id,
                minimum_gene_callable_fraction=args.minimum_gene_callable_fraction,
                maximum_database_error_rate=args.maximum_database_error_rate,
                upstream_downstream_length=args.upstream_downstream_length,
                threads=args.threads,
            )
        )
    except WorkflowError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(f"Completed comparative annotation core: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

