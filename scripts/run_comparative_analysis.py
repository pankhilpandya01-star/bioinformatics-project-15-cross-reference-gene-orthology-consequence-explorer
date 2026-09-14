"""Run the authentic Project 15 gene-level comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

from comparative_gene_explorer.comparison import (
    ComparisonConfig,
    run_comparative_analysis,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-dir", required=True, type=Path)
    parser.add_argument("--orthology-dir", required=True, type=Path)
    parser.add_argument("--project14-review-candidates", required=True, type=Path)
    parser.add_argument("--project14-annotation-effects", required=True, type=Path)
    parser.add_argument("--project13-reference-concordance", required=True, type=Path)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--minimum-gene-callable-fraction", type=float, default=0.90)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output = run_comparative_analysis(
        ComparisonConfig(
            annotation_directory=args.annotation_dir,
            orthology_directory=args.orthology_dir,
            project14_review_candidates=args.project14_review_candidates,
            project14_annotation_effects=args.project14_annotation_effects,
            project13_reference_concordance=args.project13_reference_concordance,
            output_directory=args.output_dir,
            baseline_accession=args.baseline_accession,
            minimum_gene_callable_fraction=args.minimum_gene_callable_fraction,
        )
    )
    print(f"Completed authentic comparison: {output}")


if __name__ == "__main__":
    main()
