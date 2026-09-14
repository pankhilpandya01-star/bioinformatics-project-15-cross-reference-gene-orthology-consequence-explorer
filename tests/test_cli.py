from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from comparative_gene_explorer.cli import build_parser, main


class CliTests(TestCase):
    def test_documented_defaults(self) -> None:
        args = build_parser().parse_args(
            ["--reference-manifest", "manifest.csv", "--output-dir", "output"]
        )
        self.assertEqual(args.baseline_accession, "GCF_000005845.2")
        self.assertEqual(args.sample_id, "SRR13921545")
        self.assertEqual(args.minimum_gene_callable_fraction, 0.90)
        self.assertEqual(args.maximum_database_error_rate, 0.02)
        self.assertEqual(args.upstream_downstream_length, 0)
        self.assertEqual(args.threads, 1)

    @patch("comparative_gene_explorer.cli.run_workflow")
    def test_cli_passes_public_arguments(self, run_workflow) -> None:
        run_workflow.return_value = Path("output")
        status = main(
            [
                "--reference-manifest",
                "manifest.csv",
                "--baseline-accession",
                "CTRL.1",
                "--sample-id",
                "sample",
                "--minimum-gene-callable-fraction",
                "0.8",
                "--maximum-database-error-rate",
                "0.01",
                "--upstream-downstream-length",
                "5",
                "--threads",
                "2",
                "--output-dir",
                "output",
            ]
        )
        self.assertEqual(status, 0)
        config = run_workflow.call_args.args[0]
        self.assertEqual(config.baseline_accession, "CTRL.1")
        self.assertEqual(config.minimum_gene_callable_fraction, 0.8)
        self.assertEqual(config.threads, 2)

