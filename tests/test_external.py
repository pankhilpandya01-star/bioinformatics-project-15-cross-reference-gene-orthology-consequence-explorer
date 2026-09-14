import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from comparative_gene_explorer.external import ExternalToolError, run_command


class ExternalCommandTests(TestCase):
    def test_failure_message_sanitizes_replaced_paths(self) -> None:
        with TemporaryDirectory() as directory:
            sensitive_path = str(Path(directory) / "private-input")
            with self.assertRaises(ExternalToolError) as caught:
                run_command(
                    [
                        sys.executable,
                        "-c",
                        f"import sys; print({sensitive_path!r}, file=sys.stderr); raise SystemExit(7)",
                    ],
                    cwd=Path(directory),
                    replacements=((sensitive_path, "<INPUT>"),),
                    reject_markers=(),
                )
            self.assertNotIn(sensitive_path, str(caught.exception))
            self.assertIn("<INPUT>", str(caught.exception))

    def test_logged_output_is_sanitized(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            sensitive_path = str(root / "private-input")
            log = root / "command"
            run_command(
                [sys.executable, "-c", f"print({sensitive_path!r})"],
                cwd=root,
                stdout_path=log,
                replacements=((sensitive_path, "<INPUT>"),),
                reject_markers=(),
            )
            self.assertEqual(log.read_text(encoding="utf-8").strip(), "<INPUT>")
