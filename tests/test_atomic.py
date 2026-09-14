from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from comparative_gene_explorer.atomic import atomic_output_directory


class AtomicPublicationTests(TestCase):
    def test_failure_removes_staging_and_leaves_destination_absent(self) -> None:
        with TemporaryDirectory() as directory:
            destination = Path(directory) / "result"
            with self.assertRaisesRegex(RuntimeError, "controlled failure"):
                with atomic_output_directory(destination) as staging:
                    (staging / "partial.txt").write_text("partial", encoding="utf-8")
                    raise RuntimeError("controlled failure")
            self.assertFalse(destination.exists())
            self.assertEqual(list(Path(directory).iterdir()), [])

