from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from scripts.run_source_review import ReviewError, validate_gff


class GffCoordinateValidationTests(TestCase):
    def test_accepts_ncbi_circular_origin_representation(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "annotation.gff3"
            path.write_text(
                "##gff-version 3\n"
                "chr\tRefSeq\tregion\t1\t100\t.\t+\t.\tID=chr;Is_circular=true\n"
                "chr\tRefSeq\tgene\t90\t120\t.\t+\t.\tID=gene-crossing-origin\n",
                encoding="utf-8",
            )

            counts, circular_contigs, origin_spanning = validate_gff(
                path, {"chr": "A" * 100}
            )

            self.assertEqual(counts["gene"], 1)
            self.assertEqual(circular_contigs, {"chr"})
            self.assertEqual(origin_spanning, 1)

    def test_rejects_out_of_bounds_linear_feature(self) -> None:
        with TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "annotation.gff3"
            path.write_text(
                "##gff-version 3\n"
                "chr\tRefSeq\tregion\t1\t100\t.\t+\t.\tID=chr\n"
                "chr\tRefSeq\tgene\t90\t120\t.\t+\t.\tID=invalid-gene\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ReviewError, "non-circular"):
                validate_gff(path, {"chr": "A" * 100})

