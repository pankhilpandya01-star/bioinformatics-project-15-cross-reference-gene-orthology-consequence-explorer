from unittest import TestCase

from comparative_gene_explorer.callability import callable_overlap_bases, interval_bases
from comparative_gene_explorer.validation import merge_intervals


class CallableIntervalTests(TestCase):
    def test_counts_half_open_intersections_without_double_counting(self) -> None:
        query = ((0, 10), (90, 100))
        callable_intervals = ((5, 15), (95, 100))

        self.assertEqual(callable_overlap_bases(query, callable_intervals), 10)
        self.assertEqual(interval_bases(query), 20)

    def test_merges_overlapping_and_adjacent_cds_segments(self) -> None:
        self.assertEqual(
            merge_intervals([(10, 20), (0, 5), (5, 12), (30, 40)]),
            ((0, 20), (30, 40)),
        )

