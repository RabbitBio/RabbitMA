#!/usr/bin/env python3
"""Regression tests for cross-outer-k local input reuse analysis."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'benchmarks/cross_outer_local_reuse.py'
SPEC = importlib.util.spec_from_file_location('cross_outer_local_reuse',
                                               str(SCRIPT))
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


HEADER = ('low\thigh\treads_low\treads_high\tpotential\tactual\treads\t'
          'endpoint\trounds\n')


class CrossOuterLocalReuseTest(unittest.TestCase):
    def write(self, root, name, rows):
        path = root / name
        path.write_text(HEADER + ''.join('\t'.join(map(str, row)) + '\n'
                                         for row in rows))
        return path

    def test_full_and_read_only_prefix_work(self):
        with tempfile.TemporaryDirectory(prefix='local-reuse-') as name:
            root = Path(name)
            previous = self.write(root, 'local_k21_to_k29.tsv', [
                ('a', '1', '10', '1', 200, 100, 5, 200, 2),
                ('b', '2', '20', '2', 400, 300, 8, 200, 3),
                ('c', '3', '30', '3', 100, 50, 3, 200, 1),
            ])
            current = self.write(root, 'local_k29_to_k39.tsv', [
                ('a', '1', '10', '1', 300, 150, 5, 200, 3),
                ('d', '4', '20', '2', 500, 300, 8, 200, 3),
                ('e', '5', '40', '4', 600, 400, 9, 200, 4),
            ])
            rows = ANALYZER.analyze([
                ANALYZER.load_round(previous),
                ANALYZER.load_round(current),
            ])
            self.assertEqual(len(rows), 1)
            result = rows[0]
            self.assertEqual(result['transition'], '29->39')
            self.assertEqual(result['total'], {
                'tasks': 3, 'actual': 850, 'potential': 1400, 'rounds': 10})
            self.assertEqual(result['full_input']['matched_tasks'], 1)
            self.assertEqual(result['full_input']['repeated_prefix_actual'],
                             100)
            self.assertEqual(result['full_input']['repeated_prefix_rounds'],
                             2)
            self.assertEqual(result['read_set']['matched_tasks'], 2)
            self.assertEqual(result['read_set']['repeated_prefix_actual'],
                             400)
            self.assertEqual(result['read_set']['repeated_prefix_rounds'], 5)
            aggregate = ANALYZER.aggregate(rows)
            self.assertEqual(aggregate['total']['actual'], 850)
            self.assertEqual(
                aggregate['read_set']['repeated_prefix_actual'], 400)

    def test_duplicate_fingerprints_match_expensive_tasks_first(self):
        old = (
            ANALYZER.Record((1, 1), (2, 2), 10, 10, 1, 1, 1),
            ANALYZER.Record((1, 1), (2, 2), 90, 90, 1, 1, 9),
        )
        new = (ANALYZER.Record((3, 3), (2, 2), 80, 80, 1, 1, 8),)
        stats = ANALYZER.reuse_stats(old, new, 'reads')
        self.assertEqual(stats['matched_tasks'], 1)
        self.assertEqual(stats['repeated_prefix_rounds'], 8)
        self.assertEqual(stats['repeated_prefix_actual'], 80)


if __name__ == '__main__':
    unittest.main()
