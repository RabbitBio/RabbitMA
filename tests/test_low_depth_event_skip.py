#!/usr/bin/env python3
"""Exercise exact skipping of no-op local-low-depth threshold passes."""

import os
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest


CORE = Path(os.environ.get(
    'MEGAHIT_TEST_CORE', Path(__file__).resolve().parents[1] / 'build/megahit_core'))


class LowDepthEventSkipTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='megahit-low-depth-')
        self.root = Path(self.temporary.name)
        rng = random.Random(20260920)
        dna = lambda size: ''.join(rng.choice('ACGT') for _ in range(size))
        prefix, suffix = dna(120), dna(120)
        paths = [(prefix + dna(60) + suffix, 100),
                 (prefix + dna(60) + suffix, 10)]
        contigs = self.root / 'input.fa'
        contigs.write_text(''.join(
            f'>k29_{i} flag=0 multi={depth} len={len(sequence)}\n{sequence}\n'
            for i, (sequence, depth) in enumerate(paths)))
        Path(str(contigs) + '.info').write_text('2 600\n')
        for name in ('bubble.fa', 'local.fa'):
            path = self.root / name
            path.write_text('')
            Path(str(path) + '.info').write_text('0 0\n')

    def tearDown(self):
        self.temporary.cleanup()

    def run_core(self, command, enable_event_skip=False):
        env = os.environ.copy()
        env['MEGAHIT_PROFILE_LOW_DEPTH'] = '1'
        if enable_event_skip:
            env['MEGAHIT_EXPERIMENTAL_LOW_DEPTH_EVENT_SKIP'] = '1'
        result = subprocess.run([str(CORE)] + list(map(str, command)),
                                capture_output=True, text=True, env=env,
                                timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stderr

    def test_noop_thresholds_are_skipped_without_changing_outputs(self):
        graph = self.root / 'graph'
        self.run_core([
            'seq2sdbg', '-k', 39, '--kmer_from', 29,
            '--contig', self.root / 'input.fa',
            '--bubble', self.root / 'bubble.fa',
            '--local_contig', self.root / 'local.fa',
            '--host_mem', 1000000000, '--mem_flag', 0,
            '-t', 4, '-o', graph])

        logs = {}
        for name, enabled in (('baseline', False), ('optimized', True)):
            output = self.root / name
            logs[name] = self.run_core([
                'assemble', '-s', graph, '-o', output, '-t', 4,
                '--min_depth', 2, '--prune_level', 1,
                '--bubble_level', 0, '--disconnect_ratio', 0,
                '--max_tip_len', 200], enable_event_skip=enabled)

        pattern = re.compile(
            r'Low-depth iteration \d+: threshold=[\d.]+, removed=(\d+)')
        baseline_removed = [int(x) for x in pattern.findall(logs['baseline'])]
        optimized_removed = [int(x) for x in pattern.findall(logs['optimized'])]
        self.assertEqual(len(baseline_removed), 19)
        self.assertEqual(len(optimized_removed), 3)
        self.assertEqual(sum(baseline_removed), 1)
        self.assertEqual(sum(optimized_removed), 1)
        self.assertIn('skipped=16 no-op thresholds, next=10.1089',
                      logs['optimized'])

        for suffix in ('.contigs.fa', '.contigs.fa.mgb', '.addi.fa',
                       '.addi.fa.mgb', '.bubble_seq.fa',
                       '.bubble_seq.fa.mgb', '.final.contigs.fa'):
            baseline = Path(str(self.root / 'baseline') + suffix)
            optimized = Path(str(self.root / 'optimized') + suffix)
            self.assertEqual(baseline.read_bytes(), optimized.read_bytes(), suffix)


if __name__ == '__main__':
    unittest.main()
