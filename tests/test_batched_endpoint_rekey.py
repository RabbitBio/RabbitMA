#!/usr/bin/env python3
"""Endpoint rekey batching must preserve weak-link graph/contig semantics."""
import json
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('MEGAHIT_TEST_CORE', ROOT / 'build/megahit_core')).resolve()
sys.path.insert(0, str(ROOT / 'tools'))
from compare_contigs import compare


class BatchedEndpointRekeyTest(unittest.TestCase):
    def run_core(self, args, batched=False, delta=True, packed=False,
                 disable_packed=False):
        env = {k: v for k, v in os.environ.items() if not k.startswith('MEGAHIT_')}
        env['MEGAHIT_PROFILE_GRAPH_REFRESH'] = '1'
        env['MEGAHIT_EXPERIMENTAL_PACKED_LEGACY_ORDER'] = '1'
        if batched:
            env['MEGAHIT_EXPERIMENTAL_BATCHED_ENDPOINT_REKEY'] = '1'
        if not delta:
            env['MEGAHIT_DISABLE_DELTA_REFRESH'] = '1'
        if packed:
            env['MEGAHIT_FORCE_PACKED_ENDPOINT_KEYS'] = '1'
        if disable_packed:
            env['MEGAHIT_DISABLE_PACKED_ENDPOINT_KEYS'] = '1'
        result = subprocess.run([str(CORE)] + list(map(str, args)), env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def test_weak_branch_trimming_and_merge(self):
        rng = random.Random(61004)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        items = []
        # Long low-depth alternatives have both ends trimmed; hanging branches
        # have one. Circles and independent paths exercise unchanged owners.
        for index in range(48):
            sequence = dna(2000)
            items.append((sequence, 60))
            items.append((sequence[300:339] + dna(240) + sequence[1300:1339], 1))
            items.append((sequence[700:739] + dna(360), 1))
            tip_base = rng.choice([base for base in 'ACGT' if base != sequence[1039]])
            items.append((sequence[1000:1039] + tip_base, 1))
            if index % 8 == 0:
                circle = dna(700)
                items.append((circle + circle[:39], 25))
        with tempfile.TemporaryDirectory(prefix='rabbitma-rekey-') as directory:
            directory = Path(directory)
            contigs = directory / 'input.fa'
            contigs.write_text(''.join('>k21_%d flag=0 multi=%d len=%d\n%s\n' %
                                      (i, depth, len(sequence), sequence)
                                      for i, (sequence, depth) in enumerate(items)))
            Path(str(contigs)+'.info').write_text('%d %d\n' %
                       (len(items), sum(len(sequence) for sequence, _ in items)))
            bubble = directory / 'empty_bubble.fa'
            bubble.write_text('')
            Path(str(bubble)+'.info').write_text('0 0\n')
            graph = directory / 'graph'
            self.run_core(['seq2sdbg', '-k', 39, '--kmer_from', 21, '--contig', contigs,
                           '--bubble', bubble, '--host_mem', 1000000000,
                           '--mem_flag', 0, '-t', 4, '-o', graph])
            baseline = directory / 'reference'
            base = ['assemble', '-s', graph, '--min_depth', 1, '--prune_level', 0,
                    '--bubble_level', 0, '--max_tip_len', 0, '--cleaning_rounds', 2,
                    '--disconnect_ratio', 0.1, '--min_standalone', 0]
            reference_log = self.run_core(base + ['-t', 1, '-o', baseline])
            self.assertIn('sharded sparse edge-to-unitig map', reference_log)
            reference_counts = re.findall(r'Number unitigs disconnected: (\d+)', reference_log)
            self.assertGreater(int(reference_counts[0]), 100)
            exercised = []
            for threads, delta, packed in [(1, True, False), (4, True, False),
                                          (8, False, False), (1, True, True),
                                          (4, True, True), (8, False, True)]:
                with self.subTest(threads=threads, delta=delta, packed=packed):
                    output = directory / ('batched_%d_%d_%d' % (threads, delta, packed))
                    log = self.run_core(base + ['-t', threads, '-o', output], True, delta, packed)
                    if packed:
                        self.assertIn('Endpoint hash slot payload: 12 bytes', log)
                    moves = list(map(int, re.findall(r'Batched endpoint rekey: endpoints=(\d+)', log)))
                    self.assertTrue(moves, log)
                    self.assertGreater(sum(moves), 100)
                    self.assertGreater(int(reference_counts[0]), sum(moves),
                                       'Fixture must also delete short weak unitigs')
                    self.assertEqual(re.findall(r'Number unitigs disconnected: (\d+)', log), reference_counts)
                    # prune_level=0 isolates weak-link work and emits no addi file.
                    self.assertFalse(Path(str(baseline)+'.addi.fa.mgb').exists())
                    self.assertFalse(Path(str(output)+'.addi.fa.mgb').exists())
                    for suffix in ['contigs.fa.mgb', 'bubble_seq.fa.mgb']:
                        result = compare(str(baseline)+'.'+suffix, str(output)+'.'+suffix,
                                         mode='circular', metadata=True)
                        self.assertTrue(result['equal'], json.dumps(result))
                    exercised.append((threads, delta, packed, sum(moves)))
            # The rollback switch must restore 16-byte wide slots even when
            # the fixture forces wide addressing on this small graph.
            output = directory / 'wide_rollback'
            log = self.run_core(base + ['-t', 4, '-o', output],
                                batched=True, packed=True, disable_packed=True)
            self.assertIn('Endpoint hash slot payload: 16 bytes', log)
            self.assertEqual(re.findall(r'Number unitigs disconnected: (\d+)', log),
                             reference_counts)
            for suffix in ['contigs.fa.mgb', 'bubble_seq.fa.mgb']:
                result = compare(str(baseline)+'.'+suffix, str(output)+'.'+suffix,
                                 mode='circular', metadata=True)
                self.assertTrue(result['equal'], json.dumps(result))
            print('Endpoint rekeys exercised:', exercised)


if __name__ == '__main__':
    unittest.main()
