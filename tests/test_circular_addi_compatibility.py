#!/usr/bin/env python3
"""Circular change-only evidence must preserve v0.1.0 graph semantics."""
import hashlib
import json
import os
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest

from test_stable_contigs import records

ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('MEGAHIT_TEST_CORE', ROOT / 'build/megahit_core')).resolve()
FIXTURE = json.loads((ROOT / 'test_data/circular_addi_v010.json').read_text())


class CircularAddiCompatibilityTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='megahit-circular-addi-')
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write_contigs(self, name, items, k, packed=False):
        path = self.root / name
        Path(str(path) + '.info').write_text(
            '%d %d\n' % (len(items), sum(len(sequence) for sequence, _, _ in items)))
        if not packed:
            path.write_text(''.join(
                '>k%d_%d flag=%d multi=%s len=%d\n%s\n' %
                (k, i, flag, multi, len(sequence), sequence)
                for i, (sequence, flag, multi) in enumerate(items)))
        else:
            path.touch()
            with Path(str(path) + '.mgb').open('wb') as stream:
                stream.write(struct.pack('<8sII', b'MGCTG01\0', 1, 0))
                for i, (sequence, flag, multi) in enumerate(items):
                    stream.write(struct.pack('<IIifq', len(sequence), k, flag, multi, i))
                    words = [0] * ((len(sequence) + 15) // 16)
                    for offset, base in enumerate(sequence):
                        words[offset // 16] |= 'ACGT'.index(base) << (30 - 2 * (offset % 16))
                    stream.write(struct.pack('<' + 'I' * len(words), *words))
        return path

    def run_core(self, command, extra_env=None):
        env = os.environ.copy()
        env.update(extra_env or {})
        result = subprocess.run([str(CORE)] + list(map(str, command)),
                                capture_output=True, text=True, timeout=30, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)

    def check_loop_input(self, packed, primary):
        fixture = FIXTURE['extension']
        items = [(fixture['background'], 0, 5)]
        if primary:
            items.append((fixture['loop'], 3, 1))
        contigs = self.write_contigs('primary.fa', items, fixture['k_from'], packed)
        bubble = self.write_contigs('bubble.fa', [], fixture['k_from'], packed)
        addi = self.write_contigs('addi.fa', [(fixture['loop'], 3, 1)],
                                 fixture['k_from'], packed)
        expected = sorted(tuple(record) for record in fixture[
            'expected_primary_records' if primary else 'expected_addi_records'])
        for threads in (1, 4):
            with self.subTest(threads=threads, packed=packed, primary=primary):
                graph = self.root / ('graph%d' % threads)
                output = self.root / ('out%d' % threads)
                command = ['seq2sdbg', '-k', fixture['k_to'], '--kmer_from', fixture['k_from'],
                           '--contig', contigs, '--bubble', bubble, '--host_mem', 1000000000,
                           '--mem_flag', 0, '-t', threads, '-o', graph]
                if not primary:
                    command += ['--addi_contig', addi]
                self.run_core(command)
                self.run_core(['assemble', '-s', graph, '-o', output, '-t', threads,
                               '--min_depth', 1, '--prune_level', 0, '--max_tip_len', 0,
                               '--cleaning_rounds', 0, '--bubble_level', 0, '--is_final_round'])
                self.assertEqual(records(Path(str(output) + '.contigs.fa')), expected)

    def check_changed_loop(self, packed):
        fixture = FIXTURE['changed_loop']
        primary = self.write_contigs('primary.fa', [(fixture['circle'], 0, 20),
                                                    (fixture['branch'], 0, 1)],
                                     fixture['k_from'], packed)
        bubble = self.write_contigs('bubble.fa', [], fixture['k_from'], packed)
        expected = [tuple(record) for record in fixture['expected_addi_records']]
        for threads in (1, 4):
            with self.subTest(threads=threads, packed=packed):
                graph = self.root / ('graph%d' % threads)
                output = self.root / ('out%d' % threads)
                self.run_core(['seq2sdbg', '-k', fixture['k_to'], '--kmer_from', fixture['k_from'],
                               '--contig', primary, '--bubble', bubble, '--host_mem', 1000000000,
                               '--mem_flag', 0, '-t', threads, '-o', graph])
                self.run_core(['assemble', '-s', graph, '-o', output, '-t', threads,
                               '--min_depth', 2, '--prune_level', 1,
                               '--max_tip_len', fixture['max_tip_len'],
                               '--cleaning_rounds', 0, '--bubble_level', 0])
                # Pruning the low-depth branch creates a changed loop. Its
                # origin is evidence for the next k, not output formatting.
                self.assertEqual(records(Path(str(output) + '.addi.fa')), expected)

    def test_addi_fasta_preserves_finite_windows(self):
        self.check_loop_input(packed=False, primary=False)

    def test_addi_packed_preserves_finite_windows(self):
        self.check_loop_input(packed=True, primary=False)

    def test_primary_fasta_still_extends_loops(self):
        self.check_loop_input(packed=False, primary=True)

    def test_primary_packed_still_extends_loops(self):
        self.check_loop_input(packed=True, primary=True)

    def test_changed_fasta_loop_preserves_legacy_origin(self):
        self.check_changed_loop(packed=False)

    def test_changed_packed_loop_preserves_legacy_origin(self):
        self.check_changed_loop(packed=True)

    def test_parallel_merges_match_v010_serial_records(self):
        fixture = FIXTURE['parallel_merge']
        rng = random.Random(fixture['seed'])
        items = []
        for _ in range(fixture['circles']):
            ring = ''.join(rng.choices('ACGT', k=fixture['ring_length']))
            k = fixture['k_to']
            items.append((ring + ring[:k], 0, 20))
            for position, length, depth in fixture['branches']:
                first = rng.choice([base for base in 'ACGT' if base != ring[position + k]])
                branch = (ring[position:position + k] + first +
                          ''.join(rng.choices('ACGT', k=length - 1)))
                items.append((branch, 0, depth))
        primary = self.write_contigs('parallel.fa', items, fixture['k_from'])
        bubble = self.write_contigs('parallel-bubble.fa', [], fixture['k_from'])
        graph = self.root / 'parallel-graph'
        self.run_core(['seq2sdbg', '-k', fixture['k_to'], '--kmer_from', fixture['k_from'],
                       '--contig', primary, '--bubble', bubble, '--host_mem', 1000000000,
                       '--mem_flag', 0, '-t', 1, '-o', graph])
        for run, (threads, edge_first) in enumerate(((1, False), (4, False),
                                                    (4, False), (4, True))):
            with self.subTest(threads=threads, edge_first=edge_first, run=run):
                output = self.root / ('parallel-%d' % run)
                env = {'MEGAHIT_DISABLE_UNITIG_FIRST_TIPS': '1'} if edge_first else {}
                self.run_core(['assemble', '-s', graph, '-o', output, '-t', threads,
                               '--min_depth', 2, '--prune_level', 1, '--max_tip_len', 98,
                               '--cleaning_rounds', 0, '--bubble_level', 0], env)
                result = records(Path(str(output) + '.addi.fa'))
                digest = hashlib.sha256(json.dumps(result, separators=(',', ':')).encode()).hexdigest()
                self.assertEqual(len(result), fixture['expected_addi_count'])
                self.assertEqual(digest, fixture['expected_addi_sha256'])


if __name__ == '__main__':
    unittest.main()
