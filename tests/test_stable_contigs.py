#!/usr/bin/env python3
"""Compare complete graph/assembly results with certified component reuse.

MEGAHIT_TEST_CORE selects a CPU variant. Record order and numeric IDs are
parallel scheduling artifacts; sequence strand, flags and coverage must agree.
"""

import os
from pathlib import Path
import random
import re
import struct
import subprocess
import tempfile
import unittest

CORE = Path(os.environ.get(
    'MEGAHIT_TEST_CORE', Path(__file__).resolve().parents[1] / 'build/megahit_core'))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def records(path, canonical=False):
    packed = Path(str(path) + '.mgb')
    output = []
    if packed.exists():
        with packed.open('rb') as source:
            assert source.read(16) == struct.pack('<8sII', b'MGCTG01\0', 1, 0)
            while header := source.read(24):
                length, k, flag, multi, _ = struct.unpack('<IIifq', header)
                words = struct.unpack('<' + 'I' * ((length + 15) // 16),
                                      source.read(4 * ((length + 15) // 16)))
                sequence = ''.join('ACGT'[(words[i // 16] >> (30 - 2 * (i % 16))) & 3]
                                   for i in range(length))
                output.append((sequence, k, flag, multi))
    elif path.exists():
        lines = path.read_text().splitlines()
        for i in range(0, len(lines), 2):
            header = re.fullmatch(r'>k(\d+)_\d+ flag=(\d+) multi=([\d.]+) len=(\d+)', lines[i])
            assert header, lines[i]
            k, flag, multi, length = header.groups()
            assert len(lines[i+1]) == int(length)
            output.append((lines[i+1], int(k), int(flag), float(multi)))
    if canonical:
        output = [(min(s, reverse_complement(s)), k, f, m) for s, k, f, m in output]
    return sorted(output)


class StableContigsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-carry-')
        self.root = Path(self.temp.name)
        self.rng = random.Random(82473)
        self.main = self.dna(1200)

    def tearDown(self):
        self.temp.cleanup()

    def dna(self, size):
        return ''.join(self.rng.choices('ACGT', k=size))

    def fasta(self, name, items):
        path = self.root / name
        path.write_text(''.join(f'>k29_{i} flag={flag} multi={multi} len={len(seq)}\n{seq}\n'
                                for i, (seq, flag, multi) in enumerate(items)))
        Path(str(path)+'.info').write_text(f'{len(items)} {sum(len(x[0]) for x in items)}\n')
        return path

    def run_core(self, command, exact_overlap=False):
        env = os.environ.copy()
        if exact_overlap:
            env['MEGAHIT_EXACT_STABLE_OVERLAPS'] = '1'
        run = subprocess.run([str(CORE)] + list(map(str, command)),
                             capture_output=True, text=True, timeout=60, env=env)
        self.assertEqual(run.returncode, 0, run.stderr)
        return run.stderr

    def compare(self, primary=None, additions=(), k=39, asm=(), expect_carry=True,
                background=True, rerun_without_carry=False):
        primary = primary or [(self.main, 1, 5.49)]
        # Keep a separate noncandidate component so the carry plan need not
        # take its all-isolated/empty-graph fallback.
        if background:
            primary = primary + [(self.dna(180), 0, 3)]
        contig = self.fasta('input.fa', primary)
        bubble = self.fasta('bubble.fa', [])
        additional = self.fasta('local.fa', additions)
        for label in ('reference', 'candidate'):
            graph = self.root / (label + '.graph')
            command = ['seq2sdbg', '-k', k, '--kmer_from', 29, '--contig', contig,
                       '--bubble', bubble, '--local_contig', additional,
                       '--host_mem', 1000000000, '--mem_flag', 0, '-t', 4, '-o', graph]
            if label == 'candidate':
                command.append('--carry_stable_contigs')
            log = self.run_core(command)
            if label == 'candidate':
                self.assertEqual('Carried isolated paths:' in log, expect_carry, log)
            self.run_core(['assemble', '-s', graph, '-o', self.root / label, '-t', 4,
                           '--min_depth', 2, '--prune_level', 2] + list(asm))
        for suffix in ('.contigs.fa', '.addi.fa', '.bubble_seq.fa', '.final.contigs.fa'):
            ref = records(self.root / ('reference' + suffix))
            new = records(self.root / ('candidate' + suffix))
            self.assertEqual(ref, new, suffix)
        if rerun_without_carry:
            self.run_core(command[:-1])
            self.assertFalse(Path(str(graph) + '.isolated.fa.mgb').exists())

    def test_isolated_and_contained_inputs(self):
        self.compare(additions=[(self.main[100:400], 0, 1),
                                (reverse_complement(self.main[230:620]), 0, 5)])

    def test_higher_depth_prevents_reuse(self):
        self.compare(additions=[(self.main[100:400], 0, 6)], expect_carry=False)

    def test_attachment_prevents_reuse(self):
        self.compare(additions=[(self.main[-39:] + self.dna(100), 0, 1)], expect_carry=False)

    def test_repeated_nodes_prevent_reuse(self):
        repeat = self.dna(300)
        self.compare(primary=[(repeat * 4, 1, 5)], expect_carry=False)

    def test_two_candidates_share_nodes(self):
        other = self.dna(300) + self.main[400:500] + self.dna(400)
        self.compare(primary=[(self.main, 1, 5), (other, 1, 3)], expect_carry=False)

    def test_custom_tip_threshold(self):
        self.compare(asm=['--max_tip_len', 2000])

    def test_global_depth_pruning(self):
        self.compare(asm=['--prune_level', 3, '--min_depth', 6])

    def test_inferred_depth(self):
        self.compare(asm=['--min_depth', 0])

    def test_final_standalone_output(self):
        self.compare(asm=['--is_final_round', '--output_standalone', '--min_standalone', 500])

    def test_larger_k_and_rounded_coverage(self):
        self.compare(primary=[(self.main, 1, 5.5)], k=99)

    def test_all_isolated_falls_back(self):
        self.compare(background=False, expect_carry=False)

    def test_small_k_falls_back(self):
        self.compare(k=29, expect_carry=False)

    def test_palindromic_path_is_not_carried(self):
        half = self.dna(300)
        self.compare(primary=[(half + reverse_complement(half), 1, 5)], expect_carry=False)

    def test_graph_rebuild_discards_stale_carry(self):
        self.compare(rerun_without_carry=True)

    def test_short_seed_match_is_not_a_graph_connection(self):
        seed = self.dna(31)
        first = self.main[:89] + 'A' + seed + 'A' + self.main[122:]
        other = self.dna(179) + 'C' + seed + 'C' + self.dna(700)
        # Both seeds lie on the 9-base sampling grid for k=39, but neither
        # adjacent base agrees. The exact mode must retain both candidates.
        original = self.run_core
        self.run_core = lambda command: original(command, exact_overlap=True)
        self.compare(primary=[(first, 1, 5), (other, 1, 3)])


if __name__ == '__main__':
    unittest.main()
