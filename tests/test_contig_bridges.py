#!/usr/bin/env python3
"""Exercise path compression against the complete de Bruijn graph backend."""
from pathlib import Path
import unittest
import test_stable_contigs as fixtures

records = fixtures.records
reverse_complement = fixtures.reverse_complement


class ContigBridgesTest(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.StableContigsTest()
        self.fixture.setUp()
        self.root = self.fixture.root
        self.main = self.fixture.main
        self.dna = self.fixture.dna

    def tearDown(self):
        self.fixture.tearDown()

    def compare(self, primary=None, additions=(), k=39, asm=(), expect=True, threads=4,
                mem_flag=0, rebuild=False):
        primary = primary or [(self.main, 0, 5)]
        contig = self.fixture.fasta('bridge-input.fa', primary + [(self.dna(180), 0, 3)])
        bubble = self.fixture.fasta('bridge-bubble.fa', [])
        additional = self.fixture.fasta('bridge-local.fa', additions)
        for label in ('reference', 'candidate'):
            graph = self.root / (label + '.graph')
            command = ['seq2sdbg', '-k', k, '--kmer_from', 29, '--contig', contig,
                       '--bubble', bubble, '--local_contig', additional,
                       '--host_mem', 1000000000, '--mem_flag', mem_flag, '-t', threads, '-o', graph]
            if label == 'candidate':
                command.append('--bridge_contigs')
            log = self.fixture.run_core(command)
            if label == 'candidate':
                self.assertEqual('path interiors,' in log, expect, log)
            self.fixture.run_core(['assemble', '-s', graph, '-o', self.root / label, '-t', threads,
                                   '--min_depth', 2, '--prune_level', 2] + list(asm))
        for suffix in ('.contigs.fa', '.addi.fa', '.bubble_seq.fa', '.final.contigs.fa'):
            self.assertEqual(records(self.root / ('reference' + suffix)),
                             records(self.root / ('candidate' + suffix)), suffix)
        if rebuild:
            self.fixture.run_core(command[:-1])
            self.assertFalse(Path(str(graph) + '.bridges.fa.mgb').exists())
            self.assertFalse(Path(str(graph) + '.bridges.meta').exists())

    def test_linear_path(self):
        self.compare()

    def test_both_endpoint_extensions(self):
        self.compare(additions=[(self.dna(80) + self.main[:250], 0, 1),
                                (self.main[-250:] + self.dna(90), 0, 1)])

    def test_additional_input_crosses_two_paths(self):
        second = self.dna(900)
        self.compare(primary=[(self.main, 0, 5), (second, 0, 8)],
                     additions=[(self.main + second, 0, 1)])

    def test_reverse_input_crosses_two_paths(self):
        second = self.dna(900)
        self.compare(primary=[(self.main, 0, 5), (second, 0, 8)],
                     additions=[(reverse_complement(self.main + second), 0, 1)])

    def test_interior_branch_prevents_compression(self):
        self.compare(additions=[(self.main[300:600] + self.dna(100), 0, 1)], expect=False)

    def test_branch_at_first_protected_node(self):
        branch = 'ACGT'[('ACGT'.index(self.main[78]) + 1) % 4]
        self.compare(additions=[(self.dna(99) + branch + self.main[79:600], 0, 1)], expect=False)

    def test_branch_at_last_protected_node(self):
        branch = 'ACGT'[('ACGT'.index(self.main[1121]) + 1) % 4]
        self.compare(additions=[(self.main[600:1121] + branch + self.dna(99), 0, 1)], expect=False)

    def test_interior_higher_coverage_prevents_compression(self):
        self.compare(additions=[(self.main[300:600], 0, 9)], expect=False)

    def test_endpoint_higher_coverage(self):
        self.compare(additions=[(self.main[:50], 0, 20)])

    def test_loop_minimum_inside_omitted_span(self):
        sequence = self.main[:300] + 'A' * 17 + self.main[317:]
        self.compare(primary=[(sequence, 0, 5)], additions=[(sequence[-39:] + sequence[:39], 0, 1)])

    def test_loop_minimum_in_retained_closure(self):
        sequence = self.main[:10] + 'A' * 17 + self.main[27:]
        self.compare(primary=[(sequence, 0, 5)], additions=[(sequence[-39:] + sequence[:39], 0, 1)])

    def test_two_path_loop_with_different_coverage(self):
        second = self.dna(900)
        self.compare(primary=[(self.main, 0, 5), (second, 0, 9)],
                     additions=[(self.main[-39:] + second[:39], 0, 1),
                                (second[-39:] + self.main[:39], 0, 1)])

    def test_palindromic_connection_at_end(self):
        self.compare(additions=[(self.main + reverse_complement(self.main), 0, 1)])

    def test_palindromic_connection_at_start(self):
        self.compare(additions=[(reverse_complement(self.main) + self.main, 0, 1)])

    def test_global_depth_histogram(self):
        self.compare(asm=['--min_depth', 0])

    def test_long_tip_removal(self):
        self.compare(asm=['--max_tip_len', 3000])

    def test_final_standalone_output(self):
        self.compare(asm=['--is_final_round', '--output_standalone', '--min_standalone', 500])

    def test_large_k_one_thread(self):
        self.compare(k=99, threads=1)

    def test_smallest_certified_k(self):
        self.compare(k=31)

    def test_unsupported_small_k_falls_back(self):
        self.compare(k=29, expect=False)

    def test_final_large_k(self):
        self.compare(k=141, asm=['--is_final_round'])

    def test_balanced_memory_mode(self):
        self.compare(mem_flag=1)

    def test_maximum_memory_mode(self):
        self.compare(mem_flag=2)

    def test_rebuild_removes_stale_bridges(self):
        self.compare(rebuild=True)


if __name__ == '__main__':
    unittest.main()
