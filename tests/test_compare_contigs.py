#!/usr/bin/env python3
"""Verify comparison semantics and CLI exit codes without an assembler build."""

import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'tools' / 'compare_contigs.py'
SPEC = importlib.util.spec_from_file_location('compare_contigs', str(SCRIPT))
COMPARATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARATOR)


def fasta(sequence, k=3, flag=0, multi='2.000', identifier=1):
    return '>k%d_%d flag=%d multi=%s len=%d\n%s\n' % (
        k, identifier, flag, multi, len(sequence), sequence)


def circle(core, **kwargs):
    return fasta(core + core[:3], flag=2, **kwargs)


class ContigComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.left = self.root / 'left.fa'
        self.right = self.root / 'right.fa'

    def tearDown(self):
        self.temporary.cleanup()

    def check(self, left, right, *flags, expected=0):
        self.left.write_text(left)
        self.right.write_text(right)
        process = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.left), str(self.right)] + list(flags),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(process.returncode, expected, process.stdout + process.stderr)
        if expected == 2:
            self.assertFalse(process.stdout)
            self.assertTrue(process.stderr)
            return None
        result = json.loads(process.stdout)
        self.assertEqual(result['equal'], expected == 0)
        return result

    def test_ids_order_and_line_wrapping_are_ignored(self):
        self.check(fasta('AACCGT') + fasta('TATAGG', identifier=2),
                   fasta('TATAGG', identifier=99) +
                   fasta('AACCGT', identifier=88).replace('AACCGT', 'AAC\nCGT'),
                   '--literal', '--metadata')

    def test_reverse_complement_modes(self):
        left, right = fasta('AACCGT'), fasta('ACGGTT')
        self.check(left, right, '--metadata')
        self.check(left, right, '--circular-loops', '--metadata')
        self.check(left, right, '--literal', '--metadata', expected=1)

    def test_circle_rotation_and_reverse_complement(self):
        for core in ('CCGTAA', 'CGGTTA'):
            with self.subTest(core=core):
                self.check(circle('AACCGT'), circle(core), '--metadata', expected=1)
                self.check(circle('AACCGT'), circle(core),
                           '--circular-loops', '--metadata')
                self.check(circle('AACCGT'), circle(core),
                           '--literal', '--metadata', expected=1)

    def test_changed_base_is_not_circular_equivalence(self):
        self.check(circle('AACCGT'), circle('AATCGT'),
                   '--circular-loops', '--metadata', expected=1)

    def test_unmarked_linear_records_are_not_rotated(self):
        self.check(fasta('AACCGT'), fasta('CCGTAA'), '--circular-loops', expected=1)

    def test_duplicate_multiplicity_is_preserved(self):
        result = self.check(fasta('AACCGT') * 2, fasta('AACCGT'), expected=1)
        self.assertEqual(result['missing_records'], 1)
        self.assertEqual(result['added_records'], 0)

    def test_metadata_and_flags_remain_visible(self):
        self.check(circle('AACCGT'), circle('CCGTAA', multi='2.001'),
                   '--circular-loops', '--metadata', expected=1)
        self.check(circle('AACCGT'), circle('CCGTAA', multi='2.001'), '--circular-loops')
        self.check(fasta('AACCGT'), fasta('AACCGT', k=4), '--metadata', expected=1)
        self.check(fasta('AACCGT'), fasta('AACCGT', flag=1), expected=1)

    def test_empty_files_and_one_empty_side(self):
        result = self.check('', '', '--metadata')
        self.assertEqual(result['common_records'], 0)
        self.check('', fasta('AACCGT'), expected=1)

    def test_generic_fasta_and_ambiguous_base(self):
        self.check('>left\naacgn\n', '>right\nNCGTT\n')
        self.check('>left\nAACGT\n', '>right\nAACGT\n', '--metadata', expected=2)

    def test_local_fasta_metadata(self):
        left = '>lc_1_strand_0_id_2 flag=0 multi=1.500\nAACCGT\n'
        right = '>lc_9_strand_1_id_8 flag=0 multi=1.500\nACGGTT\n'
        self.check(left, right, '--metadata')
        self.check(left, right, '--literal', '--metadata', expected=1)

    def test_invalid_fasta_is_an_error(self):
        malformed = (
            'AACCGT\n', '>id\n', '>id\nAACX\n',
            fasta('AACCGT').replace('len=6', 'len=7'),
            fasta('AACCGT', multi='nan'),
            fasta('AACCGT', multi='1e100'),
            fasta('AACCGT').replace('flag=0', 'flag=0 flag=1'),
        )
        for text in malformed:
            with self.subTest(text=text):
                self.check(text, fasta('AACCGT'), expected=2)

    def test_invalid_circular_overlap_is_an_error(self):
        self.check(fasta('AACCGT', flag=2), circle('AACCGT'),
                   '--circular-loops', expected=2)

    def test_mutually_exclusive_modes_are_an_error(self):
        self.check(fasta('AACCGT'), fasta('AACCGT'),
                   '--literal', '--circular-loops', expected=2)

    def test_minimum_rotation_against_exhaustive_reference(self):
        for size in range(1, 8):
            for bases in itertools.product('AC', repeat=size):
                sequence = ''.join(bases)
                expected = min(sequence[i:] + sequence[:i] for i in range(size))
                self.assertEqual(COMPARATOR.minimum_rotation(sequence), expected)

    def packed(self, sequence, k=3, flag=0, multi=2.0):
        result = struct.pack('<8sII', b'MGCTG01\0', 1, 0)
        result += struct.pack('<IIifq', len(sequence), k, flag, multi, 999)
        padded = sequence + 'A' * ((-len(sequence)) % 16)
        for offset in range(0, len(padded), 16):
            word = 0
            for base in padded[offset:offset + 16]:
                word = (word << 2) | 'ACGT'.index(base)
            result += struct.pack('<I', word)
        return result

    def test_packed_and_fasta_metadata_match(self):
        for length in (1, 15, 16, 17, 31, 32, 33):
            with self.subTest(length=length):
                sequence = ('ACGTTGC' * 5)[:length]
                self.left.write_bytes(self.packed(sequence, multi=0.1))
                self.right.write_text(fasta(sequence, multi='0.1'))
                self.assertTrue(COMPARATOR.compare(self.left, self.right,
                                                   'literal', True)['equal'])

    def test_packed_circle_and_fasta_rotation_match(self):
        self.left.write_bytes(self.packed('AACCGTAAC', flag=2))
        self.right.write_text(circle('CCGTAA'))
        self.assertTrue(COMPARATOR.compare(self.left, self.right, 'circular', True)['equal'])

    def test_truncated_packed_input_is_an_error(self):
        self.left = self.root / 'left.fa.mgb'
        self.left.write_bytes(self.packed('AACCGT')[:-1])
        self.right.write_text(fasta('AACCGT'))
        process = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.left), str(self.right)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(process.returncode, 2)
        self.assertIn(b'truncated packed sequence', process.stderr)

    def test_missing_input_is_an_error(self):
        process = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.left), str(self.right)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(process.returncode, 2)
        self.assertFalse(process.stdout)

    def test_script_runs_after_copying_to_another_directory(self):
        relocated = self.root / 'copied_compare.py'
        shutil.copyfile(str(SCRIPT), str(relocated))
        self.left.write_text(circle('AACCGT'))
        self.right.write_text(circle('CCGTAA'))
        process = subprocess.run(
            [sys.executable, str(relocated), 'left.fa', 'right.fa',
             '--circular-loops', '--metadata'], cwd=str(self.root),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertTrue(json.loads(process.stdout)['equal'])


if __name__ == '__main__':
    unittest.main()
