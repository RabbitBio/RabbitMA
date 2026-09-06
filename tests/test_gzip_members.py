#!/usr/bin/env python3
"""Exercise member speculation against the existing whole-stream parser.

Run with MEGAHIT_TEST_CORE=/path/to/megahit_core python3 -m unittest
discover -s tests -p test_gzip_members.py. Stored gzip blocks keep fixtures
above the normal discovery threshold without requiring large read sets.
"""

import gzip
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest


CORE = Path(os.environ.get(
    'MEGAHIT_TEST_CORE', Path(__file__).resolve().parents[1] / 'build/megahit_core'))
HEADER = struct.Struct('<QIIII6Q')
CHUNK = struct.Struct('<6Q')


def anchors(path):
    data = path.read_bytes()
    header = HEADER.unpack_from(data)
    assert header[0] == 0x5241504F53563131 and header[2] == HEADER.size
    assert header[9] == HEADER.size + header[8]
    assert len(data) == header[9] + CHUNK.size * header[10]
    return header, data[HEADER.size:header[9]], list(
        CHUNK.iter_unpack(data[header[9]:]))


class GzipMembersTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Distinct N/ambiguous/short/wrapped reads exercise packing as well as
        # read order. Long headers make boundary fixtures cheap to assemble.
        cls.record = b'@r ' + b'x' * 500 + b'\n' + b'ACGT' * 37 + b'NN\n+\n' + b'I' * 150 + b'\n'
        cls.padding = cls.record * 16000
        cls.special = b'@a\r\nNNacgt\r\nNACGTNN\r\n+\r\n!!!!!!\r\n!!!!!!!\r\n@b\nNNN\n+\n!!!\n'

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-members-')
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_build(self, config, label, threads=8, member=True, sidecar=True):
        prefix = self.root / label
        env = os.environ.copy()
        env['MEGAHIT_DISABLE_RAPIDGZIP'] = '1'
        env.pop('MEGAHIT_DISABLE_GZIP_MEMBER_PARALLEL', None)
        if not member:
            env['MEGAHIT_DISABLE_GZIP_MEMBER_PARALLEL'] = '1'
        command = [str(CORE), 'buildlib', str(config), str(prefix), str(threads)]
        if sidecar:
            command += ['19', '40']
        run = subprocess.run(command, env=env, capture_output=True, text=True, timeout=90)
        self.assertEqual(run.returncode, 0, run.stderr)
        return prefix, run.stderr

    def compare(self, parts, kind='se', threads=8, accepted=True, sidecar=True):
        source = self.root / 'input.fq.gz'
        source.write_bytes(b''.join(gzip.compress(p, compresslevel=0, mtime=0) for p in parts))
        config = self.root / 'reads.lib'
        config.write_text('one biological library\n' + kind + ' ' + str(source) + '\n')
        reference, _ = self.run_build(config, 'reference', member=False, sidecar=sidecar)
        candidate, log = self.run_build(config, 'candidate', threads=threads, sidecar=sidecar)
        self.assertIn('Discovered ', log, log)
        for suffix in ('.bin', '.lib_info'):
            self.assertEqual(Path(str(reference) + suffix).read_bytes(),
                             Path(str(candidate) + suffix).read_bytes(), suffix + '\n' + log)
        self.assertEqual('Gzip-member parallel path rejected' not in log, accepted, log)
        if accepted and sidecar:
            ref_header, ref_payload, _ = anchors(Path(str(reference) + '.bin.ridx.positions'))
            new_header, new_payload, chunks = anchors(Path(str(candidate) + '.bin.ridx.positions'))
            self.assertEqual(ref_header[:10], new_header[:10])
            self.assertEqual(ref_payload, new_payload)
            rows = [tuple(map(int, line.split())) for line in
                    Path(str(candidate) + '.bin.ridx.chunks').read_text().splitlines()[1:]]
            self.assertEqual([c[:4] for c in chunks], [r[:4] for r in rows])
            self.assertEqual(chunks[0][4], 0)
            self.assertEqual(chunks[-1][5], len(new_payload))
            self.assertTrue(all(a[5] == b[4] for a, b in zip(chunks, chunks[1:])))
        self.assertFalse(list(self.root.glob('candidate*.member.*')), log)

    def test_complete_wrapped_records(self):
        self.compare([self.padding + self.special, self.special + self.padding])

    def test_two_threads_still_generate_positions(self):
        self.compare([self.padding, self.padding + self.special], threads=2)

    def test_more_members_than_workers(self):
        self.compare([self.padding] * 4, threads=3)

    def test_split_sequence_falls_back_as_one_stream(self):
        self.compare([self.padding + b'@split\nACGT',
                      b'TGCA\n+\n!!!!!!!!\n' + self.padding], accepted=False, sidecar=False)

    def test_split_quality_falls_back_as_one_stream(self):
        self.compare([self.padding + b'@split\nACGT\n+\n!!',
                      b'!!\n' + self.padding], accepted=False, sidecar=False)

    def test_missing_member_newline_falls_back(self):
        self.compare([self.padding + b'@split\nACGT\n+\n!!!!',
                      b'\n' + self.padding], accepted=False)

    def test_pair_can_cross_a_member(self):
        self.compare([self.padding + self.record, self.record + self.padding],
                     kind='interleaved')

    def test_fasta_split_falls_back_as_one_stream(self):
        fasta = b'>padding ' + b'x' * len(self.padding) + b'\nACGT\n'
        self.compare([fasta + b'>split\nACGT', b'TGCA\n' + fasta], accepted=False, sidecar=False)

    def test_empty_member(self):
        self.compare([self.padding, b'', self.padding])

    def test_false_header_in_stored_payload_is_rejected(self):
        signature = gzip.compress(b'', compresslevel=0, mtime=0)[:10]
        record = b'@false\n' + b'A' * len(signature) + b'\n+\n' + signature + b'\n'
        self.compare([self.padding + record, self.padding], accepted=False)


if __name__ == '__main__':
    unittest.main()
