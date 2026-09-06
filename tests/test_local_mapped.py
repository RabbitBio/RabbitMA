#!/usr/bin/env python3
"""Check mapped endpoint gathering and cached minimizer gates against local."""
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

from test_stable_contigs import reverse_complement

CORE = Path(os.environ.get(
    'MEGAHIT_TEST_CORE', Path(__file__).resolve().parents[1] / 'build/megahit_core'))


def records(path):
    lines = path.read_text().splitlines()
    result = []
    for i in range(0, len(lines), 2):
        header = re.fullmatch(r'>lc_\d+_strand_\d+_id_\d+ flag=(\d+) multi=([\d.]+)', lines[i])
        assert header, lines[i]
        result.append((lines[i + 1], int(header[1]), float(header[2])))
    return sorted(result)


class LocalMappedTest(unittest.TestCase):
    def compare(self, library_kind, cache_rebuild=False):
        rng = random.Random(49141)
        dna = lambda n: ''.join(rng.choices('ACGT', k=n))
        with tempfile.TemporaryDirectory(prefix='megahit-local-mapped-') as name:
            root = Path(name)
            genomes = [dna(5000) for _ in range(3)]
            contigs = [s[1000:4000] for s in genomes]
            contig = root / 'contigs.fa'
            contig.write_text(''.join(
                f'>k39_{i} flag=0 multi=10 len={len(s)}\n{s}\n'
                for i, s in enumerate(contigs)))
            Path(str(contig) + '.info').write_text('3 9000\n')
            single = []
            paired = []
            lengths = (49, 50, 75, 99, 127, 128, 129, 150, 173, 193)
            for genome in genomes:
                for start in range(800, 4200, 2):
                    for copy in range(3):
                        length = lengths[(start + copy) % len(lengths)]
                        read = genome[start:start + length]
                        single.append(read if copy != 1 else reverse_complement(read))
                        first = genome[start:start + 150 + copy]
                        second = reverse_complement(genome[start + 220:start + 370 + copy])
                        paired.extend((first, second))
            # Force multiple packed chunks, an odd single-end library boundary,
            # and a longest source read which is never selected for assembly.
            single += [dna(512)] + [dna(67) for _ in range(18000)]
            self.assertEqual(len(single) % 2, 1)

            def fasta(path, reads):
                path.write_text(''.join(f'>r{i}\n{s}\n' for i, s in enumerate(reads)))

            fasta(root / 'single.fa', single)
            fasta(root / 'paired.fa', paired)
            config = root / 'config.lib'
            entries = []
            if library_kind in ('single', 'mixed'):
                entries += ['single reads', 'se ' + str(root / 'single.fa')]
            if library_kind in ('paired', 'mixed'):
                entries += ['paired reads', 'interleaved ' + str(root / 'paired.fa')]
            config.write_text('\n'.join(entries) + '\n')
            library = root / 'reads.lib'

            def run(command, flags=(), extra_env=None):
                env = {k: v for k, v in os.environ.items() if not k.startswith('MEGAHIT_')}
                env.update({flag: '1' for flag in flags})
                env.update(extra_env or {})
                result = subprocess.run([str(CORE)] + list(map(str, command)),
                                        env=env, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
                return result.stderr

            run(['buildlib', config, library, 4])
            mapped = 'MEGAHIT_EXPERIMENTAL_LOCAL_MAPPED_ASSEMBLY'
            gate = 'MEGAHIT_EXPERIMENTAL_LOCAL_MINIMIZER_GATE'
            validate = 'MEGAHIT_VALIDATE_LOCAL_MINIMIZER_GATE'
            outputs = {}
            for mode, flags in (('reference', ()), ('mapped', (mapped,)),
                                ('combined', (mapped, gate, validate))):
                output = root / (mode + '.fa')
                log = run(['local', '-c', contig, '-l', library, '-t', 4,
                           '-o', output, '--kmax', 59], flags)
                outputs[mode] = records(output)
                if mode == 'combined':
                    self.assertIn('bytes, built,', log)
                if mode != 'reference':
                    self.assertIn('Mapped endpoint reads:', log)
            self.assertTrue(outputs['reference'], 'fixture must produce local contigs')
            self.assertEqual(outputs['reference'], outputs['mapped'])
            self.assertEqual(outputs['reference'], outputs['combined'])
            if cache_rebuild:
                for action in ('reuse', 'source_changed', 'truncated'):
                    if action == 'source_changed':
                        # Rebuild the binary file and its matching chunk index;
                        # the old local positions must fail source identity.
                        run(['buildlib', config, library, 4])
                    elif action == 'truncated':
                        Path(str(library) + '.bin.local_seed_positions').write_bytes(b'incomplete')
                    output = root / (action + '.fa')
                    log = run(['local', '-c', contig, '-l', library, '-t', 4,
                               '-o', output, '--kmax', 59], (mapped, gate, validate))
                    self.assertIn('bytes, reused,' if action == 'reuse' else 'bytes, built,', log)
                    self.assertEqual(outputs['reference'], records(output))
                for action, extra in (
                    ('cache_unavailable', {'MEGAHIT_LOCAL_SEED_POSITIONS': str(root / 'absent' / 'positions')}),
                    ('small_budget', {'MEGAHIT_MEMORY_BUDGET_PER_JOB': '1024'}),
                ):
                    output = root / (action + '.fa')
                    log = run(['local', '-c', contig, '-l', library, '-t', 4,
                               '-o', output, '--kmax', 59], (mapped, gate, validate), extra)
                    self.assertEqual(outputs['reference'], records(output))
                    if action == 'cache_unavailable':
                        self.assertIn('Local seed position cache unavailable:', log)
                    else:
                        self.assertNotIn('Mapped endpoint reads:', log)

    def test_variable_length_single_end(self):
        self.compare('single')

    def test_paired_chunk_boundaries(self):
        self.compare('paired')

    def test_mixed_libraries_and_cache_identity(self):
        self.compare('mixed', cache_rebuild=True)


if __name__ == '__main__':
    unittest.main()
