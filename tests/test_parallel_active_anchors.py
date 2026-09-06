#!/usr/bin/env python3
"""Exact edge/context regression for parallel active-flank preparation."""
from collections import Counter
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

CORE = Path(os.environ.get(
    'MEGAHIT_TEST_CORE', Path(__file__).resolve().parents[1] / 'build/megahit_core'))


def edges(prefix):
    meta = dict(line.split() for line in Path(str(prefix) + '.edges.info').read_text().splitlines())
    width = int(meta['words_per_edge']) * 4
    output = Counter()
    for shard in range(int(meta['num_files'])):
        data = Path(str(prefix) + f'.edges.{shard}').read_bytes()
        assert len(data) % width == 0
        output.update(data[i:i + width] for i in range(0, len(data), width))
    assert sum(output.values()) == int(meta['num_edges'])
    return output


class ParallelActiveAnchorsTest(unittest.TestCase):
    def test_full_build_and_replay_across_key_widths(self):
        rng = random.Random(1959)
        dna = lambda n: ''.join(rng.choices('ACGT', k=n))
        with tempfile.TemporaryDirectory(prefix='megahit-active-anchors-') as name:
            root = Path(name)
            genomes = [dna(1800) for _ in range(30)]
            reads = []
            for genome in genomes:
                reads += [genome[start:start + 150] for start in range(450, 1150, 3)]
            fasta = root / 'reads.fa'
            fasta.write_text(''.join(f'>r{i}\n{s}\n' for i, s in enumerate(reads)))
            config = root / 'config.lib'
            config.write_text('reads\nse ' + str(fasta) + '\n')
            lib = root / 'reads.lib'

            def run(command, parallel=False, direct=False):
                env = {k: v for k, v in os.environ.items() if not k.startswith('MEGAHIT_')}
                if parallel:
                    env['MEGAHIT_EXPERIMENTAL_PARALLEL_ACTIVE_ANCHORS'] = '1'
                if direct:
                    env['MEGAHIT_EXPERIMENTAL_DIRECT_ITERATE_WINDOWS'] = '1'
                result = subprocess.run([str(CORE)] + list(map(str, command)),
                                        env=env, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
                return result.stderr

            run(['buildlib', config, lib, 4, 19, 40])
            bubble = root / 'bubble.fa'
            bubble.write_text('')
            Path(str(bubble) + '.info').write_text('0 0\n')
            for k, step in ((39, 10), (63, 10), (119, 10), (139, 2)):
                with self.subTest(k=k):
                    contig = root / f'k{k}.fa'
                    sequences = [g[600:1000] for g in genomes]
                    sequences += [g[600 + step:1000 + step] for g in genomes]
                    # Repeated endpoint keys exercise exact query count/mask reduction.
                    sequences += sequences[:5]
                    contig.write_text(''.join(
                        f'>k{k}_{i} flag=0 multi=5 len={len(s)}\n{s}\n'
                        for i, s in enumerate(sequences)))
                    Path(str(contig) + '.info').write_text(f'{len(sequences)} {sum(map(len,sequences))}\n')
                    outputs = []
                    for mode in ('serial', 'parallel', 'direct'):
                        prefix = root / f'{k}-{mode}'
                        common = ['-r', str(lib) + '.bin', '-a', 19, '-w', 40,
                                  '-c', contig, '-b', bubble, '-k', k, '-s', step,
                                  '-t', 4, '-m', 1000000000]
                        run(['read-index', '-o', prefix, '-e', str(prefix) + '-build'] + common,
                            mode != 'serial', mode == 'direct')
                        run(['read-index', '--index_prefix', prefix, '-e', str(prefix) + '-replay'] + common,
                            mode != 'serial', mode == 'direct')
                        built = edges(str(prefix) + '-build')
                        replayed = edges(str(prefix) + '-replay')
                        self.assertTrue(built)
                        self.assertEqual(built, replayed)
                        outputs.append(replayed)
                    for output in outputs[1:]:
                        self.assertEqual(outputs[0], output)


if __name__ == '__main__':
    unittest.main()
