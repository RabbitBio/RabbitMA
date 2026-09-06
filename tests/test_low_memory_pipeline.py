#!/usr/bin/env python3
"""Small memory budgets must preserve complete assembly records."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CORE = Path(os.environ.get('MEGAHIT_TEST_CORE', ROOT / 'build/megahit_core')).resolve()


def records(path):
    result = []
    header, sequence = None, []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if header is not None:
                result.append((header, ''.join(sequence)))
            identifier, attributes = line[1:].split(' ', 1)
            header = (identifier.split('_', 1)[0], attributes)
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        result.append((header, ''.join(sequence)))
    return sorted(result)


class LowMemoryPipelineTests(unittest.TestCase):
    def test_streamed_and_resident_budgets_preserve_exact_final_records(self):
        with tempfile.TemporaryDirectory(prefix='megahit-low-memory-') as name:
            temporary = Path(name)
            binary = temporary / 'bin'
            binary.mkdir()
            shutil.copy2(ROOT / 'src/megahit', binary / 'megahit')
            for executable in ('megahit_core', 'megahit_core_popcnt',
                               'megahit_core_no_hw_accel', 'megahit_toolkit'):
                (binary / executable).symlink_to(CORE)
            data = ROOT / 'test_data'
            environment = {key: value for key, value in os.environ.items()
                           if not key.startswith(('MEGAHIT_', 'OMP_', 'GOMP_', 'KMP_'))}
            reference = None
            for budget in (2 * 2**30, int(0.9 * 2**30), 512 * 2**20, 128 * 2**20):
                with self.subTest(budget=budget):
                    output = temporary / str(budget)
                    command = [sys.executable, str(binary / 'megahit'),
                               '--12', str(data / 'r1.il.fa.gz') + ',' + str(data / 'r2.il.fa.bz2'),
                               '-1', str(data / 'r3_1.fa'), '-2', str(data / 'r3_2.fa'),
                               '-r', str(data / 'r4.fa') + ',' + str(data / 'loop.fa'),
                               '--k-list', '21,29,39,59', '-t', '2', '-m', str(budget),
                               '--numa-node', 'off', '-o', str(output)]
                    run = subprocess.run(command, env=environment, capture_output=True,
                                         text=True, timeout=60)
                    detail = run.stderr
                    if (output / 'log').exists():
                        detail += (output / 'log').read_text()
                    self.assertEqual(run.returncode, 0, detail)
                    actual = records(output / 'final.contigs.fa')
                    self.assertTrue(actual)
                    if reference is None:
                        reference = actual
                    else:
                        self.assertEqual(actual, reference)


if __name__ == '__main__':
    unittest.main()
