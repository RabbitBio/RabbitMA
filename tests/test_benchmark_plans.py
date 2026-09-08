#!/usr/bin/env python3
"""Check benchmark CPU allocation without starting any assembly workload."""

import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from contextlib import redirect_stdout
from types import SimpleNamespace


PATH = Path(__file__).resolve().parents[1] / 'benchmarks' / 'common.py'
SPEC = importlib.util.spec_from_file_location('benchmark_common', str(PATH))
COMMON = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMMON)


class BenchmarkPlanTests(unittest.TestCase):
    def setUp(self):
        self.nodes = {0: list(range(32)), 1: list(range(32, 64))}

    def test_requested_thread_splits(self):
        expected = {2: [32, 32], 3: [22, 21, 21], 4: [16] * 4,
                    5: [13, 13, 13, 13, 12], 6: [11, 11, 11, 11, 10, 10]}
        for count, split in expected.items():
            self.assertEqual(COMMON.divide_threads(64, count), split)

    def test_every_layout_uses_64_disjoint_cores(self):
        for count in range(2, 7):
            threads = COMMON.divide_threads(64, count)
            for placement in ('numa', 'spread'):
                with self.subTest(count=count, placement=placement):
                    masks = COMMON.allocate(self.nodes, threads, placement)
                    self.assertEqual(list(map(len, masks)), threads)
                    self.assertEqual(sorted(cpu for mask in masks for cpu in mask), list(range(64)))
                    crossed = sum(any(cpu < 32 for cpu in mask) and
                                  any(cpu >= 32 for cpu in mask) for mask in masks)
                    self.assertEqual(crossed, count if placement == 'spread' else int(count in (3, 5)))

    def test_noncontiguous_cpu_ids_are_preserved(self):
        nodes = {2: [8, 10, 12, 14], 7: [32, 34, 36, 38]}
        for placement in ('numa', 'spread'):
            masks = COMMON.allocate(nodes, [3, 3, 2], placement)
            self.assertEqual(sorted(cpu for mask in masks for cpu in mask),
                             sorted(nodes[2] + nodes[7]))

    def test_invalid_requests_are_rejected(self):
        for count in (0, 65):
            with self.assertRaises(ValueError):
                COMMON.divide_threads(64, count)
        with self.assertRaises(ValueError):
            COMMON.allocate(self.nodes, [33, 32], 'numa')

    def test_memory_is_proportional_and_arguments_are_not_shell_strings(self):
        config = {'memory_bytes': 192 * 2**30, 'total_cpus': 64,
                  'input': '/input with spaces/sample.fq.gz', 'nodes': self.nodes,
                  'programs': {'rabbitma': {'launcher': '/package with spaces/bin/megahit'}},
                  'k_list': [39, 59, 79, 99, 119, 141]}
        record = COMMON.job(config, 'rabbitma', list(range(24)), '/output with spaces/run')
        self.assertEqual(record['memory_bytes'], 72 * 2**30)
        self.assertEqual(record['numa_nodes'], [0])
        self.assertIn(config['input'], record['command'])
        self.assertIn('/output with spaces/run', record['command'])

    def test_dry_runs_never_prepare_results_or_launch_assemblers(self):
        config = {'input': '/input/sample.fq.gz', 'output_root': '/never-created',
                  'memory_bytes': 192 * 2**30, 'total_cpus': 64, 'single_threads': 24,
                  'nodes': self.nodes, 'repetitions': 2,
                  'concurrency_levels': [2, 3, 4, 5, 6], 'placements': ['numa', 'spread'],
                  'independent_inputs': True, 'k_list': [39, 59, 79, 99, 119, 141],
                  'programs': {'rabbitma': {'launcher': '/never-executed/new'},
                               'original': {'launcher': '/never-executed/old'}}}
        for name, count in [('compare_release_24', 4), ('compare_concurrency', 20)]:
            spec = importlib.util.spec_from_file_location(name, str(PATH.parent / (name + '.py')))
            module = importlib.util.module_from_spec(spec)
            with mock.patch.dict('sys.modules', {'common': COMMON}):
                spec.loader.exec_module(module)
            output = io.StringIO()
            with mock.patch.object(module, 'parse_arguments', return_value=SimpleNamespace(
                    config='/config.json', dry_run=True)), \
                    mock.patch.object(module, 'load_config', return_value=config), \
                    mock.patch.object(module, 'prepare') as prepare, \
                    mock.patch.object(module, 'run_group') as run, redirect_stdout(output):
                module.main()
            self.assertEqual(len(json.loads(output.getvalue())['plans']), count)
            prepare.assert_not_called()
            run.assert_not_called()

    def test_manifest_records_provenance_without_running_binaries(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            launcher, input_path = root / 'megahit', root / 'input.fq'
            launcher.write_text('not an executable')
            input_path.write_text('small test input')
            config = {'output_root': str(root / 'results'), 'input': str(input_path),
                      'memory_bytes': 1024,
                      'programs': {'rabbitma': {'launcher': str(launcher)}}}
            with mock.patch.object(COMMON.subprocess, 'Popen') as launch:
                directory = COMMON.prepare(config, 'fixture', [])
            manifest = json.loads((directory / 'manifest.json').read_text())
            self.assertEqual(manifest['input_sha256'], COMMON.sha256(input_path))
            self.assertIn('machine', manifest['host'])
            launch.assert_not_called()
            with self.assertRaises(ValueError):
                COMMON.prepare(config, 'fixture', [])

    def test_wrapper_stops_active_phase_on_interruption(self):
        spec = importlib.util.spec_from_file_location(
            'benchmark_runner', str(PATH.parent / 'run_experiments.py'))
        module = importlib.util.module_from_spec(spec)
        with mock.patch.dict('sys.modules', {'common': COMMON}):
            spec.loader.exec_module(module)
        process = mock.Mock()
        process.wait.side_effect = [KeyboardInterrupt, 0]
        process.poll.return_value = None
        with mock.patch.object(module.subprocess, 'Popen', return_value=process) as launch:
            with self.assertRaises(KeyboardInterrupt):
                module.run_phase(['never-executed'])
        launch.assert_called_once_with(['never-executed'], start_new_session=True)
        process.terminate.assert_called_once_with()
        self.assertEqual(process.wait.call_count, 2)


if __name__ == '__main__':
    unittest.main()
