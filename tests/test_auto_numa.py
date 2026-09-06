#!/usr/bin/env python3
"""Resource ownership and inheritance tests; no assembly algorithms are mocked."""
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


sys.dont_write_bytecode = True
DRIVER_PATH = Path(__file__).resolve().parents[1] / 'src/megahit'
LOADER = importlib.machinery.SourceFileLoader('tested_numa_driver', str(DRIVER_PATH))
SPEC = importlib.util.spec_from_loader(LOADER.name, LOADER)
DRIVER = importlib.util.module_from_spec(SPEC)
LOADER.exec_module(DRIVER)


@unittest.skipUnless(sys.platform.startswith('linux'), 'Linux NUMA reservations')
class ReservationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = self.temporary.name
        self.reservations = []
        self.original_opt = DRIVER.opt
        DRIVER.opt = DRIVER.Options()
        self.original_reservation = DRIVER._numa_reservation
        DRIVER._numa_reservation = None
        self.nodes = {0: {(0, 0): {0, 8}, (0, 1): {1, 9}},
                      1: {(1, 0): {2, 10}, (1, 1): {3, 11}}}
        self.capacities = {0: 128, 1: 128}

    def tearDown(self):
        if DRIVER._numa_reservation is not None:
            DRIVER._numa_reservation.close()
        for reservation in self.reservations:
            reservation.close()
        DRIVER._numa_reservation = self.original_reservation
        DRIVER.opt = self.original_opt
        self.temporary.cleanup()

    def acquire(self, threads=2, **kwargs):
        result = DRIVER.NumaReservation.acquire(
            kwargs.pop('nodes', self.nodes), threads, self.capacities,
            directory=self.directory, **kwargs)
        if result is not None:
            self.reservations.append(result)
        return result

    def test_two_jobs_get_separate_nodes_and_keep_smt_siblings_together(self):
        a, b = self.acquire(), self.acquire()
        self.assertEqual((a.node, b.node), (0, 1))
        self.assertEqual(a.cpus, {0, 1, 8, 9})
        self.assertEqual(b.cpus, {2, 3, 10, 11})
        self.assertFalse(a.cpus & b.cpus)
        self.assertIsNone(self.acquire())

    def test_two_36_thread_jobs_on_72_physical_cores(self):
        nodes = {node: {(node, core): {node * 36 + core}
                        for core in range(36)} for node in (0, 1)}
        a = self.acquire(threads=36, nodes=nodes)
        b = self.acquire(threads=36, nodes=nodes)
        self.assertEqual({a.node, b.node}, {0, 1})
        self.assertEqual(len(a.cpus), 36)
        self.assertEqual(len(b.cpus), 36)
        self.assertFalse(a.cpus & b.cpus)
        self.assertEqual(a.cpus | b.cpus, set(range(72)))
        self.assertIsNone(self.acquire(threads=36, nodes=nodes))

    def test_72_logical_cpus_do_not_imply_72_physical_cores(self):
        nodes = {node: {(node, core): {node * 18 + core,
                                      node * 18 + core + 36}
                        for core in range(18)} for node in (0, 1)}
        self.assertIsNone(self.acquire(threads=36, nodes=nodes))
        self.assertIsNotNone(self.acquire(threads=18, nodes=nodes))
        self.assertIsNotNone(self.acquire(threads=18, nodes=nodes))

    def test_smaller_jobs_fill_unused_domains_before_sharing(self):
        jobs = [self.acquire(threads=1) for _ in range(4)]
        self.assertEqual([job.node for job in jobs], [0, 1, 0, 1])
        self.assertEqual(sum(len(job.cpus) for job in jobs),
                         len(set().union(*(job.cpus for job in jobs))))
        self.assertIsNone(self.acquire(threads=1))

    def test_release_allows_reuse(self):
        a, b = self.acquire(), self.acquire()
        a.close()
        self.assertEqual(self.acquire().node, 0)
        self.assertIsNone(self.acquire())

    def test_job_that_does_not_fit_is_not_partially_reserved(self):
        self.assertIsNone(self.acquire(threads=3))
        self.assertIsNotNone(self.acquire())
        self.assertIsNotNone(self.acquire())

    def test_memory_requirement_is_respected(self):
        self.capacities[0] = 32
        self.assertEqual(self.acquire(minimum_memory=64).node, 1)
        self.assertIsNone(self.acquire(minimum_memory=64))

    def test_lease_only_uses_granted_cpus(self):
        nodes = {0: {(0, 1): {9}}, 1: {(1, 0): {2}}}
        self.assertIsNone(self.acquire(threads=2, nodes=nodes))
        self.assertEqual(self.acquire(threads=1, nodes=nodes).cpus, {9})

    def test_nonprivate_directory_is_rejected(self):
        os.chmod(self.directory, 0o777)
        try:
            with self.assertRaises(OSError):
                self.acquire()
        finally:
            os.chmod(self.directory, 0o700)

    def test_lock_symlink_is_rejected(self):
        target = Path(self.directory) / 'unrelated'
        target.write_text('unchanged')
        (Path(self.directory) / 'allocation.lock').symlink_to(target)
        with self.assertRaises(OSError):
            self.acquire()
        self.assertEqual(target.read_text(), 'unchanged')

    def test_child_keeps_reservation_after_launcher_releases_it(self):
        reservation = self.acquire(nodes={0: self.nodes[0]})
        DRIVER._numa_reservation = reservation
        process = DRIVER.launch_process(
            [sys.executable, '-c', "import sys; print('ready', flush=True); sys.stdin.read()"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        try:
            self.assertEqual(process.stdout.readline().strip(), 'ready')
            reservation.close()
            DRIVER._numa_reservation = None
            self.assertIsNone(self.acquire(nodes={0: self.nodes[0]}))
            process.communicate('', timeout=5)
            self.assertIsNotNone(self.acquire(nodes={0: self.nodes[0]}))
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()

    def test_simultaneous_processes_do_not_choose_same_node(self):
        code = """
import json, runpy, sys
ns = runpy.run_path(sys.argv[1], run_name='reservation_test')
sys.stdin.readline()
nodes = {0: {(0, 0): {0}, (0, 1): {1}}, 1: {(1, 0): {2}, (1, 1): {3}}}
r = ns['NumaReservation'].acquire(nodes, 2, {0: 128, 1: 128}, directory=sys.argv[2])
print(json.dumps({'node': r.node, 'cpus': sorted(r.cpus)}), flush=True)
sys.stdin.read()
"""
        jobs = [subprocess.Popen([sys.executable, '-c', code, str(DRIVER_PATH), self.directory],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
                for _ in range(2)]
        try:
            for job in jobs:
                job.stdin.write('go\n')
                job.stdin.flush()
            results = [json.loads(job.stdout.readline()) for job in jobs]
            self.assertEqual({result['node'] for result in results}, {0, 1})
            self.assertFalse(set(results[0]['cpus']) & set(results[1]['cpus']))
        finally:
            for job in jobs:
                job.communicate('', timeout=5)

    def test_external_cpu_configuration_is_preserved(self):
        DRIVER.opt.num_cpu_threads = 2
        for environment in ({'OMP_PLACES': '{0},{2}'}, {'OMP_PROC_BIND': 'false'},
                            {'GOMP_CPU_AFFINITY': '0 2'}, {'MEGAHIT_NUMA_NODE': '1'}):
            with mock.patch.dict(os.environ, environment, clear=True), \
                    mock.patch.object(DRIVER, 'current_memory_policy') as policy:
                DRIVER.configure_automatic_numa_runtime()
                policy.assert_not_called()
                self.assertIsNone(DRIVER._numa_reservation)

    def test_inherited_memory_policy_is_preserved(self):
        DRIVER.opt.num_cpu_threads = 2
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(DRIVER, 'current_memory_policy', return_value=3), \
                mock.patch.object(DRIVER.NumaReservation, 'acquire') as acquire:
            DRIVER.configure_automatic_numa_runtime()
            acquire.assert_not_called()

    def test_default_all_thread_job_is_not_reduced_to_one_node(self):
        self.assertEqual(DRIVER.opt.num_cpu_threads, 0)
        with mock.patch.object(DRIVER, 'current_memory_policy') as policy:
            DRIVER.configure_automatic_numa_runtime()
            policy.assert_not_called()

    def test_explicit_sharing_or_disable_flag_skips_automatic_placement(self):
        DRIVER.opt.num_cpu_threads = 2
        for jobs_per_node, environment in ((2, {}),
                (1, {'MEGAHIT_DISABLE_AUTO_NUMA': '1'})):
            DRIVER.opt.jobs_per_node = jobs_per_node
            with mock.patch.dict(os.environ, environment, clear=True), \
                    mock.patch.object(DRIVER, 'current_memory_policy') as policy:
                DRIVER.configure_automatic_numa_runtime()
                policy.assert_not_called()
                self.assertIsNone(DRIVER._numa_reservation)

    def test_explicit_off_keeps_original_memory_domain(self):
        DRIVER.opt.numa_node = 'off'
        DRIVER.opt.memory = 0.5
        with mock.patch.object(DRIVER, 'configure_automatic_numa_runtime') as automatic, \
                mock.patch.object(DRIVER, 'detect_available_mem', return_value=256) as memory:
            DRIVER.configure_numa_runtime()
            automatic.assert_not_called()
            self.assertEqual(DRIVER.opt.host_mem, 128)
            memory.assert_called_once_with(None)

    def test_fractional_memory_uses_cpu_share_but_absolute_budget_is_unchanged(self):
        DRIVER._numa_reservation = self.acquire()
        DRIVER.opt.memory = 0.75
        self.assertEqual(DRIVER.opt.host_mem, 96)
        DRIVER.opt.memory = 100
        self.assertEqual(DRIVER.opt.host_mem, 100)

    def exercise_failure(self, affinity_failure=None):
        DRIVER.opt.num_cpu_threads = 2
        reservation = self.acquire()
        cpus = {0, 1, 2, 3}
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(DRIVER, 'current_memory_policy', return_value=0), \
                mock.patch.object(os, 'sched_getaffinity', return_value=cpus), \
                mock.patch.object(os, 'sched_setaffinity', side_effect=affinity_failure) as set_affinity, \
                mock.patch.object(DRIVER, 'discover_numa_cpu_sets', return_value={0: {0, 1}}), \
                mock.patch.object(DRIVER, 'physical_core_cpu_sets', return_value=self.nodes[0]), \
                mock.patch('builtins.open', mock.mock_open(read_data='0-1')), \
                mock.patch.object(DRIVER, 'detect_numa_node_memory', return_value=128), \
                mock.patch.object(DRIVER, 'detect_available_mem', return_value=256), \
                mock.patch.object(DRIVER.NumaReservation, 'acquire', return_value=reservation), \
                mock.patch.object(DRIVER, 'bind_process_memory_to_numa_node', side_effect=OSError('denied')):
            DRIVER.configure_automatic_numa_runtime()
        return reservation, cpus, set_affinity

    def test_automatic_policy_failure_restores_cpu_mask_and_releases_lease(self):
        reservation, cpus, set_affinity = self.exercise_failure()
        self.assertEqual(set_affinity.call_args_list,
                         [mock.call(0, reservation.cpus), mock.call(0, cpus)])
        self.assertEqual(reservation.fds, [])
        self.assertIsNone(DRIVER._numa_reservation)

    def test_denied_cpu_binding_falls_back_without_attempting_another_binding(self):
        reservation, cpus, set_affinity = self.exercise_failure(OSError('CPU binding denied'))
        set_affinity.assert_called_once_with(0, reservation.cpus)
        self.assertEqual(reservation.fds, [])
        self.assertIsNone(DRIVER._numa_reservation)

    def test_failed_rollback_keeps_cpu_reservation(self):
        reservation, cpus, set_affinity = self.exercise_failure([None, OSError('cpuset changed')])
        self.assertIs(DRIVER._numa_reservation, reservation)
        self.assertTrue(reservation.fds)
        self.assertEqual(self.acquire().node, 1)
        self.assertIsNone(self.acquire())


if __name__ == '__main__':
    unittest.main()
