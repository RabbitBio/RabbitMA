#!/usr/bin/env python3

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest


sys.dont_write_bytecode = True
REPOSITORY = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DRIVER_PATH = os.path.join(REPOSITORY, 'src', 'megahit')
SPEC = importlib.util.spec_from_file_location('rabbitma_driver', DRIVER_PATH)
if SPEC is None or SPEC.loader is None:
    # SourceFileLoader is required because the driver intentionally has no
    # .py suffix in installed packages.
    from importlib.machinery import SourceFileLoader
    SPEC = importlib.util.spec_from_loader(
        'rabbitma_driver', SourceFileLoader('rabbitma_driver', DRIVER_PATH))
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


class CgroupMemoryTest(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix='rabbitma-cgroup-test-')

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def write(self, relative_path, value):
        path = os.path.join(self.temp_dir, relative_path)
        directory = os.path.dirname(path)
        if not os.path.isdir(directory):
            os.makedirs(directory)
        with open(path, 'w') as output:
            output.write(value)
        return path

    def test_cgroup_v2_uses_tightest_parent_limit(self):
        cgroup_file = self.write('proc/self/cgroup', '0::/jobs/sample\n')
        unified = os.path.join(self.temp_dir, 'sys/fs/cgroup')
        self.write('sys/fs/cgroup/memory.max', 'max\n')
        self.write('sys/fs/cgroup/jobs/memory.max', '2147483648\n')
        self.write('sys/fs/cgroup/jobs/sample/memory.max', '1073741824\n')

        self.assertEqual(
            DRIVER.detect_cgroup_memory_limit(
                cgroup_file, unified, os.path.join(unified, 'memory')),
            1073741824)

    def test_cgroup_v1_ignores_unlimited_sentinel(self):
        cgroup_file = self.write(
            'proc/self/cgroup',
            '7:cpu,cpuacct:/jobs/sample\n6:memory:/jobs/sample\n')
        legacy = os.path.join(self.temp_dir, 'sys/fs/cgroup/memory')
        self.write(
            'sys/fs/cgroup/memory/memory.limit_in_bytes',
            '9223372036854771712\n')
        self.write(
            'sys/fs/cgroup/memory/jobs/memory.limit_in_bytes',
            '4294967296\n')
        self.write(
            'sys/fs/cgroup/memory/jobs/sample/memory.limit_in_bytes',
            '8589934592\n')

        self.assertEqual(
            DRIVER.detect_cgroup_memory_limit(
                cgroup_file,
                os.path.join(self.temp_dir, 'sys/fs/cgroup'),
                legacy),
            4294967296)

    def test_unlimited_or_missing_controller_returns_none(self):
        cgroup_file = self.write('proc/self/cgroup', '7:cpu:/jobs/sample\n')
        self.assertIsNone(DRIVER.detect_cgroup_memory_limit(
            cgroup_file,
            os.path.join(self.temp_dir, 'sys/fs/cgroup'),
            os.path.join(self.temp_dir, 'sys/fs/cgroup/memory')))

    def test_cgroup_path_cannot_escape_mount_root(self):
        cgroup_file = self.write('proc/self/cgroup', '0::/../../outside\n')
        outside = self.write('outside/memory.max', '1\n')
        self.assertTrue(os.path.isfile(outside))
        self.assertIsNone(DRIVER.detect_cgroup_memory_limit(
            cgroup_file,
            os.path.join(self.temp_dir, 'sys/fs/cgroup'),
            os.path.join(self.temp_dir, 'sys/fs/cgroup/memory')))


class SharedNodeMemoryPolicyTest(unittest.TestCase):

    def test_exclusive_policy_preserves_budget(self):
        self.assertEqual(DRIVER.resolve_job_memory_budget(900, 1), 900)

    def test_shared_policy_divides_default_budget(self):
        self.assertEqual(DRIVER.resolve_job_memory_budget(901, 3), 300)

    def test_explicit_per_job_budget_overrides_equal_share(self):
        self.assertEqual(
            DRIVER.resolve_job_memory_budget(900, 3, 400), 400)

    def test_explicit_budget_cannot_widen_memory_limit(self):
        self.assertEqual(
            DRIVER.resolve_job_memory_budget(900, 3, 1200), 900)


class NumaPolicyTest(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix='rabbitma-numa-test-')
        self.node_root = os.path.join(self.temp_dir, 'sys/devices/system/node')
        self.status_path = os.path.join(self.temp_dir, 'proc/self/status')

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def write(self, path, value):
        directory = os.path.dirname(path)
        if not os.path.isdir(directory):
            os.makedirs(directory)
        with open(path, 'w') as output:
            output.write(value)

    def add_node(self, node, cpus, memory_kb=1024):
        directory = os.path.join(self.node_root, 'node%d' % node)
        self.write(os.path.join(directory, 'cpulist'), cpus + '\n')
        self.write(os.path.join(directory, 'meminfo'),
                   'Node %d MemTotal: %d kB\n' % (node, memory_kb))

    def test_linux_index_list_round_trip(self):
        parsed = DRIVER.parse_linux_index_list('0-3,8,10-11\n')
        self.assertEqual(parsed, {0, 1, 2, 3, 8, 10, 11})
        self.assertEqual(DRIVER.format_linux_index_list(parsed),
                         '0-3,8,10-11')

    def test_discovery_intersects_cpu_and_memory_allowances(self):
        self.add_node(0, '0-3')
        self.add_node(2, '4-7')
        self.write(self.status_path,
                   'Name:\tpython\nMems_allowed_list:\t0,2\n')
        self.assertEqual(
            DRIVER.discover_numa_cpu_sets(
                self.node_root, self.status_path, {2, 3, 4, 9}),
            {0: {2, 3}, 2: {4}})

        self.write(self.status_path,
                   'Name:\tpython\nMems_allowed_list:\t2\n')
        self.assertEqual(
            DRIVER.discover_numa_cpu_sets(
                self.node_root, self.status_path, {2, 3, 4, 9}),
            {2: {4}})

    def test_detect_numa_node_memory(self):
        self.add_node(3, '8-9', memory_kb=12345)
        self.assertEqual(
            DRIVER.detect_numa_node_memory(3, self.node_root),
            12345 * 1024)

    def test_auto_node_requires_one_scheduler_visible_domain(self):
        self.assertEqual(DRIVER.resolve_numa_node('auto', {3: {8, 9}}), 3)
        with self.assertRaises(ValueError):
            DRIVER.resolve_numa_node('auto', {0: {0}, 1: {1}})


if __name__ == '__main__':
    unittest.main()
