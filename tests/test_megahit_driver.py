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


if __name__ == '__main__':
    unittest.main()
