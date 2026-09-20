#!/usr/bin/env python3
"""Regression tests for structured phase-profile parsing."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'benchmarks/profile_phase_breakdown.py'
SPEC = importlib.util.spec_from_file_location('profile_phase_breakdown',
                                               str(SCRIPT))
PROFILE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROFILE)


class ProfilePhaseBreakdownTest(unittest.TestCase):
    def test_local_mapping_profile(self):
        line = (
            'Local-mapping profile: k=39 next_k=59 mapped_stream=1 '
            'candidate_source=index mapper_build_time=1.25 read_lib_time=0.5 '
            'insert_size_time=2.0 endpoint_filter_time=0.25 '
            'mapping_time=3.5 collation_time=0.125 compaction_time=0.75 '
            'endpoint_assembly_time=4.0 library_reads=2000000 '
            'insert_reads_visited=524288 insert_try_map_attempts=524288 '
            'mapping_reads_visited=2000000 candidate_reads=125000 '
            'candidate_pairs=80000 endpoint_gate_attempts=150000 '
            'endpoint_selected_pairs=50000 endpoint_selected_reads=2500 '
            'try_map_attempts=102500 aligned_reads=90000 '
            'added_mappings=120000 retained_mappings=110000\n')
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8') as stream:
            stream.write(line)
            stream.flush()
            per_k, rounds, read_index, local = PROFILE.parse_log(stream.name)

        self.assertEqual(per_k, {})
        self.assertEqual(rounds, {})
        self.assertEqual(read_index, [])
        self.assertEqual(len(local), 1)
        self.assertEqual(local[0]['k'], 39)
        self.assertEqual(local[0]['candidate_source'], 'index')
        self.assertEqual(local[0]['mapping_reads_visited'], 2000000)
        self.assertEqual(local[0]['retained_mappings'], 110000)
        self.assertEqual(local[0]['mapping_time'], 3.5)
        report = PROFILE.markdown(per_k, rounds, read_index, local)
        self.assertIn('| 39 → 59 | 1.2500 |', report)
        self.assertIn('| 39 → 59 | index | 2,000,000 |', report)


if __name__ == '__main__':
    unittest.main()
