#!/usr/bin/env python3
"""Experiment 2: share a fixed physical-CPU budget across 2 through 6 RabbitMA jobs."""

import json
from pathlib import Path
import shutil

from common import (allocate, compare_outputs, divide_threads, job, load_config,
                    parse_arguments, prepare, report, run_group, sha256, write_json)


def main():
    args = parse_arguments(__doc__)
    config = load_config(args.config)
    root = Path(config['output_root']) / 'concurrency'
    maximum = max(config['concurrency_levels'])
    inputs = [root / 'inputs' / ('sample_copy_%d.fq.gz' % index)
              if config['independent_inputs'] else Path(config['input'])
              for index in range(maximum)]
    configurations = [(count, placement) for count in config['concurrency_levels']
                      for placement in config['placements']]
    plans = []
    for repetition in range(config['repetitions']):
        order = configurations if repetition % 2 == 0 else list(reversed(configurations))
        for count, placement in order:
            threads = divide_threads(config['total_cpus'], count)
            masks = allocate(config['nodes'], threads, placement)
            label = '%dj_%s_rep%d' % (count, placement, repetition + 1)
            plans.append({'label': label, 'jobs': [
                job(config, 'rabbitma', mask, root / label / ('assembly%d' % index), inputs[index])
                for index, mask in enumerate(masks)]})
    if args.dry_run:
        print(json.dumps({'experiment': 'concurrency',
                          'independent_inputs': config['independent_inputs'],
                          'plans': plans}, indent=2))
        return
    directory = prepare(config, 'concurrency', plans)
    if config['independent_inputs']:
        (directory / 'inputs').mkdir()
        expected = sha256(config['input'])
        for path in inputs:
            print('Preparing independent input ' + str(path), flush=True)
            shutil.copyfile(config['input'], str(path))
            if sha256(path) != expected:
                raise RuntimeError('input copy checksum mismatch: ' + str(path))
    results, validation = [], {}
    reference = None
    for plan in plans:
        result = run_group(config, directory, plan['label'], plan['jobs'])
        results.append(result)
        report(directory, results)
        if reference is None:
            reference = Path(result['jobs'][0]['output']) / 'final.contigs.fa'
        validation[plan['label']] = {record['output']: compare_outputs(
            reference, Path(record['output']) / 'final.contigs.fa') for record in result['jobs']}
        write_json(directory / 'validation.json', validation)
    print('Experiment 2 complete: ' + str(directory / 'report.md'), flush=True)


if __name__ == '__main__':
    main()
