#!/usr/bin/env python3
"""Experiment 1: original MEGAHIT and the pinned RabbitMA release at 24 threads."""

import json
from pathlib import Path

from common import (compare_outputs, job, load_config, parse_arguments,
                    prepare, report, run_group, write_json)


def main():
    args = parse_arguments(__doc__)
    config = load_config(args.config)
    threads = config['single_threads']
    nodes = [cpus for _, cpus in sorted(config['nodes'].items()) if len(cpus) >= threads]
    if not nodes:
        raise ValueError('the single-job comparison needs %d cores in one NUMA node' % threads)
    cpus = nodes[0][:threads]
    root = Path(config['output_root']) / 'single24'
    plans = []
    for repetition in range(config['repetitions']):
        order = ['original', 'rabbitma'] if repetition % 2 == 0 else ['rabbitma', 'original']
        for program in order:
            label = '%s_rep%d' % (program, repetition + 1)
            plans.append({'label': label, 'jobs': [job(config, program, cpus, root / label / 'assembly')]})
    if args.dry_run:
        print(json.dumps({'experiment': 'single24', 'plans': plans}, indent=2))
        return
    directory = prepare(config, 'single24', plans)
    results = []
    for plan in plans:
        results.append(run_group(config, directory, plan['label'], plan['jobs']))
        report(directory, results)
    reference = Path(next(result['jobs'][0]['output'] for result in results
                          if result['jobs'][0]['program'] == 'original')) / 'final.contigs.fa'
    validation = {result['label']: compare_outputs(
        reference, Path(result['jobs'][0]['output']) / 'final.contigs.fa') for result in results}
    write_json(directory / 'validation.json', validation)
    print('Experiment 1 complete: ' + str(directory / 'report.md'), flush=True)


if __name__ == '__main__':
    main()
