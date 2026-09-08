"""Shared planning and measurement for the two opt-in CAMI3 experiments."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import time


REPOSITORY = Path(__file__).resolve().parents[1]


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def index_list(text):
    result = set()
    for part in text.strip().split(','):
        if not part:
            continue
        limits = part.split('-')
        begin, end = int(limits[0]), int(limits[-1])
        result.update(range(begin, end + 1))
    return result


def topology(total_cpus):
    allowed = os.sched_getaffinity(0)
    allowed &= index_list(Path('/sys/devices/system/cpu/online').read_text())
    memory_nodes = None
    for line in Path('/proc/self/status').read_text().splitlines():
        if line.startswith('Mems_allowed_list:'):
            memory_nodes = index_list(line.split(':', 1)[1])
    cores = {}
    for cpu in sorted(allowed):
        root = Path('/sys/devices/system/cpu/cpu%d/topology' % cpu)
        key = (int((root / 'physical_package_id').read_text()),
               int((root / 'core_id').read_text()))
        cores.setdefault(key, cpu)  # One logical CPU per physical core.
    physical = set(cores.values())
    nodes = {}
    for path in sorted(Path('/sys/devices/system/node').glob('node[0-9]*')):
        node = int(path.name[4:])
        if memory_nodes is not None and node not in memory_nodes:
            continue
        cpus = sorted(index_list((path / 'cpulist').read_text()) & physical)
        if cpus:
            nodes[node] = cpus
    if sum(map(len, nodes.values())) < total_cpus:
        raise ValueError('need %d allowed, online physical cores; discovered %d' %
                         (total_cpus, sum(map(len, nodes.values()))))
    # Balance the selected budget across the available NUMA domains.
    selected = {node: [] for node in nodes}
    while sum(map(len, selected.values())) < total_cpus:
        for node, cpus in sorted(nodes.items()):
            index = len(selected[node])
            if index < len(cpus):
                selected[node].append(cpus[index])
            if sum(map(len, selected.values())) == total_cpus:
                break
    return {node: cpus for node, cpus in selected.items() if cpus}


def divide_threads(total, jobs):
    if jobs < 1 or jobs > total:
        raise ValueError('job count must be between 1 and the total CPU budget')
    quotient, remainder = divmod(total, jobs)
    return [quotient + (index < remainder) for index in range(jobs)]


def interleaved_cpus(nodes):
    result = []
    for index in range(max(map(len, nodes.values()))):
        for _, cpus in sorted(nodes.items()):
            if index < len(cpus):
                result.append(cpus[index])
    return result


def allocate(nodes, threads, placement):
    if sum(threads) > sum(map(len, nodes.values())):
        raise ValueError('thread requests exceed the selected physical cores')
    if placement == 'spread':
        pool = interleaved_cpus(nodes)
        masks, offset = [], 0
        for count in threads:
            masks.append(sorted(pool[offset:offset + count]))
            offset += count
        return masks
    if placement != 'numa':
        raise ValueError('unknown placement: ' + placement)
    free = {node: list(cpus) for node, cpus in nodes.items()}
    masks = []
    for count in threads:
        fitting = [node for node in free if len(free[node]) >= count]
        if fitting:
            node = min(fitting, key=lambda value: (-len(free[value]), value))
            mask, free[node] = free[node][:count], free[node][count:]
        else:
            mask = []
            while len(mask) < count:
                for node in sorted(free):
                    if free[node]:
                        mask.append(free[node].pop(0))
                    if len(mask) == count:
                        break
        masks.append(sorted(mask))
    return masks


def interrupt(signum, frame):
    # Let finally blocks stop this experiment's children on SIGTERM too.
    raise KeyboardInterrupt('received signal %d' % signum)


def parse_arguments(description):
    signal.signal(signal.SIGTERM, interrupt)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', required=True, help='JSON experiment configuration')
    parser.add_argument('--dry-run', action='store_true',
                        help='print the plan only; do not write results or start assemblers')
    return parser.parse_args()


def load_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    for key in ('input', 'output_root'):
        config[key] = str((path.parent / config[key]).resolve())
    for program in config['programs'].values():
        program['launcher'] = str((path.parent / program['launcher']).resolve())
        if not Path(program['launcher']).is_file():
            raise ValueError('missing launcher: ' + program['launcher'])
    if not Path(config['input']).is_file():
        raise ValueError('missing input: ' + config['input'])
    config.setdefault('total_cpus', 64)
    config.setdefault('single_threads', 24)
    config.setdefault('repetitions', 2)
    config.setdefault('concurrency_levels', [2, 3, 4, 5, 6])
    config.setdefault('placements', ['numa', 'spread'])
    config.setdefault('k_list', [39, 49, 59, 69, 79, 89, 99, 109, 119, 129, 139, 141])
    config.setdefault('memory_budget_gib', 192)
    config.setdefault('warm_input', True)
    config.setdefault('independent_inputs', True)
    if config['repetitions'] < 1 or config['single_threads'] > config['total_cpus']:
        raise ValueError('invalid repetition count or single-job thread budget')
    if config['single_threads'] < 1 or config['memory_budget_gib'] <= 0:
        raise ValueError('thread and memory budgets must be positive')
    if not config['placements'] or any(value not in ('numa', 'spread')
                                       for value in config['placements']):
        raise ValueError('placements must contain numa and/or spread')
    config['nodes'] = topology(config['total_cpus'])
    config['memory_bytes'] = int(config['memory_budget_gib'] * 2**30)
    config['input_bytes'] = Path(config['input']).stat().st_size
    config['config_path'] = str(path)
    return config


def job(config, program, cpus, output, input_path=None):
    threads = len(cpus)
    memory = config['memory_bytes'] * threads // config['total_cpus']
    input_path = str(input_path or config['input'])
    command = [sys.executable, config['programs'][program]['launcher'],
               '--12', input_path, '-t', str(threads), '-m', str(memory),
               '--k-list', ','.join(map(str, config['k_list'])),
               '--min-count', '2', '--min-contig-len', '200', '-o', str(output)]
    return {'program': program, 'cpus': cpus, 'threads': threads,
            'memory_bytes': memory, 'input': input_path, 'output': str(output),
            'command': command,
            'numa_nodes': [node for node, available in sorted(config['nodes'].items())
                           if set(cpus) & set(available)]}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def prepare(config, phase, plans):
    directory = Path(config['output_root']) / phase
    if directory.exists():
        raise ValueError('result phase already exists: ' + str(directory))
    available = next(int(line.split()[1]) * 1024
                     for line in Path('/proc/meminfo').read_text().splitlines()
                     if line.startswith('MemAvailable:'))
    if config['memory_bytes'] > available:
        raise ValueError('configured memory budget exceeds current MemAvailable; '
                         'choose a smaller budget or wait for the machine to be free')
    if not Path('/usr/bin/time').is_file():
        raise ValueError('GNU /usr/bin/time is required')
    programs = {}
    for name, description in config['programs'].items():
        launcher = Path(description['launcher']).resolve()
        paths = [launcher] + sorted(launcher.parent.glob('megahit_core*'))
        hashes = {str(path): sha256(path) for path in paths if path.is_file()}
        expected = description.get('expected_launcher_sha256')
        if expected and hashes[str(launcher)] != expected:
            raise ValueError('launcher changed since preparation: ' + name)
        for filename, expected in description.get('expected_core_sha256', {}).items():
            if hashes.get(str(launcher.parent / filename)) != expected:
                raise ValueError('core changed since preparation: %s/%s' % (name, filename))
        programs[name] = dict(description, sha256=hashes)
    manifest = {
        'created': utc_now(), 'config': config, 'programs': programs,
        'input_sha256': sha256(config['input']), 'plans': plans,
        'python': sys.version,
        'host': dict(zip(('sysname', 'nodename', 'release', 'version', 'machine'), os.uname())),
        'openmp': {'OMP_PLACES': 'cores', 'OMP_PROC_BIND': 'close'},
        'note': 'Input preparation and output comparison are outside assembly timing.'}
    directory.mkdir(parents=True, exist_ok=False)
    write_json(directory / 'manifest.json', manifest)
    return directory


def warm(paths):
    for path in sorted(set(paths)):
        with Path(path).open('rb') as stream:
            while stream.read(16 * 2**20):
                pass


def descendants(pid):
    pending, seen = [pid], set()
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            children = Path('/proc/%d/task/%d/children' % (current, current)).read_text()
            pending.extend(map(int, children.split()))
        except OSError:
            pass
    return seen


def snapshot(pid, allowed):
    rss = 0
    masks = Counter()
    violations = []
    for child in descendants(pid):
        root = Path('/proc/%d' % child)
        try:
            for line in (root / 'status').read_text().splitlines():
                if line.startswith('VmRSS:'):
                    rss += int(line.split()[1])
            if not (root / 'comm').read_text().startswith('megahit_core'):
                continue
            for status in (root / 'task').glob('*/status'):
                try:
                    for line in status.read_text().splitlines():
                        if line.startswith('Cpus_allowed_list:'):
                            mask = line.split(':', 1)[1].strip()
                            masks[mask] += 1
                            if not index_list(mask) <= set(allowed):
                                violations.append({'tid': status.parent.name, 'cpus': mask})
                except OSError:
                    pass
        except OSError:
            pass
    return {'rss_kib': rss, 'thread_masks': dict(masks), 'violations': violations}


def time_fields(path):
    fields = {}
    names = {'User time (seconds)': 'user_seconds',
             'System time (seconds)': 'system_seconds',
             'Maximum resident set size (kbytes)': 'max_process_rss_kib'}
    for line in Path(path).read_text().splitlines():
        key, separator, value = line.strip().partition(': ')
        if separator and key in names:
            fields[names[key]] = float(value)
    return fields


def run_group(config, directory, label, specs):
    directory = Path(directory) / label
    directory.mkdir()
    if config['warm_input']:
        print('Warming inputs for ' + label, flush=True)
        warm([spec['input'] for spec in specs])
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(('OMP_', 'GOMP_', 'KMP_', 'MEGAHIT_'))
                   and key not in ('LD_LIBRARY_PATH', 'LD_PRELOAD')}
    environment.update(OMP_PLACES='cores', OMP_PROC_BIND='close', LC_ALL='C')
    running, finished = [], []
    start = time.monotonic()
    peak_rss = 0
    violations = []
    print('START %s: %s threads' % (label, [spec['threads'] for spec in specs]), flush=True)
    try:
        with (directory / 'samples.jsonl').open('w') as samples:
            for index, spec in enumerate(specs):
                console = (directory / ('job%d.console.log' % index)).open('w')
                resource = directory / ('job%d.time.txt' % index)
                command = ['/usr/bin/time', '-v', '-o', str(resource)] + spec['command']
                try:
                    process = subprocess.Popen(
                        command, stdout=console, stderr=subprocess.STDOUT, env=environment,
                        start_new_session=True,
                        preexec_fn=lambda cpus=spec['cpus']: os.sched_setaffinity(0, cpus))
                except BaseException:
                    console.close()
                    raise
                running.append({'process': process, 'console': console, 'spec': spec,
                                'resource': resource, 'started': time.monotonic(), 'index': index})
            next_sample = 0
            while running:
                now = time.monotonic()
                if now >= next_sample:
                    observations = {str(item['index']): snapshot(
                        item['process'].pid, item['spec']['cpus']) for item in running}
                    peak_rss = max(peak_rss, sum(value['rss_kib'] for value in observations.values()))
                    for value in observations.values():
                        violations.extend(value['violations'])
                    samples.write(json.dumps({'seconds': now - start, 'jobs': observations}) + '\n')
                    samples.flush()
                    next_sample = now + 10
                for item in list(running):
                    status = item['process'].poll()
                    if status is None:
                        continue
                    item['console'].close()
                    record = dict(item['spec'], returncode=status,
                                  seconds=time.monotonic() - item['started'])
                    record.update(time_fields(item['resource']))
                    finished.append(record)
                    running.remove(item)
                    write_json(directory / 'completed_jobs.json', finished)
                    if status != 0:
                        raise RuntimeError('%s job %d failed; inspect its console log' %
                                           (label, item['index']))
                    if not (Path(record['output']) / 'final.contigs.fa').is_file():
                        raise RuntimeError('successful exit without final.contigs.fa')
                if running:
                    time.sleep(0.1)
    finally:
        for item in running:
            if item['process'].poll() is None:
                try:
                    os.killpg(item['process'].pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        for item in running:
            try:
                item['process'].wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(item['process'].pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                item['process'].wait()
            item['console'].close()
    seconds = time.monotonic() - start
    result = {'label': label, 'seconds': seconds, 'jobs': finished,
              'samples_per_hour': len(finished) * 3600 / seconds,
              'node_hours_per_sample': seconds / (3600 * len(finished)),
              'sampled_sum_rss_kib': peak_rss, 'affinity_violations': violations}
    write_json(directory / 'result.json', result)
    if violations:
        raise RuntimeError('worker affinity escaped its assigned CPUs; inspect samples.jsonl')
    print('DONE %s: %.2f s, %.3f samples/hour' %
          (label, seconds, result['samples_per_hour']), flush=True)
    return result


def compare_outputs(left, right):
    spec = importlib.util.spec_from_file_location(
        'benchmark_comparator', str(REPOSITORY / 'tools/compare_contigs.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = {}
    for mode, metadata in [('literal', True), ('circular', True), ('circular', False)]:
        key = '%s_metadata_%s' % (mode, metadata)
        try:
            result[key] = module.compare(left, right, mode, metadata)
        except (OSError, ValueError) as error:
            result[key] = {'error': str(error)}
    return result


def report(directory, results):
    write_json(Path(directory) / 'results.json', results)
    rows = [{'label': result['label'], 'jobs': len(result['jobs']),
             'seconds': result['seconds'], 'samples_per_hour': result['samples_per_hour'],
             'node_hours_per_sample': result['node_hours_per_sample'],
             'sampled_sum_rss_gib': result['sampled_sum_rss_kib'] / 2**20}
            for result in results]
    with (Path(directory) / 'timings.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    groups = {}
    for row in rows:
        key = row['label'].split('_rep', 1)[0]
        groups.setdefault(key, []).append(row)
    lines = ['# Measured assembly throughput', '',
             '| Configuration | Repetitions | Median seconds | Median samples/hour | Node-hours/sample | Sampled sum RSS GiB |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for name, values in sorted(groups.items()):
        lines.append('| %s | %d | %.2f | %.3f | %.6f | %.3f |' % (
            name, len(values), statistics.median(value['seconds'] for value in values),
            statistics.median(value['samples_per_hour'] for value in values),
            statistics.median(value['node_hours_per_sample'] for value in values),
            max(value['sampled_sum_rss_gib'] for value in values)))
    if 'original' in groups and 'rabbitma' in groups:
        baseline = statistics.median(value['seconds'] for value in groups['original'])
        release = statistics.median(value['seconds'] for value in groups['rabbitma'])
        lines += ['', 'RabbitMA speedup relative to original: %.3fx; '
                  'elapsed-time change: %+.2f%% (median timings).' %
                  (baseline / release, (release / baseline - 1) * 100)]
    concurrent = {name: statistics.median(value['samples_per_hour'] for value in values)
                  for name, values in groups.items() if 'j_' in name}
    if concurrent:
        best = max(concurrent, key=concurrent.get)
        lines += ['', 'Highest median throughput among completed configurations: '
                  '%s, %.3f samples/hour.' % (best, concurrent[best])]
        for name in sorted(concurrent):
            if name.endswith('_numa'):
                spread = name[:-5] + '_spread'
                if spread in concurrent:
                    lines += ['', '%s / %s throughput ratio: %.3fx.' %
                              (name, spread, concurrent[name] / concurrent[spread])]
    lines += ['', 'Each job processes one complete copy of the configured input. '
              'Input warming/copying and output comparisons are outside the timings. '
              'Aggregate RSS is sampled every 10 seconds and can count shared pages '
              'more than once; per-process GNU time maxima are stored separately.', '',
              'This is finite-batch throughput, not a measured multi-thousand-sample '
              'steady-state queue. At a constant node price per hour, cost per sample '
              'is that price multiplied by node_hours_per_sample. '
              'See validation.json for sequence and metadata differences.', '']
    (Path(directory) / 'report.md').write_text('\n'.join(lines))
