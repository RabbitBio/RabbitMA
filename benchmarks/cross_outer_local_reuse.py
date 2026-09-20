#!/usr/bin/env python3
"""Measure local-assembly input reuse between adjacent outer-k rounds."""

import argparse
import csv
from collections import defaultdict, namedtuple
import json
from pathlib import Path
import re


NAME = re.compile(r'k(?P<from_k>\d+)_to_k(?P<to_k>\d+)')


Record = namedtuple('Record', [
    'full', 'reads', 'potential', 'actual', 'read_count', 'endpoint_length',
    'rounds'])
Round = namedtuple('Round', ['from_k', 'to_k', 'path', 'records'])


def load_round(path):
    match = NAME.search(path.name)
    if not match:
        raise ValueError(f'{path}: filename must contain k<from>_to_k<to>')
    records = []
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        expected = {'low', 'high', 'reads_low', 'reads_high', 'potential',
                    'actual', 'reads', 'endpoint', 'rounds'}
        if set(reader.fieldnames or []) != expected:
            raise ValueError(f'{path}: unexpected fingerprint columns')
        for line, row in enumerate(reader, 2):
            try:
                records.append(Record(
                    full=(int(row['low'], 16), int(row['high'], 16)),
                    reads=(int(row['reads_low'], 16),
                           int(row['reads_high'], 16)),
                    potential=int(row['potential']), actual=int(row['actual']),
                    read_count=int(row['reads']),
                    endpoint_length=int(row['endpoint']),
                    rounds=int(row['rounds'])))
            except ValueError as error:
                raise ValueError(f'{path}:{line}: invalid value: {error}')
    return Round(int(match.group('from_k')), int(match.group('to_k')),
                 path, tuple(records))


def matched_pairs(previous, current, attribute):
    old = defaultdict(list)
    new = defaultdict(list)
    for record in previous:
        old[getattr(record, attribute)].append(record)
    for record in current:
        new[getattr(record, attribute)].append(record)
    pairs = []
    for key in old.keys() & new.keys():
        # Pair the most expensive duplicates first. This gives the maximum
        # reusable prefix for a multiset key when its multiplicity changes.
        left = sorted(old[key], key=lambda r: (r.rounds, r.actual),
                      reverse=True)
        right = sorted(new[key], key=lambda r: (r.rounds, r.actual),
                       reverse=True)
        pairs.extend(zip(left, right))
    return pairs


def reuse_stats(previous, current, attribute):
    pairs = matched_pairs(previous, current, attribute)
    repeated_rounds = 0
    repeated_actual = 0
    for old, new in pairs:
        repeated_rounds += min(old.rounds, new.rounds)
        # Equal ordered reads perform identical read-graph work for every
        # shared inner-k prefix. The record with fewer completed rounds is
        # therefore the exact prefix-work counter.
        if old.rounds < new.rounds:
            repeated_actual += old.actual
        elif new.rounds < old.rounds:
            repeated_actual += new.actual
        else:
            if old.actual != new.actual:
                raise ValueError(
                    f'equal-round {attribute} fingerprints have unequal '
                    f'actual work ({old.actual} != {new.actual})')
            repeated_actual += old.actual
    return {
        'matched_tasks': len(pairs),
        'matched_current_actual': sum(new.actual for _, new in pairs),
        'matched_current_potential': sum(new.potential for _, new in pairs),
        'repeated_prefix_actual': repeated_actual,
        'repeated_prefix_rounds': repeated_rounds,
    }


def ratio(value, total):
    return 0.0 if total == 0 else value / total


def analyze(rounds):
    output = []
    for previous, current in zip(rounds, rounds[1:]):
        if previous.to_k != current.from_k:
            raise ValueError(
                f'non-adjacent rounds: {previous.from_k}->{previous.to_k} '
                f'and {current.from_k}->{current.to_k}')
        total = {
            'tasks': len(current.records),
            'actual': sum(record.actual for record in current.records),
            'potential': sum(record.potential for record in current.records),
            'rounds': sum(record.rounds for record in current.records),
        }
        full = reuse_stats(previous.records, current.records, 'full')
        reads = reuse_stats(previous.records, current.records, 'reads')
        output.append({
            'previous_transition': f'{previous.from_k}->{previous.to_k}',
            'transition': f'{current.from_k}->{current.to_k}',
            'source': str(previous.path),
            'fingerprints': str(current.path),
            'total': total,
            'full_input': full,
            'read_set': reads,
        })
    return output


def percent(value, total):
    return f'{ratio(value, total) * 100:.2f}%'


def aggregate(rows):
    result = {
        'previous_transition': 'all',
        'transition': 'ALL',
        'source': '',
        'fingerprints': '',
        'total': {key: 0 for key in ('tasks', 'actual', 'potential', 'rounds')},
        'full_input': {key: 0 for key in (
            'matched_tasks', 'matched_current_actual',
            'matched_current_potential', 'repeated_prefix_actual',
            'repeated_prefix_rounds')},
        'read_set': {key: 0 for key in (
            'matched_tasks', 'matched_current_actual',
            'matched_current_potential', 'repeated_prefix_actual',
            'repeated_prefix_rounds')},
    }
    for row in rows:
        for section in ('total', 'full_input', 'read_set'):
            for key, value in row[section].items():
                result[section][key] += value
    return result


def markdown(rows):
    lines = [
        '| transition | tasks | full-input reuse | read-set reuse | '
        'full matched actual | read matched actual | full repeated actual | '
        'read repeated actual | full repeated rounds | read repeated rounds |',
        '| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for row in rows + [aggregate(rows)]:
        total = row['total']
        full = row['full_input']
        reads = row['read_set']
        lines.append(
            '| {transition} | {tasks:,} | {full_tasks} | {read_tasks} | '
            '{full_matched} | {read_matched} | {full_actual} | '
            '{read_actual} | {full_rounds} | {read_rounds} |'.format(
                transition=row['transition'], tasks=total['tasks'],
                full_tasks=percent(full['matched_tasks'], total['tasks']),
                read_tasks=percent(reads['matched_tasks'], total['tasks']),
                full_matched=percent(full['matched_current_actual'],
                                     total['actual']),
                read_matched=percent(reads['matched_current_actual'],
                                     total['actual']),
                full_actual=percent(full['repeated_prefix_actual'],
                                    total['actual']),
                read_actual=percent(reads['repeated_prefix_actual'],
                                    total['actual']),
                full_rounds=percent(full['repeated_prefix_rounds'],
                                    total['rounds']),
                read_rounds=percent(reads['repeated_prefix_rounds'],
                                    total['rounds'])))
    return '\n'.join(lines) + '\n'


def collect_paths(arguments):
    paths = []
    for value in arguments:
        path = Path(value)
        if path.is_dir():
            paths.extend(path.glob('local_k*_to_k*.tsv'))
        else:
            paths.append(path)
    rounds = [load_round(path) for path in paths]
    rounds.sort(key=lambda item: (item.from_k, item.to_k))
    if len(rounds) < 2:
        raise ValueError('at least two fingerprint rounds are required')
    transitions = [(item.from_k, item.to_k) for item in rounds]
    if len(set(transitions)) != len(transitions):
        raise ValueError('duplicate transition fingerprint files')
    return rounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fingerprints', nargs='+',
                        help='fingerprint TSV files or containing directories')
    parser.add_argument('--json', action='store_true', help='emit JSON')
    args = parser.parse_args()
    rows = analyze(collect_paths(args.fingerprints))
    if args.json:
        print(json.dumps({'comparisons': rows, 'aggregate': aggregate(rows)},
                         indent=2, sort_keys=True))
    else:
        print(markdown(rows), end='')


if __name__ == '__main__':
    main()
