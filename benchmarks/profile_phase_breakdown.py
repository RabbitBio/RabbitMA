#!/usr/bin/env python3
"""Summarize per-k RabbitMA phases and graph-cleaning rounds from a run log."""

import argparse
from collections import defaultdict
import json
import re


FLOAT = r"([0-9]+(?:\.[0-9]+)?)"


def empty_k():
    return {
        'count': 0.0,
        'build_graph': 0.0,
        'assemble': 0.0,
        'cleaning': 0.0,
        'cleaning_tips': 0.0,
        'cleaning_bubbles': 0.0,
        'cleaning_complex_bubbles': 0.0,
        'cleaning_disconnect': 0.0,
        'cleaning_low_depth': 0.0,
        'refresh': 0.0,
        'local_assembly': 0.0,
        'iterate': 0.0,
        'iterate_index': 0.0,
        'low_depth': 0.0,
        'input_edges': None,
        'input_sequences': None,
        'input_bases': None,
        'graph_edges': None,
        'unitigs': None,
        'stable_candidate_bases': None,
        'stable_carried_contigs': None,
        'stable_carried_bases': None,
        'stable_omitted_windows': None,
    }


def empty_round(number):
    return {
        'round': number,
        'active_begin': None,
        'active_end': None,
        'changed': None,
        'tips': 0,
        'bubbles': 0,
        'complex_bubbles': 0,
        'disconnected': 0,
        'low_depth': 0,
        'time': 0.0,
        'tips_time': 0.0,
        'bubbles_time': 0.0,
        'complex_bubbles_time': 0.0,
        'disconnected_time': 0.0,
        'low_depth_time': 0.0,
    }


def match_int(pattern, line):
    match = re.search(pattern, line)
    return int(match.group(1)) if match else None


def parse_log(path):
    per_k = defaultdict(empty_k)
    rounds = defaultdict(dict)
    read_index = []
    pending_phase = None
    build_k = None
    assemble_k = None
    cleaning_round = None

    with open(path, encoding='utf-8', errors='replace') as stream:
        for line in stream:
            marker = 'Read-index profile: '
            if marker in line:
                fields = {}
                for key, value in re.findall(r'(\w+)=([^\s\'\"]+)',
                                             line.split(marker, 1)[1]):
                    if re.match(r'^\d+$', value):
                        fields[key] = int(value)
                    elif re.match(r'^\d+\.\d+$', value):
                        fields[key] = float(value)
                    else:
                        fields[key] = value
                read_index.append(fields)
            match = re.search(r'Extract solid .* for k = (\d+)', line)
            if match:
                pending_phase = (int(match.group(1)), 'count')
                continue
            match = re.search(r'Build graph for k = (\d+)', line)
            if match:
                build_k = int(match.group(1))
                pending_phase = (build_k, 'build_graph')
                continue
            match = re.search(r'Assemble contigs from SdBG for k = (\d+)', line)
            if match:
                assemble_k = int(match.group(1))
                cleaning_round = None
                pending_phase = (assemble_k, 'assemble')
                continue
            match = re.search(
                r'Build reusable exact read index(?: and fuse k = (\d+) to \d+)?',
                line)
            if match:
                source_k = int(match.group(1)) if match.group(1) else min(per_k)
                pending_phase = (source_k, 'iterate_index')
                continue
            match = re.search(r'(?:Query reusable read index|Extract iterative edges) '
                              r'from k = (\d+) to \d+', line)
            if match:
                pending_phase = (int(match.group(1)), 'iterate')
                continue
            match = re.search(r'Local assembly for k = (\d+)', line)
            if match:
                pending_phase = (int(match.group(1)), 'local_assembly')
                continue

            match = re.search(r'Real: ' + FLOAT, line)
            if match and pending_phase:
                k, phase = pending_phase
                elapsed = float(match.group(1))
                per_k[k][phase] += elapsed
                if phase == 'iterate_index':
                    per_k[k]['iterate'] += elapsed
                pending_phase = None
                continue

            if build_k is not None:
                value = match_int(r'Number edges: (\d+)', line)
                if value is not None:
                    per_k[build_k]['input_edges'] = value
                match = re.search(
                    r'Finally, sizeof seq_package: \d+/(\d+)/(\d+)', line)
                if match:
                    per_k[build_k]['input_sequences'] = int(match.group(1))
                    per_k[build_k]['input_bases'] = int(match.group(2))
                value = match_int(r'Total number of edges: (\d+)', line)
                if value is not None:
                    per_k[build_k]['graph_edges'] = value
                match = re.search(
                    r'Stable-contig certificate: (\d+) candidate bases, '
                    r'(\d+) carried contigs / (\d+) bases, .* '
                    r'(\d+) omitted oriented windows', line)
                if match:
                    row = per_k[build_k]
                    row['stable_candidate_bases'] = int(match.group(1))
                    row['stable_carried_contigs'] = int(match.group(2))
                    row['stable_carried_bases'] = int(match.group(3))
                    row['stable_omitted_windows'] = int(match.group(4))

            if assemble_k is None:
                continue
            match = re.search(r'(?:Delta)?Refresh profile: (.*)', line)
            if match:
                per_k[assemble_k]['refresh'] += sum(
                    float(value) for value in re.findall(
                        r'(?:disconnect|invalidate|merge|loops|compact|remap)='
                        r'([0-9.]+)', match.group(1)))
            value = match_int(r'Number of Edges: (\d+); K value:', line)
            if value is not None:
                per_k[assemble_k]['graph_edges'] = value
            match = re.search(r'unitig graph size: (\d+), time for building: ' +
                              FLOAT, line)
            if match:
                per_k[assemble_k]['unitigs'] = int(match.group(1))

            match = re.search(r'Graph cleaning round (\d+)', line)
            if match:
                cleaning_round = int(match.group(1))
                rounds[assemble_k][cleaning_round] = empty_round(cleaning_round)
                continue
            if cleaning_round is not None:
                row = rounds[assemble_k][cleaning_round]
                match = re.search(
                    r'Cleaning round input: k=\d+, round=\d+, active_unitigs=(\d+)',
                    line)
                if match:
                    row['active_begin'] = int(match.group(1))
                patterns = (
                    ('tips', r'Tips removed: (\d+), time: ' + FLOAT),
                    ('bubbles', r'Number of bubbles removed: (\d+), '
                                r'Time elapsed\(sec\): ' + FLOAT),
                    ('complex_bubbles', r'Number of complex bubbles removed: (\d+), '
                                        r'Time elapsed\(sec\): ' + FLOAT),
                    ('disconnected', r'Number unitigs disconnected: (\d+), time: ' +
                                     FLOAT),
                    ('low_depth', r'Unitigs removed in (?:\(more-\))?excessive '
                                  r'pruning: (\d+), time: ' + FLOAT),
                )
                for name, pattern in patterns:
                    match = re.search(pattern, line)
                    if match:
                        row[name] = int(match.group(1))
                        elapsed = float(match.group(2))
                        row[name + '_time'] = elapsed
                        row['time'] += elapsed
                        per_k[assemble_k]['cleaning'] += elapsed
                        per_k_name = ('cleaning_disconnect' if name == 'disconnected'
                                      else 'cleaning_' + name)
                        per_k[assemble_k][per_k_name] += elapsed
                        break
                match = re.search(
                    r'Cleaning round result: k=\d+, round=\d+, '
                    r'active_unitigs=(\d+), changed=(\d+)', line)
                if match:
                    row['active_end'] = int(match.group(1))
                    row['changed'] = int(match.group(2))

            match = re.search(
                r'Number of local low depth unitigs removed: \d+, complex bubbles '
                r'removed: \d+, time: ' + FLOAT, line)
            if match:
                per_k[assemble_k]['low_depth'] = float(match.group(1))

    # Old logs lack the structured changed total; all component counters are
    # still present and can be summed exactly.
    for by_round in rounds.values():
        for row in by_round.values():
            if row['changed'] is None:
                row['changed'] = sum(row[name] for name in
                                     ('tips', 'bubbles', 'complex_bubbles',
                                      'disconnected', 'low_depth'))
    return (dict(per_k),
            {k: list(rounds[k][number] for number in sorted(rounds[k]))
             for k in sorted(rounds)},
            read_index)


def number(value):
    return '-' if value is None else format(value, ',')


def seconds(value):
    return '%.4f' % value


def markdown(per_k, rounds, read_index):
    lines = [
        '| k | graph build | assemble | cleaning | refresh | local assembly | iterate | low-depth |',
        '| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for k in sorted(per_k):
        row = per_k[k]
        graph_build = row['count'] + row['build_graph']
        lines.append('| %d | %s | %s | %s | %s | %s | %s | %s |' % (
            k, seconds(graph_build), seconds(row['assemble']),
            seconds(row['cleaning']), seconds(row['refresh']),
            seconds(row['local_assembly']), seconds(row['iterate']),
            seconds(row['low_depth'])))

    lines += ['', '| k | input edges | input sequences | input bases | graph edges | unitigs |',
              '| ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in sorted(per_k):
        row = per_k[k]
        lines.append('| %d | %s | %s | %s | %s | %s |' % (
            k, number(row['input_edges']), number(row['input_sequences']),
            number(row['input_bases']), number(row['graph_edges']),
            number(row['unitigs'])))

    lines += ['', '| k | tips | bubbles | complex bubbles | disconnect | local low-depth |',
              '| ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in sorted(per_k):
        row = per_k[k]
        lines.append('| %d | %s | %s | %s | %s | %s |' % (
            k, seconds(row['cleaning_tips']),
            seconds(row['cleaning_bubbles']),
            seconds(row['cleaning_complex_bubbles']),
            seconds(row['cleaning_disconnect']),
            seconds(row['cleaning_low_depth'])))

    if any(per_k[k]['stable_candidate_bases'] is not None for k in per_k):
        lines += ['', '| k | stable candidates (bp) | carried paths | carried bases | omitted windows |',
                  '| ---: | ---: | ---: | ---: | ---: |']
        for k in sorted(per_k):
            row = per_k[k]
            lines.append('| %d | %s | %s | %s | %s |' % (
                k, number(row['stable_candidate_bases']),
                number(row['stable_carried_contigs']),
                number(row['stable_carried_bases']),
                number(row['stable_omitted_windows'])))

    lines += ['', '| k | round | active begin | active end | changed | time |',
              '| ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in sorted(rounds):
        for row in rounds[k]:
            lines.append('| %d | %d | %s | %s | %s | %s |' % (
                k, row['round'], number(row['active_begin']),
                number(row['active_end']), number(row['changed']),
                seconds(row['time'])))
    if read_index:
        lines += [
            '',
            '| mode | k → next | status | reason | time | indexed occurrences | matched | candidates | replayed | replay reads | replay bytes | fallback reads | required / budget | scan bytes avoided |',
            '| :--- | :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
        ]
        for row in read_index:
            elapsed = row.get('index_build_time', row.get('query_time', 0.0))
            lines.append(
                '| %s | %s → %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s / %s | %s |' % (
                    row.get('mode', '-'), row.get('k', '-'),
                    row.get('next_k', '-'), row.get('status', '-'),
                    row.get('reason', '-'), seconds(elapsed),
                    number(row.get('indexed_occurrences')),
                    number(row.get('matched_occurrences')),
                    number(row.get('candidate_upper_bound')),
                    number(row.get('replayed_occurrences')),
                    number(row.get('replay_reads')),
                    number(row.get('replay_bytes')),
                    number(row.get('fallback_reads')),
                    number(row.get('required_budget_bytes')),
                    number(row.get('budget_bytes')),
                    number(row.get('full_scan_bytes_avoided'))))
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('log', help='RabbitMA output log')
    parser.add_argument('--json', action='store_true', help='emit JSON')
    args = parser.parse_args()
    per_k, rounds, read_index = parse_log(args.log)
    if args.json:
        print(json.dumps({'per_k': per_k, 'cleaning_rounds': rounds,
                          'read_index': read_index},
                         indent=2, sort_keys=True))
    else:
        print(markdown(per_k, rounds, read_index), end='')


if __name__ == '__main__':
    main()
