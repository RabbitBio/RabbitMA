#!/usr/bin/env python3
"""Compare MEGAHIT contig multisets using Python 3.6+ and its standard library."""

import argparse
from collections import Counter, namedtuple
import io
import json
import math
from pathlib import Path
import re
import struct
import sys


MAGIC = b'MGCTG01\0'
FILE_HEADER = struct.Struct('<8sII')
RECORD_HEADER = struct.Struct('<IIifq')
Record = namedtuple('Record', 'sequence k flag coverage')
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
BYTE_BASES = tuple(''.join('ACGT'[(value >> shift) & 3]
                           for shift in (6, 4, 2, 0))
                   for value in range(256))


def coverage_key(value):
    """Match the float32 precision stored in MEGAHIT packed records."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('coverage must be finite')
    return struct.pack('<f', value)


def fasta_header(header):
    fields = header.split()
    if not fields:
        raise ValueError('empty FASTA header')
    match = re.match(r'k([0-9]+)(?:_|$)', fields[0])
    k = int(match.group(1)) if match else None
    if fields[0].startswith('lc_'):
        k = 0  # Local-assembly headers do not carry a k value.
    metadata = {}
    for field in fields[1:]:
        name, separator, value = field.partition('=')
        if separator and name in ('flag', 'multi', 'len'):
            if name in metadata:
                raise ValueError('duplicate FASTA field: ' + name)
            metadata[name] = value
    flag = int(metadata['flag']) if 'flag' in metadata else None
    coverage = coverage_key(metadata['multi']) if 'multi' in metadata else None
    length = int(metadata['len']) if 'len' in metadata else None
    if (flag is not None and flag < 0) or (length is not None and length < 0):
        raise ValueError('negative FASTA flag or length')
    return k, flag, coverage, length


def fasta_record(header, lines):
    k, flag, coverage, length = fasta_header(header)
    sequence = ''.join(lines)
    if not sequence:
        raise ValueError('empty FASTA record: ' + header)
    if length is not None and len(sequence) != length:
        raise ValueError('FASTA length disagrees with len=: ' + header)
    return Record(sequence, k, flag, coverage)


def read_fasta(stream):
    header, lines = None, []
    for number, line in enumerate(stream, 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield fasta_record(header, lines)
            header, lines = line[1:], []
        else:
            if header is None:
                raise ValueError('sequence before FASTA header at line %d' % number)
            sequence = ''.join(line.split()).upper()
            if any(base not in 'ACGTN' for base in sequence):
                raise ValueError('expected A/C/G/T/N at line %d' % number)
            lines.append(sequence)
    if header is not None:
        yield fasta_record(header, lines)


def read_packed(stream):
    header = stream.read(FILE_HEADER.size)
    if len(header) != FILE_HEADER.size:
        raise ValueError('truncated packed file header')
    magic, version, _ = FILE_HEADER.unpack(header)
    if magic != MAGIC or version != 1:
        raise ValueError('expected little-endian MGCTG01 version 1')
    while True:
        raw = stream.read(RECORD_HEADER.size)
        if not raw:
            return
        if len(raw) != RECORD_HEADER.size:
            raise ValueError('truncated packed record header')
        length, k, flag, multi, _ = RECORD_HEADER.unpack(raw)
        if not length or flag < 0:
            raise ValueError('invalid packed record length or flag')
        byte_count = ((length + 15) // 16) * 4
        payload = stream.read(byte_count)
        if len(payload) != byte_count:
            raise ValueError('truncated packed sequence')
        sequence = ''.join(
            BYTE_BASES[(word >> 24) & 255] + BYTE_BASES[(word >> 16) & 255] +
            BYTE_BASES[(word >> 8) & 255] + BYTE_BASES[word & 255]
            for (word,) in struct.iter_unpack('<I', payload))[:length]
        yield Record(sequence, k, flag, coverage_key(multi))


def read_records(path):
    with Path(path).open('rb') as stream:
        magic = stream.read(len(MAGIC))
        stream.seek(0)
        if magic == MAGIC or str(path).endswith('.mgb'):
            yield from read_packed(stream)
        else:
            with io.TextIOWrapper(stream, encoding='ascii') as text:
                yield from read_fasta(text)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def minimum_rotation(sequence):
    """Find the lexicographically smallest circular rotation in linear time."""
    size = len(sequence)
    doubled = sequence + sequence
    first, second, matched = 0, 1, 0
    while first < size and second < size and matched < size:
        left, right = doubled[first + matched], doubled[second + matched]
        if left == right:
            matched += 1
            continue
        if left > right:
            first += matched + 1
            if first <= second:
                first = second + 1
        else:
            second += matched + 1
            if second <= first:
                second = first + 1
        matched = 0
    start = min(first, second)
    return doubled[start:start + size]


def sequence_key(record, mode):
    sequence = record.sequence
    if mode == 'literal':
        return sequence
    if mode == 'circular' and (record.flag or 0) & 2:
        k = record.k
        if k is None or k <= 0 or len(sequence) <= k:
            raise ValueError('circular normalization requires a positive k and length > k')
        if sequence[:k] != sequence[-k:]:
            raise ValueError('circular record does not repeat its first k bases at the end')
        sequence = sequence[:-k]
        return min(minimum_rotation(sequence),
                   minimum_rotation(reverse_complement(sequence)))
    return min(sequence, reverse_complement(sequence))


def collect(path, mode, metadata):
    records = Counter()
    total = bases = circular = 0
    try:
        for record in read_records(path):
            total += 1
            bases += len(record.sequence)
            circular += bool((record.flag or 0) & 2)
            if metadata and any(value is None for value in
                                (record.k, record.flag, record.coverage)):
                raise ValueError('--metadata requires MEGAHIT k, flag and multi fields')
            key = (sequence_key(record, mode), record.flag or 0)
            if metadata:
                key += (record.k, record.coverage)
            records[key] += 1
    except (ValueError, UnicodeError, OverflowError, struct.error) as error:
        raise ValueError('%s (record %d): %s' % (path, total or 1, error)) from error
    return records, total, bases, circular


def compare(left_path, right_path, mode='strand', metadata=False):
    left, left_total, left_bases, left_circular = collect(left_path, mode, metadata)
    right, right_total, right_bases, right_circular = collect(right_path, mode, metadata)
    missing = sum((left - right).values())
    added = sum((right - left).values())
    return {
        'equal': missing == 0 and added == 0,
        'mode': mode,
        'metadata': metadata,
        'left_records': left_total,
        'right_records': right_total,
        'left_bases': left_bases,
        'right_bases': right_bases,
        'left_circular_records': left_circular,
        'right_circular_records': right_circular,
        'common_records': left_total - missing,
        'missing_records': missing,
        'added_records': added,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Compare contig multisets, ignoring IDs and record order. '
                    'Default: allow reverse complements and preserve parsed flags. '
                    'Uses full sequences, not sequence fingerprints.',
        epilog='Exit codes: 0 equal; 1 different; 2 invalid input or arguments. '
               'Inputs: plain FASTA or little-endian MGCTG01 version 1 .mgb files.')
    parser.add_argument('left', help='reference contig file')
    parser.add_argument('right', help='candidate contig file')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--literal', action='store_true',
                       help='preserve sequence direction and circular origin')
    modes.add_argument('--circular-loops', action='store_true',
                       help='also normalize origins of flagged loops after validating '
                            'and removing their k-base terminal overlap')
    parser.add_argument('--metadata', action='store_true',
                        help='also compare k and float32 coverage; require MEGAHIT metadata')
    args = parser.parse_args(argv)
    mode = 'literal' if args.literal else 'circular' if args.circular_loops else 'strand'
    try:
        result = compare(args.left, args.right, mode, args.metadata)
    except (OSError, ValueError, MemoryError) as error:
        print('compare_contigs: %s' % error, file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result['equal'] else 1


if __name__ == '__main__':
    sys.exit(main())
