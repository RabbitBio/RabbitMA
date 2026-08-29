# Benchmarks and validation

## CAMI III reference measurement

The performance-development reference used all 20 short-read samples from the
CAMI III Toy Longitudinal Human Gut dataset (about 100 Gbp, paired-end 150 bp),
128 threads, and `--k-min 39`.

| Version | Wall time | Peak RSS | Speedup |
| --- | ---: | ---: | ---: |
| Official MEGAHIT v1.2.9 | 4:12:32 | 78.94 GiB | 1.00x |
| RabbitMA development checkpoint | 1:07:26 | 68.22 GiB | 3.74x |

These figures are historical measurements from the same machine and input.
The RabbitMA source in this initial repository also includes later bounded-
memory and correctness changes. The complete 100 Gbp workload has not been
retimed after all of those changes, so the table is a development reference,
not a newly reproduced benchmark for this exact commit.

The workload shape can be reproduced with a command of this form:

```bash
rabbitma \
  -1 sample01_R1.fastq.gz,...,sample20_R1.fastq.gz \
  -2 sample01_R2.fastq.gz,...,sample20_R2.fastq.gz \
  -t 128 --k-min 39 -o cami3_coassembly
```

Performance depends on read composition, k-mer schedule, storage, compiler,
memory topology, and CPU. RabbitMA discovers the available CPU/NUMA topology at
run time and does not encode the CAMI input size or a particular socket layout.

## Correctness work included in this source

The optimized count path emits independently bucketed shards rather than one
globally sorted edge array. Its mercy-edge lookup therefore uses an exact
compact index whose directory is at least as fine as the physical 8-base
bucket prefix. Candidate queries perform full exact comparisons; no
approximate sequence evidence is introduced.

The source also preserves deterministic upstream traversal semantics and the
handling of self-reverse-complement and circular unitigs across k rounds.

Six simulated datasets were run with identical assembly parameters against
official MEGAHIT v1.2.9. The datasets exercise:

- unique sequence;
- strain bubbles;
- repeats;
- uneven and error-prone coverage;
- circular sequence; and
- single-end tips.

Across 100 final and per-k artifact comparisons, 95 matched exactly after the
same canonicalization and five were circular-origin-only matches. There were no
missing artifacts or substantive sequence/multiplicity failures.

## Main performance changes

- Parallel multi-library FASTA/FASTQ ingestion with optional libdeflate and a
  bundled adaptive rapidgzip path for large individual gzip files.
- One-pass segmented exact counting with compact records and streamed edge
  output.
- Streamed and bucket-packed `seq2sdbg` construction with compact locators and
  bounded live working sets.
- Cached SDBG topology, parallel unitig traversal, and compact unitig
  adjacency/endpoint lookup.
- Sharded local-assembly mapping/index structures and compact mapping-result
  collection.
- Sharded iterative flank indexes and parallel edge collection.
- Runtime CPU/NUMA discovery and affinity-aware placement without assuming a
  fixed socket or core count.
