# Benchmarks and validation

## CAMI III reference measurement

The reference workload used all 20 short-read samples from the CAMI III Toy
Longitudinal Human Gut dataset (about 100 Gbp, paired-end 150 bp). The reads
were combined into one 92.13 GB (85.80 GiB) interleaved gzip file to exercise
the single-large-input path. Both versions used 128 threads and `--k-min 39` on
the same machine.

| Version | Wall time | Peak RSS | Speedup | Peak-RSS reduction |
| --- | ---: | ---: | ---: | ---: |
| Official MEGAHIT v1.2.9 | 4:12:32 | 78.94 GiB | 1.00x | - |
| RabbitMA bounded-memory run | 57:54.82 | 19.28 GiB | 4.36x | 75.6% |

RabbitMA reduced peak resident memory by 59.66 GiB, using about one quarter
(`1/4.1`) of the official version's peak RSS while also reducing wall time.
The 19.28 GiB figure is an observed peak, not a hard-coded memory cap. It comes
from general bounded-working-set policies: count records and graph buckets are
streamed, `seq2sdbg` limits simultaneously live partitions, completed mapping
chunks are released promptly, and the exact mercy-edge index avoids retaining
a second full edge collection. These policies derive their budgets from the
run-time memory allowance and workload histograms rather than the CAMI input
size or a fixed machine topology.

This table describes the earlier bounded-memory reference build. The current
source also includes later reusable read-index, graph and local-assembly
optimizations, plus automatic NUMA placement for cooperating jobs. The table
is not a timing of this updated revision. Full-pipeline timing of the exact
public commit and its resource configuration remains necessary for a current
release comparison; component tests or timings with experimental switches
should not be presented as default end-to-end performance.

The measured workload can be reproduced with a command of this form:

```bash
megahit \
  --12 cami3_20samples_merged_interleaved.fq.gz \
  -t 128 --k-min 39 -o cami3_coassembly
```

Performance depends on read composition, k-mer schedule, storage, compiler,
memory topology, and CPU. RabbitMA discovers the available CPU/NUMA topology at
run time and does not encode the CAMI input size or a particular socket layout.

## v0.2.0 NUMA experiment and synchronization checks

Two concurrent jobs each assembled the same merged 20-sample CAMI III input
on a machine with two NUMA domains and 32 physical cores per domain. Each job
used `-t 32`. The comparison held the validated computation core fixed and
changed the launcher's resource placement and memory-budget policy.

| Placement policy | Job A (seconds) | Job B (seconds) | Pair makespan (seconds) |
| --- | ---: | ---: | ---: |
| Previous placement | 4282.94 | 4323.84 | 4324.07 |
| Automatic per-domain reservations | 2353.85 | 2340.24 | 2354.42 |

Node throughput improved **1.8366x** and pair completion time decreased
45.55%. All 240 final/intermediate output checks passed. This includes the
benefits of removing overlapping CPU placement, preferring local memory and
adjusting per-job memory budgets; it does not isolate remote-memory latency
or measure a universal speedup over v0.1.0. Performance on a single exclusive
job or a different topology can differ.

After synchronizing that computation implementation into RabbitMA, all three
CPU variants passed 87 regression tests each (261 runs). Eight complete toy
pipelines and a full merged CAMI III k39 assembly comparison passed 142 exact
output checks. These checks retained complete sequence orientation, k, flags
and floating coverage and ignored only record IDs and order, including
records stored in packed intermediate sidecars. The k39 paired ABBA means
were 74.0976 seconds for the validated core and 74.0274 seconds for the
synchronized core, within run-to-run variation. These are synchronization
checks, not new full-pipeline measurements of the public release binary.

Experimental `MEGAHIT_EXPERIMENTAL_*` paths remain opt-in; timings using them
must be labeled with their enabled flags. NUMA graph replicas and native/LTO
builds that did not improve measured performance are not enabled by default.

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
same canonicalization and five represented the same circular edge sets from a
different origin. The v0.1.0 audit additionally covered four real single-end
datasets, 1/2/3/high-thread single-end local assembly, mixed library types,
ambiguous bases, variable read lengths, paths containing spaces and quotes,
and all three CPU-dispatch binaries. Every non-circular final contig matched
the official normalized sequence, flag, and multiplicity multiset; circular
contigs matched the same complete graph-edge sets.

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
