# CAMI3 release and concurrency experiments

These scripts prepare and run two separate experiments. **They do not run on
import, and `--dry-run` never launches an assembler, copies input data, or
creates a result directory.** Run the real experiments only when the selected
machine resources are available.

Copy `config.example.json` and supply the input, a new output directory, and
the two Python launcher paths. Use an actual RabbitMA release package and
record its tag/source commit in `label`. The scripts record SHA-256 hashes
of the launchers, adjacent computation cores, and input. Optional
`expected_launcher_sha256` and `expected_core_sha256` fields in each program
entry reject binaries changed since preparation. No software is downloaded
or rebuilt by these scripts.

The input is one complete interleaved paired-end FASTQ file, for example a
CAMI3 single sample of about 4.6 GB compressed. Both experiments use the same
explicit k list, minimum count, minimum contig length, and proportional
memory-budget policy. `memory_budget_gib` is the total budget for all 64 CPUs;
a 24-thread job receives 24/64 of it. Set this value for the target machine.
It is an assembler budget, not an enforced cgroup RAM limit.

Requirements: Linux with visible NUMA/core topology, Python 3.6+, GNU
`/usr/bin/time`, the runtime dependencies of both assemblers, and enough
permitted online physical cores. One logical CPU per physical core is used;
SMT siblings are not counted as additional physical cores. No root privileges
are required. Existing scheduler CPU and memory-node allowances are respected.

Preview both experiments:

```bash
python3 benchmarks/run_experiments.py --config /path/to/config.json --dry-run
```

When the machine is available, run both with one command:

```bash
python3 benchmarks/run_experiments.py --config /path/to/config.json
```

The first experiment finishes before the second begins. Running the two
experiments concurrently would confound their measurements. Each output
phase must be new; existing results are never deleted or overwritten.
To run just one experiment, call `compare_release_24.py` or
`compare_concurrency.py` with the same `--config` argument.

## Experiment 1: both versions at 24 threads

`compare_release_24.py` binds both versions to the same 24 physical cores in
one NUMA domain. Jobs run individually. Two repetitions per version use the
order original, RabbitMA, RabbitMA, original. More repetitions continue
alternating the order. Both use `OMP_PLACES=cores` and `OMP_PROC_BIND=close`;
inherited MEGAHIT/OpenMP tuning and library-path overrides are removed from
the child environment. RabbitMA retains its normal resource-aware behavior
within this CPU mask. This compares complete pipelines, including reading,
decompression and writing, not just one assembly kernel.

## Experiment 2: divide 64 physical cores among concurrent jobs

`compare_concurrency.py` tests the RabbitMA release at every configured job
count. All tasks in a batch start together, each processing a complete input.
CPU sets never overlap and sum to 64 physical cores:

| Concurrent jobs | Threads per job |
| ---: | --- |
| 2 | 32, 32 |
| 3 | 22, 21, 21 |
| 4 | 16, 16, 16, 16 |
| 5 | 13, 13, 13, 13, 12 |
| 6 | 11, 11, 11, 11, 10, 10 |

Each count has two placement controls. `numa` fills a single available NUMA
domain when possible and splits a job only when necessary. On two 32-core
domains, the 3- and 5-job configurations require one split job, while the
2-, 4-, and 6-job configurations can keep every job within one domain.
`spread` distributes each job across the selected domains, still using
disjoint CPU sets. This allows a placement comparison at the same job count.
The second repetition reverses configuration order.

By default the script makes six independently checksummed copies of the
same sample, with distinct files and independent assembly outputs. This
controls sample complexity and avoids all jobs sharing one input-file page
cache. It needs about 28 GB for six 4.6 GB inputs, in addition to assembly
workspaces. Copies and results are retained. Set `independent_inputs` to
false to reuse the same input pathname; report that shared-cache condition
when interpreting results. Set `placements` to `["numa"]` to omit the
additional spread controls.

## Timing and validation

With `warm_input: true`, each batch's input files are read before timing.
This is a warm-input experiment; input copying/warming and result comparison
are excluded from measured assembly time. No global cache dropping or
changes to unrelated processes are performed. Normal assembler temporary
cleanup is retained.

Each phase writes `manifest.json`, `results.json`, `timings.csv`, `report.md`,
`validation.json`, and per-job console/GNU-time logs. CPU affinity and summed
process-tree RSS are sampled every ten seconds in `samples.jsonl`. Sampled
aggregate RSS can miss short peaks and count shared pages more than once;
GNU time's per-process RSS maxima are reported separately. Affinity escapes
or assembler failures stop the experiment; only its own child processes are
terminated on interruption or failure.

The report calculates the release/original speedup, ranks measured concurrency
throughput, and compares NUMA and spread throughput at each job count.

Throughput is completed samples per node-hour. `node_hours_per_sample`
allows comparison at a constant node-hour price. These are finite batches,
not a measured multi-thousand-job queue; sample diversity, storage contention
and batch-tail idle time still affect production throughput. No configuration
is assumed to win before measurement.

Final contigs are compared outside the timed interval in three modes:
literal with metadata, circular/strand equivalence with metadata, and
circular/strand equivalence without k/coverage comparison. Flags and duplicate
counts remain part of every comparison. All outcomes, including input errors,
are saved; timing completion does not imply output equivalence. These scripts
do not compare every intermediate file. See the repository README for the
standalone comparator and its stricter intermediate-evidence use case.
