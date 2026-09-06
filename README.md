# RabbitMA

> **Build and CLI compatibility:** Building RabbitMA from source requires a
> **C++17-capable C++ compiler** with OpenMP support. The public command is
> **`megahit`**, with the same command-line syntax, options, and output layout
> as official MEGAHIT v1.2.9, so existing MEGAHIT commands run unchanged.
> CMake itself may be as old as 2.8.12. Users of the prebuilt package do not
> need CMake or a C++ compiler.

RabbitMA is a high-performance, memory-aware metagenome assembler derived
from [MEGAHIT v1.2.9](https://github.com/voutcn/megahit). It retains the
MEGAHIT assembly model, command-line options, and output layout while reducing
the repeated data movement, graph traversal, synchronization, and serial I/O
costs that dominate large multi-sample and single-file workloads.

The supported command is `megahit`, so existing MEGAHIT pipelines do not need
to change their command lines.

## Highlights

- Parallel multi-library input and adaptive intra-file gzip decompression. The
  rapidgzip implementation is bundled; no run-time download is required.
- One-pass segmented exact counting with compact records and streamed edge
  output instead of repeated full-read scans.
- Streamed, bucket-packed succinct de Bruijn graph construction with bounded
  working sets.
- Cached graph topology, parallel unitig traversal, and compact unitig
  adjacency and endpoint lookup.
- Sharded local-assembly mapping and iterative flank indexes.
- A reusable exact read-occurrence index across k iterations, with the first
  iterative join fused into index construction and full-scan fallback when
  the index cannot be used.
- Packed intermediate contigs and reusable read-position information to
  reduce repeated decoding and temporary I/O.
- Runtime CPU and NUMA discovery without fixed socket, core-count, or
  CAMI-specific thresholds.
- Automatic NUMA placement for concurrent jobs from the same Linux user,
  with physical-core reservations and local-preferred memory.

## Installation

### Prebuilt Linux package (recommended)

The quickest way to use RabbitMA is to download the prebuilt
`RabbitMA-v0.2.1-linux-x86_64.tar.gz` package from the
[v0.2.1 release](https://github.com/RabbitBio/RabbitMA/releases/tag/v0.2.1).
It requires Linux x86_64 with glibc 2.17 or newer, Python 3.6 or newer, gzip,
and bzip2;
CMake and a compiler are not needed. The package exposes `megahit` as its only
public command and includes the internal CPU core variants, test data, and
required non-glibc runtime libraries.

```bash
wget https://github.com/RabbitBio/RabbitMA/releases/download/v0.2.1/RabbitMA-v0.2.1-linux-x86_64.tar.gz
tar -xzf RabbitMA-v0.2.1-linux-x86_64.tar.gz
cd RabbitMA-v0.2.1-linux-x86_64
./megahit --test -t 4
```

The package includes BMI2/POPCNT, POPCNT-only, and portable core binaries. The
Python driver selects a supported variant at run time.

The binary package retains the CentOS 7 / glibc 2.17 baseline and includes
`libnuma` for automatic NUMA placement. There is no need to upgrade glibc or
set `LD_LIBRARY_PATH`. Python 3.6 or newer must be available as `python3`;
CentOS 7's default Python 2 alone is insufficient. The v0.2.1 launcher also
detects cgroup v1/v2 memory limits automatically, including inside containers
and scheduler jobs.

### Build from source

Building RabbitMA requires a C++17-capable compiler with OpenMP support, CMake
2.8.12 or newer, zlib, Python 3.6 or newer, gzip, and bzip2. C++17 enables the
bundled parallel gzip reader used by the advertised high-performance
configuration. If libdeflate is installed, RabbitMA detects and uses it
automatically.

```bash
git clone https://github.com/RabbitBio/RabbitMA.git
cd RabbitMA
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -- -j4
./megahit --test -t 4
```

Optional installation:

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install .
cmake --build . --target install
```

## Usage

The commands below assume the current directory is the extracted prebuilt
package. A source build provides the same command as `build/megahit`.

One paired-end library:

```bash
./megahit -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 32 -o assembly
```

Multiple paired-end libraries:

```bash
./megahit \
  -1 sample1_R1.fastq.gz,sample2_R1.fastq.gz \
  -2 sample1_R2.fastq.gz,sample2_R2.fastq.gz \
  -t 64 --k-min 39 -o coassembly
```

Interleaved and single-end inputs remain compatible with MEGAHIT:

```bash
./megahit --12 interleaved.fastq.gz -o assembly
./megahit -r reads.fastq.gz -o assembly
```

Run `./megahit --help` for the complete option list. Final contigs are written
to `OUT_DIR/final.contigs.fa`.

### Automatic NUMA placement

RabbitMA v0.2.1 automatically coordinates concurrent jobs from the same
Linux user when each job specifies `-t` and fits within one available NUMA
domain. No extra `numactl` command or NUMA option is needed. This behavior is
included in both the source and the v0.2.1 binary package. Upgrade older
binary packages to use it.

For example, on a machine with **72 available physical cores, split into two
NUMA domains of 36 cores each**, launch two jobs with separate output paths:

```bash
./megahit --12 sample_A.fq.gz -t 36 -o sample_A.out &
./megahit --12 sample_B.fq.gz -t 36 -o sample_B.out &
wait
```

With both domains available, one job reserves the first domain and the other
reserves the second. Both may also read the same input file. Keep
`--jobs-per-node` at its default of 1 for this mode: each job owns a separate
NUMA resource domain. The log records `Automatic NUMA placement: node ...`
with the reserved physical-core count and CPU IDs.

Reservations use host-local locks shared by participating launchers under the
same user account and `/dev/shm` mount. Simultaneous starts are coordinated,
and child processes retain the reservation until they exit. The launcher
respects inherited CPU and memory allowances and keeps SMT siblings in the
same reservation. **72 logical CPUs are not necessarily 72 physical cores.**
If the scheduler grants only 18 physical cores on each domain, a 36-thread
job cannot be moved onto 36 cores of one domain.

Automatic placement uses local-preferred memory, allowing remote allocation
when local memory is exhausted. Fractional `-m` budgets are based on the
reserved physical-core share of the NUMA domain's capacity, bounded by the
cgroup memory limit. This is a placement preference and budget, not a hard
reservation of RAM or isolation from unrelated programs.

If no single-domain reservation fits, `libnuma` or memory-policy permission
is unavailable, or explicit CPU/memory settings take precedence, the job keeps
its inherited placement and reports the reason when available. An omitted
`-t` retains the usual all-visible-CPU default. `--jobs-per-node` greater than
1 selects the explicit sharing policy below and disables automatic placement.
Use `--numa-node off` (or `MEGAHIT_DISABLE_AUTO_NUMA=1`) to disable the default
coordination. Explicit CPU IDs in `OMP_PLACES`, `GOMP_CPU_AFFINITY`,
`KMP_AFFINITY`, `OMP_PROC_BIND=false`, and inherited non-default memory policies
are respected.

### Explicit NUMA-local jobs

To select a specific NUMA domain and require strict local-memory allocation:

```bash
./megahit --12 sample.fq.gz --numa-node 0 -t 32 -o sample.out
```

`--numa-node` intersects that node with the CPU and memory masks granted by the
scheduler, binds the launcher and every worker to the resulting CPUs, and
installs a strict local-memory policy before each core process allocates its
working data. This explicit Linux mode requires the `libnuma` runtime,
which is included in the prebuilt package and must be installed separately
when building from source. If `-t` is omitted, RabbitMA uses the discovered
physical-core count rather than silently enabling SMT. A fractional `-m` is
computed from the selected NUMA node's capacity and the cgroup limit, not from
the whole machine. An unavailable or disallowed node is an error.

Explicit `--numa-node N` bypasses automatic reservations, so assign different
node IDs to different jobs or let the scheduler separate them. Keep
`--jobs-per-node 1` when each job owns its selected NUMA domain. Strict binding
means the job can run out of memory even while another NUMA node has free RAM;
use it when the job's peak working set fits in the selected domain.

When a scheduler has already restricted each task's CPU mask to exactly one
NUMA domain, `--numa-node auto` discovers that sole domain and applies strict
binding. This explicit option fails if the task can still see zero or multiple
domains. Omitting `--numa-node` uses the automatic coordination described above.

### Several jobs on one node

When several RabbitMA processes really share the same last-level caches and
node memory allocation, tell every process how many jobs are co-located:

```bash
./megahit --12 sample.fq.gz -t 16 --jobs-per-node 4 -o sample.out
```

This does not silently change `-t`. It gives each job an equal share of the
`-m` memory budget and lets bounded-workset kernels size their active data from
the CPU affinity, discovered cache topology, operating-system page size, and
the explicit sharing count. For unequal inputs, an exact byte cap can override
the equal memory split with `--memory-budget-per-job BYTES`.

Leave `--jobs-per-node` at 1 when the scheduler already gives each process an
isolated CPU/cache allocation and a per-job memory cgroup. A cpuset alone is
not memory isolation: jobs pinned to disjoint cores can still share LLC or
DRAM bandwidth. For paid clusters, choose `jobs × threads` by measuring
completed samples per node-hour; the fastest isolated `-t` is not necessarily
the most cost-efficient concurrent configuration.

## Compatibility and validation

The current source includes the computational optimizations validated in the
full CAMI3 development runs, together with automatic NUMA placement. Default
adaptive optimizations use the existing command-line interface. Additional
`MEGAHIT_EXPERIMENTAL_*` paths remain opt-in; their enabled benchmark timings
are not a claim about default performance on every input or machine. The
per-NUMA SDBG replica prototype that failed to show a useful gain is excluded.

For source regression tests, use Python 3.8 or newer after building:

```bash
python3 -B -m unittest discover -s tests -p 'test_*.py' -v
cmake --build build --target check_stream_stores check_local_minimizer_gate
```

`MEGAHIT_TEST_CORE=/path/to/megahit_core_popcnt` (or the portable variant)
selects another core for the same regression suite. These developer tests
cover complete sequence orientation and metadata, exact counting, compressed
input boundaries, local mapping, index build/replay, and graph reuse. The
application launcher itself continues to support Python 3.6 or newer.

RabbitMA is intended to preserve MEGAHIT v1.2.9 assembly semantics. Structured
simulations cover unique sequence, strain bubbles, repeats, uneven/error-prone
coverage, circular sequence, single-end tips, ambiguous bases, variable read
lengths, paired-end, interleaved, gzip, and bzip2 inputs. The v0.1.0
compatibility audit also compared four real single-end datasets and exercised
1-, 2-, 3-, and high-thread local assembly. All non-circular final contigs
matched the official normalized sequence, flag, and multiplicity multisets;
circular contigs matched the same complete graph-edge sets.

Raw FASTA MD5 values are not a reliable semantic comparison for parallel
assemblers: record order, reverse-complement orientation, and the chosen origin
of a circular contig can change without changing the assembled graph. Official
MEGAHIT can select different circular origins at different thread counts; its
seed-weighted circular depth can consequently differ slightly as well.
Validation should compare normalized sequence, flag, and multiplicity
multisets for non-circular contigs, and exact graph-edge sets for circular
contigs.

See [BENCHMARKS.md](BENCHMARKS.md) for the CAMI III reference measurement and
its scope.

## Attribution

RabbitMA is a derivative work, not an official MEGAHIT release. Please retain
the original MEGAHIT attribution and cite the MEGAHIT papers when using this
software:

- Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an ultra-fast
  single-node solution for large and complex metagenomics assembly via
  succinct de Bruijn graph. *Bioinformatics* (2015).
- Li D, Luo R, Liu C-M, et al. MEGAHIT v1.0: A fast and scalable metagenome
  assembler driven by advanced methodologies and community practices.
  *Methods* (2016).

Additional provenance is recorded in [NOTICE](NOTICE). RabbitMA is distributed
under the GNU General Public License v3 or later; see [LICENSE](LICENSE).
