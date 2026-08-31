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
- Runtime CPU and NUMA discovery without fixed socket, core-count, or
  CAMI-specific thresholds.

## Installation

### Prebuilt Linux package (recommended)

The quickest way to use RabbitMA is to download the prebuilt
`RabbitMA-v0.1.0-linux-x86_64.tar.gz` package from the
[v0.1.0 release](https://github.com/RabbitBio/RabbitMA/releases/tag/v0.1.0).
It requires Linux x86_64 with glibc 2.17 or newer, Python 3.6 or newer, gzip,
and bzip2;
CMake and a compiler are not needed. The package exposes `megahit` as its only
public command and includes the internal CPU core variants, test data, and
required non-glibc runtime libraries.

```bash
wget https://github.com/RabbitBio/RabbitMA/releases/download/v0.1.0/RabbitMA-v0.1.0-linux-x86_64.tar.gz
tar -xzf RabbitMA-v0.1.0-linux-x86_64.tar.gz
cd RabbitMA-v0.1.0-linux-x86_64
./megahit --test -t 4
```

The package includes BMI2/POPCNT, POPCNT-only, and portable core binaries. The
Python driver selects a supported variant at run time.

The v0.1.0 prebuilt launcher predates automatic cgroup memory detection. When
running that package inside Docker, Kubernetes, or a scheduler job with a hard
memory limit, pass the job limit explicitly with `-m`. The current source tree
detects cgroup v1/v2 limits automatically; this fix will be included in the
next binary release.

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

## Compatibility and validation

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
