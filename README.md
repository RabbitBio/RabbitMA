# RabbitMA

> **Build requirement:** RabbitMA's supported high-performance build requires
> a **C++17-capable C++ compiler** with OpenMP support. CMake itself may be as
> old as 2.8.12; the CMake version does not remove the C++17 compiler
> requirement.

RabbitMA is a high-performance, memory-aware metagenome assembler derived
from [MEGAHIT v1.2.9](https://github.com/voutcn/megahit). It retains the
MEGAHIT assembly model, command-line options, and output layout while reducing
the repeated data movement, graph traversal, synchronization, and serial I/O
costs that dominate large multi-sample and single-file workloads.

The primary command is `rabbitma`. A compatible `megahit` command is built and
installed as well, so existing pipelines do not need to change immediately.

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

## Build

RabbitMA requires a C++17-capable compiler with OpenMP support, CMake 2.8.12 or
newer, zlib, Python 3, gzip, and bzip2. C++17 enables the bundled parallel gzip
reader used by the advertised high-performance configuration. If libdeflate is
installed, RabbitMA detects and uses it automatically.

```bash
git clone https://github.com/RabbitBio/RabbitMA.git
cd RabbitMA
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -- -j
./rabbitma --test -t 4
```

Optional installation:

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install .
cmake --build . --target install
```

### Prebuilt Linux package

Users who do not have a suitable CMake/compiler toolchain can download
`RabbitMA-v0.1.0-linux-x86_64.tar.gz` from the
[v0.1.0 release](https://github.com/RabbitBio/RabbitMA/releases/tag/v0.1.0).
It requires Linux x86_64 with glibc 2.28 or newer, Python 3, gzip, and bzip2;
CMake and a compiler are not needed. The package contains both `rabbitma` and
the compatible `megahit` command, all three CPU core variants, test data, and
the required non-glibc runtime libraries.

```bash
sha256sum -c RabbitMA-v0.1.0-linux-x86_64.tar.gz.sha256
tar -xzf RabbitMA-v0.1.0-linux-x86_64.tar.gz
export PATH="$PWD/RabbitMA-v0.1.0-linux-x86_64/bin:$PATH"
rabbitma --test -t 4
```

The build produces BMI2/POPCNT, POPCNT-only, and portable core binaries. The
Python driver selects a supported variant at run time.

## Usage

One paired-end library:

```bash
rabbitma -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 32 -o assembly
```

Multiple paired-end libraries:

```bash
rabbitma \
  -1 sample1_R1.fastq.gz,sample2_R1.fastq.gz \
  -2 sample1_R2.fastq.gz,sample2_R2.fastq.gz \
  -t 64 --k-min 39 -o coassembly
```

Interleaved and single-end inputs remain compatible with MEGAHIT:

```bash
rabbitma --12 interleaved.fastq.gz -o assembly
rabbitma -r reads.fastq.gz -o assembly
```

Run `rabbitma --help` for the complete option list. Final contigs are written
to `OUT_DIR/final.contigs.fa`.

## Compatibility and validation

RabbitMA is intended to preserve MEGAHIT v1.2.9 assembly semantics. Six
simulated datasets covering unique sequence, strain bubbles, repeats,
uneven/error-prone coverage, circular sequence, and single-end tips were
compared with official MEGAHIT. Across 100 final and per-k checks, 95 artifacts
were byte/canonical-sequence exact and the remaining five differed only by the
rotation origin of circular contigs; no substantive sequence or multiplicity
difference was found.

Raw FASTA MD5 values are not a reliable semantic comparison for parallel
assemblers: record order, reverse-complement orientation, and the chosen origin
of a circular contig can change without changing the assembled sequence.
Validation should compare normalized sequence and multiplicity multisets, with
circular sequences compared modulo rotation.

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
