# Bundled rapidgzip

MEGAHIT vendors the header-only decompression core from rapidgzip v0.16.0 so
that parallel decompression of a single gzip input does not require a separate
runtime executable or a network connection while configuring the build.

- Upstream: https://github.com/mxmlnkn/rapidgzip
- rapidgzip tag: `rapidgzip-v0.16.0`
- rapidgzip commit: `d2350e9c9ba54398cd64e45bfc8c631beec017f0`
- librapidarchive commit: `1221a30bb548b305a69e5715f2bc348ba37ac243`
- Included upstream directories: `core`, `filereader`, `huffman`,
  `indexed_bzip2`, and `rapidgzip`

The upstream sources are available under either the MIT license or the Apache
License 2.0.  Both license texts are kept next to this file.  MEGAHIT compiles
the C++17-only templates behind a small C++11-compatible adapter; the rest of
MEGAHIT retains its original language level.

MEGAHIT assigns decoder threads from one process-wide `-t` budget.  For a
single or small number of large FASTQ inputs, decompressed bytes are cut only
at complete record boundaries, encoded by a bounded parallel parser/packer
pipeline, and committed in original batch order.  The raw in-flight queue has
a fixed byte cap, and wrapped FASTQ records are supported.  FASTA and unusual
FASTQ layouts automatically replay through the existing kseq path from a
private part file, preserving compatibility and output semantics.
