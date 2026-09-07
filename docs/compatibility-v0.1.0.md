# Output compatibility with RabbitMA v0.1.0

These fixes are included in the consolidated v0.2.0 source and binary package.
This publication replaces the earlier v0.2.0 and v0.2.1 packages. The public
`megahit` command and its options are unchanged.

## Reference and comparison

The reference for circular change-only evidence is the published v0.1.0 core's
single-worker block builder, including its 32 interleaved paths. v0.1.0's
multi-worker execution can choose different path owners and circular origins
on repeated runs of the same graph. Those intermediate differences can also
change later graphs and final contigs, so an arbitrary old parallel run is
not a unique reference result.

Comparisons preserve the exact sequence, including strand and circular origin,
k, flags, and coverage at the stored precision. Numeric contig IDs and record
ordering are excluded. Equal contig counts or total lengths alone are not a
successful comparison.

## Corrections

In v0.1.0, `k*.addi.fa` records are finite sequence evidence for the next k.
Rotating a circular additional contig, or extending its overlap to the next k,
changes the windows contributed by that record. The current source preserves
the legacy finite windows. Primary circular contigs still extend their
overlap as before. Both FASTA and packed intermediates use these rules.

Parallel unitig construction now restores the reference owner, stored strand,
and linear-merge priority from path lengths and endpoints. It simulates only
the reference builder's admission and completion events; the expensive graph
walks remain parallel. The metadata is reconstructed after initial tip removal
so both initial-pruning strategies use the same reference boundary. Subsequent
linear merges select their owner before either worker changes the path.

The compatibility tests include frozen v0.1.0 records and a fixture containing
64 circles with short initial tips and low-depth branches. They exercise
multiple worker counts and both initial-pruning strategies. Source CI and
portable-package CI run the suite against all three CPU variants.

## Validation scope

The pre-release source audit used exact record comparisons, including sequence
orientation and circular cut, rather than only assembly statistics:

| Comparison | Matching files |
| --- | ---: |
| CAMI3 subset of 2,000,000 read pairs, complete new 32-thread pipeline versus published v0.1.0 single-thread pipeline | 60 / 60 |
| Same subset, two concurrent new 32-thread jobs versus an existing v0.1.0 32-thread run | 120 / 120 |
| Same subset, the two new concurrent jobs versus each other | 60 / 60 |
| Full CAMI3 sample_0, new 32-thread versus 64-thread pipeline | 60 / 60 |
| Identical full-sample k109 graph, new repeated 32-thread, 64-thread and edge-first paths versus published v0.1.0 single-thread assembly | 16 / 16 |
| Twelve mixed-length synthetic graphs at k=21/49/99, new default and edge-first paths versus v0.1.0 single-thread assembly | 72 / 72 |

The subset reference produced 65,442 contigs and 62,340,941 bp. The full sample
produced 154,269 contigs and 221,229,065 bp with both new thread counts. The
full sample was not run through the entire old single-thread pipeline; its
direct old single-thread comparison covered the k109 graph that exposed the
unstable circular cuts. These results do not imply equality with every old
parallel execution.

All 95 regression tests passed on each of the three native CPU variants
(285 tests), including seven new circular-evidence tests per variant. Portable
release builds also run the full suite against every packaged core. The audit
ran concurrent verification workloads and is not a new controlled performance
benchmark.

## A baseline failure is not an output comparison

On the tested CAMI3 paired FASTQ input, the published v0.1.0 binary crashes in
`Read2SdbgS1::Initialize()` under `--presets meta-sensitive`. Its compact read
length index is incompatible with its old `BuildIndex()` implementation.
Later RabbitMA versions already contain the read-index fix. A successful run
of the current version does not establish equality with a baseline run that
produced no assembly.
