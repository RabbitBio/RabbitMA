#ifndef MEGAHIT_STARTUP_AFFINITY_H
#define MEGAHIT_STARTUP_AFFINITY_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Runtime NUMA information restricted to the CPUs and memory nodes available
// to this process (taskset/cpuset/cgroup restrictions are honoured).  Node IDs
// are deliberately kept separate from dense domain indices: Linux node IDs
// need not be contiguous, and callers should never bake a socket count or a
// particular CPU numbering scheme into their algorithms.
struct NumaTopology {
  std::vector<int> node_ids;
  std::vector<unsigned> cpu_counts;
  std::vector<unsigned> physical_core_counts;
  std::vector<uint64_t> last_level_cache_bytes;
  std::vector<int> cpu_to_domain;
  // Sum of unique last-level cache instances intersecting the process CPU
  // affinity.  This is kept separately from the per-NUMA fallback above:
  // neither "one LLC per socket" nor "one LLC per NUMA node" is portable.
  uint64_t unique_last_level_cache_bytes{0};

  size_t domain_count() const { return node_ids.empty() ? 1 : node_ids.size(); }
  uint64_t total_last_level_cache_bytes() const;
  unsigned total_physical_core_count() const;
};

// Resource sharing is an explicit scheduling fact, not something a process
// can infer reliably from instantaneous load.  The launcher publishes these
// values for every core subprocess.  Direct megahit_core users can set the
// same MEGAHIT_JOBS_PER_NODE and MEGAHIT_MEMORY_BUDGET_PER_JOB environment
// variables.  With the defaults below every policy decision is identical to
// the exclusive path.
struct RuntimeResourcePolicy {
  unsigned jobs_per_node{1};
  uint64_t memory_budget_per_job{0};

  bool shared_node() const { return jobs_per_node > 1u; }
  uint64_t last_level_cache_budget_bytes() const;
};

/**
 * OpenMP thread binding (OMP_PROC_BIND/OMP_PLACES) narrows the affinity mask
 * of every OpenMP thread, including the master.  A std::thread spawned from a
 * bound thread inherits that narrow mask and may end up serialized on a
 * single core.  The process's original mask is captured by a global
 * constructor before OpenMP initializes; helper threads call this to widen
 * themselves back to it.
 */
void ResetThreadAffinityToStartupMask();

// If the launcher supplied MEGAHIT_NUMA_NODE, install a strict MPOL_BIND
// policy for all future private allocations in this core process. CPU
// affinity is narrowed by the launcher and inherited across exec. Returning
// false is fatal for explicit NUMA mode: silently continuing would violate
// the user's resource-isolation request. With no requested node this is a
// no-op that returns true.
bool ConfigureProcessNumaMemoryPolicy(std::string *error_message);

// The returned object is discovered once and remains valid for the process
// lifetime.  On non-Linux systems, or when topology files are unavailable, it
// describes one portable fallback domain.
const NumaTopology &GetNumaTopology();

// Parsed once from the launcher-provided resource contract.  Invalid or
// absent values conservatively select the existing exclusive policy.
const RuntimeResourcePolicy &GetRuntimeResourcePolicy();

// Dense domain index of the calling worker.  Returns zero for the portable
// fallback domain or when the current CPU cannot be mapped.
unsigned CurrentNumaDomain();

// Apply a stable interleaved NUMA policy to the page-aligned interior of a
// shared allocation.  This is intended for large arrays that are written and
// subsequently consumed by different worker partitions; thread-local arrays
// should continue to rely on first touch.  Returns true only when a policy was
// installed across at least two allowed memory nodes.
bool InterleaveMemoryPages(void *address, size_t bytes);

// Install a future-allocation policy for one runtime-discovered NUMA domain.
// Callers normally follow this with DiscardMemoryPages before first touch so
// pages recycled by malloc cannot retain an old placement.  No data is moved.
bool BindMemoryPagesToNumaDomain(void *address, size_t bytes,
                                 unsigned domain);

// Ask the kernel to back a large, long-lived random-access region with
// transparent huge pages when available.  This is only advice and therefore
// preserves the portable fallback and exact program semantics.
bool AdviseHugePages(void *address, size_t bytes);

// Drop the page-aligned interior of an anonymous allocation.  A later write
// receives a fresh zero page under the allocation's current NUMA policy.
bool DiscardMemoryPages(void *address, size_t bytes);

#endif  // MEGAHIT_STARTUP_AFFINITY_H
