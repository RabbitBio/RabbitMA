// Read-only graph replay: exact transforms and reverse lookup against the
// uncached SDBG, including positions above UINT32_MAX on a retained graph.
#include "sdbg/sdbg.h"
#include "assembly/sdbg_pruning.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Probe {
  uint64_t edge, backward, forward, next, previous, reverse;
};

static uint64_t Mix(uint64_t x) {
  x ^= x >> 30u; x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27u; x *= UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31u);
}

static void Require(bool value, const char *message) {
  if (!value) { std::fprintf(stderr, "%s\n", message); std::exit(2); }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s graph-prefix [threads=8] [probes=100000]\n", argv[0]);
    return 2;
  }
  omp_set_num_threads(argc > 2 ? std::atoi(argv[2]) : 8);
  const size_t count = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 100000;
  SDBG graph;
  double start = omp_get_wtime();
  graph.LoadFromFile(argv[1]);
  std::printf("load_s=%.6f edges=%llu\n", omp_get_wtime() - start,
              (unsigned long long)graph.size());
  Require(graph.size() != 0 && count != 0, "empty graph/probe set");
  std::vector<Probe> probes(count);
  const uint64_t boundary = std::numeric_limits<uint32_t>::max();
  uint64_t wide_probes = 0;
  start = omp_get_wtime();
#pragma omp parallel for reduction(+ : wide_probes) schedule(static)
  for (size_t i = 0; i < count; ++i) {
    uint64_t edge = Mix(i) % graph.size();
    // Always include the address-width boundary and both ends of the graph.
    if (i < 1024 && graph.size() > boundary + 512) edge = boundary + i - 512;
    else if (i < 2048) edge = (i - 1024 + graph.size()) % graph.size();
    while (!graph.IsValidEdge(edge)) edge = (edge + 1) % graph.size();
    Probe &p = probes[i];
    p.edge = edge;
    p.backward = graph.Backward(edge);
    p.forward = graph.Forward(edge);
    p.next = graph.NextSimplePathEdge(edge);
    p.previous = graph.PrevSimplePathEdge(edge);
    p.reverse = graph.EdgeReverseComplement(edge);
    wide_probes += edge >= boundary;
  }
  std::printf("reference_s=%.6f probes=%zu wide_probes=%llu\n",
              omp_get_wtime() - start, count, (unsigned long long)wide_probes);

  setenv("MEGAHIT_EXPERIMENTAL_LARGE_GRAPH_CACHE", "1", 1);
  setenv("MEGAHIT_FORCE_LARGE_GRAPH_CACHE", "1", 1);
  const uint64_t budget = uint64_t{200} << 30u;
  start = omp_get_wtime();
  Require(graph.BuildTopologyCache(budget), "no topology/lookup cache built");
  std::printf("topology_build_s=%.6f bytes=%llu\n", omp_get_wtime() - start,
              (unsigned long long)graph.TopologyCacheBytes());
  auto verify = [&]() {
    uint64_t errors = 0;
#pragma omp parallel for reduction(+ : errors) schedule(static)
    for (size_t i = 0; i < probes.size(); ++i) {
      const auto &p = probes[i];
      errors += graph.Backward(p.edge) != p.backward;
      errors += graph.Forward(p.edge) != p.forward;
      errors += graph.CachedOrLiveNextSimplePathEdge(p.edge) != p.next;
      errors += graph.PrevSimplePathEdge(p.edge) != p.previous;
      errors += graph.EdgeReverseComplement(p.edge) != p.reverse;
      uint64_t previous = SDBG::kNullID;
      if (graph.TryCachedPrevSimplePathEdge(p.edge, &previous))
        errors += previous != p.previous;
    }
    Require(errors == 0, "cache differs from uncached graph");
  };
  start = omp_get_wtime();
  verify();
  std::printf("topology_verify_s=%.6f\n", omp_get_wtime() - start);
  start = omp_get_wtime();
  Require(graph.BuildSimplePathSnapshot(budget), "no shared snapshot built");
  Require(graph.SimplePathSnapshotUsesWideBases(), "wide-base path not exercised");
  std::printf("snapshot_build_s=%.6f bytes=%llu inverse_bytes=%llu overflow=%llu\n",
              omp_get_wtime() - start,
              (unsigned long long)graph.SimplePathSnapshotBytes(),
              (unsigned long long)graph.InverseSimplePathSnapshotBytes(),
              (unsigned long long)graph.SimplePathSnapshotOverflowCount());
  start = omp_get_wtime();
  verify();
  std::printf("snapshot_verify_s=%.6f\n", omp_get_wtime() - start);
  // Invalidate sparse targets and check that cached links fall back correctly.
  for (size_t i = 0; i < std::min<size_t>(probes.size(), 128); ++i) {
    if (probes[i].next != SDBG::kNullID) graph.SetInvalidEdge(probes[i].next);
  }
  uint64_t errors = 0;
#pragma omp parallel for reduction(+ : errors) schedule(static)
  for (size_t i = 0; i < probes.size(); ++i) {
    const uint64_t edge = probes[i].edge;
    if (graph.IsValidEdge(edge)) {
      errors += graph.CachedOrLiveNextSimplePathEdge(edge) != graph.NextSimplePathEdge(edge);
      uint64_t previous = SDBG::kNullID;
      if (graph.TryCachedPrevSimplePathEdge(edge, &previous))
        errors += previous != graph.PrevSimplePathEdge(edge);
    }
  }
  Require(errors == 0, "stale snapshot target not invalidated");
  // A graph reload must discard the directory and both snapshot directions.
  graph.LoadFromFile(argv[1]);
  Require(graph.TopologyCacheBytes() == 0 && !graph.HasSimplePathSnapshot() &&
          !graph.HasInverseSimplePathSnapshot(), "cache survived graph reload");
  Require(!graph.BuildTopologyCache(1) && !graph.BuildSimplePathSnapshot(1),
          "cache did not honor tiny memory budget");
  Require(graph.TopologyCacheBytes() == 0 && !graph.HasSimplePathSnapshot(),
          "failed optional allocation left a published cache");
  // Dense/word pruning used to retain the inverse array through all later
  // assembly phases. Force the non-directed path: the directed-frontier
  // trimmer has its own earlier release and would mask this lifecycle leak.
  unsetenv("MEGAHIT_DISABLE_INVERSE_SIMPLE_SNAPSHOT");
  setenv("MEGAHIT_DISABLE_TIP_FRONTIER", "1", 1);
  Require(graph.BuildTopologyCache(budget) && graph.BuildSimplePathSnapshot(budget),
          "no cache built for pruning lifecycle check");
  Require(graph.HasInverseSimplePathSnapshot(), "inverse cache not exercised");
  sdbg_pruning::RemoveTips(graph, 2);
  Require(!graph.HasInverseSimplePathSnapshot(), "inverse snapshot survived tip pruning");
  Require(graph.HasSimplePathSnapshot(), "reusable forward snapshot was discarded");
  std::puts("Exact SDBG cache, reverse-complement, invalidation, reload and pruning lifecycle checks passed");
}
