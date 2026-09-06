#ifndef MEGAHIT_ASSEMBLY_CONTIG_BRIDGES_H_
#define MEGAHIT_ASSEMBLY_CONTIG_BRIDGES_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include "sequence/io/contig/contig_reader.h"
#include "sdbg/sdbg.h"

// Exact sequence spans omitted from the SDBG. The two retained stubs have
// unique interior endpoints; only those endpoint pairs acquire a new link.
class ContigBridges {
 public:
  struct Link { uint64_t target; uint32_t sequence; bool reverse; };
  struct Minimum {
    uint64_t physical{SDBG::kNullID};
    const Link *gap{nullptr};
  };

  void Load(const std::string &prefix, SDBG &graph) {
    const std::string path = prefix + ".bridges.fa";
    struct stat st {};
    if (::stat(packed_contig::Path(path).c_str(), &st) != 0) return;
    std::ifstream metadata(prefix + ".bridges.meta");
    unsigned version = 0, k = 0;
    uint64_t count = 0;
    metadata >> version >> k >> margin_ >> count;
    if (!metadata || version != 1u || k != graph.k()) xfatal("Invalid bridge metadata\n");
    k_ = k;
    ContigReader reader(path);
    reader.ReadAllWithMultiplicity(&sequences_, &multi_, false);
    if (sequences_.seq_count() != count || count >= std::numeric_limits<uint32_t>::max())
      xfatal("Invalid bridge sequence count\n");
    minima_.resize(count);
    std::vector<std::array<uint64_t, 4>> endpoints(count);
#pragma omp parallel for schedule(dynamic, 16)
    for (int64_t i = 0; i < static_cast<int64_t>(count); ++i) {
      const auto sequence = sequences_.GetSeqView(i);
      if (uint64_t(sequence.length()) <= uint64_t{2} * margin_ + k_ + 2u)
        xfatal("Invalid bridge span\n");
      const uint64_t from = FindEdge(graph, sequence, margin_);
      const uint64_t to = FindEdge(graph, sequence, sequence.length() - margin_ - k_ - 1u);
      const uint64_t reverse_from = graph.EdgeReverseComplement(to);
      const uint64_t reverse_to = graph.EdgeReverseComplement(from);
      if (from == SDBG::kNullID || to == SDBG::kNullID ||
          reverse_from == SDBG::kNullID || reverse_to == SDBG::kNullID)
        xfatal("Missing bridge endpoint\n");
      endpoints[i] = {{from, to, reverse_from, reverse_to}};

    }
    next_.reserve(count * 2u);
    previous_.reserve(count * 2u);
    for (uint32_t i = 0; i < count; ++i) {
      for (unsigned strand = 0; strand < 2; ++strand) {
        const uint64_t from = endpoints[i][2u * strand], to = endpoints[i][2u * strand + 1u];
        uint64_t outgoing[4];
        if (graph.OutgoingEdges(from, outgoing) != 0 ||
            graph.PrevSimplePathEdge(to) != SDBG::kNullID ||
            !next_.emplace(from, Link{to, i, strand != 0}).second ||
            !previous_.emplace(to, from).second) xfatal("Conflicting bridge endpoints\n");
      }
    }
  }

  bool empty() const { return next_.empty(); }
  uint64_t size() const { return sequences_.seq_count(); }
  unsigned margin() const { return margin_; }
  const Link *Next(uint64_t edge) const {
    const auto found = next_.find(edge);
    return found == next_.end() ? nullptr : &found->second;
  }
  uint64_t Previous(uint64_t edge) const {
    const auto found = previous_.find(edge);
    return found == previous_.end() ? SDBG::kNullID : found->second;
  }
  uint32_t Length(const Link &link) const {
    return sequences_.GetSeqView(link.sequence).length() - 2u * margin_ - k_ - 2u;
  }
  uint64_t Depth(const Link &link) const { return uint64_t(Length(link)) * multi_[link.sequence]; }
  void AddHistogram(std::vector<uint64_t> *counts) const {
    counts->resize(kMaxMul + 1u, 0);
    for (uint32_t i = 0; i < sequences_.seq_count(); ++i)
      (*counts)[multi_[i]] += 2u * uint64_t(sequences_.GetSeqView(i).length() - 2u * margin_ - k_ - 2u);
  }
  void WriteGap(const Link &link, char *output) const {
    const auto sequence = sequences_.GetSeqView(link.sequence);
    for (uint32_t i = 0; i < Length(link); ++i)
      output[i] = "ACGT"[Base(sequence, margin_ + k_ + 1u + i, link.reverse)];
  }

  void IncludeGapMinimum(Minimum *minimum, const Link *gap) const {
    if (!minimum->gap || CompareSequenceEdges(
        gap->sequence, MinimumPosition(*gap), gap->reverse,
        minimum->gap->sequence, MinimumPosition(*minimum->gap),
        minimum->gap->reverse) < 0) minimum->gap = gap;
  }
  Minimum NormalizeMinimum(Minimum minimum, const SDBG &graph) const {
    if (minimum.gap && minimum.physical != SDBG::kNullID) {
      Minimum physical; physical.physical = minimum.physical;
      Minimum gap; gap.gap = minimum.gap;
      if (!MinimumLess(gap, physical, graph)) minimum.gap = nullptr;
    }
    return minimum;
  }
  bool MinimumLess(const Minimum &a, const Minimum &b, const SDBG &graph) const {
    if (!a.gap && !b.gap) return a.physical < b.physical;
    std::vector<uint8_t> left(k_ + 1u), right(k_ + 1u);
    MinimumLabel(a, graph, left.data()); MinimumLabel(b, graph, right.data());
    for (unsigned i = k_; i > 0; --i) if (left[i-1] != right[i-1]) return left[i-1] < right[i-1];
    return left[k_] < right[k_];
  }
  uint64_t MinimumCoverage(const Minimum &minimum, const SDBG &graph) const {
    return minimum.gap ? multi_[minimum.gap->sequence] : graph.EdgeMultiplicity(minimum.physical);
  }
  // The original raw-loop constructor starts at the successor of the
  // globally smallest edge and recognizes a palindrome when that successor
  // is also the reverse complement of the seed edge.
  std::string PalindromeOrigin(const Minimum &minimum, const SDBG &graph) const {
    std::vector<uint8_t> seed(k_ + 1u), next(k_ + 1u);
    MinimumLabel(minimum, graph, seed.data());
    if (minimum.gap) {
      SequenceLabel(*minimum.gap, MinimumPosition(*minimum.gap) + 1u, next.data());
    } else if (const Link *gap = Next(minimum.physical)) {
      SequenceLabel(*gap, margin_ + 1u, next.data());
    } else {
      const uint64_t edge = graph.NextSimplePathEdge(minimum.physical);
      if (edge == SDBG::kNullID) xfatal("Broken bridge loop\n");
      PhysicalLabel(graph, edge, next.data());
    }
    for (unsigned i = 0; i <= k_; ++i) if (next[i] != 5u - seed[k_ - i]) return {};
    std::string origin(k_ + 1u, 'A');
    for (unsigned i = 0; i <= k_; ++i) origin[i] = "ACGT"[next[i] - 1u];
    return origin;
  }

 private:
  uint32_t MinimumPosition(const Link &link) const {
    auto *slot = &minima_[link.sequence][link.reverse];
    const uint32_t cached = __atomic_load_n(slot, __ATOMIC_RELAXED);
    if (cached) return cached;
    const auto sequence = sequences_.GetSeqView(link.sequence);
    uint32_t best = margin_ + 1u;
    const uint32_t end = sequence.length() - margin_ - k_ - 1u;
    for (uint32_t position = best + 1u; position < end; ++position) {
      int difference = 0;
      for (unsigned i = k_; i > 0 && !difference; --i)
        difference = int(Base(sequence, position + i - 1u, link.reverse)) -
                     int(Base(sequence, best + i - 1u, link.reverse));
      if (!difference) difference = int(Base(sequence, position + k_, link.reverse)) -
                                   int(Base(sequence, best + k_, link.reverse));
      if (difference < 0) best = position;
    }
    __atomic_store_n(slot, best, __ATOMIC_RELAXED);
    return best;
  }
  static unsigned Base(const SeqPackage::SeqView &sequence, uint32_t position, bool reverse) {
    return reverse ? sequence.base_at(sequence.length() - 1u - position) ^ 3u : sequence.base_at(position);
  }
  int CompareSequenceEdges(uint32_t a, uint32_t pa, bool ra,
                           uint32_t b, uint32_t pb, bool rb) const {
    const auto left = sequences_.GetSeqView(a), right = sequences_.GetSeqView(b);
    for (unsigned i = k_; i > 0; --i) {
      const int difference = int(Base(left, pa + i - 1u, ra)) - int(Base(right, pb + i - 1u, rb));
      if (difference) return difference;
    }
    return int(Base(left, pa + k_, ra)) - int(Base(right, pb + k_, rb));
  }
  void SequenceLabel(const Link &link, uint32_t position, uint8_t *label) const {
    const auto sequence = sequences_.GetSeqView(link.sequence);
    for (unsigned i = 0; i <= k_; ++i) label[i] = Base(sequence, position + i, link.reverse) + 1u;
  }
  static void PhysicalLabel(const SDBG &graph, uint64_t edge, uint8_t *label) {
    graph.GetLabel(edge, label);
    const unsigned w = graph.GetW(edge);
    label[graph.k()] = w > 4u ? w - 4u : w;
  }
  void MinimumLabel(const Minimum &minimum, const SDBG &graph, uint8_t *label) const {
    if (minimum.gap) SequenceLabel(*minimum.gap,
        MinimumPosition(*minimum.gap), label);
    else PhysicalLabel(graph, minimum.physical, label);
  }
  static uint64_t FindEdge(const SDBG &graph, const SeqPackage::SeqView &sequence, unsigned position) {
    std::vector<uint8_t> label(graph.k() + 1u);
    for (unsigned i = 0; i <= graph.k(); ++i) label[i] = sequence.base_at(position + i) + 1u;
    uint64_t edge = graph.IndexBinarySearch(label.data());
    if (edge == SDBG::kNullID) xfatal("Missing bridge node\n");
    do {
      const unsigned w = graph.GetW(edge);
      if ((w == label[graph.k()] || w == label[graph.k()] + 4u) && graph.IsValidEdge(edge)) return edge;
      --edge;
    } while (edge != SDBG::kNullID && !graph.IsLastOrTip(edge));
    xfatal("Missing bridge edge\n");
    return SDBG::kNullID;
  }

  SeqPackage sequences_;
  std::vector<mul_t> multi_;
  mutable std::vector<std::array<uint32_t, 2>> minima_;
  std::unordered_map<uint64_t, Link> next_;
  std::unordered_map<uint64_t, uint64_t> previous_;
  unsigned k_{0}, margin_{0};
};
#endif
