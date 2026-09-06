/**
 * @file hash_graph.h
 * @brief HashGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_HASH_GRAPH_H_

#define __GRAPH_HASH_GRAPH_H_

#include <algorithm>
#include <array>
#include <deque>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "bit_operation.h"
#include "idba/contig_info.h"
#include "idba/hash_graph_vertex.h"
#include "idba/contig_graph_vertex.h"
#include "idba/hash_graph_table.h"
#include "idba/kmer.h"
#include "idba/sequence.h"
#include "utils/histgram.h"

class IdbaKmer;
class Sequence;

struct HashGraphPackedRead {
  const uint64_t *forward_words;
  uint32_t length;
};

struct HashGraphBulkStats {
  double generate_seconds{0};
  double key_sort_seconds{0};
  double reduce_seconds{0};
  double ordinal_sort_seconds{0};
  double replay_seconds{0};
  uint64_t occurrences{0};
  uint64_t aggregates{0};
};

/**
 * @brief It is a hash table based de Bruijn graph implementation.
 */
class HashGraph {
 public:
  typedef HashGraphVertexTable vertex_table_type;

  explicit HashGraph(uint32_t kmer_size = 0) { set_kmer_size(kmer_size); }
  ~HashGraph() {}

  HashGraphVertexAdaptor FindVertexAdaptor(const IdbaKmer &kmer) {
    IdbaKmer key = kmer.unique_format();
    auto p = vertex_table_.find(key);
    return ((p != vertex_table_.end())
                ? HashGraphVertexAdaptor(&*p, kmer != key)
                : HashGraphVertexAdaptor(NULL));
  }

  int64_t InsertKmers(const Sequence &seq);
  int64_t InsertPackedKmers(const uint64_t *forward_words,
                            uint32_t sequence_length);
  int64_t InsertPackedKmersIndexed(const uint64_t *forward_words,
                                   uint32_t sequence_length,
                                   uint32_t *forward_groups,
                                   uint32_t *reverse_groups);
  int64_t InsertPackedUncountKmers(const uint64_t *forward_words,
                                   uint32_t sequence_length);
  int64_t InsertPackedKmersBulk(const HashGraphPackedRead *reads,
                                size_t num_reads,
                                HashGraphBulkStats *stats = nullptr);
  int64_t InsertUncountKmers(const Sequence &seq);


  // Internal entry point used by the exact sort/reduce builder.  Calls arrive
  // in historical first-occurrence order, so the compact table recreates the
  // same bucket chains while doing only one formal insertion per unique key.
  void InsertAggregate(const IdbaKmer &key, uint32_t count,
                       uint8_t in_edges, uint8_t out_edges,
                       uint32_t *vertex_index);

  // Add another occurrence to a vertex whose exact oriented string was
  // already resolved by a structural builder.  The cached index is stable
  // across rehashes, so this avoids rematerializing and reprobeing the DNA key
  // while preserving the historical pending-rehash boundary of an ordinary
  // occurrence.
  void AccumulateResolvedOccurrence(uint32_t vertex_index,
                                    uint8_t in_edges,
                                    uint8_t out_edges);

  // Replay an already observed oriented transition and count the destination
  // occurrence. Returns false when this base has not yet been resolved; the
  // caller then performs its exact cold-path lookup and ObserveTransition().
  bool ReplayResolvedTransition(uint32_t from_code, uint8_t base,
                                uint32_t *to_code);

  void SetAggregateNeighbor(uint32_t vertex_index, bool strand,
                            uint32_t neighbor_index, bool neighbor_strand);

  // Replays an observed read transition using compact oriented vertex IDs.
  // This is used by exact aggregate builders after counts/edges have already
  // been reduced; it changes only the same unique-neighbor cache maintained
  // by occurrence-wise insertion.
  void ObserveTransition(uint32_t from_code, uint32_t to_code);

  // Aggregate builders replay only unique transitions.  Normalize the cache
  // slots whose final oriented degree is non-unique to the same sentinel the
  // occurrence-wise insertion path installs when its second edge appears.
  void FinalizeAggregateTransitions();

  void FinishBulkOccurrences(bool had_occurrence_after_last_unique) {
    vertex_table_.finish_bulk_occurrences(
        had_occurrence_after_last_unique);
  }

  // Validation hook for structural builders. It compares only semantic graph
  // state, deliberately ignoring table layout and transition caches.
  bool SameGraphState(const HashGraph &other, std::string *difference) const;
  bool SameTraversalState(const HashGraph &other,
                          std::string *difference) const;
  bool SameTransitionCache(const HashGraph &other,
                           std::string *difference) const;
  size_t DebugBucketCount() const { return vertex_table_.bucket_count(); }
  void DebugResetEmptyBucketCount(size_t count) {
    vertex_table_.reset_empty_bucket_count(count);
  }

  void ClearStatus() {
    ClearStatusFunc func;
    vertex_table_.for_each(func);
  }

  int64_t Assemble(std::vector<ContigGraphVertex> &unitigs);

  // Derive the unitig graph's explicit adjacency while the compact k-mer
  // table is still alive. This avoids rebuilding an endpoint hash map from
  // materialized unitig sequences. Returns false for exceptional topologies
  // that must retain the historical ContigGraph initialization path.
  bool BuildUnitigAdjacency(std::vector<ContigGraphVertex> &unitigs,
                            std::vector<uint32_t> &neighbor_codes,
                            uint64_t *num_edges);

  void reserve(uint64_t capacity) { vertex_table_.reserve(capacity); }

  uint32_t kmer_size() const { return kmer_size_; }
  void set_kmer_size(uint32_t kmer_size) {
    kmer_size_ = kmer_size;
    packed_hash_fixed_accumulator_ =
        IdbaKmer::PackedHashFixedAccumulator(kmer_size);
  }

  int coverage_percentile(double percentile) {
    ++coverage_epoch_;
    if (coverage_epoch_ == 0u) {
      std::fill(coverage_epochs_.begin(), coverage_epochs_.end(), 0u);
      coverage_epoch_ = 1u;
    }
    CoveragePercentileFunc func(coverage_bins_, coverage_epochs_,
                                coverage_epoch_);
    vertex_table_.for_each_value(func);
    return func.percentile(percentile);
  }

  void swap(HashGraph &hash_graph) {
    if (this != &hash_graph) {
      vertex_table_.swap(hash_graph.vertex_table_);
      assembled_endpoint_codes_.swap(
          hash_graph.assembled_endpoint_codes_);
      assemble_arm_bases_[0].swap(hash_graph.assemble_arm_bases_[0]);
      assemble_arm_bases_[1].swap(hash_graph.assemble_arm_bases_[1]);
      coverage_bins_.swap(hash_graph.coverage_bins_);
      coverage_epochs_.swap(hash_graph.coverage_epochs_);
      std::swap(coverage_epoch_, hash_graph.coverage_epoch_);
      std::swap(kmer_size_, hash_graph.kmer_size_);
      std::swap(packed_hash_fixed_accumulator_,
                hash_graph.packed_hash_fixed_accumulator_);
      branch_transition_descriptors_.swap(
          hash_graph.branch_transition_descriptors_);
      branch_transition_blocks_.swap(hash_graph.branch_transition_blocks_);
      std::swap(branch_transition_epoch_,
                hash_graph.branch_transition_epoch_);
      std::swap(branch_transition_hits_, hash_graph.branch_transition_hits_);
      std::swap(branch_transition_misses_,
                hash_graph.branch_transition_misses_);
    }
  }

  uint64_t num_vertices() const { return vertex_table_.size(); }
  void clear() {
    vertex_table_.clear();
    ResetBranchTransitionCache();
  }
  void reset_for_endpoint() {
    vertex_table_.reset_for_endpoint();
    ResetBranchTransitionCache();
    branch_transition_hits_ = 0;
    branch_transition_misses_ = 0;
  }
  uint64_t DebugBranchTransitionHits() const {
    return branch_transition_hits_;
  }
  uint64_t DebugBranchTransitionMisses() const {
    return branch_transition_misses_;
  }

 private:
  template <unsigned Words, bool EmitGroups>
  int64_t InsertPackedKmersImpl(const uint64_t *forward_words,
                                uint32_t sequence_length,
                                bool increment_count,
                                uint32_t *forward_groups,
                                uint32_t *reverse_groups);

#if __cplusplus >= 201103L
  HashGraph(const HashGraph &) = delete;
  const HashGraph &operator=(const HashGraph &) = delete;
#else
  HashGraph(const HashGraph &);
  const HashGraph &operator=(const HashGraph &);
#endif

  bool GetNextVertexAdaptor(const HashGraphVertexAdaptor &current,
                            HashGraphVertexAdaptor &next) {
    if (current.out_edges().size() != 1) return false;

    const uint32_t cached = current.next_neighbor();
    if (cached < HashGraphVertex::kUncacheableNeighbor) {
      next.set_vertex(&vertex_table_.value_at(cached >> 1u), cached & 1u);
    } else {
      IdbaKmer kmer = current.kmer();
      kmer.ShiftAppend(bit_operation::BitToIndex(current.out_edges()));
      next = FindVertexAdaptor(kmer);
    }

    return ((kmer_size_ & 1u) != 0u ||
            !next.vertex().kmer().IsPalindrome()) &&
           next.in_edges().size() == 1;
  }

  bool IsLoop(const IdbaKmer &begin_kmer,
              const HashGraphVertexAdaptor &next) {
    return begin_kmer == next.kmer();
  }

  bool IsPalindromeLoop(const HashGraphVertexAdaptor &current,
                        const HashGraphVertexAdaptor &next) {
    HashGraphVertexAdaptor reverse_next = next;
    reverse_next.ReverseComplement();
    return current == reverse_next;
  }

  void ResetBranchTransitionCache() {
    ++branch_transition_epoch_;
    if (branch_transition_epoch_ == 0u) {
      std::fill(branch_transition_descriptors_.begin(),
                branch_transition_descriptors_.end(), uint64_t{0});
      branch_transition_epoch_ = 1u;
    }
    branch_transition_blocks_.clear();
  }

  bool FindBranchTransition(uint32_t vertex_index, bool is_reverse,
                            uint8_t base, uint32_t *neighbor_code) const {
    if (vertex_index >= branch_transition_descriptors_.size()) return false;
    const uint64_t descriptor =
        branch_transition_descriptors_[vertex_index];
    if (static_cast<uint32_t>(descriptor >> 32u) !=
        branch_transition_epoch_) {
      return false;
    }
    const uint32_t block_index = static_cast<uint32_t>(descriptor);
    const uint32_t code = branch_transition_blocks_[block_index]
        [(is_reverse ? 4u : 0u) + (base & 3u)];
    if (code >= HashGraphVertex::kUncacheableNeighbor) return false;
    *neighbor_code = code;
    return true;
  }

  void RememberBranchTransition(uint32_t vertex_index, bool is_reverse,
                                uint8_t base, uint32_t neighbor_code) {
    if (neighbor_code >= HashGraphVertex::kUncacheableNeighbor) return;
    if (vertex_index >= branch_transition_descriptors_.size()) {
      branch_transition_descriptors_.resize(vertex_index + 1u, uint64_t{0});
    }
    uint64_t &descriptor = branch_transition_descriptors_[vertex_index];
    uint32_t block_index;
    if (static_cast<uint32_t>(descriptor >> 32u) !=
        branch_transition_epoch_) {
      if (branch_transition_blocks_.size() >= UINT32_MAX) return;
      block_index = static_cast<uint32_t>(branch_transition_blocks_.size());
      std::array<uint32_t, 8> block;
      block.fill(HashGraphVertex::kNoCachedNeighbor);
      branch_transition_blocks_.push_back(block);
      descriptor = (uint64_t(branch_transition_epoch_) << 32u) | block_index;
    } else {
      block_index = static_cast<uint32_t>(descriptor);
    }
    uint32_t &slot = branch_transition_blocks_[block_index]
        [(is_reverse ? 4u : 0u) + (base & 3u)];
    if (slot == HashGraphVertex::kNoCachedNeighbor) {
      slot = neighbor_code;
    } else if (slot != neighbor_code) {
      slot = HashGraphVertex::kUncacheableNeighbor;
    }
  }

  class ClearStatusFunc {
   public:
    ClearStatusFunc() {}

    void operator()(HashGraphVertex &vertex) { vertex.status().clear(); }
  };

  class AssembleFunc {
   public:
    AssembleFunc(HashGraph *hash_graph,
                 std::vector<ContigGraphVertex> *unitigs,
                 std::vector<uint32_t> *endpoint_codes)
        : hash_graph_(hash_graph),
          unitigs_(unitigs),
          endpoint_codes_(endpoint_codes) {}
    ~AssembleFunc() {}

    void operator()(HashGraphVertex &vertex);

   private:
    HashGraph *hash_graph_;
    std::vector<ContigGraphVertex> *unitigs_;
    std::vector<uint32_t> *endpoint_codes_;
  };

  class CoveragePercentileFunc {
   public:
    CoveragePercentileFunc(std::vector<size_t> &bins,
                           std::vector<uint32_t> &epochs, uint32_t epoch)
        : bins_(bins), epochs_(epochs), epoch_(epoch) {}

    void operator()(HashGraphVertex &vertex) {
      const size_t coverage = static_cast<size_t>(vertex.count());
      if (coverage >= bins_.size()) {
        bins_.resize(coverage + 1u, 0);
        epochs_.resize(coverage + 1u, 0);
      }
      if (epochs_[coverage] != epoch_) {
        epochs_[coverage] = epoch_;
        bins_[coverage] = 0;
      }
      ++bins_[coverage];
      ++size_;
      max_coverage_ = std::max(max_coverage_, coverage);
    }

    int percentile(double p) const {
      if (size_ == 0) {
        return 0;
      }
      const size_t rank = static_cast<size_t>(size_ * p);
      size_t cumulative = 0;
      for (size_t coverage = 0; coverage <= max_coverage_; ++coverage) {
        if (epochs_[coverage] == epoch_) cumulative += bins_[coverage];
        if (cumulative > rank) {
          return static_cast<int>(coverage);
        }
      }
      return 0;
    }

   private:
    std::vector<size_t> &bins_;
    std::vector<uint32_t> &epochs_;
    uint32_t epoch_;
    size_t size_{0};
    size_t max_coverage_{0};
  };

  vertex_table_type vertex_table_;
  uint32_t kmer_size_;
  uint64_t packed_hash_fixed_accumulator_{0};
  std::vector<uint64_t> branch_transition_descriptors_;
  std::vector<std::array<uint32_t, 8>> branch_transition_blocks_;
  uint32_t branch_transition_epoch_{1};
  uint64_t branch_transition_hits_{0};
  uint64_t branch_transition_misses_{0};
  std::vector<uint32_t> assembled_endpoint_codes_;
  // Per-worker HashGraph instances survive every endpoint task. Keep the two
  // unitig arm tapes here rather than in the short-lived AssembleFunc, so a
  // graph round only rewinds their sizes and never repeats allocator work.
  std::vector<uint8_t> assemble_arm_bases_[2];
  std::vector<size_t> coverage_bins_;
  std::vector<uint32_t> coverage_epochs_;
  uint32_t coverage_epoch_{0};
};

namespace std {
inline void swap(HashGraph &x, HashGraph &y) { x.swap(y); }
}  // namespace std

#endif
