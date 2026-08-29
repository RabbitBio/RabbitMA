//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_H
#define MEGAHIT_SDBG_H

#include "sdbg_def.h"
#include "sdbg_raw_content.h"

#include <omp.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <vector>
#include "kmlib/kmbitvector.h"
#include "kmlib/kmrns.h"

/**
 * Succicent De Bruijn graph
 */
class SDBG {
 public:
  static const uint64_t kNullID = static_cast<uint64_t>(-1);
  SDBG() = default;
  ~SDBG() = default;

  void LoadFromFile(const char *dbg_name) {
    // A loader object can be reused for another graph.  Topology transforms
    // belong to the previous content and must never survive that reload.
    compact_topology_cache_.reset();
    forward_cache_.reset();
    backward_cache_.reset();
    compact_backward_blocks_.clear();
    compact_backward_high_.reset();
    compact_backward_low_.reset();
    compact_backward_high_words_ = 0;
    compact_backward_low_words_ = 0;
    compact_backward_ready_ = false;
    forward_cache_borrowed_ = false;
    simple_path_codes_.reset();
    simple_path_bases_.reset();
    simple_path_edge_count_ = 0;
    simple_path_overflow_count_ = 0;
    reverse_lookup_samples_.clear();
    reverse_lookup_sample_offsets_.clear();
    reverse_lookup_prefix_bases_ = 0;
    reverse_lookup_key_bases_ = 0;
    LoadSdbgRawContent(&content_, dbg_name);
    k_ = content_.meta.k();
    rs_is_tip_.from_packed_array(content_.tip.data(),
                                 content_.meta.item_count());
    rs_w_.from_packed_array(content_.w.data(), content_.meta.item_count());
    rs_last_.from_packed_array(content_.last.data(),
                               content_.meta.item_count());
    invalid_ = kmlib::AtomicBitVector<uint64_t>(
        content_.tip.data(), content_.tip.data() + content_.tip.word_count());
    // An empty prefix bucket must be represented by an empty interval.  The
    // former value-initialized {0, 0} range caused pointless comparisons
    // against edge zero and also made auxiliary indexes sample empty buckets.
    prefix_look_up_.assign(content_.meta.bucket_count(), {1, 0});
    std::fill(f_, f_ + kAlphabetSize + 2, 0);
    f_[0] = -1;
    for (auto it = content_.meta.begin_bucket();
         it != content_.meta.end_bucket() && it->bucket_id != it->kNullID;
         ++it) {
      f_[it->bucket_id / (content_.meta.bucket_count() / kAlphabetSize) + 2] +=
          it->num_items;
      prefix_look_up_[it->bucket_id].first = it->accumulate_item_count;
      prefix_look_up_[it->bucket_id].second =
          it->accumulate_item_count + it->num_items - 1;
    }
    for (unsigned i = 2; i < kAlphabetSize + 2; ++i) {
      f_[i] += f_[i - 1];
    }

    for (unsigned i = 1; i < kAlphabetSize + 2; ++i) {
      rank_f_[i] = rs_last_.rank(f_[i] - 1);
    }

#pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < content_.meta.item_count(); ++i) {
      if (GetW(i) == 0) {
        SetInvalidEdge(i);
      }
    }
  }

  uint64_t size() const { return content_.meta.item_count(); }

  uint32_t k() const { return k_; }

  uint64_t TipSentinelCount() const { return content_.meta.tip_count(); }

  uint8_t GetW(uint64_t x) const { return content_.w[x]; }

  bool IsLast(uint64_t x) const { return content_.last[x]; }

  bool IsLastOrTip(uint64_t x) const {
    return ((content_.last.data()[x / 64] | content_.tip.data()[x / 64]) >>
            (x % 64)) &
           1u;
  }

  int64_t GetLastIndex(uint64_t x) const { return rs_last_.succ(x); }

  uint8_t LastCharOf(uint64_t x) const {
    for (uint8_t i = 1; i < kAlphabetSize + 2; ++i) {
      if (f_[i] > int64_t(x)) {
        return i - 1;
      }
    }
    return kAlphabetSize + 2;
  }

  // Forward() and Backward() are pure topology transforms, but their generic
  // rank/select implementations touch several sparse sampling tables.  The
  // assembly pipeline invokes both transforms repeatedly for almost every
  // edge (tip detection, unitig construction, reverse complements and graph
  // cleaning).  When the caller's memory budget permits, materialize the two
  // 32-bit transforms once and turn those repeated random walks into a single
  // dense load.  The cache is optional and never changes graph state.
  bool BuildTopologyCache(uint64_t host_mem) {
    host_memory_budget_ = host_mem;
    if (compact_topology_cache_ || (forward_cache_ && backward_cache_) ||
        compact_backward_ready_) {
      BuildReverseLookupSamples(host_mem);
      return true;
    }
    const uint64_t num_edges = size();
    const uint64_t max_cached_id =
        std::numeric_limits<uint32_t>::max();
    if (num_edges == 0 || num_edges >= max_cached_id ||
        num_edges > std::numeric_limits<uint64_t>::max() /
                        (2 * sizeof(uint32_t))) {
      return false;
    }
    if (std::getenv("MEGAHIT_DISABLE_SDBG_TOPOLOGY_CACHE") != nullptr) {
      BuildReverseLookupSamples(host_mem);
      return !reverse_lookup_samples_.empty();
    }

    const bool use_compact_cache =
        num_edges < kCompactTopologyNull &&
        std::getenv("MEGAHIT_DISABLE_COMPACT_SDBG_TOPOLOGY_CACHE") == nullptr;
    // A wide cache spends eight bytes per SDBG row and is the source of the
    // 40+ GiB assembly peak on large graphs.  It is retained solely for
    // explicit diagnostics; the normal path uses the bounded reverse-label
    // directory plus compact simple-neighbor blocks below.
    if (!use_compact_cache &&
        std::getenv("MEGAHIT_ENABLE_WIDE_SDBG_TOPOLOGY_CACHE") == nullptr) {
      BuildReverseLookupSamples(host_mem);
      return !reverse_lookup_samples_.empty();
    }
    const uint64_t cache_bytes =
        num_edges * (use_compact_cache ? kCompactTopologyBytesPerEdge
                                       : 2 * sizeof(uint32_t));
    // Keep this opportunistic acceleration below one eighth of the explicit
    // MEGAHIT memory budget.  The remaining budget is needed by graph content,
    // the transient simple-neighbour arrays and downstream unitig structures.
    if (host_mem == 0 || cache_bytes > host_mem / 8) {
      if (!BuildCompactBackwardTopologyCache(host_mem)) {
        BuildReverseLookupSamples(host_mem);
        return !reverse_lookup_samples_.empty();
      }
      BuildReverseLookupSamples(host_mem);
      return true;
    }

    try {
      if (use_compact_cache) {
        compact_topology_cache_.reset(
            new uint8_t[num_edges * kCompactTopologyBytesPerEdge]);
      } else {
        forward_cache_.reset(new uint32_t[num_edges]);
        backward_cache_.reset(new uint32_t[num_edges]);
      }
    } catch (const std::bad_alloc &) {
      compact_topology_cache_.reset();
      forward_cache_.reset();
      backward_cache_.reset();
      if (!BuildCompactBackwardTopologyCache(host_mem)) {
        BuildReverseLookupSamples(host_mem);
        return !reverse_lookup_samples_.empty();
      }
      BuildReverseLookupSamples(host_mem);
      return true;
    }

    if (std::getenv("MEGAHIT_DISABLE_LINEAR_TOPOLOGY_BUILD") == nullptr) {
      const int num_threads = omp_get_max_threads();

      // Every worker's W ranks are monotone within its static edge interval.
      // Keep one select(last) cursor per character: initialize each cursor
      // with one sampled select, then advance by streaming packed last bits.
      // This has the same O(E) work as materializing a dense select array but
      // removes the former 4-byte-per-edge scratch allocation and its extra
      // write/read traffic.
      BuildForwardTopologyCacheLinear(num_edges, num_threads);

      // Backward() is symmetric.  Last ranks advance monotonically over the
      // worker's edge interval, so one select(W,c) cursor per character turns
      // all remaining lookups into sequential scans of the packed W stream.
      std::array<uint64_t, kAlphabetSize + 2> w_select_offset{};
      for (uint8_t c = 0; c <= kAlphabetSize; ++c) {
        // rank(size - 1) returns the precomputed full-word total.  For c=0
        // that total also counts zero padding in the final packed word, which
        // is correct for the generic representation but not for a dense
        // edge-sized select array.  Compute the exact logical-length total.
        const uint64_t exact_count =
            num_edges == 1
                ? static_cast<uint64_t>(GetW(0) == c)
                : static_cast<uint64_t>(rs_w_.rank(c, num_edges - 2)) +
                      static_cast<uint64_t>(GetW(num_edges - 1) == c);
        w_select_offset[c + 1] =
            w_select_offset[c] + exact_count;
      }

#pragma omp parallel num_threads(num_threads)
      {
        const uint64_t tid = static_cast<uint64_t>(omp_get_thread_num());
        const uint64_t team = static_cast<uint64_t>(omp_get_num_threads());
        const uint64_t begin = num_edges * tid / team;
        const uint64_t end = num_edges * (tid + 1) / team;
        int64_t last_rank =
            begin == 0 ? 0 : rs_last_.rank(static_cast<int64_t>(begin - 1));
        const int64_t kUninitializedRank =
            std::numeric_limits<int64_t>::min();
        std::array<int64_t, kAlphabetSize + 1> cursor_rank;
        cursor_rank.fill(kUninitializedRank);
        std::array<uint64_t, kAlphabetSize + 1> cursor_edge{};

        auto select_w = [&](uint8_t c, int64_t target_rank) -> uint64_t {
          const uint64_t char_count =
              w_select_offset[c + 1] - w_select_offset[c];
          (void)char_count;
          assert(target_rank >= 0 &&
                 static_cast<uint64_t>(target_rank) < char_count);
          if (cursor_rank[c] == kUninitializedRank) {
            cursor_edge[c] = rs_w_.select(c, target_rank);
            cursor_rank[c] = target_rank;
            return cursor_edge[c];
          }
          assert(target_rank >= cursor_rank[c]);
          while (cursor_rank[c] < target_rank) {
            do {
              ++cursor_edge[c];
            } while (cursor_edge[c] < num_edges &&
                     GetW(cursor_edge[c]) != c);
            assert(cursor_edge[c] < num_edges);
            ++cursor_rank[c];
          }
          return cursor_edge[c];
        };

        for (uint64_t edge_id = begin; edge_id < end; ++edge_id) {
          const uint8_t c = LastCharOf(edge_id);
          const int64_t select_rank = last_rank - rank_f_[c];
          const uint64_t char_count =
              w_select_offset[c + 1] - w_select_offset[c];
          uint64_t backward;
          if (select_rank >= 0 &&
              static_cast<uint64_t>(select_rank) < char_count) {
            backward = select_w(c, select_rank);
          } else if (select_rank >= 0 &&
                     static_cast<uint64_t>(select_rank) == char_count) {
            backward = num_edges;
          } else {
            backward = kNullID;
          }
          SetCachedBackward(edge_id, backward);
          last_rank += IsLast(edge_id);
        }
      }
    } else {
#pragma omp parallel for schedule(static)
      for (uint64_t edge_id = 0; edge_id < num_edges; ++edge_id) {
        const uint64_t forward = ForwardUncached(edge_id);
        const uint64_t backward = BackwardUncached(edge_id);
        SetCachedForward(edge_id, forward);
        SetCachedBackward(edge_id, backward);
      }
    }
    BuildReverseLookupSamples(host_mem);
    return true;
  }

  uint64_t TopologyCacheBytes() const {
    const uint64_t lookup_bytes =
        reverse_lookup_samples_.size() * sizeof(reverse_lookup_samples_[0]) +
        reverse_lookup_sample_offsets_.size() *
            sizeof(reverse_lookup_sample_offsets_[0]);
    if (compact_topology_cache_) {
      return size() * kCompactTopologyBytesPerEdge + lookup_bytes;
    }
    if (forward_cache_) {
      return size() * 2 * sizeof(uint32_t) + lookup_bytes;
    }
    return CompactBackwardTopologyBytes() + lookup_bytes;
  }

  uint64_t HostMemoryBudget() const { return host_memory_budget_; }

  // Keep the simple-path topology in the same compact representation used by
  // UnitigGraph.  Initial SDBG tip walking and later unitig construction both
  // ask the identical question, so building this snapshot once avoids a
  // second all-edge rank/select pass.  A non-zero cached link remains exact
  // after pruning while both endpoints are valid (deletion cannot introduce a
  // competing edge).  Zero, overflow and stale entries deliberately fall
  // back to the live graph because pruning can expose a new simple link.
  static const uint64_t kSimplePathBlockEdges = 256;

  uint64_t SimplePathSnapshotBytes() const {
    if (!simple_path_codes_) return 0;
    const uint64_t blocks =
        (simple_path_edge_count_ + kSimplePathBlockEdges - 1u) /
        kSimplePathBlockEdges;
    return simple_path_edge_count_ +
           blocks * (kAlphabetSize + 1u) * sizeof(uint32_t);
  }

  uint64_t SimplePathSnapshotOverflowCount() const {
    return simple_path_overflow_count_;
  }

  bool HasSimplePathSnapshot() const {
    return simple_path_codes_ && simple_path_bases_ &&
           simple_path_edge_count_ == size();
  }

  bool BuildSimplePathSnapshot(uint64_t host_mem) {
    if (HasSimplePathSnapshot()) return true;
    if (std::getenv("MEGAHIT_DISABLE_SHARED_SIMPLE_PATH_SNAPSHOT") !=
        nullptr) {
      return false;
    }

    const uint64_t num_edges = size();
    if (num_edges == 0 ||
        num_edges >= std::numeric_limits<uint32_t>::max()) {
      return false;
    }
    const uint64_t num_blocks =
        (num_edges + kSimplePathBlockEdges - 1u) /
        kSimplePathBlockEdges;
    if (num_blocks > std::numeric_limits<uint64_t>::max() /
                         ((kAlphabetSize + 1u) * sizeof(uint32_t))) {
      return false;
    }
    const uint64_t bytes =
        num_edges +
        num_blocks * (kAlphabetSize + 1u) * sizeof(uint32_t);
    // Match the portable auto-selection contract used by UnitigGraph: an
    // optional acceleration may consume at most a quarter of the caller's
    // declared memory budget.  No socket count or machine-specific constant
    // enters the decision.
    if (host_mem != 0 && bytes > host_mem / 4u) {
      return false;
    }

    std::unique_ptr<uint8_t[]> codes;
    std::unique_ptr<uint32_t[]> bases;
    try {
      codes.reset(new uint8_t[num_edges]);
      bases.reset(new uint32_t[num_blocks * (kAlphabetSize + 1u)]);
    } catch (const std::bad_alloc &) {
      return false;
    }

    uint64_t overflow_count = 0;
#pragma omp parallel reduction(+ : overflow_count)
    {
      std::array<uint32_t, kSimplePathBlockEdges> local_next;
#pragma omp for schedule(static)
      for (uint64_t block_id = 0; block_id < num_blocks; ++block_id) {
        uint32_t *block_bases =
            bases.get() + block_id * (kAlphabetSize + 1u);
        std::fill(block_bases, block_bases + kAlphabetSize + 1u,
                  std::numeric_limits<uint32_t>::max());
        const uint64_t edge_begin = block_id * kSimplePathBlockEdges;
        const uint64_t edge_end =
            std::min<uint64_t>(num_edges,
                               edge_begin + kSimplePathBlockEdges);

        for (uint64_t edge = edge_begin; edge < edge_end; ++edge) {
          uint32_t next_id = std::numeric_limits<uint32_t>::max();
          if (IsValidEdge(edge)) {
            const uint64_t next = NextSimplePathEdge(edge);
            if (next != kNullID) {
              next_id = static_cast<uint32_t>(next);
              uint8_t c = GetW(edge);
              if (c > kAlphabetSize) c -= kAlphabetSize;
              block_bases[c] = std::min(block_bases[c], next_id);
            }
          }
          local_next[edge - edge_begin] = next_id;
        }

        for (uint64_t edge = edge_begin; edge < edge_end; ++edge) {
          const uint32_t next_id = local_next[edge - edge_begin];
          if (next_id == std::numeric_limits<uint32_t>::max()) {
            codes[edge] = 0;
            continue;
          }
          uint8_t c = GetW(edge);
          if (c > kAlphabetSize) c -= kAlphabetSize;
          const uint32_t delta = next_id - block_bases[c];
          if (delta < 254u) {
            codes[edge] = static_cast<uint8_t>(delta + 1u);
          } else {
            codes[edge] = 255u;
            ++overflow_count;
          }
        }
      }
    }

    simple_path_codes_ = std::move(codes);
    simple_path_bases_ = std::move(bases);
    simple_path_edge_count_ = num_edges;
    simple_path_overflow_count_ = overflow_count;
    return true;
  }

  // Return true only when the compact snapshot itself supplied an exact,
  // still-live simple successor.  Callers must use the live graph on false.
  bool TryCachedNextSimplePathEdge(uint64_t edge, uint64_t *next) const {
    if (!HasSimplePathSnapshot() || edge >= simple_path_edge_count_) {
      return false;
    }
    const uint8_t code = simple_path_codes_[edge];
    if (code == 0 || code == 255u) return false;
    uint8_t c = GetW(edge);
    if (c > kAlphabetSize) c -= kAlphabetSize;
    const uint32_t base = simple_path_bases_[
        (edge / kSimplePathBlockEdges) * (kAlphabetSize + 1u) + c];
    if (base == std::numeric_limits<uint32_t>::max()) return false;
    const uint64_t cached =
        static_cast<uint64_t>(base) + static_cast<uint64_t>(code - 1u);
    if (cached >= simple_path_edge_count_ || !IsValidEdge(cached)) {
      return false;
    }
    *next = cached;
    return true;
  }

  uint64_t CachedOrLiveNextSimplePathEdge(uint64_t edge) const {
    uint64_t next = kNullID;
    return TryCachedNextSimplePathEdge(edge, &next)
               ? next
               : NextSimplePathEdge(edge);
  }

  bool HasReusableForwardTopologyCache() const {
    return forward_cache_ && backward_cache_ && !forward_cache_borrowed_;
  }

  // Unitig construction needs a dense simple-next snapshot for one phase.
  // A wide Forward cache has the same element type and is not otherwise read
  // by that phase after each entry has been converted.  Let the graph builder
  // reuse it in place, avoiding a second edge-sized allocation.  Backward and
  // reverse-lookup acceleration remain live throughout the conversion.
  uint32_t *BorrowForwardCacheForSimpleNeighbors() {
    if (!forward_cache_ || !backward_cache_ || forward_cache_borrowed_ ||
        std::getenv("MEGAHIT_DISABLE_TOPOLOGY_SCRATCH_REUSE") != nullptr) {
      return nullptr;
    }
    forward_cache_borrowed_ = true;
    return forward_cache_.get();
  }

  uint64_t ConvertBorrowedForwardEntryToSimpleNeighbor(uint64_t edge_id) {
    assert(forward_cache_borrowed_ && forward_cache_ &&
           edge_id < size());
    const uint32_t encoded_forward = forward_cache_[edge_id];
    const uint64_t first_outgoing =
        encoded_forward == std::numeric_limits<uint32_t>::max()
            ? kNullID
            : static_cast<uint64_t>(encoded_forward);
    uint64_t next_edge = 0;
    uint64_t simple_next = kNullID;
    if (ComputeOutgoingsFromFirst<kFlagWriteOut | kFlagMustEq1>(
            edge_id, first_outgoing, &next_edge) == 1 &&
        UniquePrevEdge(next_edge) != kNullID) {
      simple_next = next_edge;
    }
    forward_cache_[edge_id] =
        simple_next < std::numeric_limits<uint32_t>::max()
            ? static_cast<uint32_t>(simple_next)
            : std::numeric_limits<uint32_t>::max();
    return simple_next;
  }

  void RestoreBorrowedForwardTopologyCache() {
    if (!forward_cache_borrowed_) {
      return;
    }
    BuildForwardTopologyCacheLinear(size(), omp_get_max_threads());
    forward_cache_borrowed_ = false;
  }

  bool IsValidEdge(uint64_t edge_id) const { return !invalid_.at(edge_id); }

  bool IsTip(uint64_t edge_id) const { return content_.tip[edge_id]; }

  void SetValidEdge(uint64_t edge_id) { invalid_.unset(edge_id); }

  void SetInvalidEdge(uint64_t edge_id) { invalid_.set(edge_id); }

  mul_t EdgeMultiplicity(uint64_t edge_id) const {
    if (!content_.full_mul.empty()) {
      return content_.full_mul[edge_id];
    }
    if (content_.small_mul[edge_id] != kSmallMulSentinel) {
      return content_.small_mul[edge_id];
    } else {
      const uint64_t block_edges = SdbgRawContent::kLargeMulRankBlockEdges;
      const uint64_t block = edge_id / block_edges;
      const uint64_t first_word = block * (block_edges / 64u);
      const uint64_t edge_word = edge_id >> 6u;
      uint64_t rank = content_.large_mul_rank[block];
      for (uint64_t word = first_word; word < edge_word; ++word) {
        rank += static_cast<uint64_t>(
            __builtin_popcountll(content_.large_mul_bits[word]));
      }
      const unsigned bit = static_cast<unsigned>(edge_id & 63u);
      if (bit != 0u) {
        rank += static_cast<uint64_t>(__builtin_popcountll(
            content_.large_mul_bits[edge_word] &
            ((uint64_t{1} << bit) - 1u)));
      }
      assert(rank < content_.large_mul_values.size());
      return content_.large_mul_values[rank];
    }
  }

  void PrefetchMultiplicity(uint64_t edge_id) const {
    if (!content_.full_mul.empty()) {
      __builtin_prefetch(content_.full_mul.data() + edge_id, 0, 1);
    } else {
      __builtin_prefetch(content_.small_mul.data() + edge_id, 0, 1);
    }
  }

  uint64_t Forward(uint64_t edge_id) const {  // the last edge edge_id points to
    if (compact_topology_cache_) {
      const uint32_t cached = PackedForward(edge_id);
      return cached == kCompactTopologyNull ? kNullID
                                             : static_cast<uint64_t>(cached);
    }
    if (forward_cache_ && !forward_cache_borrowed_) {
      const uint32_t cached = forward_cache_[edge_id];
      return cached == std::numeric_limits<uint32_t>::max()
                 ? kNullID
                 : static_cast<uint64_t>(cached);
    }
    return ForwardUncached(edge_id);
  }

  uint64_t Backward(
      uint64_t edge_id) const {  // the first edge points to edge_id
    if (compact_topology_cache_) {
      const uint32_t cached = PackedBackward(edge_id);
      return cached == kCompactTopologyNull ? kNullID
                                             : static_cast<uint64_t>(cached);
    }
    if (backward_cache_) {
      const uint32_t cached = backward_cache_[edge_id];
      return cached == std::numeric_limits<uint32_t>::max()
                 ? kNullID
                 : static_cast<uint64_t>(cached);
    }
    if (compact_backward_ready_) {
      return CompactCachedBackward(edge_id);
    }
    return BackwardUncached(edge_id);
  }

 private:
  uint64_t ForwardUncached(uint64_t edge_id) const {
    uint8_t a = GetW(edge_id);
    if (a > kAlphabetSize) {
      a -= kAlphabetSize;
    }
    int64_t count_a = rs_w_.rank(a, edge_id);
    return rs_last_.select(rank_f_[a] + count_a - 1);
  }

  uint64_t BackwardUncached(uint64_t edge_id) const {
    uint8_t a = LastCharOf(edge_id);
    int64_t count_a = rs_last_.rank(edge_id - 1) - rank_f_[a];
    return rs_w_.select(a, count_a);
  }

  const label_word_t *TipLabelStartPtr(uint64_t edge_id) const {
    return content_.tip_lables.data() +
           content_.meta.words_per_tip_label() * (rs_is_tip_.rank(edge_id) - 1);
  }
  static uint8_t CharAtTipLabel(const label_word_t *label_start_ptr,
                                unsigned offset) {
    return kmlib::CompactVector<kBitsPerChar, label_word_t>::at(label_start_ptr,
                                                                offset) +
           1;
  }

 public:
  /**
   * Find the index given a sequence
   * @param seq the sequence encoded in [0-3]
   * @return the index of the sequence in the graph, kNullID if not exists
   */
  uint64_t IndexBinarySearch(const uint8_t *seq) const {
    uint64_t prefix = 0;
    for (uint64_t i = 0; (1u << (i * kBitsPerChar)) < prefix_look_up_.size();
         ++i) {
      prefix = prefix * kAlphabetSize + seq[k_ - 1 - i] - 1;
    }
    auto l = prefix_look_up_[prefix].first;
    auto r = prefix_look_up_[prefix].second;

    // A compact, exact second-level index narrows large metadata buckets to
    // one or a few fixed-size blocks.  Equal sampled keys deliberately widen
    // the candidate interval; the full comparison loop below remains the
    // source of truth, including for repetitive labels and tip nodes.
    if (!reverse_lookup_samples_.empty() && l <= r &&
        prefix < reverse_lookup_sample_offsets_.size() - 1u) {
      uint32_t lookup_key = 0;
      for (unsigned i = 0; i < reverse_lookup_key_bases_; ++i) {
        lookup_key = (lookup_key << kBitsPerChar) |
                     static_cast<uint32_t>(
                         seq[k_ - 1u - reverse_lookup_prefix_bases_ - i] - 1u);
      }
      const uint64_t sample_begin =
          reverse_lookup_sample_offsets_[prefix];
      const uint64_t sample_end =
          reverse_lookup_sample_offsets_[prefix + 1u];
      if (sample_begin < sample_end) {
        const auto begin_it = reverse_lookup_samples_.begin() + sample_begin;
        const auto end_it = reverse_lookup_samples_.begin() + sample_end;
        const auto lower_it =
            std::lower_bound(begin_it, end_it, lookup_key);
        const auto upper_it =
            std::upper_bound(lower_it, end_it, lookup_key);
        const uint64_t lower_sample = lower_it - begin_it;
        const uint64_t upper_sample = upper_it - begin_it;
        const uint64_t sample_count = sample_end - sample_begin;
        const int64_t bucket_begin = l;
        if (lower_it != upper_it) {
          const uint64_t first_block =
              lower_sample == 0 ? 0 : lower_sample - 1u;
          l = bucket_begin + static_cast<int64_t>(
                                   first_block * kReverseLookupSampleStride);
          if (upper_sample < sample_count) {
            r = std::min<int64_t>(
                r, bucket_begin + static_cast<int64_t>(
                                      upper_sample *
                                          kReverseLookupSampleStride) -
                       1);
          }
        } else {
          // The first sample is the smallest key in this bucket.  A smaller
          // query cannot occur anywhere in the bucket.
          if (lower_sample == 0) {
            return kNullID;
          }
          l = bucket_begin + static_cast<int64_t>(
                                   (lower_sample - 1u) *
                                   kReverseLookupSampleStride);
          if (lower_sample < sample_count) {
            r = std::min<int64_t>(
                r, bucket_begin + static_cast<int64_t>(
                                      lower_sample *
                                          kReverseLookupSampleStride) -
                       1);
          }
        }
      }
    }

    while (l <= r) {
      int cmp = 0;
      uint64_t mid = (l + r) / 2;
      uint64_t y = mid;

      for (int i = k_ - 1; i >= 0; --i) {
        if (IsTip(y)) {
          const label_word_t *tip_label = TipLabelStartPtr(y);
          for (int j = 0; j < i; ++j) {
            auto c = CharAtTipLabel(tip_label, j);
            if (c < seq[i - j]) {
              cmp = -1;
              break;
            } else if (c > seq[i - j]) {
              cmp = 1;
              break;
            }
          }

          if (cmp == 0) {
            if (IsTip(mid)) {
              cmp = -1;
            } else {
              auto c = CharAtTipLabel(tip_label, i);
              if (c < seq[0]) {
                cmp = -1;
                break;
              } else if (c > seq[0]) {
                cmp = 1;
                break;
              }
            }
          }
          break;
        }

        y = Backward(y);
        uint8_t c = GetW(y);

        if (c < seq[i]) {
          cmp = -1;
          break;
        } else if (c > seq[i]) {
          cmp = 1;
          break;
        }
      }

      if (cmp == 0) {
        return GetLastIndex(mid);
      } else if (cmp > 0) {
        r = mid - 1;
      } else {
        l = mid + 1;
      }
    }
    return kNullID;
  }
  /**
   * Fetch the label of an edge from the graph
   * @param id the index in the graph
   * @param seq the label will be written to this address
   * @return the length of label (always k)
   */
  uint32_t GetLabel(uint64_t id, uint8_t *seq) const {
    uint64_t x = id;
    for (int i = k_ - 1; i >= 0; --i) {
      if (IsTip(x)) {
        const label_word_t *tip_label = TipLabelStartPtr(x);
        for (int j = 0; j <= i; ++j) {
          seq[i - j] = CharAtTipLabel(tip_label, j);
        }
        break;
      }
      x = Backward(x);
      seq[i] = GetW(x);
      if (seq[i] > kAlphabetSize) {
        seq[i] -= kAlphabetSize;
      }
    }
    return k_;
  }

 private:
  static const uint8_t kFlagWriteOut = 0x1;
  static const uint8_t kFlagMustEq0 = 0x2;
  static const uint8_t kFlagMustEq1 = 0x4;
  /**
   * An internal function to collect incoming edges & in-degrees
   * @tparam flag
   * @param edge_id
   * @param incomings the incoming edges will be written in the address if
   * kFlagWriteOut is set
   * @return in degree of the edge; -1 if edge or flag invalid
   */
  template <uint8_t flag = 0>
  int ComputeIncomings(uint64_t edge_id, uint64_t *incomings) const {
    if (!IsValidEdge(edge_id)) {
      return -1;
    }

    uint64_t first_income = Backward(edge_id);
    uint8_t c = GetW(first_income);
    unsigned count_ones = IsLastOrTip(first_income);
    int indegree = IsValidEdge(first_income);

    if (flag & kFlagMustEq0) {
      if (indegree) return -1;
    }

    if (flag & kFlagWriteOut) {
      if (indegree > 0) {
        incomings[0] = first_income;
      }
    }

    for (uint64_t y = first_income + 1;
         count_ones < kAlphabetSize + 1 && y < size(); ++y) {
      count_ones += IsLastOrTip(y);
      uint8_t cur_char = GetW(y);

      if (cur_char == c) {
        break;
      } else if (cur_char == c + kAlphabetSize && IsValidEdge(y)) {
        if (flag & kFlagMustEq0) {
          return -1;
        } else if (flag & kFlagMustEq1) {
          if (indegree == 1) return -1;
        }
        if (flag & kFlagWriteOut) {
          assert(incomings != nullptr);
          incomings[indegree] = y;
        }
        ++indegree;
      }
    }
    return indegree;
  }
  /**
   * An internal function to collect outgoing edges & in-degrees
   * @tparam flag
   * @param edge_id
   * @param outgoings the outgoing edges will be written in the address if
   * kFlagWriteOut is set
   * @return out degree of the edge; -1 if edge or flag invalid
   */
  template <uint8_t flag = 0>
  int ComputeOutgoingsFromFirst(uint64_t edge_id, uint64_t next_edge,
                                uint64_t *outgoings) const {
    if (!IsValidEdge(edge_id)) {
      return -1;
    }
    uint64_t outdegree = 0;
    do {
      if (IsValidEdge(next_edge)) {
        if (flag & kFlagMustEq0) {
          return -1;
        } else if (flag & kFlagMustEq1) {
          if (outdegree == 1) return -1;
        }
        if (flag & kFlagWriteOut) {
          assert(outgoings != nullptr);
          outgoings[outdegree] = next_edge;
        }
        ++outdegree;
      }
      --next_edge;
    } while (next_edge != kNullID && !IsLastOrTip(next_edge));

    return outdegree;
  }

  template <uint8_t flag = 0>
  int ComputeOutgoings(uint64_t edge_id, uint64_t *outgoings) const {
    return ComputeOutgoingsFromFirst<flag>(edge_id, Forward(edge_id),
                                           outgoings);
  }

 public:
  /**
   * the in-degree of a node/edge
   * @param edge_id
   * @return the in-degree. -1 if id invalid.
   */
  int EdgeIndegree(uint64_t edge_id) const {
    return ComputeIncomings(edge_id, nullptr);
  }
  /**
   * the out-degree of an edge
   * @param edge_id
   * @return the out-degree. -1 if id invalid
   */
  int EdgeOutdegree(uint64_t edge_id) const {
    return ComputeOutgoings(edge_id, nullptr);
  }
  /**
   * get all incoming edges of an edge
   * @param edge_id
   * @param incomings all incoming edges' id will be written here
   * @return in-degree
   */
  int IncomingEdges(uint64_t edge_id, uint64_t *incomings) const {
    return ComputeIncomings<kFlagWriteOut>(edge_id, incomings);
  }
  /**
   * get all outgoing edges of an edge
   * @param edge_id
   * @param outgoings all outgoing edges' id will be written here
   * @return out-degree
   */
  int OutgoingEdges(uint64_t edge_id, uint64_t *outgoings) const {
    return ComputeOutgoings<kFlagWriteOut>(edge_id, outgoings);
  }
  /**
   * Query zero/unique/branch degree while writing the sole neighbor when it
   * exists.  A return value of 0 means zero degree, 1 means unique, and -1
   * means invalid or at least two neighbors.  Unlike IncomingEdges and
   * OutgoingEdges, branch detection stops as soon as the second valid edge is
   * found, which is useful to path walkers that never need the full list.
   */
  int UniqueIncomingEdge(uint64_t edge_id, uint64_t *incoming) const {
    return ComputeIncomings<kFlagWriteOut | kFlagMustEq1>(edge_id, incoming);
  }
  int UniqueOutgoingEdge(uint64_t edge_id, uint64_t *outgoing) const {
    return ComputeOutgoings<kFlagWriteOut | kFlagMustEq1>(edge_id, outgoing);
  }
  /**
   * a more efficient way to judge whether an edge's in-degree is 0
   * @param edge_id
   * @return true if the edge's in-degree is 0
   */
  bool EdgeIndegreeZero(uint64_t edge_id) const {
    return ComputeIncomings<kFlagMustEq0>(edge_id, nullptr) == 0;
  }
  /**
   * A more efficient way to judge whether an edge's out-degree is 0
   * @param edge_id
   * @return true if the edge's out-degree is 0
   */
  bool EdgeOutdegreeZero(uint64_t edge_id) const {
    return ComputeOutgoings<kFlagMustEq0>(edge_id, nullptr) == 0;
  }

  /**
   * Initialize the exact set of valid source/sink-adjacent edges used by tip
   * pruning without issuing two degree queries for every row in the BOSS
   * table.  seq2sdbg materializes a $->v sentinel for every source and a
   * v->$ sentinel (W==0) for every sink.  Traversing only those sentinels
   * enumerates precisely the valid zero-indegree/zero-outdegree edges.  All
   * other bits stay ignored.  Hashing, graph ordering and pruning barriers are
   * unchanged; this is only a cheaper way to construct the initial frontier.
   */
  uint64_t InitializeTipEndpointMask(AtomicBitVector *ignored) const {
    assert(ignored != nullptr && ignored->size() == size());
    const uint64_t words = ignored->word_count();
#pragma omp parallel for schedule(static)
    for (uint64_t word = 0; word < words; ++word) {
      ignored->store_word(word, ~AtomicBitVector::word_type(0));
    }

    uint64_t endpoints = 0;
#pragma omp parallel for schedule(static) reduction(+ : endpoints)
    for (uint64_t sentinel = 0; sentinel < size(); ++sentinel) {
      if (IsTip(sentinel)) {
        uint64_t outgoing = Forward(sentinel);
        while (outgoing != kNullID && outgoing < size()) {
          if (IsValidEdge(outgoing)) {
            endpoints += ignored->try_unset(outgoing);
          }
          if (outgoing == 0) break;
          --outgoing;
          if (IsLastOrTip(outgoing)) break;
        }
      }

      if (GetW(sentinel) == 0) {
        const uint64_t first_incoming = Backward(sentinel);
        if (first_incoming == kNullID || first_incoming >= size()) continue;
        const uint8_t first_char = GetW(first_incoming);
        unsigned boundaries = IsLastOrTip(first_incoming);
        if (IsValidEdge(first_incoming)) {
          endpoints += ignored->try_unset(first_incoming);
        }
        for (uint64_t incoming = first_incoming + 1u;
             boundaries < kAlphabetSize + 1u && incoming < size();
             ++incoming) {
          boundaries += IsLastOrTip(incoming);
          const uint8_t current_char = GetW(incoming);
          if (current_char == first_char) break;
          if (current_char == first_char + kAlphabetSize &&
              IsValidEdge(incoming)) {
            endpoints += ignored->try_unset(incoming);
          }
        }
      }
    }
    return endpoints;
  }
  /**
   * @param edge_id
   * @return if the edge has only one outgoing edge, return that one; otherwise
   * -1
   */
  uint64_t UniqueNextEdge(uint64_t edge_id) const {
    uint64_t ret = 0;
    if (UniqueOutgoingEdge(edge_id, &ret) == 1) {
      return ret;
    } else {
      return kNullID;
    }
  }
  /**
   * @param edge_id
   * @return if the edge has only one incoming edge, return that one; otherwise
   * -1
   */
  uint64_t UniquePrevEdge(uint64_t edge_id) const {
    uint64_t ret = 0;
    if (UniqueIncomingEdge(edge_id, &ret) == 1) {
      return ret;
    } else {
      return kNullID;
    }
  }

  /**
   * @param edge_id
   * @return if the edge has only one incoming edge which is on a simple path,
   * return that one;
   * otherwise -1
   */
  uint64_t PrevSimplePathEdge(uint64_t edge_id) const {
    uint64_t prev_edge = UniquePrevEdge(edge_id);
    if (prev_edge != kNullID && UniqueNextEdge(prev_edge) != kNullID) {
      return prev_edge;
    } else {
      return kNullID;
    }
  }
  /**
   * @param edge_id
   * @return if the edge has only one outgoing edge which is on a simple path,
   * return that one;
   * otherwise -1
   */
  uint64_t NextSimplePathEdge(uint64_t edge_id) const {
    uint64_t next_edge = UniqueNextEdge(edge_id);
    if (next_edge != kNullID && UniquePrevEdge(next_edge) != kNullID) {
      assert(next_edge < size());
      return next_edge;
    } else {
      return kNullID;
    }
  }

  // A BOSS/SDBG node has at most one edge per alphabet symbol.  Therefore a
  // simple successor is one of the few consecutive entries ending at
  // Forward(edge_id).  Return that small offset together with the exact edge
  // ID so unitig construction can retain four bits per edge instead of a
  // dense 32-bit successor table.  The offset is only a representation of the
  // already-computed topology; validity and degree tests are identical to
  // NextSimplePathEdge().
  uint64_t NextSimplePathEdgeWithOffset(uint64_t edge_id,
                                        uint8_t *forward_offset) const {
    const uint64_t first_outgoing = Forward(edge_id);
    uint64_t next_edge = 0;
    if (ComputeOutgoingsFromFirst<kFlagWriteOut | kFlagMustEq1>(
            edge_id, first_outgoing, &next_edge) == 1 &&
        UniquePrevEdge(next_edge) != kNullID) {
      assert(next_edge < size());
      assert(first_outgoing >= next_edge);
      const uint64_t offset = first_outgoing - next_edge;
      // Four bits reserve zero for "no simple successor" and encode offsets
      // 0..14 as values 1..15.  The DNA alphabet needs at most four entries,
      // but retain the wider assertion for format robustness.
      assert(offset < 15u);
      *forward_offset = static_cast<uint8_t>(offset);
      return next_edge;
    }
    *forward_offset = 0;
    return kNullID;
  }
  /**
   * @param edge_id
   * @return the index of the reverse-complement edge of the input edge
   */
  uint64_t EdgeReverseComplement(uint64_t edge_id) const {
    if (!IsValidEdge(edge_id)) {
      return kNullID;
    }

    uint8_t seq[kMaxK + 1];
    GetLabel(edge_id, seq);
    seq[k_] = GetW(edge_id);

    if (seq[k_] > kAlphabetSize) {
      seq[k_] -= kAlphabetSize;
    }

    for (int i = 0, j = k_; i < j; ++i, --j) {
      std::swap(seq[i], seq[j]);
    }
    for (unsigned i = 0; i < k_ + 1; ++i) {
      seq[i] = kAlphabetSize + 1 - seq[i];
    }

    uint64_t rev_node = IndexBinarySearch(seq);
    if (rev_node == kNullID) return kNullID;
    do {
      uint8_t edge_label = GetW(rev_node);
      if (edge_label == seq[k_] || edge_label - kAlphabetSize == seq[k_]) {
        assert(rev_node < size());
        return rev_node;
      }
      --rev_node;
    } while (rev_node != kNullID && !IsLastOrTip(rev_node));

    return kNullID;
  }

  // Resolve several independent reverse-complement endpoints in lockstep.
  // The scalar algorithm is preserved exactly, but its dependent Backward()
  // walks are interleaved so a worker can keep multiple memory requests in
  // flight.  This is especially useful with the exact succinct topology,
  // where one endpoint by itself exposes little memory-level parallelism.
  void EdgeReverseComplementBatch(const uint64_t *edge_ids,
                                  uint64_t *results, size_t count) const {
    static const size_t kMaxBatch = 64;
    assert(count <= kMaxBatch);
    if (count == 0) {
      return;
    }
    if (count == 1) {
      results[0] = EdgeReverseComplement(edge_ids[0]);
      return;
    }

    std::array<std::array<uint8_t, kMaxK + 1>, kMaxBatch> sequences{};
    std::array<uint64_t, kMaxBatch> label_edge{};
    std::array<uint8_t, kMaxBatch> label_done{};
    std::array<uint8_t, kMaxBatch> valid{};

    for (size_t q = 0; q < count; ++q) {
      results[q] = kNullID;
      if (!IsValidEdge(edge_ids[q])) {
        continue;
      }
      valid[q] = 1;
      label_edge[q] = edge_ids[q];
      uint8_t edge_char = GetW(edge_ids[q]);
      if (edge_char > kAlphabetSize) {
        edge_char -= kAlphabetSize;
      }
      sequences[q][k_] = edge_char;
    }

    // This is GetLabel() with the outer query dimension interleaved.
    for (int i = static_cast<int>(k_) - 1; i >= 0; --i) {
      for (size_t q = 0; q < count; ++q) {
        if (!valid[q] || label_done[q]) {
          continue;
        }
        const uint64_t x = label_edge[q];
        if (IsTip(x)) {
          const label_word_t *tip_label = TipLabelStartPtr(x);
          for (int j = 0; j <= i; ++j) {
            sequences[q][i - j] = CharAtTipLabel(tip_label, j);
          }
          label_done[q] = 1;
          continue;
        }
        label_edge[q] = Backward(x);
        uint8_t c = GetW(label_edge[q]);
        if (c > kAlphabetSize) {
          c -= kAlphabetSize;
        }
        sequences[q][i] = c;
      }
    }

    for (size_t q = 0; q < count; ++q) {
      if (!valid[q]) {
        continue;
      }
      for (unsigned i = 0, j = k_; i < j; ++i, --j) {
        std::swap(sequences[q][i], sequences[q][j]);
      }
      for (unsigned i = 0; i < k_ + 1u; ++i) {
        sequences[q][i] = kAlphabetSize + 1u - sequences[q][i];
      }
    }

    std::array<int64_t, kMaxBatch> left{};
    std::array<int64_t, kMaxBatch> right{};
    std::array<uint64_t, kMaxBatch> mid{};
    std::array<uint64_t, kMaxBatch> walk_edge{};
    std::array<int, kMaxBatch> walk_pos{};
    std::array<int, kMaxBatch> compare{};
    std::array<uint8_t, kMaxBatch> searching{};
    std::array<uint8_t, kMaxBatch> comparing{};

    for (size_t q = 0; q < count; ++q) {
      if (!valid[q]) {
        continue;
      }
      const uint8_t *seq = sequences[q].data();
      uint64_t prefix = 0;
      for (uint64_t i = 0;
           (uint64_t{1} << (i * kBitsPerChar)) < prefix_look_up_.size();
           ++i) {
        prefix = prefix * kAlphabetSize + seq[k_ - 1u - i] - 1u;
      }
      int64_t l = prefix_look_up_[prefix].first;
      int64_t r = prefix_look_up_[prefix].second;

      if (!reverse_lookup_samples_.empty() && l <= r &&
          prefix < reverse_lookup_sample_offsets_.size() - 1u) {
        uint32_t lookup_key = 0;
        for (unsigned i = 0; i < reverse_lookup_key_bases_; ++i) {
          lookup_key = (lookup_key << kBitsPerChar) |
                       static_cast<uint32_t>(
                           seq[k_ - 1u - reverse_lookup_prefix_bases_ - i] -
                           1u);
        }
        const uint64_t sample_begin =
            reverse_lookup_sample_offsets_[prefix];
        const uint64_t sample_end =
            reverse_lookup_sample_offsets_[prefix + 1u];
        if (sample_begin < sample_end) {
          const auto begin_it = reverse_lookup_samples_.begin() + sample_begin;
          const auto end_it = reverse_lookup_samples_.begin() + sample_end;
          const auto lower_it =
              std::lower_bound(begin_it, end_it, lookup_key);
          const auto upper_it = std::upper_bound(lower_it, end_it, lookup_key);
          const uint64_t lower_sample = lower_it - begin_it;
          const uint64_t upper_sample = upper_it - begin_it;
          const uint64_t sample_count = sample_end - sample_begin;
          const int64_t bucket_begin = l;
          if (lower_it != upper_it) {
            const uint64_t first_block =
                lower_sample == 0 ? 0 : lower_sample - 1u;
            l = bucket_begin + static_cast<int64_t>(
                                   first_block * kReverseLookupSampleStride);
            if (upper_sample < sample_count) {
              r = std::min<int64_t>(
                  r, bucket_begin + static_cast<int64_t>(
                                        upper_sample *
                                            kReverseLookupSampleStride) -
                         1);
            }
          } else {
            if (lower_sample == 0) {
              continue;
            }
            l = bucket_begin + static_cast<int64_t>(
                                   (lower_sample - 1u) *
                                   kReverseLookupSampleStride);
            if (lower_sample < sample_count) {
              r = std::min<int64_t>(
                  r, bucket_begin + static_cast<int64_t>(
                                        lower_sample *
                                            kReverseLookupSampleStride) -
                         1);
            }
          }
        }
      }
      left[q] = l;
      right[q] = r;
      searching[q] = l <= r;
    }

    for (;;) {
      size_t active_searches = 0;
      for (size_t q = 0; q < count; ++q) {
        if (!searching[q]) {
          continue;
        }
        if (left[q] > right[q]) {
          searching[q] = 0;
          continue;
        }
        mid[q] = static_cast<uint64_t>((left[q] + right[q]) / 2);
        walk_edge[q] = mid[q];
        walk_pos[q] = static_cast<int>(k_) - 1;
        compare[q] = 0;
        comparing[q] = 1;
        ++active_searches;
      }
      if (active_searches == 0) {
        break;
      }

      for (;;) {
        size_t active_comparisons = 0;
        for (size_t q = 0; q < count; ++q) {
          if (!comparing[q]) {
            continue;
          }
          ++active_comparisons;
          const int i = walk_pos[q];
          uint64_t y = walk_edge[q];
          const uint8_t *seq = sequences[q].data();
          if (IsTip(y)) {
            const label_word_t *tip_label = TipLabelStartPtr(y);
            for (int j = 0; j < i; ++j) {
              const uint8_t c = CharAtTipLabel(tip_label, j);
              if (c < seq[i - j]) {
                compare[q] = -1;
                break;
              }
              if (c > seq[i - j]) {
                compare[q] = 1;
                break;
              }
            }
            if (compare[q] == 0) {
              if (IsTip(mid[q])) {
                compare[q] = -1;
              } else {
                const uint8_t c = CharAtTipLabel(tip_label, i);
                if (c < seq[0]) {
                  compare[q] = -1;
                } else if (c > seq[0]) {
                  compare[q] = 1;
                }
              }
            }
            comparing[q] = 0;
            continue;
          }

          y = Backward(y);
          walk_edge[q] = y;
          const uint8_t c = GetW(y);
          if (c < seq[i]) {
            compare[q] = -1;
            comparing[q] = 0;
          } else if (c > seq[i]) {
            compare[q] = 1;
            comparing[q] = 0;
          } else if (i == 0) {
            comparing[q] = 0;
          } else {
            walk_pos[q] = i - 1;
          }
        }
        if (active_comparisons == 0) {
          break;
        }
      }

      for (size_t q = 0; q < count; ++q) {
        if (!searching[q]) {
          continue;
        }
        if (compare[q] == 0) {
          results[q] = GetLastIndex(mid[q]);
          searching[q] = 0;
        } else if (compare[q] > 0) {
          right[q] = static_cast<int64_t>(mid[q]) - 1;
        } else {
          left[q] = static_cast<int64_t>(mid[q]) + 1;
        }
      }
    }

    for (size_t q = 0; q < count; ++q) {
      uint64_t rev_node = results[q];
      if (rev_node == kNullID) {
        continue;
      }
      results[q] = kNullID;
      do {
        const uint8_t edge_label = GetW(rev_node);
        if (edge_label == sequences[q][k_] ||
            edge_label - kAlphabetSize == sequences[q][k_]) {
          results[q] = rev_node;
          break;
        }
        --rev_node;
      } while (rev_node != kNullID && !IsLastOrTip(rev_node));
    }
  }

  /**
   * free multiplicity of all edges to reduce memory
   * WARNING: use this with cautions
   * After that EdgeMultiplicity() is invalid
   */
  void FreeMultiplicity() {
    content_.small_mul = std::vector<small_mul_t>();
    content_.large_mul_bits = std::vector<uint64_t>();
    content_.large_mul_rank = std::vector<uint64_t>();
    content_.large_mul_values = std::vector<mul_t>();
    content_.full_mul = std::vector<mul_t>();
  }

 private:
  static const uint32_t kCompactTopologyNull = (uint32_t{1} << 28u) - 1u;
  static const uint64_t kCompactTopologyBytesPerEdge = 7;
  static const uint32_t kCompactBackwardBlockSize = 512;
  static const uint64_t kReverseLookupSampleStride = 128;
  static const unsigned kMaxReverseLookupKeyBases =
      sizeof(uint32_t) * 8u / kBitsPerChar;

  void BuildForwardTopologyCacheLinear(uint64_t num_edges, int num_threads) {
#pragma omp parallel num_threads(num_threads)
    {
      const uint64_t tid = static_cast<uint64_t>(omp_get_thread_num());
      const uint64_t team = static_cast<uint64_t>(omp_get_num_threads());
      const uint64_t begin = num_edges * tid / team;
      const uint64_t end = num_edges * (tid + 1) / team;
      std::array<int64_t, kAlphabetSize + 1> w_rank{};
      for (uint8_t c = 0; c <= kAlphabetSize; ++c) {
        w_rank[c] = begin == 0
                        ? 0
                        : rs_w_.rank(c, static_cast<int64_t>(begin - 1));
      }
      const int64_t total_last =
          num_edges == 1
              ? static_cast<int64_t>(IsLast(0))
              : rs_last_.rank(static_cast<int64_t>(num_edges - 2)) +
                    static_cast<int64_t>(IsLast(num_edges - 1));
      const int64_t kUninitializedRank =
          std::numeric_limits<int64_t>::min();
      std::array<int64_t, kAlphabetSize + 1> cursor_rank;
      cursor_rank.fill(kUninitializedRank);
      std::array<uint64_t, kAlphabetSize + 1> cursor_edge{};

      auto select_last = [&](uint8_t c, int64_t target_rank) -> uint64_t {
        assert(target_rank >= 0 && target_rank < total_last);
        if (cursor_rank[c] == kUninitializedRank) {
          cursor_edge[c] = rs_last_.select(target_rank);
          cursor_rank[c] = target_rank;
          return cursor_edge[c];
        }
        assert(target_rank >= cursor_rank[c]);
        while (cursor_rank[c] < target_rank) {
          do {
            ++cursor_edge[c];
          } while (cursor_edge[c] < num_edges && !IsLast(cursor_edge[c]));
          assert(cursor_edge[c] < num_edges);
          ++cursor_rank[c];
        }
        return cursor_edge[c];
      };

      for (uint64_t edge_id = begin; edge_id < end; ++edge_id) {
        const uint8_t w = GetW(edge_id);
        const uint8_t c = w > kAlphabetSize ? w - kAlphabetSize : w;
        if (w == c) {
          ++w_rank[c];
        }
        const int64_t select_rank = rank_f_[c] + w_rank[c] - 1;
        uint64_t forward;
        if (edge_id + 1 == num_edges && w == 0) {
          // Preserve the historical all-edge API exactly: rank(0,size-1)
          // includes zero padding in the final packed word.  W==0 edges are
          // invalid for assembly, but callers may still query them.
          forward = ForwardUncached(edge_id);
        } else if (select_rank >= 0 && select_rank < total_last) {
          forward = select_last(c, select_rank);
        } else if (select_rank == total_last) {
          forward = num_edges;
        } else {
          forward = kNullID;
        }
        SetCachedForward(edge_id, forward);
      }
    }
  }

  // Backward(edge) is monotone inside each first-character interval.  Encode
  // that monotone sequence in independent blocks: low bits are bit-packed,
  // while the high parts are represented in unary as
  //   (value >> low_bits) - block_base + item_index.
  // One select in the small block-local unary bitmap recovers the high part.
  // Blocks are word-aligned so construction and queries are parallel without
  // atomics.  This keeps exact O(1) lookup but uses roughly four to five bits
  // per edge for typical BOSS character intervals instead of 32 bits.
  struct CompactBackwardGroup {
    uint64_t edge_begin{0};
    uint64_t edge_end{0};
    uint64_t valid_begin{0};
    uint64_t char_count{0};
    uint32_t block_begin{0};
    uint8_t low_bits{0};
  };

  struct CompactBackwardBlock {
    uint32_t base_high{0};
    uint32_t high_word_offset{0};
    uint32_t low_word_offset{0};
  };

  uint64_t CompactBackwardTopologyBytes() const {
    return compact_backward_blocks_.size() *
               sizeof(compact_backward_blocks_[0]) +
           compact_backward_high_words_ * sizeof(uint64_t) +
           compact_backward_low_words_ * sizeof(uint64_t);
  }

  unsigned CompactBackwardGroupForBlock(uint32_t block_id) const {
    unsigned c = 1;
    while (c < kAlphabetSize &&
           block_id >= compact_backward_groups_[c + 1u].block_begin) {
      ++c;
    }
    return c;
  }

  bool BuildCompactBackwardTopologyCache(uint64_t host_mem) {
    if (compact_backward_ready_) {
      return true;
    }
    const uint64_t num_edges = size();
    if (num_edges == 0 ||
        num_edges >= static_cast<uint64_t>(
                         std::numeric_limits<uint32_t>::max()) ||
        host_mem == 0 ||
        std::getenv("MEGAHIT_ENABLE_COMPACT_BACKWARD_TOPOLOGY_CACHE") ==
            nullptr) {
      return false;
    }

    uint64_t total_blocks = 0;
    const uint64_t universe = num_edges + 1u;
    for (unsigned c = 0; c <= kAlphabetSize; ++c) {
      CompactBackwardGroup &group = compact_backward_groups_[c];
      const int64_t signed_begin = c == 0 ? 0 : f_[c];
      const int64_t signed_end = c == kAlphabetSize ? f_[c + 1u]
                                                    : f_[c + 1u];
      group.edge_begin = static_cast<uint64_t>(std::max<int64_t>(
          0, std::min<int64_t>(signed_begin, num_edges)));
      group.edge_end = static_cast<uint64_t>(std::max<int64_t>(
          0, std::min<int64_t>(signed_end, num_edges)));
      if (group.edge_end < group.edge_begin) {
        group.edge_end = group.edge_begin;
      }

      group.valid_begin = group.edge_begin;
      if (group.edge_begin < group.edge_end && rank_f_[c] > 0) {
        const int64_t last_before_first_valid =
            rs_last_.select(rank_f_[c] - 1);
        if (last_before_first_valid >= 0) {
          group.valid_begin = std::max<uint64_t>(
              group.valid_begin,
              static_cast<uint64_t>(last_before_first_valid) + 1u);
        }
      }
      group.valid_begin = std::min(group.valid_begin, group.edge_end);
      group.char_count =
          group.edge_begin < group.edge_end
              ? static_cast<uint64_t>(rs_w_.rank(c, num_edges - 1u))
              : 0;
      group.block_begin = static_cast<uint32_t>(total_blocks);

      const uint64_t count = group.edge_end - group.valid_begin;
      unsigned low_bits = 0;
      if (count != 0) {
        const uint64_t ratio = std::max<uint64_t>(1u, universe / count);
        while (low_bits < 31u &&
               (uint64_t{1} << (low_bits + 1u)) <= ratio) {
          ++low_bits;
        }
      }
      group.low_bits = static_cast<uint8_t>(low_bits);
      total_blocks +=
          (count + kCompactBackwardBlockSize - 1u) /
          kCompactBackwardBlockSize;
      if (total_blocks > std::numeric_limits<uint32_t>::max()) {
        return false;
      }
    }
    // Sentinel block boundary makes block->character routing branch-light.
    compact_backward_groups_[kAlphabetSize + 1u].block_begin =
        static_cast<uint32_t>(total_blocks);

    try {
      compact_backward_blocks_.resize(static_cast<size_t>(total_blocks));
    } catch (const std::bad_alloc &) {
      compact_backward_blocks_.clear();
      return false;
    }

    // Each block is independent.  Two exact endpoint queries determine its
    // private high/low storage sizes without materializing an edge-sized
    // scratch array.
#pragma omp parallel for schedule(static)
    for (int64_t signed_block = 0;
         signed_block < static_cast<int64_t>(total_blocks); ++signed_block) {
      const uint32_t block_id = static_cast<uint32_t>(signed_block);
      const unsigned c = CompactBackwardGroupForBlock(block_id);
      const CompactBackwardGroup &group = compact_backward_groups_[c];
      const uint64_t local_block = block_id - group.block_begin;
      const uint64_t edge_begin =
          group.valid_begin + local_block * kCompactBackwardBlockSize;
      const uint64_t edge_end = std::min<uint64_t>(
          group.edge_end, edge_begin + kCompactBackwardBlockSize);
      const uint64_t count = edge_end - edge_begin;
      const uint64_t first = BackwardUncached(edge_begin);
      const uint64_t last = BackwardUncached(edge_end - 1u);
      assert(first != kNullID && last != kNullID && last >= first &&
             last <= num_edges);
      const unsigned low_bits = group.low_bits;
      const uint64_t base_high = first >> low_bits;
      const uint64_t last_high = last >> low_bits;
      assert(base_high <= std::numeric_limits<uint32_t>::max());
      CompactBackwardBlock &block = compact_backward_blocks_[block_id];
      block.base_high = static_cast<uint32_t>(base_high);
      // Temporarily hold word counts; the serial prefix sum below replaces
      // them with final offsets.
      block.high_word_offset = static_cast<uint32_t>(
          (count + last_high - base_high + 63u) / 64u);
      block.low_word_offset = static_cast<uint32_t>(
          (count * low_bits + 63u) / 64u);
    }

    uint64_t high_words = 0;
    uint64_t low_words = 0;
    for (CompactBackwardBlock &block : compact_backward_blocks_) {
      const uint32_t high_count = block.high_word_offset;
      const uint32_t low_count = block.low_word_offset;
      if (high_words > std::numeric_limits<uint32_t>::max() - high_count ||
          low_words > std::numeric_limits<uint32_t>::max() - low_count) {
        compact_backward_blocks_.clear();
        return false;
      }
      block.high_word_offset = static_cast<uint32_t>(high_words);
      block.low_word_offset = static_cast<uint32_t>(low_words);
      high_words += high_count;
      low_words += low_count;
    }

    const uint64_t cache_bytes =
        compact_backward_blocks_.size() *
            sizeof(compact_backward_blocks_[0]) +
        (high_words + low_words) * sizeof(uint64_t);
    // Unlike the wide two-direction cache, this representation is useful
    // under a bounded-memory run.  Still retain three quarters of the caller's
    // budget for the SDBG, unitigs and output structures.
    if (cache_bytes > host_mem / 4u) {
      compact_backward_blocks_.clear();
      return false;
    }

    try {
      compact_backward_high_.reset(new uint64_t[high_words]);
      compact_backward_low_.reset(new uint64_t[low_words]);
    } catch (const std::bad_alloc &) {
      compact_backward_blocks_.clear();
      compact_backward_high_.reset();
      compact_backward_low_.reset();
      return false;
    }
    compact_backward_high_words_ = high_words;
    compact_backward_low_words_ = low_words;

#pragma omp parallel for schedule(static)
    for (int64_t signed_block = 0;
         signed_block < static_cast<int64_t>(total_blocks); ++signed_block) {
      const uint32_t block_id = static_cast<uint32_t>(signed_block);
      const unsigned c = CompactBackwardGroupForBlock(block_id);
      const CompactBackwardGroup &group = compact_backward_groups_[c];
      const uint64_t local_block = block_id - group.block_begin;
      const uint64_t edge_begin =
          group.valid_begin + local_block * kCompactBackwardBlockSize;
      const uint64_t edge_end = std::min<uint64_t>(
          group.edge_end, edge_begin + kCompactBackwardBlockSize);
      const CompactBackwardBlock &block = compact_backward_blocks_[block_id];
      const uint32_t next_high_offset =
          block_id + 1u < total_blocks
              ? compact_backward_blocks_[block_id + 1u].high_word_offset
              : static_cast<uint32_t>(high_words);
      const uint32_t next_low_offset =
          block_id + 1u < total_blocks
              ? compact_backward_blocks_[block_id + 1u].low_word_offset
              : static_cast<uint32_t>(low_words);
      uint64_t *high = compact_backward_high_.get() + block.high_word_offset;
      uint64_t *low = compact_backward_low_.get() + block.low_word_offset;
      std::fill(high, high + next_high_offset - block.high_word_offset, 0);
      std::fill(low, low + next_low_offset - block.low_word_offset, 0);

      int64_t target_rank =
          (edge_begin == 0
               ? 0
               : rs_last_.rank(static_cast<int64_t>(edge_begin - 1u))) -
          rank_f_[c];
      assert(target_rank >= 0 &&
             static_cast<uint64_t>(target_rank) <= group.char_count);
      uint64_t cursor =
          static_cast<uint64_t>(target_rank) == group.char_count
              ? num_edges
              : static_cast<uint64_t>(rs_w_.select(c, target_rank));
      const unsigned low_bits = group.low_bits;
      const uint64_t low_mask =
          low_bits == 0 ? 0 : (uint64_t{1} << low_bits) - 1u;

      for (uint64_t i = 0, edge = edge_begin; edge < edge_end;
           ++i, ++edge) {
        const uint64_t value = cursor;
        const uint64_t high_position =
            (value >> low_bits) - block.base_high + i;
        high[high_position >> 6u] |= uint64_t{1} << (high_position & 63u);
        if (low_bits != 0) {
          const uint64_t bit_position = i * low_bits;
          const uint64_t low_value = value & low_mask;
          low[bit_position >> 6u] |= low_value << (bit_position & 63u);
          if ((bit_position & 63u) + low_bits > 64u) {
            low[(bit_position >> 6u) + 1u] |=
                low_value >> (64u - (bit_position & 63u));
          }
        }

        if (IsLast(edge)) {
          ++target_rank;
          assert(static_cast<uint64_t>(target_rank) <= group.char_count);
          if (static_cast<uint64_t>(target_rank) == group.char_count) {
            cursor = num_edges;
          } else {
            do {
              ++cursor;
            } while (cursor < num_edges && GetW(cursor) != c);
            assert(cursor < num_edges);
          }
        }
      }
    }

    compact_backward_ready_ = true;
    return true;
  }

  uint64_t CompactCachedBackward(uint64_t edge_id) const {
    const unsigned c = LastCharOf(edge_id);
    assert(c <= kAlphabetSize);
    const CompactBackwardGroup &group = compact_backward_groups_[c];
    if (edge_id < group.valid_begin) {
      return kNullID;
    }
    assert(edge_id < group.edge_end);
    const uint64_t local = edge_id - group.valid_begin;
    const uint32_t block_id = static_cast<uint32_t>(
        group.block_begin + local / kCompactBackwardBlockSize);
    const uint32_t in_block =
        static_cast<uint32_t>(local % kCompactBackwardBlockSize);
    const CompactBackwardBlock &block = compact_backward_blocks_[block_id];

    const uint64_t *high =
        compact_backward_high_.get() + block.high_word_offset;
    uint32_t remaining = in_block;
    uint64_t high_word_index = 0;
    uint64_t word = high[0];
    for (;;) {
      const uint32_t ones = static_cast<uint32_t>(__builtin_popcountll(word));
      if (remaining < ones) {
        break;
      }
      remaining -= ones;
      word = high[++high_word_index];
    }
    while (remaining-- != 0u) {
      word &= word - 1u;
    }
    const uint64_t selected_high_position =
        high_word_index * 64u + static_cast<uint64_t>(__builtin_ctzll(word));
    const uint64_t high_value =
        static_cast<uint64_t>(block.base_high) +
        selected_high_position - in_block;

    uint64_t low_value = 0;
    const unsigned low_bits = group.low_bits;
    if (low_bits != 0) {
      const uint64_t bit_position =
          static_cast<uint64_t>(in_block) * low_bits;
      const uint64_t *low =
          compact_backward_low_.get() + block.low_word_offset;
      low_value = low[bit_position >> 6u] >> (bit_position & 63u);
      if ((bit_position & 63u) + low_bits > 64u) {
        low_value |= low[(bit_position >> 6u) + 1u]
                     << (64u - (bit_position & 63u));
      }
      low_value &= (uint64_t{1} << low_bits) - 1u;
    }
    const uint64_t value = (high_value << low_bits) | low_value;
    assert(value <= size());
    return value;
  }

  void BuildReverseLookupSamples(uint64_t host_mem) {
    if (!reverse_lookup_samples_.empty() || prefix_look_up_.empty() ||
        std::getenv("MEGAHIT_DISABLE_SDBG_REVERSE_LOOKUP_SAMPLES") !=
            nullptr) {
      return;
    }

    unsigned prefix_bases = 0;
    uint64_t prefix_capacity = 1;
    while (prefix_capacity < prefix_look_up_.size()) {
      if (prefix_capacity >
          std::numeric_limits<uint64_t>::max() / kAlphabetSize) {
        return;
      }
      prefix_capacity *= kAlphabetSize;
      ++prefix_bases;
    }
    if (prefix_bases >= k_) {
      return;
    }
    const unsigned key_bases =
        std::min<unsigned>(kMaxReverseLookupKeyBases, k_ - prefix_bases);

    std::vector<uint64_t> offsets(prefix_look_up_.size() + 1u, 0);
    for (size_t bucket = 0; bucket < prefix_look_up_.size(); ++bucket) {
      const int64_t first = prefix_look_up_[bucket].first;
      const int64_t last = prefix_look_up_[bucket].second;
      uint64_t count = 0;
      if (first >= 0 && first <= last) {
        count = (static_cast<uint64_t>(last - first) +
                 kReverseLookupSampleStride) /
                kReverseLookupSampleStride;
      }
      if (offsets[bucket] >
          std::numeric_limits<uint64_t>::max() - count) {
        return;
      }
      offsets[bucket + 1u] = offsets[bucket] + count;
    }

    const uint64_t num_samples = offsets.back();
    if (num_samples >
        std::numeric_limits<size_t>::max() / sizeof(uint32_t)) {
      return;
    }
    const uint64_t sample_bytes =
        num_samples * sizeof(uint32_t) +
        offsets.size() * sizeof(offsets[0]);
    const uint64_t transform_bytes =
        compact_topology_cache_
            ? size() * kCompactTopologyBytesPerEdge
            : (forward_cache_ ? size() * 2u * sizeof(uint32_t)
                              : CompactBackwardTopologyBytes());
    if (host_mem != 0) {
      if (transform_bytes > host_mem / 8u ||
          sample_bytes > host_mem / 8u - transform_bytes) {
        return;
      }
    }

    std::vector<uint32_t> samples;
    try {
      samples.resize(static_cast<size_t>(num_samples));
    } catch (const std::bad_alloc &) {
      return;
    }

    const bool full_sample_labels =
        std::getenv("MEGAHIT_DISABLE_PARTIAL_REVERSE_SAMPLE_LABELS") !=
        nullptr;
    const bool validate_partial_sample_labels =
        std::getenv("MEGAHIT_VALIDATE_PARTIAL_REVERSE_SAMPLE_LABELS") !=
        nullptr;

#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t bucket = 0;
         bucket < static_cast<int64_t>(prefix_look_up_.size()); ++bucket) {
      const int64_t first = prefix_look_up_[bucket].first;
      const uint64_t out_begin = offsets[bucket];
      const uint64_t out_end = offsets[bucket + 1u];
      uint8_t label[kMaxK];
      for (uint64_t out = out_begin; out < out_end; ++out) {
        const uint64_t sample_in_bucket = out - out_begin;
        const uint64_t edge_id =
            static_cast<uint64_t>(first) +
            sample_in_bucket * kReverseLookupSampleStride;
        if (full_sample_labels) {
          GetLabel(edge_id, label);
        } else {
          // The lookup directory stores only the bases immediately preceding
          // the already-known prefix.  GetLabel used to reconstruct all k
          // bases for every sample even though the remaining prefix and the
          // earlier suffix were discarded.  Stop the identical Backward walk
          // at the lowest base that contributes to the key.  Tip-label copies
          // use the same offsets as GetLabel, restricted to that suffix.
          const int lowest = static_cast<int>(
              k_ - prefix_bases - key_bases);
          uint64_t x = edge_id;
          for (int i = static_cast<int>(k_) - 1; i >= lowest; --i) {
            if (IsTip(x)) {
              const label_word_t *tip_label = TipLabelStartPtr(x);
              for (int j = 0; j <= i - lowest; ++j) {
                label[i - j] = CharAtTipLabel(tip_label, j);
              }
              break;
            }
            x = Backward(x);
            uint8_t c = GetW(x);
            if (c > kAlphabetSize) c -= kAlphabetSize;
            label[i] = c;
          }
        }
        uint32_t key = 0;
        for (unsigned i = 0; i < key_bases; ++i) {
          key = (key << kBitsPerChar) |
                static_cast<uint32_t>(
                    label[k_ - 1u - prefix_bases - i] - 1u);
        }
        if (validate_partial_sample_labels && !full_sample_labels) {
          uint8_t full_label[kMaxK];
          GetLabel(edge_id, full_label);
          uint32_t reference = 0;
          for (unsigned i = 0; i < key_bases; ++i) {
            reference =
                (reference << kBitsPerChar) |
                static_cast<uint32_t>(
                    full_label[k_ - 1u - prefix_bases - i] - 1u);
          }
          if (reference != key) std::abort();
        }
        samples[out] = key;
      }
#ifndef NDEBUG
      assert(std::is_sorted(samples.begin() + out_begin,
                            samples.begin() + out_end));
#endif
    }

    reverse_lookup_prefix_bases_ = prefix_bases;
    reverse_lookup_key_bases_ = key_bases;
    reverse_lookup_samples_.swap(samples);
    reverse_lookup_sample_offsets_.swap(offsets);
  }

  void SetCachedForward(uint64_t edge_id, uint64_t value) {
    if (compact_topology_cache_) {
      const uint32_t cached =
          value < kCompactTopologyNull ? static_cast<uint32_t>(value)
                                       : kCompactTopologyNull;
      uint8_t *record = compact_topology_cache_.get() +
                        edge_id * kCompactTopologyBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      std::memcpy(record, &cached, sizeof(cached));
#else
      record[0] = static_cast<uint8_t>(cached);
      record[1] = static_cast<uint8_t>(cached >> 8u);
      record[2] = static_cast<uint8_t>(cached >> 16u);
      // Backward occupies the high nibble and is filled in the next phase.
      record[3] = static_cast<uint8_t>((cached >> 24u) & 0x0Fu);
#endif
    } else {
      forward_cache_[edge_id] =
          value < std::numeric_limits<uint32_t>::max()
              ? static_cast<uint32_t>(value)
              : std::numeric_limits<uint32_t>::max();
    }
  }

  void SetCachedBackward(uint64_t edge_id, uint64_t value) {
    if (compact_topology_cache_) {
      const uint32_t cached =
          value < kCompactTopologyNull ? static_cast<uint32_t>(value)
                                       : kCompactTopologyNull;
      uint8_t *record = compact_topology_cache_.get() +
                        edge_id * kCompactTopologyBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      const uint32_t encoded =
          (cached << 4u) | static_cast<uint32_t>(record[3] & 0x0Fu);
      std::memcpy(record + 3, &encoded, sizeof(encoded));
#else
      record[3] = static_cast<uint8_t>((record[3] & 0x0Fu) |
                                       ((cached & 0x0Fu) << 4u));
      record[4] = static_cast<uint8_t>(cached >> 4u);
      record[5] = static_cast<uint8_t>(cached >> 12u);
      record[6] = static_cast<uint8_t>(cached >> 20u);
#endif
    } else {
      backward_cache_[edge_id] =
          value < std::numeric_limits<uint32_t>::max()
              ? static_cast<uint32_t>(value)
              : std::numeric_limits<uint32_t>::max();
    }
  }

  uint32_t PackedForward(uint64_t edge_id) const {
    const uint8_t *record = compact_topology_cache_.get() +
                            edge_id * kCompactTopologyBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    uint32_t packed;
    std::memcpy(&packed, record, sizeof(packed));
    return packed & kCompactTopologyNull;
#else
    return static_cast<uint32_t>(record[0]) |
           (static_cast<uint32_t>(record[1]) << 8u) |
           (static_cast<uint32_t>(record[2]) << 16u) |
           (static_cast<uint32_t>(record[3] & 0x0Fu) << 24u);
#endif
  }

  uint32_t PackedBackward(uint64_t edge_id) const {
    const uint8_t *record = compact_topology_cache_.get() +
                            edge_id * kCompactTopologyBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    uint32_t packed;
    std::memcpy(&packed, record + 3, sizeof(packed));
    return packed >> 4u;
#else
    return static_cast<uint32_t>(record[3] >> 4u) |
           (static_cast<uint32_t>(record[4]) << 4u) |
           (static_cast<uint32_t>(record[5]) << 12u) |
           (static_cast<uint32_t>(record[6]) << 20u);
#endif
  }

  uint32_t k_{};
  SdbgRawContent content_;
  kmlib::AtomicBitVector<uint64_t> invalid_;
  std::vector<std::pair<int64_t, int64_t>> prefix_look_up_;
  int64_t f_[kAlphabetSize + 2]{};
  int64_t rank_f_[kAlphabetSize + 2]{};  // = rs_last_.Rank(f_[i] - 1)
  kmlib::RankAndSelect<kAlphabetSize, kWAlphabetSize> rs_w_;
  kmlib::RankAndSelect<1, 2> rs_last_;
  kmlib::RankAndSelect<1, 2, kmlib::rnsmode::kRankOnly> rs_is_tip_;
  std::unique_ptr<uint8_t[]> compact_topology_cache_;
  std::unique_ptr<uint32_t[]> forward_cache_;
  std::unique_ptr<uint32_t[]> backward_cache_;
  std::array<CompactBackwardGroup, kAlphabetSize + 2>
      compact_backward_groups_{};
  std::vector<CompactBackwardBlock> compact_backward_blocks_;
  std::unique_ptr<uint64_t[]> compact_backward_high_;
  std::unique_ptr<uint64_t[]> compact_backward_low_;
  uint64_t compact_backward_high_words_{0};
  uint64_t compact_backward_low_words_{0};
  bool compact_backward_ready_{false};
  uint64_t host_memory_budget_{0};
  bool forward_cache_borrowed_{false};
  std::unique_ptr<uint8_t[]> simple_path_codes_;
  std::unique_ptr<uint32_t[]> simple_path_bases_;
  uint64_t simple_path_edge_count_{0};
  uint64_t simple_path_overflow_count_{0};
  std::vector<uint32_t> reverse_lookup_samples_;
  std::vector<uint64_t> reverse_lookup_sample_offsets_;
  unsigned reverse_lookup_prefix_bases_{0};
  unsigned reverse_lookup_key_bases_{0};
};

#endif  // MEGAHIT_SDBG_H
