//
// Created by vout on 7/13/19.
//

#ifndef MEGAHIT_LOCALASM_MAPPING_RESULT_H
#define MEGAHIT_LOCALASM_MAPPING_RESULT_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

#include "hash_mapper.h"
#include "utils/startup_affinity.h"

/**
 * Collect local-assembly mappings without a contig lock and a heap-allocated
 * deque block on the mapping hot path.
 *
 * Workers first append to cache-local vectors.  Finalize groups those entries
 * by `(contig, end)`, writes one compact CSR value array, and sorts each final
 * range by the historical 64-bit encoding.  GetMappingResults therefore sees
 * exactly the same ordered values as the old per-contig deque implementation,
 * while mapping itself performs no shared writes.
 */
class MappingResultCollector {
 public:
  class MappingRange {
   public:
    typedef std::vector<uint64_t>::const_iterator const_iterator;

    MappingRange(const_iterator begin, const_iterator end)
        : begin_(begin), end_(end) {}

    const_iterator begin() const { return begin_; }
    const_iterator end() const { return end_; }
    size_t size() const { return static_cast<size_t>(end_ - begin_); }

   private:
    const_iterator begin_;
    const_iterator end_;
  };

  explicit MappingResultCollector(size_t n_ref)
      : n_ref_(n_ref), thread_entries_(omp_get_max_threads()) {
    if (n_ref_ >
        (static_cast<size_t>(std::numeric_limits<uint32_t>::max()) + 1u) /
            2u) {
      throw std::length_error(
          "local mapping target index exceeds compact 32-bit range");
    }
    for (auto &buffer : thread_entries_) {
      buffer.entries.reserve(4096);
    }
  }

  unsigned AddSingle(const MappingRecord &rec, int32_t contig_len,
                     int32_t read_len, int32_t local_range) {
    if (rec.contig_to < local_range && rec.query_from != 0 &&
        rec.query_to == read_len - 1) {
      Add(rec.contig_id, 0,
          EncodeMappingRead(rec.contig_to, 0, rec.mismatch, rec.strand,
                            rec.query_id));
      return 1;
    }
    if (rec.contig_from + local_range >= contig_len &&
        rec.query_to < read_len - 1 && rec.query_from == 0) {
      Add(rec.contig_id, 1,
          EncodeMappingRead(contig_len - 1 - rec.contig_from, 0,
                            rec.mismatch, rec.strand, rec.query_id));
      return 1;
    }
    return 0;
  }

  unsigned AddMate(const MappingRecord &rec1, const MappingRecord &rec2,
                   int32_t contig_len, uint64_t mate_id,
                   int32_t local_range) {
    if (rec2.valid && rec2.contig_id == rec1.contig_id) {
      return 0;
    }
    if (rec1.contig_to < local_range && rec1.strand == 1) {
      Add(rec1.contig_id, 0,
          EncodeMappingRead(rec1.contig_to, 1, rec1.mismatch, rec1.strand,
                            mate_id));
      return 1;
    }
    if (rec1.contig_from + local_range >= contig_len && rec1.strand == 0) {
      Add(rec1.contig_id, 1,
          EncodeMappingRead(contig_len - 1 - rec1.contig_from, 1,
                            rec1.mismatch, rec1.strand, mate_id));
      return 1;
    }
    return 0;
  }

  void Finalize() {
    assert(!finalized_);
    const size_t num_targets = n_ref_ * 2u;
    std::vector<uint64_t> target_counts(num_targets, 0);

    // Sorting only by target is sufficient for grouped counting/copying.  The
    // final per-target sort below exactly restores the old observable order.
#pragma omp parallel for schedule(static)
    for (int thread_id = 0;
         thread_id < static_cast<int>(thread_entries_.size()); ++thread_id) {
      auto &entries = thread_entries_[thread_id].entries;
      std::sort(entries.begin(), entries.end(),
                [](const MappingEntry &lhs, const MappingEntry &rhs) {
                  return lhs.target < rhs.target;
                });
      size_t begin = 0;
      while (begin < entries.size()) {
        size_t end = begin + 1;
        while (end < entries.size() &&
               entries[end].target == entries[begin].target) {
          ++end;
        }
        __atomic_fetch_add(&target_counts[entries[begin].target], end - begin,
                           __ATOMIC_RELAXED);
        begin = end;
      }
    }

    offsets_.resize(num_targets + 1u);
    offsets_[0] = 0;
    std::vector<uint32_t> active_targets;
    active_targets.reserve(num_targets / 4u);
    for (size_t target = 0; target < num_targets; ++target) {
      if (target_counts[target] >
          std::numeric_limits<uint64_t>::max() - offsets_[target]) {
        throw std::length_error("local mapping result count overflow");
      }
      offsets_[target + 1] = offsets_[target] + target_counts[target];
      if (target_counts[target] > 1) {
        active_targets.push_back(static_cast<uint32_t>(target));
      }
    }
    if (offsets_.back() > std::numeric_limits<size_t>::max()) {
      throw std::length_error("local mapping result array exceeds size_t");
    }

    values_.resize(static_cast<size_t>(offsets_.back()));
    if (!values_.empty()) {
      const size_t bytes = values_.size() * sizeof(values_[0]);
      // resize() has already faulted and zeroed these pages.  Discarding them
      // after mbind makes the parallel copy pay one kernel fault per page and
      // was substantially slower than the subsequent remote-memory traffic
      // on large collectors.  Keep the resident pages; the copy below is a
      // single streaming pass and later endpoint reads are distributed.
      AdviseHugePages(values_.data(), bytes);
    }

    std::vector<uint64_t> write_cursor(offsets_.begin(), offsets_.end() - 1);
#pragma omp parallel for schedule(static)
    for (int thread_id = 0;
         thread_id < static_cast<int>(thread_entries_.size()); ++thread_id) {
      const auto &entries = thread_entries_[thread_id].entries;
      size_t begin = 0;
      while (begin < entries.size()) {
        size_t end = begin + 1;
        while (end < entries.size() &&
               entries[end].target == entries[begin].target) {
          ++end;
        }
        const uint32_t target = entries[begin].target;
        const uint64_t output_begin = __atomic_fetch_add(
            &write_cursor[target], end - begin, __ATOMIC_RELAXED);
        for (size_t i = begin; i < end; ++i) {
          values_[static_cast<size_t>(output_begin + i - begin)] =
              entries[i].encoded;
        }
        begin = end;
      }
    }

    // Release the wider temporary entries before sorting/assembly so the
    // steady-state collector is only 8 bytes per mapping plus CSR offsets.
#pragma omp parallel for schedule(static)
    for (int thread_id = 0;
         thread_id < static_cast<int>(thread_entries_.size()); ++thread_id) {
      std::vector<MappingEntry>().swap(thread_entries_[thread_id].entries);
    }
    std::vector<ThreadEntries>().swap(thread_entries_);

#pragma omp parallel for schedule(dynamic, 256)
    for (int64_t active_id = 0;
         active_id < static_cast<int64_t>(active_targets.size()); ++active_id) {
      const uint32_t target = active_targets[active_id];
      std::sort(values_.begin() + offsets_[target],
                values_.begin() + offsets_[target + 1]);
    }
    finalized_ = true;
  }

  MappingRange GetMappingResults(int64_t contig_id, uint8_t strand) const {
    assert(finalized_);
    assert(contig_id >= 0 && static_cast<size_t>(contig_id) < n_ref_);
    assert(strand <= 1);
    const size_t target = static_cast<size_t>(contig_id) * 2u + strand;
    return MappingRange(values_.begin() + offsets_[target],
                        values_.begin() + offsets_[target + 1]);
  }

  size_t size() const { return values_.size(); }

  /**
   * Replace the original library read ids in every retained mapping with a
   * dense, order-preserving id and return the corresponding original ids.
   *
   * The mapping encoding is sorted lexicographically and read id is its least
   * significant tie breaker.  Ranking the selected ids (rather than assigning
   * ids in discovery order) therefore leaves every finalized target range in
   * exactly the same order.  Local assembly can subsequently gather from a
   * compact package containing only referenced reads without changing the
   * historical "first three reads at a position" rule.
   */
  std::vector<uint64_t> CompactReadIds(size_t num_library_reads) {
    assert(finalized_);
    constexpr uint64_t kReadIdMask = (uint64_t{1} << 44u) - 1u;
    if (num_library_reads > kReadIdMask + 1u) {
      throw std::length_error("read library exceeds local mapping id range");
    }
    if (values_.empty()) {
      return {};
    }

    const size_t num_words = (num_library_reads + 63u) / 64u;
    std::vector<uint64_t> selected_words(num_words, 0);
    const size_t selected_bytes = num_words * sizeof(uint64_t);
    if (omp_get_max_threads() > 1 && selected_bytes != 0 &&
        InterleaveMemoryPages(selected_words.data(), selected_bytes)) {
      DiscardMemoryPages(selected_words.data(), selected_bytes);
    }

#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(values_.size()); ++i) {
      const uint64_t read_id = GetReadId(values_[i]);
      if (read_id >= num_library_reads) {
        xfatal("Invalid read id {} in local mapping result ({} reads)\n",
               read_id, num_library_reads);
      }
      __atomic_fetch_or(&selected_words[read_id >> 6u],
                        uint64_t{1} << (read_id & 63u), __ATOMIC_RELAXED);
    }

    // One rank sample per eight words keeps the directory small while a
    // lookup needs at most seven independent register popcounts.
    constexpr size_t kRankWords = 8u;
    const size_t num_rank_blocks =
        (num_words + kRankWords - 1u) / kRankWords;
    std::vector<uint64_t> rank_blocks(num_rank_blocks + 1u, 0);
    uint64_t num_selected = 0;
    for (size_t block = 0; block < num_rank_blocks; ++block) {
      rank_blocks[block] = num_selected;
      const size_t end = std::min(num_words, (block + 1u) * kRankWords);
      for (size_t word = block * kRankWords; word < end; ++word) {
        num_selected += static_cast<unsigned>(
            __builtin_popcountll(selected_words[word]));
      }
    }
    rank_blocks[num_rank_blocks] = num_selected;
    if (num_selected > kReadIdMask + 1u ||
        num_selected > std::numeric_limits<size_t>::max()) {
      throw std::length_error("compact local read package is too large");
    }

    std::vector<uint64_t> selected_ids(static_cast<size_t>(num_selected));
    size_t selected_id = 0;
    for (size_t word = 0; word < num_words; ++word) {
      uint64_t bits = selected_words[word];
      while (bits != 0) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
        selected_ids[selected_id++] = (word << 6u) + bit;
        bits &= bits - 1u;
      }
    }
    assert(selected_id == selected_ids.size());

#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(values_.size()); ++i) {
      const uint64_t encoded = values_[i];
      const uint64_t read_id = GetReadId(encoded);
      const size_t word = static_cast<size_t>(read_id >> 6u);
      const size_t block = word / kRankWords;
      uint64_t dense_id = rank_blocks[block];
      for (size_t w = block * kRankWords; w < word; ++w) {
        dense_id += static_cast<unsigned>(
            __builtin_popcountll(selected_words[w]));
      }
      const unsigned bit = static_cast<unsigned>(read_id & 63u);
      const uint64_t before =
          bit == 0u ? 0u : selected_words[word] & ((uint64_t{1} << bit) - 1u);
      dense_id += static_cast<unsigned>(__builtin_popcountll(before));
      assert(dense_id < selected_ids.size());
      values_[i] = (encoded & ~kReadIdMask) | dense_id;
    }

    return selected_ids;
  }

  static uint64_t GetContigAbsPos(uint64_t encoded) {
    return encoded >> (44u + 1u + 4u);
  }

  static uint64_t GetReadId(uint64_t encoded) {
    return encoded & ((1ull << 44u) - 1);
  }

 private:
  struct MappingEntry {
    uint64_t encoded;
    uint32_t target;
  };

  struct alignas(64) ThreadEntries {
    std::vector<MappingEntry> entries;
  };

  void Add(uint32_t contig_id, uint8_t strand, uint64_t encoded) {
    assert(!finalized_);
    assert(contig_id < n_ref_);
    assert(strand <= 1);
    const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
    assert(thread_id >= 0 &&
           static_cast<size_t>(thread_id) < thread_entries_.size());
    thread_entries_[thread_id].entries.push_back(
        {encoded, static_cast<uint32_t>(contig_id * 2u + strand)});
  }

  static uint64_t EncodeMappingRead(uint32_t contig_offset, uint8_t is_mate,
                                    uint32_t mismatch, uint8_t strand,
                                    uint64_t read_id) {
    assert(contig_offset <= (1u << 14u));
    assert(strand <= 1);
    uint64_t ret = contig_offset;
    ret = (ret << 1u) | is_mate;
    ret = (ret << 4u) | (mismatch < 15 ? mismatch : 15u);
    ret = (ret << 1u) | strand;
    ret = (ret << 44u) | read_id;  // 44 bits for read id
    return ret;
  }

  size_t n_ref_;
  std::vector<ThreadEntries> thread_entries_;
  std::vector<uint64_t> offsets_;
  std::vector<uint64_t> values_;
  bool finalized_{false};
};

#endif  // MEGAHIT_LOCALASM_MAPPING_RESULT_H
