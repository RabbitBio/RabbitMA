//
// Created by vout on 11/24/18.
//

#ifndef MEGAHIT_SPANNING_KMER_COLLECTOR_H
#define MEGAHIT_SPANNING_KMER_COLLECTOR_H

#include <algorithm>
#include <omp.h>

#include "kmlib/kmsort.h"
#include "parallel_hashmap/phmap.h"
#include "sdbg/sdbg_def.h"
#include "sequence/io/edge/edge_writer.h"
#include "sequence/kmer_plus.h"
#include "utils/mutex.h"

template <class KmerType>
class KmerCollector {
 public:
  using kmer_type = KmerType;
  using kmer_plus = KmerPlus<KmerType, mul_t>;
  using hash_set = phmap::parallel_flat_hash_set<
      kmer_plus, KmerHash,
      phmap::container_internal::hash_default_eq<kmer_plus>,
      phmap::container_internal::Allocator<kmer_plus>, 12, SpinLock>;
  using shard_set = phmap::flat_hash_set<
      kmer_plus, KmerHash,
      phmap::container_internal::hash_default_eq<kmer_plus>,
      phmap::container_internal::Allocator<kmer_plus>>;

  KmerCollector(unsigned k, const std::string &out_prefix,
                bool batch_partitioned = false)
      : k_(k),
        output_prefix_(out_prefix),
        batch_partitioned_(batch_partitioned) {
    last_shift_ = k_ % 16;
    last_shift_ = (last_shift_ == 0 ? 0 : 16 - last_shift_) * 2;
    words_per_kmer_ = DivCeiling(k_ * 2 + kBitsPerMul, 32);
    buffer_.resize(words_per_kmer_);

    writer_.SetFilePrefix(out_prefix);
    writer_.SetUnordered();
    writer_.SetKmerSize(k_ - 1);
    writer_.InitFiles();

    if (batch_partitioned_) {
      batch_buffers_.resize(std::max(1, omp_get_max_threads()));
      shard_sets_.resize(kNumShards);
    }
  }

  void Insert(const KmerType &kmer, mul_t mul) {
    if (batch_partitioned_) {
      const int tid = omp_in_parallel() ? omp_get_thread_num() : 0;
      assert(tid >= 0 && tid < static_cast<int>(batch_buffers_.size()));
      batch_buffers_[tid].emplace_back(kmer, mul);
      return;
    }
    collection_.insert({kmer, mul});
  }

  const hash_set &collection() const { return collection_; }
  size_t size() const {
    return batch_partitioned_ ? unique_count_ : collection_.size();
  }
  bool empty() const { return size() == 0; }

  // Move one input batch through a stable prefix radix before touching the
  // persistent exact sets.  Every shard is owned by the same OpenMP worker on
  // each call (static scheduling), so its table has no lock and retains NUMA
  // and cache locality across read batches.
  uint64_t CommitBatch() {
    if (!batch_partitioned_) return 0u;
    const size_t num_threads = batch_buffers_.size();
    std::vector<uint64_t> counts(num_threads * kNumShards, 0);

#pragma omp parallel for schedule(static)
    for (int64_t tid = 0; tid < static_cast<int64_t>(num_threads); ++tid) {
      uint64_t *local_counts =
          counts.data() + static_cast<size_t>(tid) * kNumShards;
      for (const auto &item : batch_buffers_[tid]) {
        ++local_counts[ShardOf(item.kmer)];
      }
    }

    std::vector<uint64_t> shard_begin(kNumShards + 1u, 0);
    std::vector<uint64_t> cursors(num_threads * kNumShards, 0);
    for (size_t shard = 0; shard < kNumShards; ++shard) {
      uint64_t cursor = shard_begin[shard];
      for (size_t tid = 0; tid < num_threads; ++tid) {
        const size_t cell = tid * kNumShards + shard;
        cursors[cell] = cursor;
        cursor += counts[cell];
      }
      shard_begin[shard + 1u] = cursor;
    }
    const uint64_t batch_size = shard_begin.back();
    partitioned_batch_.resize(static_cast<size_t>(batch_size));

#pragma omp parallel for schedule(static)
    for (int64_t tid = 0; tid < static_cast<int64_t>(num_threads); ++tid) {
      uint64_t *local_cursors =
          cursors.data() + static_cast<size_t>(tid) * kNumShards;
      for (const auto &item : batch_buffers_[tid]) {
        partitioned_batch_[local_cursors[ShardOf(item.kmer)]++] = item;
      }
      batch_buffers_[tid].clear();
    }

    uint64_t added = 0;
    uint64_t batch_unique = 0;
#pragma omp parallel for schedule(static) reduction(+ : added, batch_unique)
    for (int64_t shard = 0; shard < static_cast<int64_t>(kNumShards);
         ++shard) {
      shard_set &set = shard_sets_[shard];
      const size_t before = set.size();
      auto first = partitioned_batch_.begin() + shard_begin[shard];
      auto last = partitioned_batch_.begin() + shard_begin[shard + 1u];
      // A large first-iteration batch contains tens of millions of generated
      // edges but often only a small fraction of distinct k-mers. Collapse
      // repeats in contiguous memory before probing the much larger
      // persistent table.
      kmlib::kmsort(first, last);
      last = std::unique(first, last);
      batch_unique += static_cast<uint64_t>(last - first);
      for (auto item = first; item != last; ++item) {
        set.insert(*item);
      }
      added += set.size() - before;
    }
    unique_count_ += added;
    partitioned_candidates_ += batch_size;
    partitioned_unique_candidates_ += batch_unique;
    return batch_size;
  }

  void FlushToFile() {
    std::vector<kmer_plus>().swap(partitioned_batch_);
    if (batch_partitioned_) {
      // Convert a bounded wave of prefix shards in parallel, then issue one
      // large sequential write per shard in the original prefix order.  This
      // removes per-edge ostream calls and serial base packing without
      // perturbing the prefix locality consumed by seq2sdbg.  Only one shard
      // per worker is resident, so temporary memory is O(output / 4096 *
      // workers), not O(total output).
      const int workers = static_cast<int>(batch_buffers_.size());
      std::vector<std::vector<uint32_t>> packed_shards(workers);
      for (size_t first = 0; first < kNumShards;
           first += static_cast<size_t>(workers)) {
        const int wave = static_cast<int>(std::min<size_t>(
            workers, kNumShards - first));
#pragma omp parallel for schedule(static) num_threads(workers)
        for (int lane = 0; lane < wave; ++lane) {
          const shard_set &set = shard_sets_[first + lane];
          std::vector<uint32_t> &packed = packed_shards[lane];
          packed.resize(set.size() * words_per_kmer_);
          size_t edge = 0;
          for (const auto &item : set) {
            PackEdge(item.kmer, item.aux,
                     packed.data() + edge * words_per_kmer_);
            ++edge;
          }
        }
        for (int lane = 0; lane < wave; ++lane) {
          const std::vector<uint32_t> &packed = packed_shards[lane];
          writer_.WriteUnorderedBatch(
              packed.data(), packed.size() / words_per_kmer_);
        }
      }
    } else {
      for (const auto &item : collection_) {
        PackEdge(item.kmer, item.aux, buffer_.data());
        writer_.WriteUnordered(buffer_.data());
      }
    }
  }

  uint64_t partitioned_candidates() const { return partitioned_candidates_; }
  uint64_t partitioned_unique_candidates() const {
    return partitioned_unique_candidates_;
  }

 private:
  static constexpr size_t kNumShards = 1u << 12u;

  static size_t ShardOf(const KmerType &kmer) {
    static_assert(KmerType::kBitsPerWord >= 12u,
                  "iterate radix requires at least 12 key bits");
    return static_cast<size_t>(
        kmer.data()[0] >> (KmerType::kBitsPerWord - 12u));
  }

  void PackEdge(const KmerType &kmer, mul_t mul, uint32_t *buffer) const {
    std::fill(buffer, buffer + words_per_kmer_, 0u);

    uint32_t *ptr = buffer;
    uint32_t w = 0;

    for (unsigned j = 0; j < k_; ++j) {
      w = (w << 2) | kmer.GetBase(k_ - 1 - j);
      if (j % 16 == 15) {
        *ptr = w;
        w = 0;
        ++ptr;
      }
    }

    assert(static_cast<unsigned>(ptr - buffer) < words_per_kmer_);
    *ptr = (w << last_shift_);
    assert((buffer[words_per_kmer_ - 1u] & kMaxMul) == 0);
    buffer[words_per_kmer_ - 1u] |= mul;
  }

 private:
  unsigned k_;
  std::string output_prefix_;
  hash_set collection_;
  EdgeWriter writer_;
  unsigned last_shift_;
  unsigned words_per_kmer_;
  std::vector<uint32_t> buffer_;
  bool batch_partitioned_{false};
  std::vector<std::vector<kmer_plus>> batch_buffers_;
  std::vector<kmer_plus> partitioned_batch_;
  std::vector<shard_set> shard_sets_;
  size_t unique_count_{0};
  uint64_t partitioned_candidates_{0};
  uint64_t partitioned_unique_candidates_{0};
};

#endif  // MEGAHIT_SPANNING_KMER_COLLECTOR_H
