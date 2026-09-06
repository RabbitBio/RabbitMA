/*
 * Profile the exact, reusable read-occurrence index used by iterative
 * assembly.  This first stage deliberately writes no index: it measures the
 * content-defined anchor density and posting distribution before committing
 * disk or memory to a particular representation.
 */

#include <fcntl.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <time.h>
#include <vector>

#include "iterate/contig_flank_index.h"
#include "iterate/kmer_collector.h"
#include "kmlib/kmsort.h"
#include "localasm/hash_mapper.h"
#include "parallel_hashmap/phmap.h"
#include "sequence/io/async_sequence_reader.h"
#include "sequence/io/local_candidate_index.h"
#include "sequence/io/read_anchor_positions.h"
#include "sequence/io/read_chunk_index.h"
#include "sequence/sequence_package.h"
#include "utils/mutex.h"
#include "utils/options_description.h"
#include "utils/startup_affinity.h"
#include "utils/utils.h"

namespace {

constexpr unsigned kHllBits = 14;
constexpr size_t kHllRegisters = size_t{1} << kHllBits;

uint64_t Mix64(uint64_t value) {
  value ^= value >> 30u;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27u;
  value *= UINT64_C(0x94d049bb133111eb);
  value ^= value >> 31u;
  return value;
}

uint64_t AnchorHash(uint64_t key) {
  // Salting is essential here: the SplitMix finalizer maps zero to zero.  An
  // unsalted minimizer would therefore select poly-A/zero keys whenever they
  // occur, creating artificial hot postings and badly skewed external-sort
  // buckets.  The fixed salt keeps construction and lookup deterministic.
  return ReadAnchorRankHash(key);
}

uint64_t AnchorBucketHash(uint64_t key) {
  // Bucket ownership must be independent of minimizer rank.  Reusing the
  // minimizer hash would put almost all selected minima into low-numbered
  // buckets by construction, defeating load balance.
  return Mix64(key ^ UINT64_C(0xd1b54a32d192ed03));
}

uint64_t LowBitMask(unsigned bits) {
  return bits == 0
             ? 0
             : (bits >= 64 ? std::numeric_limits<uint64_t>::max()
                           : (uint64_t{1} << bits) - 1u);
}

uint64_t PermuteReadIndexKey(uint64_t key, unsigned key_bits) {
  // A bijection in the exact 2*anchor_len-bit domain.  Unlike a fingerprint,
  // this loses no information: the high bits can be implicit in bucket
  // ownership while the remaining suffix is still collision-free.  Masked
  // odd multiplication and xor-shifts are individually invertible modulo
  // 2^key_bits; no inverse is needed because construction and lookup compare
  // the same encoded suffix.
  assert(key_bits > 0 && key_bits < 64u);
  const uint64_t mask = LowBitMask(key_bits);
  uint64_t value = (key ^ UINT64_C(0x2d358dccaa6c78a5)) & mask;
  value ^= value >> std::max(1u, key_bits / 3u);
  value = (value * UINT64_C(0x9e3779b97f4a7c15)) & mask;
  value ^= value >> std::max(1u, key_bits / 2u);
  value = (value * UINT64_C(0xbf58476d1ce4e5b9)) & mask;
  value ^= value >> std::max(1u, key_bits / 3u + 1u);
  return value & mask;
}

uint64_t ReadIndexCardinalityHash(uint64_t key) {
  // Independent from both minimizer rank and bucket ownership so HLL remains
  // unbiased even though selected minimizers have deliberately small rank
  // hashes.
  return Mix64(key ^ UINT64_C(0x94d049bb133111eb));
}

unsigned BitsNeeded(uint64_t max_value) {
  unsigned bits = 0;
  do {
    ++bits;
    max_value >>= 1u;
  } while (max_value != 0);
  return bits;
}

struct ReadIndexOptions {
  std::string read_file;
  std::string contig_file;
  std::string bubble_file;
  std::string cache_contig_file;
  std::string cache_bubble_file;
  std::string output_prefix;
  std::string index_prefix;
  std::string edge_output_prefix;
  std::string local_candidate_output;
  unsigned anchor_len{21};
  unsigned window_len{40};
  unsigned sample_bits{12};
  int kmer_k{0};
  int step{10};
  int cache_kmer_k{0};
  int cache_step{10};
  int local_endpoint_range{650};
  int local_seed_len{31};
  int local_sparsity{8};
  int local_min_contig_len{200};
  int num_threads{0};
  int max_reads{0};
  double memory_bytes{static_cast<double>(uint64_t{4} << 30u)};
};

class ActiveQuerySet {
  struct PendingQuery {
    uint64_t key;
    uint8_t offset;
  };
  static constexpr unsigned kQueryShardBits = 8u;
  static constexpr unsigned kNumQueryShards = 1u << kQueryShardBits;

 public:
  struct Batch {
    std::array<std::vector<PendingQuery>, kNumQueryShards> queries;
    std::array<std::vector<uint64_t>, kNumQueryShards> contexts;
  };

  struct QueryInfo {
    uint64_t offset_mask{0};
    uint64_t query_count{0};
  };

  void ConfigureContext(unsigned context_bits, unsigned anchor_len,
                        unsigned window_len, unsigned context_distance) {
    context_bits_ = context_bits;
    anchor_len_ = anchor_len;
    window_len_ = window_len;
    context_distance_ = context_distance;
  }

  unsigned context_bits() const { return context_bits_; }
  unsigned context_distance() const { return context_distance_; }

  void AddToBatch(Batch *batch, uint64_t key, unsigned offset,
                  uint8_t context) const {
    assert(offset < 64u);
    const unsigned shard = Shard(key);
    batch->queries[shard].push_back(
        PendingQuery{key, static_cast<uint8_t>(offset)});
    const uint8_t mask = ContextMask(offset);
    if (mask != 0u) {
      batch->contexts[shard].push_back(
          ContextHash(key, offset, context & mask, mask));
    }
  }

  void MergeBatches(std::vector<Batch> *batches) {
    // Each destination shard has one writer.  Counts and offset/context masks
    // are commutative; the existing exact sort/reduce establishes lookup order.
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned shard = 0; shard < kNumQueryShards; ++shard) {
      size_t query_count = pending_queries_[shard].size();
      size_t context_count = context_hashes_[shard].size();
      for (const Batch &batch : *batches) {
        query_count += batch.queries[shard].size();
        context_count += batch.contexts[shard].size();
      }
      auto &queries = pending_queries_[shard];
      auto &contexts = context_hashes_[shard];
      queries.reserve(query_count);
      contexts.reserve(context_count);
      for (Batch &batch : *batches) {
        queries.insert(queries.end(), batch.queries[shard].begin(),
                       batch.queries[shard].end());
        contexts.insert(contexts.end(), batch.contexts[shard].begin(),
                        batch.contexts[shard].end());
        std::vector<PendingQuery>().swap(batch.queries[shard]);
        std::vector<uint64_t>().swap(batch.contexts[shard]);
      }
    }
  }

  void Add(uint64_t key, unsigned offset, uint8_t context = 0u) {
    assert(offset < 64u);
    const unsigned shard = Shard(key);
    pending_queries_[shard].push_back(
        PendingQuery{key, static_cast<uint8_t>(offset)});
    const uint8_t mask = ContextMask(offset);
    if (mask != 0u) {
      context_hashes_[shard].push_back(
          ContextHash(key, offset, context & mask, mask));
    }
  }

  void AddConcurrent(uint64_t key, unsigned offset, uint8_t context = 0u) {
    assert(offset < 64u);
    const unsigned shard = Shard(key);
    std::lock_guard<SpinLock> lock(query_locks_[shard]);
    pending_queries_[shard].push_back(
        PendingQuery{key, static_cast<uint8_t>(offset)});
    const uint8_t mask = ContextMask(offset);
    if (mask != 0u) {
      context_hashes_[shard].push_back(
          ContextHash(key, offset, context & mask, mask));
    }
  }

  void Finalize() {
    // Active flanks are an immutable batch.  Reduce a contiguous record stream
    // once instead of issuing one random hash-table mutation per oriented
    // flank occurrence.  Sharding by an order-independent key hash exposes
    // enough parallel sort tasks without imposing any machine-sized cutoff.
    const int workers = std::max(1, omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
    for (unsigned shard = 0; shard < kNumQueryShards; ++shard) {
      std::vector<PendingQuery> &pending = pending_queries_[shard];
      std::sort(pending.begin(), pending.end(),
                [](const PendingQuery &lhs, const PendingQuery &rhs) {
                  if (lhs.key != rhs.key) return lhs.key < rhs.key;
                  return lhs.offset < rhs.offset;
                });
      std::vector<QueryEntry> &entries = query_info_[shard];
      entries.clear();
      entries.reserve(pending.size());
      for (size_t begin = 0; begin < pending.size();) {
        size_t end = begin + 1u;
        QueryInfo info;
        info.query_count = 1u;
        info.offset_mask = uint64_t{1} << pending[begin].offset;
        while (end < pending.size() &&
               pending[end].key == pending[begin].key) {
          ++info.query_count;
          info.offset_mask |= uint64_t{1} << pending[end].offset;
          ++end;
        }
        entries.push_back(QueryEntry{pending[begin].key, info});
        begin = end;
      }
      std::vector<PendingQuery>().swap(pending);
    }

    size_t context_count = 0u;
    for (const auto &hashes : context_hashes_) context_count += hashes.size();
    if (context_count != 0u) {
      size_t context_words = 1u;
      const size_t target_context_words =
          DivCeiling(context_count, size_t{8});
      while (context_words < target_context_words) {
        if (context_words > std::numeric_limits<size_t>::max() / 2u) {
          throw std::length_error("active context filter is too large");
        }
        context_words *= 2u;
      }
      context_filter_.assign(context_words, 0u);
      context_filter_mask_ = context_words - 1u;
#pragma omp parallel for schedule(static) num_threads(workers)
      for (unsigned shard = 0; shard < kNumQueryShards; ++shard) {
        for (uint64_t hash : context_hashes_[shard]) {
          __atomic_fetch_or(&context_filter_[hash & context_filter_mask_],
                            BitMask(hash), __ATOMIC_RELAXED);
        }
      }
      for (auto &hashes : context_hashes_) {
        std::vector<uint64_t>().swap(hashes);
      }
    }
  }

  bool MayMatchContext(uint64_t key, unsigned offset,
                       uint8_t context) const {
    const uint8_t mask = ContextMask(offset);
    if (mask == 0u || context_filter_.empty()) return true;
    const uint64_t hash = ContextHash(key, offset, context & mask, mask);
    const uint64_t bits = BitMask(hash);
    return (context_filter_[hash & context_filter_mask_] & bits) == bits;
  }

  uint64_t QueryOffsets(uint64_t key) const {
    if (filter_.empty()) {
      return 0;
    }
    const uint64_t hash = FilterHash(key);
    const uint64_t mask = BitMask(hash);
    if ((filter_[hash & filter_mask_] & mask) != mask) {
      return 0;
    }
    const unsigned shard = Shard(key);
    if (random_lookup_ready_) {
      const QueryMap &map = random_lookup_[shard];
      auto iterator = map.find(key);
      return iterator == map.end() ? 0 : iterator->second.offset_mask;
    }
    const std::vector<QueryEntry> &entries = query_info_[shard];
    const auto iterator = std::lower_bound(
        entries.begin(), entries.end(), key,
        [](const QueryEntry &entry, uint64_t value) {
          return entry.key < value;
        });
    return iterator == entries.end() || iterator->key != key
               ? 0u
               : iterator->info.offset_mask;
  }

  void PrepareRandomLookup() {
    const int workers = std::max(1, omp_get_max_threads());
    const size_t keys = key_count();
    if (keys != 0u && filter_.empty()) {
      size_t target_words = DivCeiling(keys, size_t{8});
      size_t num_words = 1u;
      while (num_words < target_words) {
        if (num_words > std::numeric_limits<size_t>::max() / 2u) {
          throw std::length_error("active anchor filter is too large");
        }
        num_words *= 2u;
      }
      filter_.assign(num_words, 0u);
      filter_mask_ = num_words - 1u;
#pragma omp parallel for schedule(static) num_threads(workers)
      for (unsigned shard = 0; shard < kNumQueryShards; ++shard) {
        for (const QueryEntry &entry : query_info_[shard]) {
          const uint64_t hash = FilterHash(entry.key);
          __atomic_fetch_or(&filter_[hash & filter_mask_], BitMask(hash),
                            __ATOMIC_RELAXED);
        }
      }
    }
#pragma omp parallel for schedule(static) num_threads(workers)
    for (unsigned shard = 0; shard < kNumQueryShards; ++shard) {
      QueryMap &map = random_lookup_[shard];
      map.clear();
      map.reserve(query_info_[shard].size());
      for (const QueryEntry &entry : query_info_[shard]) {
        map.emplace(entry.key, entry.info);
      }
    }
    random_lookup_ready_ = true;
  }

  size_t key_count() const {
    size_t count = 0;
    for (const auto &entries : query_info_) count += entries.size();
    return count;
  }

  uint64_t query_count() const {
    uint64_t count = 0;
    for (const auto &entries : query_info_) {
      for (const QueryEntry &entry : entries) count += entry.info.query_count;
    }
    return count;
  }

  uint64_t candidate_offset_count() const {
    uint64_t count = 0;
    for (const auto &entries : query_info_) {
      for (const QueryEntry &entry : entries) {
        count += static_cast<unsigned>(
            __builtin_popcountll(entry.info.offset_mask));
      }
    }
    return count;
  }

  size_t filter_bytes() const { return filter_.size() * sizeof(uint64_t); }

  bool Contains(uint64_t key) const {
    const std::vector<QueryEntry> &entries = query_info_[Shard(key)];
    const auto iterator = std::lower_bound(
        entries.begin(), entries.end(), key,
        [](const QueryEntry &entry, uint64_t value) {
          return entry.key < value;
        });
    return iterator != entries.end() && iterator->key == key;
  }

  template <class Visitor>
  void ForEach(const Visitor &visitor) const {
    for (const auto &entries : query_info_) {
      for (const QueryEntry &entry : entries) {
        visitor(entry.key, entry.info);
      }
    }
  }

 private:
  using QueryMap = phmap::flat_hash_map<uint64_t, QueryInfo>;
  struct QueryEntry {
    uint64_t key;
    QueryInfo info;
  };

  static unsigned Shard(uint64_t key) {
    return static_cast<unsigned>(Mix64(key) >>
                                 (64u - kQueryShardBits));
  }

  static uint64_t FilterHash(uint64_t key) {
    return Mix64(key ^ UINT64_C(0x3c6ef372fe94f82b));
  }

  static uint64_t BitMask(uint64_t hash) {
    return (uint64_t{1} << ((hash >> 32u) & 63u)) |
           (uint64_t{1} << ((hash >> 40u) & 63u)) |
           (uint64_t{1} << ((hash >> 48u) & 63u)) |
           (uint64_t{1} << ((hash >> 56u) & 63u));
  }

  uint8_t ContextMask(unsigned offset) const {
    uint8_t mask = 0u;
    if (context_bits_ >= 2u && offset >= context_distance_) mask |= 0x3u;
    if (context_bits_ >= 4u &&
        offset + anchor_len_ + context_distance_ <= window_len_) {
      mask |= 0xcu;
    }
    return mask;
  }

  static uint64_t ContextHash(uint64_t key, unsigned offset,
                              uint8_t context, uint8_t mask) {
    uint64_t value = key ^
                     (uint64_t(offset) * UINT64_C(0x9e3779b97f4a7c15));
    value ^= uint64_t(context | uint8_t(mask << 4u)) *
             UINT64_C(0xbf58476d1ce4e5b9);
    return Mix64(value ^ UINT64_C(0x243f6a8885a308d3));
  }

  std::array<std::vector<PendingQuery>, kNumQueryShards> pending_queries_;
  std::array<std::vector<QueryEntry>, kNumQueryShards> query_info_;
  std::array<QueryMap, kNumQueryShards> random_lookup_;
  bool random_lookup_ready_{false};
  std::array<SpinLock, kNumQueryShards> query_locks_;
  std::array<std::vector<uint64_t>, kNumQueryShards> context_hashes_;
  std::vector<uint64_t> filter_;
  size_t filter_mask_{0};
  std::vector<uint64_t> context_filter_;
  size_t context_filter_mask_{0};
  unsigned context_bits_{0};
  unsigned anchor_len_{0};
  unsigned window_len_{0};
  unsigned context_distance_{1};
};

using AnchorWindow = Kmer<4, uint32_t>;  // up to 64 exact bases

class ActiveWindowSet {
 public:
  void Reserve(size_t expected_windows) {
    if (expected_windows == 0u) {
      words_.clear();
      word_mask_ = 0u;
      return;
    }
    // Sixteen bits per oriented window keeps the blocked-filter false
    // positive rate well below one percent.  False positives are harmless:
    // the prepared flank index performs the authoritative exact lookup during
    // replay.  Unlike the former full-window hash set, this structure cannot
    // reject a real flank and occupies a compact, contiguous working set.
    const size_t target_words = DivCeiling(expected_windows, size_t{4});
    size_t num_words = 1u;
    while (num_words < target_words) {
      if (num_words > std::numeric_limits<size_t>::max() / 2u) {
        throw std::length_error("active window filter is too large");
      }
      num_words *= 2u;
    }
    words_.assign(num_words, 0u);
    word_mask_ = num_words - 1u;
  }

  void Add(const AnchorWindow &window) {
    const uint64_t hash = Hash(window);
    words_[hash & word_mask_] |= BitMask(hash);
    ++insertions_;
  }
  void AddConcurrent(const AnchorWindow &window) {
    const uint64_t hash = Hash(window);
    __atomic_fetch_or(&words_[hash & word_mask_], BitMask(hash),
                      __ATOMIC_RELAXED);
  }
  void SetInsertionCount(size_t count) { insertions_ = count; }
  bool Contains(const AnchorWindow &window) const {
    if (words_.empty()) return false;
    const uint64_t hash = Hash(window);
    const uint64_t mask = BitMask(hash);
    return (words_[hash & word_mask_] & mask) == mask;
  }
  size_t size() const { return insertions_; }
  size_t byte_size() const { return words_.size() * sizeof(uint64_t); }

 private:
  static uint64_t Hash(const AnchorWindow &window) {
    return Mix64(static_cast<uint64_t>(KmerHash()(window)) ^
                 UINT64_C(0xa4093822299f31d0));
  }

  static uint64_t BitMask(uint64_t hash) {
    return (uint64_t{1} << ((hash >> 32u) & 63u)) |
           (uint64_t{1} << ((hash >> 40u) & 63u)) |
           (uint64_t{1} << ((hash >> 48u) & 63u)) |
           (uint64_t{1} << ((hash >> 56u) & 63u));
  }

  std::vector<uint64_t> words_;
  size_t word_mask_{0u};
  size_t insertions_{0u};
};

using BinaryChunk = PackedReadChunk;

class MappedWords {
 public:
  explicit MappedWords(const std::string &path) : path_(path) {
    fd_ = open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
      throw std::runtime_error("cannot open packed read library: " + path);
    }

    struct stat st;
    if (fstat(fd_, &st) != 0 || st.st_size <= 0 ||
        st.st_size % static_cast<off_t>(sizeof(uint32_t)) != 0 ||
        static_cast<uint64_t>(st.st_size) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
      close(fd_);
      fd_ = -1;
      throw std::runtime_error("invalid packed read library: " + path);
    }

    byte_size_ = static_cast<size_t>(st.st_size);
    void *mapping =
        mmap(nullptr, byte_size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (mapping == MAP_FAILED) {
      close(fd_);
      fd_ = -1;
      throw std::runtime_error("cannot mmap packed read library: " + path);
    }
    words_ = static_cast<const uint32_t *>(mapping);
#if defined(POSIX_FADV_SEQUENTIAL)
    (void)posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
  }

  ~MappedWords() {
    if (words_ != nullptr) {
      munmap(const_cast<uint32_t *>(words_), byte_size_);
    }
    if (fd_ >= 0) {
      close(fd_);
    }
  }

  MappedWords(const MappedWords &) = delete;
  MappedWords &operator=(const MappedWords &) = delete;

  const uint32_t *data() const { return words_; }
  size_t word_size() const { return byte_size_ / sizeof(uint32_t); }
  size_t byte_size() const { return byte_size_; }
  const std::string &path() const { return path_; }

 private:
  std::string path_;
  int fd_{-1};
  size_t byte_size_{0};
  const uint32_t *words_{nullptr};
};

class MappedAnchorPositions {
 public:
  MappedAnchorPositions(const std::string &read_path, unsigned anchor_len,
                        unsigned window_len, uint64_t source_bytes,
                        uint64_t num_reads, uint64_t num_bases) {
    const std::string path = ReadAnchorPositionPath(read_path);
    fd_ = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    struct stat status {};
    if (fd_ < 0 || fstat(fd_, &status) != 0 || status.st_size <= 0 ||
        static_cast<uint64_t>(status.st_size) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
      Close();
      return;
    }
    bytes_ = static_cast<size_t>(status.st_size);
    void *address = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (address == MAP_FAILED) {
      mapping_ = nullptr;
      Close();
      return;
    }
    mapping_ = static_cast<const uint8_t *>(address);
    if (bytes_ < sizeof(ReadAnchorPositionHeader)) {
      Close();
      return;
    }
    std::memcpy(&header_, mapping_, sizeof(header_));
    if (header_.magic != kReadAnchorPositionMagic ||
        header_.version != kReadAnchorPositionVersion ||
        header_.header_bytes != sizeof(ReadAnchorPositionHeader) ||
        header_.anchor_len != anchor_len ||
        header_.window_len != window_len ||
        header_.source_bytes != source_bytes ||
        header_.num_reads != num_reads || header_.num_bases != num_bases ||
        header_.num_chunks == 0u ||
        header_.num_chunks >
            std::numeric_limits<size_t>::max() /
                sizeof(ReadAnchorPositionChunk) ||
        header_.payload_bytes > UINT64_MAX - header_.header_bytes ||
        header_.chunk_index_offset !=
            header_.header_bytes + header_.payload_bytes ||
        header_.chunk_index_offset > bytes_) {
      Close();
      return;
    }
    const size_t chunk_bytes =
        static_cast<size_t>(header_.num_chunks) *
        sizeof(ReadAnchorPositionChunk);
    const size_t chunk_index_offset =
        static_cast<size_t>(header_.chunk_index_offset);
    if (chunk_bytes > bytes_ - chunk_index_offset ||
        chunk_index_offset + chunk_bytes != bytes_) {
      Close();
      return;
    }
    chunk_storage_.resize(static_cast<size_t>(header_.num_chunks));
    std::memcpy(chunk_storage_.data(),
                mapping_ + chunk_index_offset, chunk_bytes);
    chunks_ = chunk_storage_.data();
    for (size_t i = 0; i < static_cast<size_t>(header_.num_chunks); ++i) {
      const ReadAnchorPositionChunk &chunk = chunks_[i];
      if (chunk.word_begin >= chunk.word_end ||
          chunk.read_begin >= chunk.read_end ||
          chunk.payload_begin > chunk.payload_end ||
          chunk.payload_end > header_.payload_bytes ||
          (i == 0u && (chunk.word_begin != 0u || chunk.read_begin != 0u ||
                      chunk.payload_begin != 0u)) ||
          (i != 0u &&
           (chunk.word_begin != chunks_[i - 1u].word_end ||
            chunk.read_begin != chunks_[i - 1u].read_end ||
            chunk.payload_begin != chunks_[i - 1u].payload_end))) {
        Close();
        return;
      }
    }
    const ReadAnchorPositionChunk &last =
        chunks_[static_cast<size_t>(header_.num_chunks) - 1u];
    if (last.word_end * sizeof(uint32_t) != source_bytes ||
        last.read_end != num_reads ||
        last.payload_end != header_.payload_bytes) {
      Close();
      return;
    }
  }

  ~MappedAnchorPositions() { Close(); }
  MappedAnchorPositions(const MappedAnchorPositions &) = delete;
  MappedAnchorPositions &operator=(const MappedAnchorPositions &) = delete;

  bool valid() const { return mapping_ != nullptr; }
  uint64_t payload_bytes() const { return header_.payload_bytes; }

  bool Range(const BinaryChunk &wanted, const uint8_t **begin,
             const uint8_t **end) const {
    if (!valid()) return false;
    size_t low = 0;
    size_t high = static_cast<size_t>(header_.num_chunks);
    while (low < high) {
      const size_t middle = low + (high - low) / 2u;
      if (chunks_[middle].word_begin < wanted.word_begin)
        low = middle + 1u;
      else
        high = middle;
    }
    if (low >= static_cast<size_t>(header_.num_chunks) ||
        chunks_[low].word_begin != wanted.word_begin ||
        chunks_[low].read_begin != wanted.read_begin) {
      return false;
    }
    size_t last = low;
    while (last < static_cast<size_t>(header_.num_chunks) &&
           chunks_[last].word_end < wanted.word_end) {
      ++last;
    }
    if (last >= static_cast<size_t>(header_.num_chunks) ||
        chunks_[last].word_end != wanted.word_end ||
        chunks_[last].read_end != wanted.read_end) {
      return false;
    }
    *begin = mapping_ + header_.header_bytes + chunks_[low].payload_begin;
    *end = mapping_ + header_.header_bytes + chunks_[last].payload_end;
    return true;
  }

 private:
  void Close() {
    chunks_ = nullptr;
    chunk_storage_.clear();
    if (mapping_ != nullptr) {
      munmap(const_cast<uint8_t *>(mapping_), bytes_);
      mapping_ = nullptr;
    }
    if (fd_ >= 0) {
      close(fd_);
      fd_ = -1;
    }
    bytes_ = 0;
  }

  int fd_{-1};
  size_t bytes_{0};
  const uint8_t *mapping_{nullptr};
  ReadAnchorPositionHeader header_;
  std::vector<ReadAnchorPositionChunk> chunk_storage_;
  const ReadAnchorPositionChunk *chunks_{nullptr};
};

class MappedBytes {
 public:
  explicit MappedBytes(const std::string &path) {
    fd_ = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) {
      throw std::runtime_error("cannot open mapped index file: " + path);
    }
    struct stat status;
    if (fstat(fd_, &status) != 0 || status.st_size < 0) {
      close(fd_);
      fd_ = -1;
      throw std::runtime_error("cannot stat mapped index file: " + path);
    }
    size_ = static_cast<size_t>(status.st_size);
    if (size_ == 0) {
      return;
    }
    void *mapping = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (mapping == MAP_FAILED) {
      close(fd_);
      fd_ = -1;
      throw std::runtime_error("cannot mmap index file: " + path);
    }
    data_ = static_cast<const uint8_t *>(mapping);
#if defined(POSIX_FADV_RANDOM)
    (void)posix_fadvise(fd_, 0, 0, POSIX_FADV_RANDOM);
#endif
  }

  ~MappedBytes() {
    if (data_ != nullptr) {
      munmap(const_cast<uint8_t *>(data_), size_);
    }
    if (fd_ >= 0) close(fd_);
  }

  MappedBytes(const MappedBytes &) = delete;
  MappedBytes &operator=(const MappedBytes &) = delete;

  const uint8_t *data() const { return data_; }
  size_t size() const { return size_; }

 private:
  int fd_{-1};
  size_t size_{0};
  const uint8_t *data_{nullptr};
};

struct alignas(64) ThreadProfile {
  uint64_t reads{0};
  uint64_t bases{0};
  uint64_t windows{0};
  uint64_t anchors{0};
  uint64_t active_postings{0};
  uint64_t candidate_pairs{0};
  std::array<uint8_t, kHllRegisters> hll{};
  phmap::flat_hash_map<uint64_t, uint64_t> sampled_postings;
  std::vector<uint32_t> queue_pos;
  std::vector<uint64_t> queue_hash;
  std::vector<uint64_t> queue_key;
  std::vector<uint8_t> candidate_seen;
};

void UpdateHll(uint64_t key, std::array<uint8_t, kHllRegisters> *registers) {
  const uint64_t hash = Mix64(key ^ UINT64_C(0x6a09e667f3bcc909));
  const size_t index = hash & (kHllRegisters - 1u);
  const uint64_t remaining = hash >> kHllBits;
  const unsigned rank =
      remaining == 0
          ? 64u - kHllBits + 1u
          : static_cast<unsigned>(__builtin_clzll(remaining) + 1u -
                                  kHllBits);
  (*registers)[index] =
      std::max((*registers)[index], static_cast<uint8_t>(rank));
}

double EstimateHll(const std::array<uint8_t, kHllRegisters> &registers) {
  double inverse_sum = 0;
  size_t zero_count = 0;
  for (uint8_t value : registers) {
    inverse_sum += std::ldexp(1.0, -static_cast<int>(value));
    zero_count += value == 0;
  }
  const double count = static_cast<double>(kHllRegisters);
  const double alpha = 0.7213 / (1.0 + 1.079 / count);
  double estimate = alpha * count * count / inverse_sum;
  if (estimate <= 2.5 * count && zero_count != 0) {
    estimate = count * std::log(count / static_cast<double>(zero_count));
  }
  return estimate;
}

uint8_t PackedBase(const uint32_t *sequence, unsigned position) {
  return PackedReadBase(sequence, position);
}

template <class Visitor>
void ForEachReadAnchor(const uint32_t *sequence, unsigned length,
                       const ReadIndexOptions &options,
                       std::vector<uint32_t> *queue_pos,
                       std::vector<uint64_t> *queue_hash,
                       std::vector<uint64_t> *queue_key,
                       const Visitor &visitor) {
  ForEachPackedReadAnchor(sequence, length, options.anchor_len,
                          options.window_len, queue_pos, queue_hash,
                          queue_key, visitor);
}

void ProfileRead(const uint32_t *sequence, unsigned length,
                 const ReadIndexOptions &options, uint64_t sample_mask,
                 const ActiveQuerySet *active_queries,
                 ThreadProfile *profile) {
  ++profile->reads;
  profile->bases += length;
  if (length < options.window_len) {
    return;
  }

  profile->windows += length - options.window_len + 1u;
  if (active_queries != nullptr) {
    profile->candidate_seen.assign(length - options.window_len + 1u, 0);
  }

  ForEachReadAnchor(
      sequence, length, options, &profile->queue_pos, &profile->queue_hash,
      &profile->queue_key, [&](uint64_t selected_key, uint32_t selected_pos) {
    ++profile->anchors;
    UpdateHll(selected_key, &profile->hll);
    const uint64_t sample_hash =
        Mix64(selected_key ^ UINT64_C(0xbb67ae8584caa73b));
    if ((sample_hash & sample_mask) == 0) {
      ++profile->sampled_postings[selected_key];
    }
    if (active_queries != nullptr) {
      const uint64_t query_offsets =
          active_queries->QueryOffsets(selected_key);
      if (query_offsets != 0) {
        ++profile->active_postings;
        uint64_t remaining_offsets = query_offsets;
        while (remaining_offsets != 0) {
          const unsigned offset =
              static_cast<unsigned>(__builtin_ctzll(remaining_offsets));
          remaining_offsets &= remaining_offsets - 1u;
          if (selected_pos >= offset) {
            const unsigned candidate_start = selected_pos - offset;
            if (candidate_start + options.window_len <= length) {
              profile->candidate_seen[candidate_start] = 1;
            }
          }
        }
      }
    }
  });
  if (active_queries != nullptr) {
    for (uint8_t seen : profile->candidate_seen) {
      profile->candidate_pairs += seen != 0;
    }
  }
}

struct SelectedAnchor {
  uint64_t key;
  unsigned offset;
  uint8_t context;
};

template <class KmerType>
SelectedAnchor SelectQueryAnchor(const KmerType &kmer,
                                 const ReadIndexOptions &options,
                                 unsigned context_bits,
                                 unsigned context_distance) {
  const uint64_t key_mask =
      (uint64_t{1} << (2u * options.anchor_len)) - 1u;
  uint64_t key = 0;
  uint64_t selected_key = 0;
  unsigned selected_offset = 0;
  uint64_t selected_hash = std::numeric_limits<uint64_t>::max();
  for (unsigned i = 0; i < options.window_len; ++i) {
    key = ((key << 2u) | kmer.GetBase(i)) & key_mask;
    if (i + 1u < options.anchor_len) {
      continue;
    }
    const uint64_t hash = AnchorHash(key);
    // Match the read-side minimizer's right-most tie breaking.
    if (hash <= selected_hash) {
      selected_hash = hash;
      selected_key = key;
      selected_offset = i + 1u - options.anchor_len;
    }
  }
  uint8_t context = 0u;
  if (context_bits >= 2u && selected_offset >= context_distance) {
    context |= kmer.GetBase(selected_offset - context_distance);
  }
  if (context_bits >= 4u &&
      selected_offset + options.anchor_len + context_distance <=
          options.window_len) {
    context |= static_cast<uint8_t>(
        kmer.GetBase(selected_offset + options.anchor_len +
                     context_distance - 1u)
        << 2u);
  }
  return SelectedAnchor{selected_key, selected_offset, context};
}

struct PreparedFlankIndexBase {
  virtual ~PreparedFlankIndexBase() = default;
};

template <class KmerType>
struct PreparedFlankIndex final : public PreparedFlankIndexBase {
  static void *operator new(std::size_t bytes) {
    void *memory = nullptr;
    if (posix_memalign(&memory, alignof(PreparedFlankIndex), bytes) != 0) {
      throw std::bad_alloc();
    }
    return memory;
  }
  static void operator delete(void *memory) noexcept { free(memory); }

  PreparedFlankIndex(unsigned k, unsigned step) : index(k, step) {}
  ContigFlankIndex<KmerType> index;
};

template <class KmerType>
bool BuildActiveQueriesForType(const ReadIndexOptions &options,
                               ActiveQuerySet *queries,
                               ActiveWindowSet *windows,
                               std::unique_ptr<PreparedFlankIndexBase>
                                   *prepared_index) {
  if (KmerType::max_size() < static_cast<unsigned>(options.kmer_k + 1)) {
    return false;
  }
  std::unique_ptr<PreparedFlankIndex<KmerType>> prepared(
      new PreparedFlankIndex<KmerType>(options.kmer_k, options.step));
  ContigFlankIndex<KmerType> &flank_index = prepared->index;
  const double feed_begin = omp_get_wtime();
  const std::string files[] = {options.contig_file, options.bubble_file};
  for (const std::string &file : files) {
    if (file.empty()) {
      continue;
    }
    AsyncContigReader reader(file);
    while (true) {
      auto &batch = reader.Next();
      if (batch.first.seq_count() == 0) {
        break;
      }
      flank_index.FeedBatchContigs(batch.first, batch.second);
    }
  }
  const double feed_end = omp_get_wtime();
  if (windows != nullptr) {
    if (flank_index.size() > std::numeric_limits<size_t>::max() / 2u) {
      throw std::length_error("active oriented window count is too large");
    }
    windows->Reserve(flank_index.size() * 2u);
  }
  const bool parallel_anchors =
      std::getenv("MEGAHIT_EXPERIMENTAL_PARALLEL_ACTIVE_ANCHORS") != nullptr;
  std::vector<ActiveQuerySet::Batch> batches(
      parallel_anchors ? std::max(1, omp_get_max_threads()) : 0);
  const auto add_anchor = [&](const SelectedAnchor &anchor) {
    if (parallel_anchors) {
      queries->AddToBatch(&batches[omp_get_thread_num()], anchor.key,
                          anchor.offset, anchor.context);
    } else {
      queries->Add(anchor.key, anchor.offset, anchor.context);
    }
  };
  const auto add_window = [&](const AnchorWindow &window) {
    if (parallel_anchors) windows->AddConcurrent(window);
    else windows->Add(window);
  };
  const auto extract = [&](const KmerType &canonical) {
    const SelectedAnchor forward =
        SelectQueryAnchor(canonical, options, queries->context_bits(),
                          queries->context_distance());
    add_anchor(forward);
    if (windows != nullptr) {
      AnchorWindow window;
      for (unsigned i = 0; i < options.window_len; ++i) {
        window.ShiftAppend(canonical.GetBase(i), options.window_len);
      }
      add_window(window);
    }
    KmerType reverse = canonical;
    reverse.ReverseComplement(options.kmer_k + 1u);
    const SelectedAnchor backward =
        SelectQueryAnchor(reverse, options, queries->context_bits(),
                          queries->context_distance());
    add_anchor(backward);
    if (windows != nullptr) {
      AnchorWindow window;
      for (unsigned i = 0; i < options.window_len; ++i) {
        window.ShiftAppend(reverse.GetBase(i), options.window_len);
      }
      add_window(window);
    }
  };
  if (parallel_anchors) {
    flank_index.ForEachCanonicalKmerParallel(extract);
    queries->MergeBatches(&batches);
    if (windows != nullptr) windows->SetInsertionCount(flank_index.size() * 2u);
  } else {
    flank_index.ForEachCanonicalKmer(extract);
  }
  const double enumerate_end = omp_get_wtime();
  queries->Finalize();
  const double reduce_end = omp_get_wtime();
  xinfo("Active query anchors: {} oriented exact flanks in {} distinct "
        "anchor keys / {} distinct (key, offset) candidates, filter {} "
        "bytes\n",
        queries->query_count(), queries->key_count(),
        queries->candidate_offset_count(), queries->filter_bytes());
  if (prepared_index != nullptr) {
    // Active-anchor extraction only needs the canonical key table.  Exact
    // replay additionally needs the membership filters normally finalized by
    // FeedReplayFlankIndex(); build them once before transferring ownership.
    flank_index.Finalize();
    *prepared_index = std::move(prepared);
  }
  const double replay_index_end = omp_get_wtime();
  xinfo("Active query phases: flank feed {.4}, anchor records {.4}, "
        "sort/reduce+context {.4}, replay filters {.4} s; window filter {} "
        "bytes\n",
        feed_end - feed_begin, enumerate_end - feed_end,
        reduce_end - enumerate_end, replay_index_end - reduce_end,
        windows == nullptr ? 0u : windows->byte_size());
  return true;
}

void BuildActiveQueries(const ReadIndexOptions &options,
                        ActiveQuerySet *queries,
                        ActiveWindowSet *windows = nullptr,
                        std::unique_ptr<PreparedFlankIndexBase>
                            *prepared_index = nullptr) {
  if (BuildActiveQueriesForType<Kmer<1, uint64_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<3, uint32_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<2, uint64_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<5, uint32_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<3, uint64_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<7, uint32_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<4, uint64_t>>(
          options, queries, windows, prepared_index)) return;
  if (BuildActiveQueriesForType<Kmer<kUint32PerKmerMaxK, uint32_t>>(
          options, queries, windows, prepared_index)) {
    return;
  }
  throw std::logic_error("k is too large for active anchor profiling");
}

void ParallelSortUint64(std::vector<uint64_t> *values,
                        unsigned significant_bytes, int requested_threads) {
  if (values->size() < 4096u || requested_threads <= 1) {
    kmlib::kmsort(values->begin(), values->end());
    return;
  }
  significant_bytes = std::max(1u, std::min(8u, significant_bytes));
  const size_t size = values->size();
  const int threads = static_cast<int>(std::min<size_t>(
      static_cast<size_t>(std::max(1, requested_threads)), size));
  constexpr unsigned kBuckets = 256u;
  const unsigned shift = (significant_bytes - 1u) * 8u;
  std::vector<size_t> histogram(static_cast<size_t>(threads) * kBuckets, 0u);
#pragma omp parallel num_threads(threads)
  {
    const int thread = omp_get_thread_num();
    const size_t begin = size * static_cast<size_t>(thread) / threads;
    const size_t end = size * static_cast<size_t>(thread + 1) / threads;
    size_t *counts = histogram.data() +
                     static_cast<size_t>(thread) * kBuckets;
    for (size_t i = begin; i < end; ++i) {
      ++counts[((*values)[i] >> shift) & 0xffu];
    }
  }
  std::array<size_t, kBuckets + 1u> bucket_offsets{};
  for (unsigned bucket = 0; bucket < kBuckets; ++bucket) {
    size_t count = 0;
    for (int thread = 0; thread < threads; ++thread) {
      const size_t index = static_cast<size_t>(thread) * kBuckets + bucket;
      const size_t thread_count = histogram[index];
      histogram[index] = bucket_offsets[bucket] + count;
      count += thread_count;
    }
    bucket_offsets[bucket + 1u] = bucket_offsets[bucket] + count;
  }
  std::vector<uint64_t> partitioned(size);
#pragma omp parallel num_threads(threads)
  {
    const int thread = omp_get_thread_num();
    const size_t begin = size * static_cast<size_t>(thread) / threads;
    const size_t end = size * static_cast<size_t>(thread + 1) / threads;
    size_t *cursors = histogram.data() +
                      static_cast<size_t>(thread) * kBuckets;
    for (size_t i = begin; i < end; ++i) {
      const uint64_t value = (*values)[i];
      const unsigned bucket = static_cast<unsigned>((value >> shift) & 0xffu);
      partitioned[cursors[bucket]++] = value;
    }
  }
  values->swap(partitioned);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
  for (unsigned bucket = 0; bucket < kBuckets; ++bucket) {
    kmlib::kmsort(values->begin() + bucket_offsets[bucket],
                  values->begin() + bucket_offsets[bucket + 1u]);
  }
}

struct CompactLocalEndpointQueries {
  struct Entry {
    uint64_t key;
    uint64_t offset_mask;
  };

  template <class Visitor>
  void ForEach(const Visitor &visitor) const {
    for (const Entry &entry : anchors) {
      ActiveQuerySet::QueryInfo info;
      info.offset_mask = entry.offset_mask;
      visitor(entry.key, info);
    }
  }

  bool Contains(const AnchorWindow &window) const {
    uint64_t forward = 0;
    for (unsigned i = 0; i < seed_len; ++i) {
      forward = (forward << 2u) | window.GetBase(i);
    }
    const unsigned padding = 64u - seed_len * 2u;
    if (padding != 0u) forward <<= padding;
    uint64_t reverse = kmlib::bit::ReverseComplement<2>(forward);
    if (padding != 0u) reverse <<= padding;
    return std::binary_search(canonical_seeds.begin(),
                              canonical_seeds.end(),
                              std::min(forward, reverse));
  }

  uint64_t candidate_offset_count() const {
    uint64_t count = 0;
    for (const Entry &entry : anchors) {
      count += static_cast<unsigned>(__builtin_popcountll(entry.offset_mask));
    }
    return count;
  }

  unsigned seed_len{0};
  std::vector<uint64_t> canonical_seeds;
  std::vector<Entry> anchors;
};

void BuildLocalEndpointQueries(const ReadIndexOptions &options,
                               CompactLocalEndpointQueries *queries) {
  if (options.local_seed_len <= 0 || options.local_seed_len > 32 ||
      options.local_sparsity <= 0 || options.local_endpoint_range < 0 ||
      options.local_min_contig_len < 0) {
    throw std::logic_error("invalid local endpoint query parameters");
  }
  if (options.window_len != static_cast<unsigned>(options.local_seed_len)) {
    throw std::logic_error(
        "compact local query needs index window == exact mapper seed");
  }

  {
    // Rebuild exactly the mapper's sampled canonical seed set, including its
    // global-repeat rejection.  The persistent read index then reverses the
    // full-read join; the normal mapper remains the final semantic verifier.
    HashMapper endpoint_mapper;
    endpoint_mapper.LoadAndBuild(
        options.contig_file, options.local_min_contig_len,
        options.local_seed_len, options.local_sparsity, false);
    queries->canonical_seeds = endpoint_mapper.CollectUniqueEndpointSeeds(
        options.local_endpoint_range);
    endpoint_mapper.ReleaseIndex();
  }

  const unsigned seed_len = static_cast<unsigned>(options.local_seed_len);
  queries->seed_len = seed_len;
  const unsigned seed_padding = 64u - seed_len * 2u;
  ParallelSortUint64(&queries->canonical_seeds, sizeof(uint64_t),
                     options.num_threads);

  unsigned offset_bits = 0u;
  const unsigned maximum_offset = options.window_len - options.anchor_len;
  while ((uint64_t{1} << offset_bits) <= maximum_offset) ++offset_bits;
  const unsigned anchor_key_bits = options.anchor_len * 2u;
  if (anchor_key_bits + offset_bits > 64u ||
      queries->canonical_seeds.size() >
          std::numeric_limits<size_t>::max() / 2u) {
    throw std::length_error("compact local anchor records are too large");
  }
  std::vector<uint64_t> packed_anchors(queries->canonical_seeds.size() * 2u);
  const uint64_t anchor_mask =
      (uint64_t{1} << anchor_key_bits) - 1u;
  const auto select_anchor = [&](uint64_t oriented) {
    uint64_t key = 0u;
    uint64_t selected_key = 0u;
    uint64_t selected_hash = std::numeric_limits<uint64_t>::max();
    unsigned selected_offset = 0u;
    for (unsigned i = 0; i < seed_len; ++i) {
      const uint8_t base = static_cast<uint8_t>(
          (oriented >> (62u - i * 2u)) & 3u);
      key = ((key << 2u) | base) & anchor_mask;
      if (i + 1u < options.anchor_len) continue;
      const uint64_t hash = AnchorHash(key);
      if (hash <= selected_hash) {
        selected_hash = hash;
        selected_key = key;
        selected_offset = i + 1u - options.anchor_len;
      }
    }
    return (selected_key << offset_bits) | selected_offset;
  };
#pragma omp parallel for schedule(static) num_threads(options.num_threads)
  for (int64_t seed_id = 0;
       seed_id < static_cast<int64_t>(queries->canonical_seeds.size());
       ++seed_id) {
    const uint64_t canonical =
        queries->canonical_seeds[static_cast<size_t>(seed_id)];
    uint64_t reverse = kmlib::bit::ReverseComplement<2>(canonical);
    if (seed_padding != 0u) reverse <<= seed_padding;
    packed_anchors[static_cast<size_t>(seed_id) * 2u] =
        select_anchor(canonical);
    packed_anchors[static_cast<size_t>(seed_id) * 2u + 1u] =
        select_anchor(reverse);
  }
  const unsigned packed_anchor_bytes =
      DivCeiling(anchor_key_bits + offset_bits, 8u);
  ParallelSortUint64(&packed_anchors, packed_anchor_bytes,
                     options.num_threads);
  queries->anchors.clear();
  queries->anchors.reserve(packed_anchors.size());
  const uint64_t offset_mask = LowBitMask(offset_bits);
  for (uint64_t packed : packed_anchors) {
    const uint64_t key = offset_bits == 0u ? packed : packed >> offset_bits;
    const unsigned offset =
        offset_bits == 0u ? 0u : static_cast<unsigned>(packed & offset_mask);
    if (queries->anchors.empty() || queries->anchors.back().key != key) {
      queries->anchors.push_back(
          CompactLocalEndpointQueries::Entry{key, uint64_t{1} << offset});
    } else {
      queries->anchors.back().offset_mask |= uint64_t{1} << offset;
    }
  }
  std::vector<uint64_t>().swap(packed_anchors);
  xinfo("Local endpoint index query: {} exact oriented {}-mer windows in {} "
        "anchor keys / {} (key, offset) candidates; compact exact set {} "
        "bytes\n",
        queries->canonical_seeds.size() * 2u, options.window_len,
        queries->anchors.size(), queries->candidate_offset_count(),
        queries->canonical_seeds.size() * sizeof(uint64_t));
}

std::vector<BinaryChunk> ParseChunks(const MappedWords &mapping,
                                     uint64_t max_reads,
                                     int num_threads,
                                     uint64_t *num_reads,
                                     uint64_t *num_bases,
                                     unsigned *max_read_len,
                                     size_t *parsed_words,
                                     size_t max_chunk_bytes =
                                         size_t{16} << 20u) {
  constexpr size_t kMinChunkBytes = size_t{1} << 20u;
  max_chunk_bytes = std::max(kMinChunkBytes, max_chunk_bytes);
  const size_t desired_chunks =
      static_cast<size_t>(std::max(1, num_threads)) * 8u;
  const size_t adaptive_bytes =
      DivCeiling(mapping.byte_size(), desired_chunks);
  const size_t target_chunk_bytes =
      std::max(kMinChunkBytes,
               std::min(max_chunk_bytes, adaptive_bytes));
  const size_t target_chunk_words =
      target_chunk_bytes / sizeof(uint32_t);
  const size_t fine_chunk_words = kMinChunkBytes / sizeof(uint32_t);

  const auto coalesce_chunks = [&](const std::vector<BinaryChunk> &fine) {
    std::vector<BinaryChunk> output;
    output.reserve(std::min<size_t>(
        fine.size(), DivCeiling(mapping.word_size(), target_chunk_words)));
    BinaryChunk accumulated;
    bool have_chunk = false;
    for (size_t i = 0; i < fine.size(); ++i) {
      const BinaryChunk &chunk = fine[i];
      if (!have_chunk) {
        accumulated = chunk;
        have_chunk = true;
      } else {
        accumulated.word_end = chunk.word_end;
        accumulated.read_end = chunk.read_end;
        accumulated.num_bases += chunk.num_bases;
        accumulated.max_read_len =
            std::max(accumulated.max_read_len, chunk.max_read_len);
      }
      if (accumulated.word_end - accumulated.word_begin >=
              target_chunk_words ||
          i + 1u == fine.size()) {
        output.push_back(accumulated);
        have_chunk = false;
      }
    }
    return output;
  };

  if (max_reads == 0) {
    PackedReadChunkIndex cached;
    if (LoadPackedReadChunkIndex(mapping.path(), &cached)) {
      *num_reads = cached.num_reads;
      *num_bases = cached.num_bases;
      *max_read_len = cached.max_read_len;
      *parsed_words = mapping.word_size();
      return coalesce_chunks(cached.chunks);
    }
  }

  const uint32_t *words = mapping.data();
  const size_t total_words = mapping.word_size();
  std::vector<BinaryChunk> fine_chunks;
  fine_chunks.reserve(DivCeiling(total_words, fine_chunk_words));

  size_t word_pos = 0;
  size_t chunk_begin = 0;
  uint64_t read_id = 0;
  uint64_t chunk_read_begin = 0;
  uint64_t bases = 0;
  uint64_t chunk_bases = 0;
  unsigned longest = 0;
  unsigned chunk_longest = 0;
  while (word_pos < total_words && (max_reads == 0 || read_id < max_reads)) {
    const uint32_t length = words[word_pos];
    if (length == 0) {
      throw std::runtime_error("packed read library contains a zero length");
    }
    const size_t sequence_words =
        DivCeiling(static_cast<size_t>(length), SeqPackage::kBasesPerWord);
    if (sequence_words > total_words - word_pos - 1u) {
      throw std::runtime_error("truncated packed read library");
    }
    word_pos += 1u + sequence_words;
    ++read_id;
    bases += length;
    chunk_bases += length;
    longest = std::max(longest, static_cast<unsigned>(length));
    chunk_longest =
        std::max(chunk_longest, static_cast<unsigned>(length));

    if (word_pos - chunk_begin >= fine_chunk_words ||
        word_pos == total_words || (max_reads != 0 && read_id == max_reads)) {
      fine_chunks.push_back({chunk_begin, word_pos, chunk_read_begin, read_id,
                             chunk_bases, chunk_longest});
      // The boundary pass only needs record headers.  Drop completed PTEs so
      // profiling a very large library does not make the complete file
      // mapping look resident before parallel workers start their scan.  The
      // file pages remain available in the kernel cache and are faulted back
      // by the worker that owns this chunk.
      DiscardMemoryPages(
          const_cast<uint32_t *>(words + chunk_begin),
          (word_pos - chunk_begin) * sizeof(uint32_t));
      chunk_begin = word_pos;
      chunk_read_begin = read_id;
      chunk_bases = 0;
      chunk_longest = 0;
    }
  }
  if (max_reads == 0 && word_pos != total_words) {
    throw std::runtime_error("packed read library has trailing data");
  }
  *num_reads = read_id;
  *num_bases = bases;
  *max_read_len = longest;
  *parsed_words = word_pos;
  if (max_reads == 0) {
    (void)PublishPackedReadChunkIndex(mapping.path(), read_id, bases, longest,
                                      fine_chunks);
  }
  return coalesce_chunks(fine_chunks);
}

uint64_t Quantile(const std::vector<uint64_t> &values, double fraction) {
  if (values.empty()) {
    return 0;
  }
  const size_t index = static_cast<size_t>(
      fraction * static_cast<double>(values.size() - 1u));
  return values[index];
}

// Read-order buckets bound the mapped-read working set during exact replay.
// Their count is selected from the requested concurrency, input extent and
// memory budget; it is deliberately independent of any particular socket or
// NUMA topology.
constexpr unsigned kMinimumReadOrderBuckets = 256u;
constexpr unsigned kMinimumReadIndexBucketBits = 8;
constexpr unsigned kMaximumReadIndexBucketBits = 24;
constexpr unsigned kBuildHllBits = 6;
constexpr unsigned kBuildHllRegisters = 1u << kBuildHllBits;

template <int KeyBytes>
struct ReadIndexKeyRadixTraits {
  static const int n_bytes = KeyBytes;
  int kth_byte(uint64_t key, int byte) const {
    return static_cast<int>((key >> (byte * 8u)) & 0xFFu);
  }
  bool operator()(uint64_t lhs, uint64_t rhs) const {
    return lhs < rhs;
  }
};

void SortReadIndexKeys(std::vector<uint64_t> *keys, unsigned key_bytes) {
  switch (key_bytes) {
    case 1:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<1>());
      return;
    case 2:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<2>());
      return;
    case 3:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<3>());
      return;
    case 4:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<4>());
      return;
    case 5:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<5>());
      return;
    case 6:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<6>());
      return;
    case 7:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<7>());
      return;
    case 8:
      kmlib::kmsort(keys->begin(), keys->end(),
                    ReadIndexKeyRadixTraits<8>());
      return;
    default:
      throw std::logic_error("unsupported read-index key width");
  }
}

double ProcessCpuSeconds() {
  struct timespec value;
  if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &value) != 0) {
    return 0.0;
  }
  return static_cast<double>(value.tv_sec) +
         static_cast<double>(value.tv_nsec) * 1e-9;
}

uint64_t ClockNanoseconds(clockid_t clock_id) {
  struct timespec value;
  if (clock_gettime(clock_id, &value) != 0) return 0;
  return static_cast<uint64_t>(value.tv_sec) * UINT64_C(1000000000) +
         static_cast<uint64_t>(value.tv_nsec);
}

class CompactMemoryLimiter {
 public:
  explicit CompactMemoryLimiter(uint64_t budget)
      : budget_(std::max<uint64_t>(1, budget)) {}

  uint64_t Acquire(uint64_t requested) {
    const uint64_t charged = std::min(requested, budget_);
    std::unique_lock<std::mutex> lock(mutex_);
    ready_.wait(lock, [&] { return charged <= budget_ - used_; });
    used_ += charged;
    peak_used_ = std::max(peak_used_, used_);
    ++active_workers_;
    peak_workers_ = std::max(peak_workers_, active_workers_);
    return charged;
  }

  void Release(uint64_t charged) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      assert(charged <= used_ && active_workers_ != 0);
      used_ -= charged;
      --active_workers_;
    }
    ready_.notify_all();
  }

  uint64_t peak_used() const { return peak_used_; }
  unsigned peak_workers() const { return peak_workers_; }

 private:
  const uint64_t budget_;
  uint64_t used_{0};
  uint64_t peak_used_{0};
  unsigned active_workers_{0};
  unsigned peak_workers_{0};
  mutable std::mutex mutex_;
  std::condition_variable ready_;
};

class CompactMemoryLease {
 public:
  CompactMemoryLease(CompactMemoryLimiter *limiter, uint64_t requested)
      : limiter_(limiter), charged_(limiter->Acquire(requested)) {}
  ~CompactMemoryLease() { limiter_->Release(charged_); }

  CompactMemoryLease(const CompactMemoryLease &) = delete;
  CompactMemoryLease &operator=(const CompactMemoryLease &) = delete;

 private:
  CompactMemoryLimiter *limiter_;
  uint64_t charged_;
};

struct ReadIndexBucketMeta {
  ReadIndexBucketMeta() = default;
  ReadIndexBucketMeta(uint64_t directory_offset_arg,
                      uint64_t directory_count_arg,
                      uint64_t posting_offset_arg,
                      uint64_t posting_count_arg,
                      unsigned file_id_arg = 0)
      : file_id(file_id_arg),
        directory_offset(directory_offset_arg),
        directory_count(directory_count_arg),
        posting_offset(posting_offset_arg),
        posting_count(posting_count_arg) {}

  unsigned file_id{0};
  uint64_t directory_offset{0};
  uint64_t directory_count{0};
  uint64_t posting_offset{0};
  uint64_t posting_count{0};
};

struct ReadIndexBuildScratch {
  std::vector<uint32_t> queue_pos;
  std::vector<uint64_t> queue_hash;
  std::vector<uint64_t> queue_key;
  std::vector<uint32_t> replay_state;
  std::vector<uint32_t> replay_candidate_positions;
  std::vector<std::vector<uint8_t>> bucket_records;
  std::vector<uint8_t> packed_records;
  std::vector<uint8_t> bucket_hll;
  std::vector<size_t> offsets;
  std::vector<size_t> cursors;
  uint64_t replay_candidate_reads{0};
  uint64_t replay_candidate_windows{0};
  uint64_t replay_aligned_reads{0};
  uint64_t replay_generated_edges{0};
  size_t buffered_record_bytes{0};
  unsigned reads_until_flush_check{64};
};

unsigned ReadIndexBucketFromHash(uint64_t hash, unsigned bucket_bits) {
  assert(bucket_bits > 0 && bucket_bits < 32u);
  return static_cast<unsigned>(hash >> (64u - bucket_bits));
}

unsigned ReadIndexBucket(uint64_t key, unsigned bucket_bits) {
  return ReadIndexBucketFromHash(AnchorBucketHash(key), bucket_bits);
}

void UpdateBuildHll(uint64_t hash, unsigned bucket,
                    std::vector<uint8_t> *registers) {
  const unsigned register_id = static_cast<unsigned>(
      hash >> (64u - kBuildHllBits));
  const uint64_t remainder = hash << kBuildHllBits;
  const unsigned rank =
      remainder == 0
          ? 64u - kBuildHllBits + 1u
          : static_cast<unsigned>(__builtin_clzll(remainder)) + 1u;
  uint8_t &stored =
      (*registers)[static_cast<size_t>(bucket) * kBuildHllRegisters +
                   register_id];
  stored = std::max<uint8_t>(stored, static_cast<uint8_t>(rank));
}

void StoreLowBytes(uint8_t *destination, uint64_t value, unsigned bytes) {
  assert(bytes <= sizeof(value));
  for (unsigned byte = 0; byte < bytes; ++byte) {
    destination[byte] = static_cast<uint8_t>(value >> (byte * 8u));
  }
}

uint64_t LoadLowBytes(const uint8_t *source, unsigned bytes) {
  assert(bytes <= sizeof(uint64_t));
  uint64_t value = 0;
  for (unsigned byte = 0; byte < bytes; ++byte) {
    value |= static_cast<uint64_t>(source[byte]) << (byte * 8u);
  }
  return value;
}

bool PwriteFully(int fd, const void *buffer, size_t bytes, uint64_t offset) {
  const uint8_t *source = static_cast<const uint8_t *>(buffer);
  while (bytes != 0) {
    const ssize_t written =
        pwrite(fd, source, bytes, static_cast<off_t>(offset));
    if (written < 0 && errno == EINTR) {
      continue;
    }
    if (written <= 0) {
      return false;
    }
    source += written;
    bytes -= static_cast<size_t>(written);
    offset += static_cast<uint64_t>(written);
  }
  return true;
}

bool PwritevFully(int fd, std::vector<struct iovec> *vectors,
                  uint64_t offset) {
  size_t first = 0;
  long maximum_iov = sysconf(_SC_IOV_MAX);
  if (maximum_iov <= 0) maximum_iov = 1024;
  while (first < vectors->size()) {
    const int count = static_cast<int>(std::min<size_t>(
        vectors->size() - first, static_cast<size_t>(maximum_iov)));
    ssize_t written =
        pwritev(fd, vectors->data() + first, count, static_cast<off_t>(offset));
    if (written < 0 && errno == EINTR) continue;
    if (written <= 0) return false;
    offset += static_cast<uint64_t>(written);
    size_t consumed = static_cast<size_t>(written);
    while (first < vectors->size() &&
           consumed >= (*vectors)[first].iov_len) {
      consumed -= (*vectors)[first].iov_len;
      ++first;
    }
    if (consumed != 0) {
      uint8_t *base = static_cast<uint8_t *>((*vectors)[first].iov_base);
      (*vectors)[first].iov_base = base + consumed;
      (*vectors)[first].iov_len -= consumed;
    }
  }
  return true;
}

bool PreadFully(int fd, void *buffer, size_t bytes, uint64_t offset) {
  uint8_t *destination = static_cast<uint8_t *>(buffer);
  while (bytes != 0) {
    const ssize_t amount =
        pread(fd, destination, bytes, static_cast<off_t>(offset));
    if (amount < 0 && errno == EINTR) {
      continue;
    }
    if (amount <= 0) {
      return false;
    }
    destination += amount;
    bytes -= static_cast<size_t>(amount);
    offset += static_cast<uint64_t>(amount);
  }
  return true;
}

struct ReadIndexMetadata {
  unsigned version{0};
  unsigned num_files{1};
  unsigned anchor_len{0};
  unsigned window_len{0};
  unsigned bucket_bits{0};
  unsigned position_bits{0};
  unsigned locator_bytes{0};
  unsigned locator_context_bits{0};
  unsigned locator_context_distance{1};
  unsigned key_bytes{0};
  unsigned bucket_begin_bytes{0};
  unsigned directory_entry_bytes{0};
  bool postings_sorted{true};
  uint64_t num_reads{0};
  uint64_t num_bases{0};
  uint64_t parsed_words{0};
  uint64_t num_keys{0};
  uint64_t num_occurrences{0};
  std::vector<ReadIndexBucketMeta> buckets;
};

unsigned ReadIndexBucketForMetadata(uint64_t key,
                                    const ReadIndexMetadata &metadata) {
  if (metadata.version == 1) {
    return ReadIndexBucket(key, metadata.bucket_bits);
  }
  const unsigned key_bits = metadata.anchor_len * 2u;
  const unsigned suffix_bits = key_bits - metadata.bucket_bits;
  return static_cast<unsigned>(
      PermuteReadIndexKey(key, key_bits) >> suffix_bits);
}

uint64_t ReadIndexDirectoryKey(uint64_t key,
                               const ReadIndexMetadata &metadata) {
  if (metadata.version == 1) return key;
  const unsigned key_bits = metadata.anchor_len * 2u;
  const unsigned suffix_bits = key_bits - metadata.bucket_bits;
  return PermuteReadIndexKey(key, key_bits) & LowBitMask(suffix_bits);
}

ReadIndexMetadata LoadReadIndexMetadata(const std::string &prefix) {
  std::ifstream input((prefix + ".ridx.info").c_str());
  if (!input) {
    throw std::runtime_error("cannot open read-index metadata");
  }
  ReadIndexMetadata metadata;
  std::string field;
  while (input >> field) {
    if (field == "bucket") {
      unsigned bucket = 0;
      ReadIndexBucketMeta entry;
      input >> bucket;
      if (metadata.version >= 3) input >> entry.file_id;
      input >> entry.directory_offset >> entry.directory_count >>
          entry.posting_offset >> entry.posting_count;
      if (metadata.bucket_bits == 0 ||
          bucket >= (1u << metadata.bucket_bits)) {
        throw std::runtime_error("invalid read-index bucket metadata");
      }
      if (metadata.buckets.empty()) {
        metadata.buckets.resize(1u << metadata.bucket_bits);
      }
      metadata.buckets[bucket] = entry;
    } else {
      uint64_t value = 0;
      input >> value;
      if (field == "version") metadata.version = value;
      else if (field == "num_files") metadata.num_files = value;
      else if (field == "anchor_len") metadata.anchor_len = value;
      else if (field == "window_len") metadata.window_len = value;
      else if (field == "bucket_bits") metadata.bucket_bits = value;
      else if (field == "position_bits") metadata.position_bits = value;
      else if (field == "locator_bytes") metadata.locator_bytes = value;
      else if (field == "locator_context_bits") {
        metadata.locator_context_bits = value;
      }
      else if (field == "locator_context_distance") {
        metadata.locator_context_distance = value;
      }
      else if (field == "key_bytes") metadata.key_bytes = value;
      else if (field == "bucket_begin_bytes") {
        metadata.bucket_begin_bytes = value;
      } else if (field == "directory_entry_bytes") {
        metadata.directory_entry_bytes = value;
      } else if (field == "postings_sorted") {
        metadata.postings_sorted = value != 0;
      } else if (field == "num_reads") metadata.num_reads = value;
      else if (field == "num_bases") metadata.num_bases = value;
      else if (field == "parsed_words") metadata.parsed_words = value;
      else if (field == "num_keys") metadata.num_keys = value;
      else if (field == "num_occurrences") metadata.num_occurrences = value;
      else throw std::runtime_error("unknown read-index metadata field");
    }
    if (!input) {
      throw std::runtime_error("truncated read-index metadata");
    }
  }
  const bool supported_version = metadata.version == 1 ||
                                 metadata.version == 2 ||
                                 metadata.version == 3 ||
                                 metadata.version == 4;
  const unsigned key_bits = metadata.anchor_len * 2u;
  const unsigned locator_value_bits =
      metadata.position_bits +
      BitsNeeded(metadata.parsed_words == 0 ? 0 : metadata.parsed_words - 1u);
  const unsigned expected_key_bytes =
      metadata.version >= 2 && key_bits > metadata.bucket_bits
          ? DivCeiling(key_bits - metadata.bucket_bits, 8u)
          : DivCeiling(key_bits, 8u);
  if (!supported_version || metadata.anchor_len == 0 ||
      metadata.anchor_len > 31 ||
      metadata.bucket_bits < kMinimumReadIndexBucketBits ||
      metadata.bucket_bits > kMaximumReadIndexBucketBits ||
      metadata.locator_bytes == 0 || metadata.locator_bytes > 8 ||
      locator_value_bits + metadata.locator_context_bits >
          metadata.locator_bytes * 8u ||
      (metadata.locator_context_bits != 0u &&
       metadata.locator_context_bits != 2u &&
       metadata.locator_context_bits != 4u) ||
      (metadata.locator_context_bits != 0u &&
       (metadata.locator_context_distance == 0u ||
        metadata.window_len <= metadata.anchor_len ||
        metadata.locator_context_distance >
            metadata.window_len - metadata.anchor_len)) ||
      (metadata.version >= 2 && key_bits <= metadata.bucket_bits) ||
      metadata.num_files == 0 ||
      metadata.key_bytes != expected_key_bytes ||
      metadata.bucket_begin_bytes == 0 ||
      metadata.bucket_begin_bytes > 8 ||
      metadata.directory_entry_bytes !=
          metadata.key_bytes + metadata.bucket_begin_bytes + 4u ||
      metadata.buckets.size() != (size_t{1} << metadata.bucket_bits)) {
    throw std::runtime_error("unsupported or invalid read-index metadata");
  }
  for (const ReadIndexBucketMeta &bucket : metadata.buckets) {
    if (bucket.file_id >= metadata.num_files) {
      throw std::runtime_error("read-index bucket refers to invalid shard");
    }
  }
  return metadata;
}

struct ActiveAnchorEntry {
  uint64_t key;
  uint64_t anchor_key;
  uint64_t offset_mask;
};

struct MatchedPostingGroup {
  unsigned file_id;
  uint64_t posting_byte_offset;
  uint64_t count;
  uint64_t anchor_key;
  uint64_t offset_mask;
};

uint64_t FileSize(int fd);

std::string ReadIndexShardPath(const std::string &prefix,
                               const char *component, unsigned file_id) {
  const std::string base = prefix + ".ridx." + component;
  return file_id == 0 ? base : base + "." + std::to_string(file_id);
}

class ReadIndexShardFiles {
 public:
  ReadIndexShardFiles(const std::string &prefix,
                      const ReadIndexMetadata &metadata,
                      bool map_directories)
      : directory_fds_(metadata.num_files, -1),
        posting_fds_(metadata.num_files, -1),
        directory_sizes_(metadata.num_files, 0),
        posting_sizes_(metadata.num_files, 0),
        mapped_directories_(metadata.num_files) {
    try {
      for (unsigned file_id = 0; file_id < metadata.num_files; ++file_id) {
        const std::string directory_path =
            ReadIndexShardPath(prefix, "keys", file_id);
        const std::string posting_path =
            ReadIndexShardPath(prefix, "postings", file_id);
        directory_fds_[file_id] =
            open(directory_path.c_str(), O_RDONLY | O_CLOEXEC);
        posting_fds_[file_id] =
            open(posting_path.c_str(), O_RDONLY | O_CLOEXEC);
        if (directory_fds_[file_id] < 0 || posting_fds_[file_id] < 0) {
          throw std::runtime_error("cannot open compact read-index shard");
        }
        directory_sizes_[file_id] = FileSize(directory_fds_[file_id]);
        posting_sizes_[file_id] = FileSize(posting_fds_[file_id]);
        total_directory_size_ += directory_sizes_[file_id];
        total_posting_size_ += posting_sizes_[file_id];
        if (map_directories) {
          mapped_directories_[file_id].reset(
              new MappedBytes(directory_path));
          if (mapped_directories_[file_id]->size() !=
              directory_sizes_[file_id]) {
            throw std::runtime_error(
                "read-index directory changed during query");
          }
        }
      }
    } catch (...) {
      Close();
      throw;
    }
  }

  ~ReadIndexShardFiles() { Close(); }
  ReadIndexShardFiles(const ReadIndexShardFiles &) = delete;
  ReadIndexShardFiles &operator=(const ReadIndexShardFiles &) = delete;

  int directory_fd(unsigned file_id) const {
    return directory_fds_[file_id];
  }
  int posting_fd(unsigned file_id) const { return posting_fds_[file_id]; }
  uint64_t directory_size(unsigned file_id) const {
    return directory_sizes_[file_id];
  }
  uint64_t posting_size(unsigned file_id) const {
    return posting_sizes_[file_id];
  }
  uint64_t total_directory_size() const { return total_directory_size_; }
  uint64_t total_posting_size() const { return total_posting_size_; }
  const uint8_t *mapped_directory(unsigned file_id) const {
    return mapped_directories_[file_id]
               ? mapped_directories_[file_id]->data()
               : nullptr;
  }
  void DiscardMappedDirectories() const {
    for (const auto &mapping : mapped_directories_) {
      if (mapping && mapping->size() != 0) {
        DiscardMemoryPages(const_cast<uint8_t *>(mapping->data()),
                           mapping->size());
      }
    }
  }

 private:
  void Close() {
    mapped_directories_.clear();
    for (int &fd : directory_fds_) {
      if (fd >= 0) {
        close(fd);
        fd = -1;
      }
    }
    for (int &fd : posting_fds_) {
      if (fd >= 0) {
        close(fd);
        fd = -1;
      }
    }
  }

  std::vector<int> directory_fds_;
  std::vector<int> posting_fds_;
  std::vector<uint64_t> directory_sizes_;
  std::vector<uint64_t> posting_sizes_;
  std::vector<std::unique_ptr<MappedBytes>> mapped_directories_;
  uint64_t total_directory_size_{0};
  uint64_t total_posting_size_{0};
};

constexpr size_t kDirectoryFenceHeaderBytes = 32u;
const char kDirectoryFenceMagic[8] = {'M', 'H', 'R', 'I', 'D', 'X', 'F', '1'};

struct DirectoryFenceIndex {
  unsigned stride{0};
  std::vector<uint8_t> keys;
  std::vector<uint64_t> bucket_offsets;
};

unsigned NativeDirectoryFenceStride(unsigned directory_entry_bytes) {
  long page_size = sysconf(_SC_PAGESIZE);
  if (page_size <= 0) page_size = 4096;
  return std::max<unsigned>(
      1u, static_cast<unsigned>(page_size) / directory_entry_bytes);
}

bool LoadDirectoryFenceIndex(const std::string &path,
                             const ReadIndexMetadata &metadata,
                             uint64_t directory_file_size,
                             DirectoryFenceIndex *index) {
  const int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (fd < 0) return false;
  const uint64_t file_size = FileSize(fd);
  std::array<uint8_t, kDirectoryFenceHeaderBytes> header{};
  const bool header_ok =
      file_size >= header.size() &&
      PreadFully(fd, header.data(), header.size(), 0) &&
      std::memcmp(header.data(), kDirectoryFenceMagic,
                  sizeof(kDirectoryFenceMagic)) == 0;
  if (!header_ok) {
    close(fd);
    return false;
  }
  const unsigned stride =
      static_cast<unsigned>(LoadLowBytes(header.data() + 8u, 4u));
  const unsigned key_bytes =
      static_cast<unsigned>(LoadLowBytes(header.data() + 12u, 4u));
  const uint64_t stored_keys = LoadLowBytes(header.data() + 16u, 8u);
  const uint64_t stored_directory_size =
      LoadLowBytes(header.data() + 24u, 8u);
  if (stride == 0 || key_bytes != metadata.key_bytes ||
      stored_keys != metadata.num_keys ||
      stored_directory_size != directory_file_size) {
    close(fd);
    return false;
  }
  index->stride = stride;
  index->bucket_offsets.assign(metadata.buckets.size() + 1u, 0);
  for (size_t bucket = 0; bucket < metadata.buckets.size(); ++bucket) {
    index->bucket_offsets[bucket + 1u] =
        index->bucket_offsets[bucket] +
        DivCeiling(metadata.buckets[bucket].directory_count,
                   static_cast<uint64_t>(stride));
  }
  const uint64_t fence_bytes =
      index->bucket_offsets.back() * metadata.key_bytes;
  if (fence_bytes > std::numeric_limits<size_t>::max() ||
      file_size != header.size() + fence_bytes) {
    close(fd);
    return false;
  }
  index->keys.resize(static_cast<size_t>(fence_bytes));
  const bool data_ok =
      index->keys.empty() ||
      PreadFully(fd, index->keys.data(), index->keys.size(), header.size());
  close(fd);
  return data_ok;
}

void PublishDirectoryFenceIndex(const std::string &path, unsigned key_bytes,
                                uint64_t num_keys,
                                uint64_t directory_file_size,
                                const DirectoryFenceIndex &index) {
  std::array<uint8_t, kDirectoryFenceHeaderBytes> header{};
  std::memcpy(header.data(), kDirectoryFenceMagic,
              sizeof(kDirectoryFenceMagic));
  StoreLowBytes(header.data() + 8u, index.stride, 4u);
  StoreLowBytes(header.data() + 12u, key_bytes, 4u);
  StoreLowBytes(header.data() + 16u, num_keys, 8u);
  StoreLowBytes(header.data() + 24u, directory_file_size, 8u);
  const std::string temporary_path =
      path + ".tmp." + std::to_string(static_cast<unsigned long long>(getpid()));
  const int output_fd = open(temporary_path.c_str(),
                             O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0664);
  if (output_fd < 0 ||
      !PwriteFully(output_fd, header.data(), header.size(), 0) ||
      (!index.keys.empty() &&
       !PwriteFully(output_fd, index.keys.data(), index.keys.size(),
                    header.size()))) {
    if (output_fd >= 0) close(output_fd);
    std::remove(temporary_path.c_str());
    throw std::runtime_error("cannot write read-index directory fences");
  }
  close(output_fd);
  if (rename(temporary_path.c_str(), path.c_str()) != 0) {
    std::remove(temporary_path.c_str());
    throw std::runtime_error("cannot publish read-index directory fences");
  }
}

DirectoryFenceIndex EnsureDirectoryFenceIndex(
    const std::string &prefix, const ReadIndexMetadata &metadata,
    const ReadIndexShardFiles &files) {
  const std::string path = prefix + ".ridx.fences";
  const uint64_t directory_file_size = files.total_directory_size();
  DirectoryFenceIndex index;
  if (LoadDirectoryFenceIndex(path, metadata, directory_file_size, &index)) {
    return index;
  }

  SimpleTimer timer;
  timer.start();
  index.stride = NativeDirectoryFenceStride(metadata.directory_entry_bytes);
  index.bucket_offsets.assign(metadata.buckets.size() + 1u, 0);
  for (size_t bucket = 0; bucket < metadata.buckets.size(); ++bucket) {
    index.bucket_offsets[bucket + 1u] =
        index.bucket_offsets[bucket] +
        DivCeiling(metadata.buckets[bucket].directory_count,
                   static_cast<uint64_t>(index.stride));
  }
  const uint64_t fence_bytes = index.bucket_offsets.back() * metadata.key_bytes;
  if (fence_bytes > std::numeric_limits<size_t>::max()) {
    throw std::length_error("read-index directory fences are too large");
  }
  index.keys.resize(static_cast<size_t>(fence_bytes));
  for (size_t bucket = 0; bucket < metadata.buckets.size(); ++bucket) {
    const ReadIndexBucketMeta &meta = metadata.buckets[bucket];
    const uint64_t shard_size = files.directory_size(meta.file_id);
    const uint64_t directory_bytes =
        meta.directory_count * metadata.directory_entry_bytes;
    if (meta.directory_offset > shard_size ||
        directory_bytes > shard_size - meta.directory_offset) {
      throw std::runtime_error("read-index directory bucket exceeds file");
    }
    std::vector<uint8_t> directory(static_cast<size_t>(directory_bytes));
    if (!directory.empty() &&
        !PreadFully(files.directory_fd(meta.file_id), directory.data(),
                    directory.size(),
                    meta.directory_offset)) {
      throw std::runtime_error("cannot build read-index directory fences");
    }
    uint8_t *destination =
        index.keys.data() + index.bucket_offsets[bucket] * metadata.key_bytes;
    for (uint64_t group = 0, fence = 0; group < meta.directory_count;
         group += index.stride, ++fence) {
      std::memcpy(destination + fence * metadata.key_bytes,
                  directory.data() + group * metadata.directory_entry_bytes,
                  metadata.key_bytes);
    }
  }

  PublishDirectoryFenceIndex(path, metadata.key_bytes, metadata.num_keys,
                             directory_file_size, index);
  timer.stop();
  xinfo("Built sparse directory fences: {} entries / {} bytes, stride {}, "
        "{.4} s\n",
        index.bucket_offsets.back(), index.keys.size(), index.stride,
        timer.elapsed());
  return index;
}

using CandidateBuckets = std::vector<std::vector<uint64_t>>;

unsigned ReadOrderBucketCount(const ReadIndexOptions &options,
                              const ReadIndexMetadata &metadata,
                              uint64_t candidate_upper_bound) {
  const uint64_t requested_memory =
      options.memory_bytes >=
              static_cast<double>(std::numeric_limits<uint64_t>::max())
          ? std::numeric_limits<uint64_t>::max()
          : std::max<uint64_t>(1u,
                               static_cast<uint64_t>(options.memory_bytes));
  const uint64_t candidate_bytes =
      candidate_upper_bound >
              std::numeric_limits<uint64_t>::max() / sizeof(uint64_t)
          ? std::numeric_limits<uint64_t>::max()
          : candidate_upper_bound * sizeof(uint64_t);
  const uint64_t remaining =
      requested_memory > candidate_bytes ? requested_memory - candidate_bytes
                                         : 1u;
  const uint64_t workers =
      static_cast<uint64_t>(std::max(1, options.num_threads));
  // Retain half of the non-candidate budget for the directory, flank index,
  // edge collector and allocator metadata.  The other half bounds all mapped
  // read ranges that may be live concurrently.
  const uint64_t bytes_per_live_bucket =
      std::max<uint64_t>(1u, remaining / (workers * 2u));
  const uint64_t read_bytes =
      metadata.parsed_words >
              std::numeric_limits<uint64_t>::max() / sizeof(uint32_t)
          ? std::numeric_limits<uint64_t>::max()
          : metadata.parsed_words * sizeof(uint32_t);
  uint64_t desired = DivCeiling(read_bytes, bytes_per_live_bucket);
  desired = std::max<uint64_t>(desired, workers * 2u);
  desired = std::max<uint64_t>(desired, kMinimumReadOrderBuckets);
  size_t buckets = 1u;
  while (buckets < desired) {
    if (buckets > std::numeric_limits<size_t>::max() / 2u) {
      throw std::length_error("read-order bucket count is too large");
    }
    buckets *= 2u;
  }
  if (buckets > std::numeric_limits<unsigned>::max()) {
    throw std::length_error("read-order bucket count exceeds ID range");
  }
  return static_cast<unsigned>(buckets);
}

void DiscardCandidateReadPages(const std::vector<uint64_t> &candidates,
                               const ReadIndexMetadata &metadata,
                               const uint32_t *read_words) {
  if (candidates.empty()) return;
  const uint64_t first_word = candidates.front() >> metadata.position_bits;
  const uint64_t last_word = candidates.back() >> metadata.position_bits;
  if (last_word >= metadata.parsed_words) return;
  const uint64_t last_end =
      last_word + 1u +
      DivCeiling(static_cast<uint64_t>(read_words[last_word]),
                 static_cast<uint64_t>(SeqPackage::kBasesPerWord));
  if (first_word >= last_end || last_end > metadata.parsed_words) return;
  DiscardMemoryPages(
      const_cast<uint32_t *>(read_words + first_word),
      static_cast<size_t>(last_end - first_word) * sizeof(uint32_t));
}

unsigned ReadWorkingSetThreads(const ReadIndexOptions &options,
                               const ReadIndexMetadata &metadata,
                               const CandidateBuckets &candidate_buckets) {
  uint64_t candidate_bytes = 0;
  for (const auto &bucket : candidate_buckets) {
    const uint64_t bytes = bucket.size() * sizeof(uint64_t);
    candidate_bytes =
        bytes > std::numeric_limits<uint64_t>::max() - candidate_bytes
            ? std::numeric_limits<uint64_t>::max()
            : candidate_bytes + bytes;
  }
  const uint64_t requested_memory =
      options.memory_bytes >=
              static_cast<double>(std::numeric_limits<uint64_t>::max())
          ? std::numeric_limits<uint64_t>::max()
          : static_cast<uint64_t>(options.memory_bytes);
  const uint64_t bucket_words =
      DivCeiling(metadata.parsed_words,
                 static_cast<uint64_t>(candidate_buckets.size()));
  const uint64_t bucket_bytes = std::max<uint64_t>(
      1u, bucket_words * static_cast<uint64_t>(sizeof(uint32_t)));
  // Candidate storage is persistent during verification/replay.  Split the
  // remainder between mapped read pages and structures not represented here:
  // active/flank tables, the edge collector, allocator metadata, and uneven
  // bucket boundaries.  This makes the user memory value a real working-set
  // budget instead of consuming all of it with only the easiest-to-count
  // arrays.
  const uint64_t remaining =
      requested_memory > candidate_bytes ? requested_memory - candidate_bytes
                                         : bucket_bytes;
  const uint64_t read_page_budget =
      std::max<uint64_t>(bucket_bytes, remaining / 2u);
  const uint64_t budget_threads =
      std::max<uint64_t>(1u, read_page_budget / bucket_bytes);
  return static_cast<unsigned>(std::max<uint64_t>(
      1u, std::min<uint64_t>(std::max(1, options.num_threads),
                             budget_threads)));
}

class MappedReadView {
 public:
  using word_type = uint32_t;
  MappedReadView(const uint32_t *sequence, unsigned length)
      : sequence_(sequence), length_(length) {}
  unsigned length() const { return length_; }
  uint8_t base_at(unsigned position) const {
    assert(position < length_);
    return PackedBase(sequence_, position);
  }
  std::pair<const uint32_t *, unsigned> raw_address() const {
    return std::make_pair(sequence_, 0u);
  }

 private:
  const uint32_t *sequence_;
  unsigned length_;
};

template <class KmerType>
void FeedReplayFlankIndex(const ReadIndexOptions &options,
                          ContigFlankIndex<KmerType> *index) {
  const std::string files[] = {options.contig_file, options.bubble_file};
  for (const std::string &file : files) {
    if (file.empty()) continue;
    AsyncContigReader reader(file);
    while (true) {
      auto &batch = reader.Next();
      if (batch.first.seq_count() == 0) break;
      index->FeedBatchContigs(batch.first, batch.second);
    }
  }
  index->Finalize();
}

template <class FlankKmerType, class NextKmerType>
bool ReplayCandidateEdgesForNextType(
    const ReadIndexOptions &options, const ReadIndexMetadata &metadata,
    const CandidateBuckets &candidate_buckets, const uint32_t *read_words,
    uint64_t window_filter_hits, unsigned working_threads,
    const ContigFlankIndex<FlankKmerType> &flank_index,
    uint64_t *aligned_reads_out) {
  if (NextKmerType::max_size() <
      static_cast<unsigned>(options.kmer_k + options.step + 1)) {
    return false;
  }
  KmerCollector<NextKmerType> collector(
      options.kmer_k + options.step + 1u, options.edge_output_prefix);
  const uint64_t position_mask =
      (uint64_t{1} << metadata.position_bits) - 1u;
  uint64_t candidate_reads = 0;
  uint64_t aligned_reads = 0;
  SimpleTimer replay_timer;
  replay_timer.start();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : candidate_reads, aligned_reads) \
    num_threads(working_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(candidate_buckets.size()); ++bucket) {
    const auto &candidates = candidate_buckets[static_cast<size_t>(bucket)];
    std::vector<uint32_t> positions;
    std::vector<uint32_t> kmer_state;
    for (size_t begin = 0; begin < candidates.size();) {
      const uint64_t word_offset =
          candidates[begin] >> metadata.position_bits;
      size_t end = begin + 1u;
      while (end < candidates.size() &&
             (candidates[end] >> metadata.position_bits) == word_offset) {
        ++end;
      }
      const unsigned read_len = read_words[word_offset];
      positions.clear();
      positions.reserve(end - begin);
      for (size_t i = begin; i < end; ++i) {
        const uint32_t start =
            static_cast<uint32_t>(candidates[i] & position_mask);
        if (static_cast<uint64_t>(start) + metadata.window_len > read_len) {
          continue;
        }
        positions.push_back(start);
      }
      if (!positions.empty()) {
        MappedReadView read(read_words + word_offset + 1u, read_len);
        aligned_reads += flank_index.FindNextKmersFromReadCandidates(
            read, positions, &collector, &kmer_state);
        ++candidate_reads;
      }
      begin = end;
    }
    DiscardCandidateReadPages(candidates, metadata, read_words);
  }
  replay_timer.stop();
  SimpleTimer flush_timer;
  flush_timer.start();
  collector.FlushToFile();
  flush_timer.stop();
  xinfo("Exact index replay: {} window-filter hits ({} bases) in {} "
        "candidate reads, "
        "{} aligned, {} iterative edges; state replay {.4}, "
        "edge output {.4} s\n",
        window_filter_hits, metadata.window_len, candidate_reads, aligned_reads,
        collector.collection().size(), replay_timer.elapsed(),
        flush_timer.elapsed());
  *aligned_reads_out = aligned_reads;
  return true;
}

template <class FlankKmerType>
bool ReplayCandidateEdgesForFlankType(
    const ReadIndexOptions &options, const ReadIndexMetadata &metadata,
    const CandidateBuckets &candidate_buckets, const uint32_t *read_words,
    uint64_t window_filter_hits, unsigned working_threads,
    PreparedFlankIndexBase *prepared_index,
    uint64_t *aligned_reads_out) {
  if (FlankKmerType::max_size() <
      static_cast<unsigned>(options.kmer_k + 1)) {
    return false;
  }
  std::unique_ptr<ContigFlankIndex<FlankKmerType>> fallback_index;
  ContigFlankIndex<FlankKmerType> *flank_index = nullptr;
  if (prepared_index != nullptr) {
    auto *typed =
        dynamic_cast<PreparedFlankIndex<FlankKmerType> *>(prepared_index);
    if (typed == nullptr) return false;
    flank_index = &typed->index;
  } else {
    fallback_index.reset(
        new ContigFlankIndex<FlankKmerType>(options.kmer_k, options.step));
    FeedReplayFlankIndex(options, fallback_index.get());
    flank_index = fallback_index.get();
  }
#define TRY_NEXT_KMER(...)                                                    \
  if (ReplayCandidateEdgesForNextType<FlankKmerType, __VA_ARGS__>(           \
          options, metadata, candidate_buckets, read_words,                   \
          window_filter_hits, working_threads, *flank_index,                 \
          aligned_reads_out)) {                                              \
    return true;                                                             \
  }
  TRY_NEXT_KMER(Kmer<1, uint64_t>)
  TRY_NEXT_KMER(Kmer<3, uint32_t>)
  TRY_NEXT_KMER(Kmer<2, uint64_t>)
  TRY_NEXT_KMER(Kmer<5, uint32_t>)
  TRY_NEXT_KMER(Kmer<3, uint64_t>)
  TRY_NEXT_KMER(Kmer<7, uint32_t>)
  TRY_NEXT_KMER(Kmer<4, uint64_t>)
  TRY_NEXT_KMER(Kmer<kUint32PerKmerMaxK, uint32_t>)
#undef TRY_NEXT_KMER
  throw std::logic_error("next k exceeds replay k-mer capacity");
}

void ReplayCandidateEdges(const ReadIndexOptions &options,
                          const ReadIndexMetadata &metadata,
                          const CandidateBuckets &candidate_buckets,
                          const uint32_t *read_words,
                          uint64_t window_filter_hits,
                          unsigned working_threads,
                          PreparedFlankIndexBase *prepared_index,
                          uint64_t *aligned_reads_out) {
#define TRY_FLANK_KMER(...)                                                   \
  if (ReplayCandidateEdgesForFlankType<__VA_ARGS__>(                         \
          options, metadata, candidate_buckets, read_words,                  \
          window_filter_hits, working_threads, prepared_index,               \
          aligned_reads_out)) {                                              \
    return;                                                                  \
  }
  TRY_FLANK_KMER(Kmer<1, uint64_t>)
  TRY_FLANK_KMER(Kmer<3, uint32_t>)
  TRY_FLANK_KMER(Kmer<2, uint64_t>)
  TRY_FLANK_KMER(Kmer<5, uint32_t>)
  TRY_FLANK_KMER(Kmer<3, uint64_t>)
  TRY_FLANK_KMER(Kmer<7, uint32_t>)
  TRY_FLANK_KMER(Kmer<4, uint64_t>)
  TRY_FLANK_KMER(Kmer<kUint32PerKmerMaxK, uint32_t>)
#undef TRY_FLANK_KMER
  throw std::logic_error("current k exceeds replay k-mer capacity");
}

// Type-erased bridge used only at the read-index orchestration boundary.
// The hot exact read kernel remains a statically dispatched template; the
// single predictable virtual call occurs once per gated read, not per base or
// k-mer.
class FusedEdgeReplayBase {
 public:
  virtual ~FusedEdgeReplayBase() = default;
  virtual bool ReplayRead(const uint32_t *sequence, unsigned length,
                          const std::vector<uint32_t> &candidate_positions,
                          std::vector<uint32_t> *state_scratch,
                          uint64_t *generated_edges) = 0;
  virtual uint64_t CommitBatch() = 0;
  virtual size_t candidate_record_bytes() const = 0;
  virtual uint64_t batch_unique_candidates() const = 0;
  virtual size_t unique_edges() const = 0;
  virtual void Flush() = 0;
};

template <class FlankKmerType, class NextKmerType>
class FusedEdgeReplay final : public FusedEdgeReplayBase {
 public:
  static void *operator new(std::size_t bytes) {
    void *memory = nullptr;
    if (posix_memalign(&memory, alignof(FusedEdgeReplay), bytes) != 0) {
      throw std::bad_alloc();
    }
    return memory;
  }
  static void operator delete(void *memory) noexcept { free(memory); }

  FusedEdgeReplay(const ReadIndexOptions &options,
                  const ContigFlankIndex<FlankKmerType> &index)
      : index_(index),
        collector_(options.kmer_k + options.step + 1u,
                   options.edge_output_prefix, true) {}

  bool ReplayRead(const uint32_t *sequence, unsigned length,
                  const std::vector<uint32_t> &candidate_positions,
                  std::vector<uint32_t> *state_scratch,
                  uint64_t *generated_edges) override {
    return index_.FindNextKmersFromRead(
        MappedReadView(sequence, length), &collector_, state_scratch,
        generated_edges, &candidate_positions);
  }

  uint64_t CommitBatch() override { return collector_.CommitBatch(); }
  size_t candidate_record_bytes() const override {
    return sizeof(typename KmerCollector<NextKmerType>::kmer_plus);
  }
  uint64_t batch_unique_candidates() const override {
    return collector_.partitioned_unique_candidates();
  }
  size_t unique_edges() const override { return collector_.size(); }
  void Flush() override { collector_.FlushToFile(); }

 private:
  const ContigFlankIndex<FlankKmerType> &index_;
  KmerCollector<NextKmerType> collector_;
};

template <class FlankKmerType>
std::unique_ptr<FusedEdgeReplayBase> MakeFusedEdgeReplayForFlankType(
    const ReadIndexOptions &options,
    PreparedFlankIndexBase *prepared_index) {
  auto *typed =
      dynamic_cast<PreparedFlankIndex<FlankKmerType> *>(prepared_index);
  if (typed == nullptr) return std::unique_ptr<FusedEdgeReplayBase>();

#define TRY_FUSED_NEXT_KMER(...)                                             \
  if (__VA_ARGS__::max_size() >=                                            \
      static_cast<unsigned>(options.kmer_k + options.step + 1)) {           \
    return std::unique_ptr<FusedEdgeReplayBase>(                            \
        new FusedEdgeReplay<FlankKmerType, __VA_ARGS__>(options,            \
                                                        typed->index));      \
  }
  TRY_FUSED_NEXT_KMER(Kmer<1, uint64_t>)
  TRY_FUSED_NEXT_KMER(Kmer<3, uint32_t>)
  TRY_FUSED_NEXT_KMER(Kmer<2, uint64_t>)
  TRY_FUSED_NEXT_KMER(Kmer<5, uint32_t>)
  TRY_FUSED_NEXT_KMER(Kmer<3, uint64_t>)
  TRY_FUSED_NEXT_KMER(Kmer<7, uint32_t>)
  TRY_FUSED_NEXT_KMER(Kmer<4, uint64_t>)
  TRY_FUSED_NEXT_KMER(Kmer<kUint32PerKmerMaxK, uint32_t>)
#undef TRY_FUSED_NEXT_KMER
  throw std::logic_error("next k exceeds fused replay k-mer capacity");
}

std::unique_ptr<FusedEdgeReplayBase> MakeFusedEdgeReplay(
    const ReadIndexOptions &options,
    PreparedFlankIndexBase *prepared_index) {
  std::unique_ptr<FusedEdgeReplayBase> replay;
#define TRY_FUSED_FLANK_KMER(...)                                            \
  replay = MakeFusedEdgeReplayForFlankType<__VA_ARGS__>(                    \
      options, prepared_index);                                             \
  if (replay) return replay;
  TRY_FUSED_FLANK_KMER(Kmer<1, uint64_t>)
  TRY_FUSED_FLANK_KMER(Kmer<3, uint32_t>)
  TRY_FUSED_FLANK_KMER(Kmer<2, uint64_t>)
  TRY_FUSED_FLANK_KMER(Kmer<5, uint32_t>)
  TRY_FUSED_FLANK_KMER(Kmer<3, uint64_t>)
  TRY_FUSED_FLANK_KMER(Kmer<7, uint32_t>)
  TRY_FUSED_FLANK_KMER(Kmer<4, uint64_t>)
  TRY_FUSED_FLANK_KMER(Kmer<kUint32PerKmerMaxK, uint32_t>)
#undef TRY_FUSED_FLANK_KMER
  throw std::logic_error("current k exceeds fused replay k-mer capacity");
}

void PublishLocalCandidates(const ReadIndexOptions &options,
                            const ReadIndexMetadata &metadata,
                            CandidateBuckets *candidate_buckets) {
  const uint64_t position_mask =
      (uint64_t{1} << metadata.position_bits) - 1u;
  std::atomic<bool> valid(true);
#pragma omp parallel for schedule(static) num_threads(options.num_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(candidate_buckets->size()); ++bucket) {
    std::vector<uint64_t> &candidates =
        (*candidate_buckets)[static_cast<size_t>(bucket)];
    size_t output = 0;
    uint64_t previous = std::numeric_limits<uint64_t>::max();
    for (uint64_t candidate : candidates) {
      const uint64_t word_offset = candidate >> metadata.position_bits;
      if ((candidate & position_mask) + metadata.window_len >
              std::numeric_limits<uint32_t>::max() ||
          word_offset >= metadata.parsed_words) {
        valid.store(false, std::memory_order_relaxed);
        continue;
      }
      if (word_offset != previous) {
        candidates[output++] = word_offset;
        previous = word_offset;
      }
    }
    candidates.resize(output);
  }
  if (!valid.load(std::memory_order_relaxed)) {
    throw std::runtime_error("invalid exact candidate while publishing local "
                             "read offsets");
  }

  std::vector<uint64_t> bucket_offsets(candidate_buckets->size() + 1u, 0u);
  for (size_t bucket = 0; bucket < candidate_buckets->size(); ++bucket) {
    const uint64_t count = (*candidate_buckets)[bucket].size();
    if (count > std::numeric_limits<uint64_t>::max() -
                    bucket_offsets[bucket]) {
      throw std::length_error("local candidate count overflow");
    }
    bucket_offsets[bucket + 1u] = bucket_offsets[bucket] + count;
  }
  const uint64_t candidate_count = bucket_offsets.back();
  if (candidate_count >
      (std::numeric_limits<uint64_t>::max() -
       sizeof(LocalCandidateFileHeader)) /
          sizeof(uint64_t)) {
    throw std::length_error("local candidate file is too large");
  }
  const uint64_t output_bytes = sizeof(LocalCandidateFileHeader) +
                                candidate_count * sizeof(uint64_t);
  if (output_bytes >
      static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
    throw std::length_error("local candidate file exceeds file-size limit");
  }

  const std::string temporary = options.local_candidate_output + ".tmp." +
                                std::to_string(static_cast<uint64_t>(getpid()));
  const int fd = open(temporary.c_str(),
                      O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0666);
  if (fd < 0 || ftruncate(fd, static_cast<off_t>(output_bytes)) != 0) {
    if (fd >= 0) close(fd);
    std::remove(temporary.c_str());
    throw std::runtime_error("cannot create local candidate output");
  }
  const LocalCandidateFileHeader header = MakeLocalCandidateFileHeader(
      metadata.parsed_words, metadata.num_reads, candidate_count);
  std::atomic<bool> write_ok(
      PwriteFully(fd, &header, sizeof(header), 0u));
#pragma omp parallel for schedule(static) num_threads(options.num_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(candidate_buckets->size()); ++bucket) {
    const std::vector<uint64_t> &candidates =
        (*candidate_buckets)[static_cast<size_t>(bucket)];
    if (candidates.empty() || !write_ok.load(std::memory_order_relaxed)) {
      continue;
    }
    const uint64_t offset = sizeof(LocalCandidateFileHeader) +
                            bucket_offsets[static_cast<size_t>(bucket)] *
                                sizeof(uint64_t);
    if (!PwriteFully(fd, candidates.data(),
                     candidates.size() * sizeof(uint64_t), offset)) {
      write_ok.store(false, std::memory_order_relaxed);
    }
  }
  const bool close_ok = close(fd) == 0;
  if (!write_ok.load(std::memory_order_relaxed) || !close_ok ||
      std::rename(temporary.c_str(), options.local_candidate_output.c_str()) !=
          0) {
    std::remove(temporary.c_str());
    throw std::runtime_error("cannot publish local candidate output");
  }
  xinfo("Published {} exact candidate reads ({} bytes) for local assembly: "
        "{}\n",
        candidate_count, output_bytes,
        options.local_candidate_output.c_str());
}

// Large active endpoint sets can expand to billions of candidate windows
// before exact verification, even though the final set contains only read
// offsets.  Publish read-order bucket waves directly instead of retaining all
// filtered candidates until the end.  Waves and buckets are monotonically
// ordered by packed-read offset, so this is byte-for-byte the same sorted,
// unique stream produced by PublishLocalCandidates().
class StreamingLocalCandidateWriter {
 public:
  StreamingLocalCandidateWriter(const ReadIndexOptions &options,
                                const ReadIndexMetadata &metadata)
      : output_(options.local_candidate_output),
        temporary_(output_ + ".tmp." +
                   std::to_string(static_cast<uint64_t>(getpid()))),
        metadata_(metadata),
        num_threads_(options.num_threads) {
    fd_ = open(temporary_.c_str(),
               O_CREAT | O_TRUNC | O_WRONLY | O_CLOEXEC, 0666);
    if (fd_ < 0) {
      throw std::runtime_error("cannot create streamed local candidate output");
    }
  }

  StreamingLocalCandidateWriter(const StreamingLocalCandidateWriter &) =
      delete;
  StreamingLocalCandidateWriter &operator=(
      const StreamingLocalCandidateWriter &) = delete;

  ~StreamingLocalCandidateWriter() {
    if (fd_ >= 0) close(fd_);
    if (!published_) std::remove(temporary_.c_str());
  }

  void AppendWave(CandidateBuckets *candidate_buckets) {
    const uint64_t position_mask =
        (uint64_t{1} << metadata_.position_bits) - 1u;
    std::atomic<bool> valid(true);
#pragma omp parallel for schedule(static) num_threads(num_threads_)
    for (int64_t bucket = 0;
         bucket < static_cast<int64_t>(candidate_buckets->size()); ++bucket) {
      std::vector<uint64_t> &candidates =
          (*candidate_buckets)[static_cast<size_t>(bucket)];
      size_t output = 0;
      uint64_t previous = std::numeric_limits<uint64_t>::max();
      for (uint64_t candidate : candidates) {
        const uint64_t word_offset = candidate >> metadata_.position_bits;
        if ((candidate & position_mask) + metadata_.window_len >
                std::numeric_limits<uint32_t>::max() ||
            word_offset >= metadata_.parsed_words) {
          valid.store(false, std::memory_order_relaxed);
          continue;
        }
        if (word_offset != previous) {
          candidates[output++] = word_offset;
          previous = word_offset;
        }
      }
      candidates.resize(output);
    }
    if (!valid.load(std::memory_order_relaxed)) {
      throw std::runtime_error(
          "invalid exact candidate while streaming local read offsets");
    }

    std::vector<uint64_t> offsets(candidate_buckets->size() + 1u, 0u);
    for (size_t bucket = 0; bucket < candidate_buckets->size(); ++bucket) {
      const uint64_t count = (*candidate_buckets)[bucket].size();
      if (count > std::numeric_limits<uint64_t>::max() - offsets[bucket]) {
        throw std::length_error("local candidate count overflow");
      }
      offsets[bucket + 1u] = offsets[bucket] + count;
    }
    if (offsets.back() >
        (std::numeric_limits<uint64_t>::max() - candidate_count_)) {
      throw std::length_error("local candidate count overflow");
    }
    std::atomic<bool> write_ok(true);
#pragma omp parallel for schedule(static) num_threads(num_threads_)
    for (int64_t bucket = 0;
         bucket < static_cast<int64_t>(candidate_buckets->size()); ++bucket) {
      const std::vector<uint64_t> &candidates =
          (*candidate_buckets)[static_cast<size_t>(bucket)];
      if (candidates.empty() || !write_ok.load(std::memory_order_relaxed)) {
        continue;
      }
      const uint64_t item_offset =
          candidate_count_ + offsets[static_cast<size_t>(bucket)];
      if (item_offset >
          (static_cast<uint64_t>(std::numeric_limits<off_t>::max()) -
           sizeof(LocalCandidateFileHeader)) /
              sizeof(uint64_t) ||
          !PwriteFully(fd_, candidates.data(),
                       candidates.size() * sizeof(uint64_t),
                       sizeof(LocalCandidateFileHeader) +
                           item_offset * sizeof(uint64_t))) {
        write_ok.store(false, std::memory_order_relaxed);
      }
    }
    if (!write_ok.load(std::memory_order_relaxed)) {
      throw std::runtime_error("cannot stream local candidate output");
    }
    candidate_count_ += offsets.back();
  }

  void Finish() {
    if (candidate_count_ >
        (static_cast<uint64_t>(std::numeric_limits<off_t>::max()) -
         sizeof(LocalCandidateFileHeader)) /
            sizeof(uint64_t)) {
      throw std::length_error("local candidate file is too large");
    }
    const uint64_t output_bytes = sizeof(LocalCandidateFileHeader) +
                                  candidate_count_ * sizeof(uint64_t);
    const LocalCandidateFileHeader header = MakeLocalCandidateFileHeader(
        metadata_.parsed_words, metadata_.num_reads, candidate_count_);
    if (ftruncate(fd_, static_cast<off_t>(output_bytes)) != 0 ||
        !PwriteFully(fd_, &header, sizeof(header), 0u) || close(fd_) != 0) {
      fd_ = -1;
      throw std::runtime_error("cannot finish streamed local candidate output");
    }
    fd_ = -1;
    if (std::rename(temporary_.c_str(), output_.c_str()) != 0) {
      throw std::runtime_error("cannot publish streamed local candidate output");
    }
    published_ = true;
    xinfo("Published {} exact candidate reads ({} bytes) for local assembly: "
          "{}\n",
          candidate_count_, output_bytes, output_.c_str());
  }

 private:
  std::string output_;
  std::string temporary_;
  const ReadIndexMetadata &metadata_;
  int num_threads_{1};
  int fd_{-1};
  uint64_t candidate_count_{0};
  bool published_{false};
};

struct StreamedLocalCandidateStats {
  uint64_t unique_candidates{0};
  uint64_t window_filter_hits{0};
  double count_seconds{0};
  double gather_seconds{0};
  double sort_seconds{0};
  double verify_seconds{0};
  unsigned verify_threads{1};
};

StreamedLocalCandidateStats QueryLocalCandidatesInReadOrderWaves(
    const ReadIndexOptions &options, const ReadIndexMetadata &metadata,
    ReadIndexShardFiles *index_files,
    const std::vector<std::vector<MatchedPostingGroup>> &groups_by_file,
    const uint32_t *read_words,
    const CompactLocalEndpointQueries &local_queries,
    uint64_t candidate_upper_bound, uint64_t requested_memory,
    unsigned num_read_order_buckets) {
  StreamedLocalCandidateStats stats;
  stats.verify_threads = static_cast<unsigned>(std::max(1, options.num_threads));
  const uint64_t position_mask =
      (uint64_t{1} << metadata.position_bits) - 1u;
  const unsigned locator_value_bits =
      metadata.position_bits +
      BitsNeeded(metadata.parsed_words == 0 ? 0 : metadata.parsed_words - 1u);
  const uint64_t locator_value_mask = LowBitMask(locator_value_bits);
  const uint64_t read_order_bucket_words = std::max<uint64_t>(
      1u, DivCeiling(metadata.parsed_words,
                     static_cast<uint64_t>(num_read_order_buckets)));
  const uint64_t candidate_wave_bytes = std::max<uint64_t>(
      uint64_t{64} << 20u,
      requested_memory > (uint64_t{256} << 20u)
          ? requested_memory * 2u / 5u
          : requested_memory / 2u);
  const uint64_t candidate_wave_items =
      std::max<uint64_t>(1u, candidate_wave_bytes / sizeof(uint64_t));
  const int worker_count = std::max(1, options.num_threads);

  if (static_cast<uint64_t>(worker_count) * num_read_order_buckets >
      std::numeric_limits<size_t>::max() / sizeof(uint64_t)) {
    throw std::length_error("streamed candidate histogram is too large");
  }
  std::vector<uint64_t> thread_counts(
      static_cast<size_t>(worker_count) * num_read_order_buckets, 0u);
  std::atomic<bool> scan_ok(true);
  SimpleTimer count_timer;
  count_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(metadata.num_files);
       ++file_id) {
    const unsigned file = static_cast<unsigned>(file_id);
    uint64_t *local_counts =
        thread_counts.data() + static_cast<size_t>(omp_get_thread_num()) *
                                   num_read_order_buckets;
    std::vector<uint8_t> posting_buffer;
    for (const MatchedPostingGroup &group : groups_by_file[file]) {
      if (!scan_ok.load(std::memory_order_relaxed)) break;
      const uint64_t posting_file_size = index_files->posting_size(file);
      const uint64_t group_bytes = group.count * metadata.locator_bytes;
      if (group.posting_byte_offset > posting_file_size ||
          group_bytes > posting_file_size - group.posting_byte_offset) {
        scan_ok.store(false, std::memory_order_relaxed);
        break;
      }
      posting_buffer.resize(static_cast<size_t>(group_bytes));
      if (!posting_buffer.empty() &&
          !PreadFully(index_files->posting_fd(file), posting_buffer.data(),
                      posting_buffer.size(), group.posting_byte_offset)) {
        scan_ok.store(false, std::memory_order_relaxed);
        break;
      }
      const uint8_t *postings = posting_buffer.data();
      for (uint64_t posting = 0; posting < group.count; ++posting) {
        const uint64_t encoded = LoadLowBytes(
            postings + posting * metadata.locator_bytes,
            metadata.locator_bytes);
        const uint64_t locator = encoded & locator_value_mask;
        const uint64_t anchor_pos = locator & position_mask;
        const uint64_t word_offset = locator >> metadata.position_bits;
        uint64_t offsets = group.offset_mask;
        while (offsets != 0u) {
          const unsigned offset =
              static_cast<unsigned>(__builtin_ctzll(offsets));
          offsets &= offsets - 1u;
          if (anchor_pos < offset || word_offset >= metadata.parsed_words) {
            continue;
          }
          const unsigned bucket = std::min(
              static_cast<unsigned>(word_offset / read_order_bucket_words),
              num_read_order_buckets - 1u);
          ++local_counts[bucket];
        }
      }
    }
  }
  if (!scan_ok.load(std::memory_order_relaxed)) {
    throw std::runtime_error("failed while counting streamed candidates");
  }
  std::vector<uint64_t> bucket_counts(num_read_order_buckets, 0u);
#pragma omp parallel for schedule(static) num_threads(options.num_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(num_read_order_buckets); ++bucket) {
    uint64_t count = 0u;
    for (int thread = 0; thread < worker_count; ++thread) {
      const uint64_t value =
          thread_counts[static_cast<size_t>(thread) * num_read_order_buckets +
                        static_cast<size_t>(bucket)];
      if (value > std::numeric_limits<uint64_t>::max() - count) {
        scan_ok.store(false, std::memory_order_relaxed);
        break;
      }
      count += value;
    }
    bucket_counts[static_cast<size_t>(bucket)] = count;
  }
  std::vector<uint64_t>().swap(thread_counts);
  if (!scan_ok.load(std::memory_order_relaxed)) {
    throw std::length_error("streamed candidate count overflow");
  }

  std::vector<std::pair<unsigned, unsigned>> waves;
  for (unsigned begin = 0; begin < num_read_order_buckets;) {
    unsigned end = begin;
    uint64_t items = 0u;
    while (end < num_read_order_buckets) {
      const uint64_t count = bucket_counts[end];
      if (count > candidate_wave_items) {
        throw std::runtime_error(
            "one read-order bucket exceeds streamed candidate budget");
      }
      if (end != begin && count > candidate_wave_items - items) break;
      items += count;
      ++end;
    }
    waves.emplace_back(begin, end);
    begin = end;
  }
  count_timer.stop();
  stats.count_seconds = count_timer.elapsed();
  xinfo("Streamed local candidate plan: {} read-order buckets in {} waves, "
        "at most {.3} GiB candidates per wave; counting {.4} s\n",
        num_read_order_buckets, waves.size(),
        static_cast<double>(candidate_wave_items * sizeof(uint64_t)) /
            static_cast<double>(uint64_t{1} << 30u),
        stats.count_seconds);

  const uint64_t workers = static_cast<uint64_t>(worker_count);
  const uint64_t requested_stage_items =
      requested_memory / (uint64_t{64} * sizeof(uint64_t));
  const size_t stage_limit = static_cast<size_t>(std::max<uint64_t>(
      1024u, std::min<uint64_t>(DivCeiling(candidate_upper_bound, workers),
                                requested_stage_items / workers)));
  StreamingLocalCandidateWriter writer(options, metadata);
  SimpleTimer gather_timer;
  SimpleTimer sort_timer;
  SimpleTimer verify_timer;
  uint64_t window_filter_hits = 0u;

  for (const auto &wave : waves) {
    const unsigned wave_begin = wave.first;
    const unsigned wave_end = wave.second;
    const unsigned wave_buckets = wave_end - wave_begin;
    CandidateBuckets candidates(wave_buckets);
    for (unsigned bucket = wave_begin; bucket < wave_end; ++bucket) {
      if (bucket_counts[bucket] > std::numeric_limits<size_t>::max()) {
        throw std::length_error("candidate bucket exceeds address space");
      }
      candidates[bucket - wave_begin].reserve(
          static_cast<size_t>(bucket_counts[bucket]));
    }
    std::vector<std::mutex> bucket_mutexes(wave_buckets);
    scan_ok.store(true, std::memory_order_relaxed);
    gather_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
    for (int file_id = 0; file_id < static_cast<int>(metadata.num_files);
         ++file_id) {
      const unsigned file = static_cast<unsigned>(file_id);
      std::vector<uint8_t> posting_buffer;
      CandidateBuckets staged(wave_buckets);
      size_t staged_items = 0u;
      const auto flush_bucket = [&](unsigned bucket) {
        auto &source = staged[bucket];
        if (source.empty()) return;
        {
          std::lock_guard<std::mutex> lock(bucket_mutexes[bucket]);
          candidates[bucket].insert(candidates[bucket].end(), source.begin(),
                                    source.end());
        }
        source.clear();
      };
      const auto flush_all = [&]() {
        for (unsigned bucket = 0; bucket < wave_buckets; ++bucket) {
          flush_bucket(bucket);
        }
        staged_items = 0u;
      };
      for (const MatchedPostingGroup &group : groups_by_file[file]) {
        if (!scan_ok.load(std::memory_order_relaxed)) break;
        const uint64_t posting_file_size = index_files->posting_size(file);
        const uint64_t group_bytes = group.count * metadata.locator_bytes;
        if (group.posting_byte_offset > posting_file_size ||
            group_bytes > posting_file_size - group.posting_byte_offset) {
          scan_ok.store(false, std::memory_order_relaxed);
          break;
        }
        posting_buffer.resize(static_cast<size_t>(group_bytes));
        if (!posting_buffer.empty() &&
            !PreadFully(index_files->posting_fd(file), posting_buffer.data(),
                        posting_buffer.size(), group.posting_byte_offset)) {
          scan_ok.store(false, std::memory_order_relaxed);
          break;
        }
        const uint8_t *postings = posting_buffer.data();
        for (uint64_t posting = 0; posting < group.count; ++posting) {
          const uint64_t encoded = LoadLowBytes(
              postings + posting * metadata.locator_bytes,
              metadata.locator_bytes);
          const uint64_t locator = encoded & locator_value_mask;
          const uint64_t anchor_pos = locator & position_mask;
          const uint64_t word_offset = locator >> metadata.position_bits;
          uint64_t offsets = group.offset_mask;
          while (offsets != 0u) {
            const unsigned offset =
                static_cast<unsigned>(__builtin_ctzll(offsets));
            offsets &= offsets - 1u;
            if (anchor_pos < offset || word_offset >= metadata.parsed_words) {
              continue;
            }
            const unsigned bucket = std::min(
                static_cast<unsigned>(word_offset / read_order_bucket_words),
                num_read_order_buckets - 1u);
            if (bucket < wave_begin || bucket >= wave_end) continue;
            const uint64_t candidate =
                (word_offset << metadata.position_bits) |
                (anchor_pos - offset);
            staged[bucket - wave_begin].push_back(candidate);
            if (++staged_items >= stage_limit) flush_all();
          }
        }
      }
      flush_all();
    }
    gather_timer.stop();
    if (!scan_ok.load(std::memory_order_relaxed)) {
      throw std::runtime_error("failed while gathering streamed candidates");
    }

    sort_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
    for (int64_t bucket = 0;
         bucket < static_cast<int64_t>(candidates.size()); ++bucket) {
      auto &values = candidates[static_cast<size_t>(bucket)];
      kmlib::kmsort(values.begin(), values.end());
      values.erase(std::unique(values.begin(), values.end()), values.end());
    }
    sort_timer.stop();
    for (const auto &bucket : candidates) {
      stats.unique_candidates += bucket.size();
    }

    verify_timer.start();
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : window_filter_hits) num_threads(options.num_threads)
    for (int64_t bucket = 0;
         bucket < static_cast<int64_t>(candidates.size()); ++bucket) {
      auto &values = candidates[static_cast<size_t>(bucket)];
      size_t output = 0u;
      for (uint64_t candidate : values) {
        const uint64_t start = candidate & position_mask;
        const uint64_t word_offset = candidate >> metadata.position_bits;
        if (start + metadata.window_len > read_words[word_offset]) continue;
        AnchorWindow window;
        window.InitFromPtr(read_words + word_offset + 1u,
                           static_cast<unsigned>(start), metadata.window_len);
        if (local_queries.Contains(window)) {
          values[output++] = candidate;
          ++window_filter_hits;
        }
      }
      DiscardCandidateReadPages(values, metadata, read_words);
      values.resize(output);
      if (values.size() < values.capacity()) {
        DiscardMemoryPages(values.data() + values.size(),
                           (values.capacity() - values.size()) *
                               sizeof(uint64_t));
      }
    }
    verify_timer.stop();
    writer.AppendWave(&candidates);
  }
  writer.Finish();
  stats.gather_seconds = gather_timer.elapsed();
  stats.sort_seconds = sort_timer.elapsed();
  stats.verify_seconds = verify_timer.elapsed();
  stats.window_filter_hits = window_filter_hits;
  return stats;
}

void ProfileIndexQuery(const ReadIndexOptions &options) {
  ReadIndexMetadata metadata = LoadReadIndexMetadata(options.index_prefix);
  if (metadata.window_len > AnchorWindow::max_size()) {
    throw std::logic_error("query-profile window exceeds exact window type");
  }
  ReadIndexOptions query_options = options;
  query_options.anchor_len = metadata.anchor_len;
  query_options.window_len = metadata.window_len;

  ActiveQuerySet active_queries;
  active_queries.ConfigureContext(metadata.locator_context_bits,
                                  metadata.anchor_len,
                                  metadata.window_len,
                                  metadata.locator_context_distance);
  ActiveWindowSet active_windows;
  CompactLocalEndpointQueries local_queries;
  std::unique_ptr<PreparedFlankIndexBase> prepared_flank_index;
  const bool local_candidate_mode =
      !options.local_candidate_output.empty();
  SimpleTimer active_timer;
  active_timer.start();
  if (local_candidate_mode) {
    BuildLocalEndpointQueries(query_options, &local_queries);
  } else {
    BuildActiveQueries(query_options, &active_queries, &active_windows,
                       &prepared_flank_index);
  }
  active_timer.stop();

  const unsigned num_index_buckets =
      static_cast<unsigned>(metadata.buckets.size());
  std::vector<std::vector<ActiveAnchorEntry>> query_buckets(
      num_index_buckets);
  const auto add_query =
      [&](uint64_t key, const ActiveQuerySet::QueryInfo &info) {
        query_buckets[ReadIndexBucketForMetadata(key, metadata)].push_back(
            ActiveAnchorEntry{ReadIndexDirectoryKey(key, metadata),
                              key,
                              info.offset_mask});
      };
  if (local_candidate_mode) {
    local_queries.ForEach(add_query);
  } else {
    active_queries.ForEach(add_query);
  }
  for (auto &bucket : query_buckets) {
    std::sort(bucket.begin(), bucket.end(),
              [](const ActiveAnchorEntry &lhs, const ActiveAnchorEntry &rhs) {
                return lhs.key < rhs.key;
              });
  }

  ReadIndexShardFiles index_files(options.index_prefix, metadata, true);
  const DirectoryFenceIndex directory_fences = EnsureDirectoryFenceIndex(
      options.index_prefix, metadata, index_files);
  MappedWords reads(options.read_file);
  struct BucketLookupResult {
    std::vector<MatchedPostingGroup> groups;
    uint64_t matched_postings{0};
    uint64_t candidate_upper_bound{0};
    uint64_t sparse_directory_blocks{0};
    uint64_t sparse_directory_bytes{0};
    uint64_t dense_directory_bytes{0};
  };
  std::vector<BucketLookupResult> lookup_results(num_index_buckets);
  std::atomic<bool> lookup_ok(true);
  SimpleTimer lookup_timer;
  lookup_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
  for (int bucket_id = 0; bucket_id < static_cast<int>(num_index_buckets);
       ++bucket_id) {
    const unsigned bucket = static_cast<unsigned>(bucket_id);
    BucketLookupResult &result = lookup_results[bucket];
    const ReadIndexBucketMeta &meta = metadata.buckets[bucket];
    const uint64_t directory_file_size =
        index_files.directory_size(meta.file_id);
    const auto &queries = query_buckets[bucket];
    const uint64_t directory_bytes =
        meta.directory_count * metadata.directory_entry_bytes;
    if (meta.directory_offset > directory_file_size ||
        directory_bytes > directory_file_size - meta.directory_offset) {
      lookup_ok.store(false, std::memory_order_relaxed);
      continue;
    }
    if (queries.empty() || meta.directory_count == 0) continue;
    result.groups.reserve(queries.size());

    const auto add_match = [&](const ActiveAnchorEntry &query,
                               const uint8_t *entry) {
      const uint64_t begin = LoadLowBytes(
          entry + metadata.key_bytes, metadata.bucket_begin_bytes);
      const uint64_t count = LoadLowBytes(
          entry + metadata.key_bytes + metadata.bucket_begin_bytes, 4u);
      const uint64_t offset_count = static_cast<unsigned>(
          __builtin_popcountll(query.offset_mask));
      if (count != 0 &&
          offset_count >
              (std::numeric_limits<uint64_t>::max() -
               result.candidate_upper_bound) /
                  count) {
        lookup_ok.store(false, std::memory_order_relaxed);
        return;
      }
      result.groups.push_back(MatchedPostingGroup{
          meta.file_id,
          meta.posting_offset + begin * metadata.locator_bytes, count,
          query.anchor_key,
          query.offset_mask});
      result.matched_postings += count;
      result.candidate_upper_bound += count * offset_count;
    };

    const long double sparse_upper_entries =
        static_cast<long double>(queries.size()) * directory_fences.stride;
    const bool use_sparse =
        sparse_upper_entries * 4.0L <
        static_cast<long double>(meta.directory_count) * 3.0L;
    if (!use_sparse) {
      std::vector<uint8_t> directory(static_cast<size_t>(directory_bytes));
      if (!directory.empty() &&
          !PreadFully(index_files.directory_fd(meta.file_id),
                      directory.data(), directory.size(),
                      meta.directory_offset)) {
        lookup_ok.store(false, std::memory_order_relaxed);
        continue;
      }
      result.dense_directory_bytes += directory_bytes;
      size_t query_index = 0;
      for (uint64_t group = 0;
           group < meta.directory_count && query_index < queries.size();
           ++group) {
        const uint8_t *entry =
            directory.data() + group * metadata.directory_entry_bytes;
        const uint64_t key = LoadLowBytes(entry, metadata.key_bytes);
        while (query_index < queries.size() &&
               queries[query_index].key < key) {
          ++query_index;
        }
        if (query_index == queries.size()) break;
        if (queries[query_index].key != key) continue;
        add_match(queries[query_index], entry);
        ++query_index;
      }
      continue;
    }

    const uint64_t fence_begin = directory_fences.bucket_offsets[bucket];
    const uint64_t fence_count =
        directory_fences.bucket_offsets[bucket + 1u] - fence_begin;
    const auto fence_key = [&](uint64_t fence) {
      return LoadLowBytes(
          directory_fences.keys.data() +
              (fence_begin + fence) * metadata.key_bytes,
          metadata.key_bytes);
    };
    size_t query_begin = 0;
    uint64_t next_fence = std::min<uint64_t>(1u, fence_count);
    while (query_begin < queries.size()) {
      while (next_fence < fence_count &&
             fence_key(next_fence) <= queries[query_begin].key) {
        ++next_fence;
      }
      const uint64_t block = next_fence == 0 ? 0 : next_fence - 1u;
      const uint64_t next_key =
          next_fence < fence_count ? fence_key(next_fence)
                                   : std::numeric_limits<uint64_t>::max();
      size_t query_end = query_begin + 1u;
      while (query_end < queries.size() &&
             (next_fence == fence_count ||
              queries[query_end].key < next_key)) {
        ++query_end;
      }

      const uint64_t group_begin = block * directory_fences.stride;
      const uint64_t group_count = std::min<uint64_t>(
          directory_fences.stride, meta.directory_count - group_begin);
      const uint64_t block_bytes =
          group_count * metadata.directory_entry_bytes;
      const uint64_t block_offset =
          meta.directory_offset +
          group_begin * metadata.directory_entry_bytes;
      const uint8_t *block_data =
          index_files.mapped_directory(meta.file_id) + block_offset;
      size_t query_index = query_begin;
      for (uint64_t group = 0;
           group < group_count && query_index < query_end; ++group) {
        const uint8_t *entry =
            block_data + group * metadata.directory_entry_bytes;
        const uint64_t key = LoadLowBytes(entry, metadata.key_bytes);
        while (query_index < query_end && queries[query_index].key < key) {
          ++query_index;
        }
        if (query_index == query_end) break;
        if (queries[query_index].key != key) continue;
        add_match(queries[query_index], entry);
        ++query_index;
      }
      ++result.sparse_directory_blocks;
      result.sparse_directory_bytes += block_bytes;
      query_begin = query_end;
    }
  }
  lookup_timer.stop();
  if (!lookup_ok.load(std::memory_order_relaxed)) {
    throw std::runtime_error("failed while looking up read-index directory");
  }
  index_files.DiscardMappedDirectories();

  uint64_t matched_keys = 0;
  uint64_t matched_postings = 0;
  uint64_t candidate_upper_bound = 0;
  uint64_t sparse_directory_blocks = 0;
  uint64_t sparse_directory_bytes = 0;
  uint64_t dense_directory_bytes = 0;
  std::vector<std::vector<MatchedPostingGroup>> groups_by_file(
      metadata.num_files);
  for (BucketLookupResult &result : lookup_results) {
    matched_keys += result.groups.size();
    matched_postings += result.matched_postings;
    if (result.candidate_upper_bound >
        std::numeric_limits<uint64_t>::max() - candidate_upper_bound) {
      throw std::length_error("candidate locator count overflow");
    }
    candidate_upper_bound += result.candidate_upper_bound;
    sparse_directory_blocks += result.sparse_directory_blocks;
    sparse_directory_bytes += result.sparse_directory_bytes;
    dense_directory_bytes += result.dense_directory_bytes;
    for (const MatchedPostingGroup &group : result.groups) {
      groups_by_file[group.file_id].push_back(group);
    }
    std::vector<MatchedPostingGroup>().swap(result.groups);
  }
  std::vector<BucketLookupResult>().swap(lookup_results);
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(metadata.num_files);
       ++file_id) {
    std::vector<MatchedPostingGroup> &groups =
        groups_by_file[static_cast<unsigned>(file_id)];
    std::sort(groups.begin(), groups.end(),
              [](const MatchedPostingGroup &lhs,
                 const MatchedPostingGroup &rhs) {
                return lhs.posting_byte_offset < rhs.posting_byte_offset;
              });
  }

  const uint64_t requested_memory =
      options.memory_bytes >=
              static_cast<double>(std::numeric_limits<uint64_t>::max())
          ? std::numeric_limits<uint64_t>::max()
          : static_cast<uint64_t>(options.memory_bytes);
  xinfo("Read-index candidate plan: {} upper-bound locators / {.3} GiB raw, "
        "memory budget {.3} GiB\n",
        candidate_upper_bound,
        static_cast<double>(candidate_upper_bound) * sizeof(uint64_t) /
            static_cast<double>(uint64_t{1} << 30u),
        static_cast<double>(requested_memory) /
            static_cast<double>(uint64_t{1} << 30u));
  const bool streamed_local_query =
      local_candidate_mode &&
      candidate_upper_bound >
          requested_memory / sizeof(uint64_t) * 3u / 4u;
  if (!local_candidate_mode && candidate_upper_bound >
      requested_memory / sizeof(uint64_t) * 3u / 4u) {
    throw std::runtime_error(
        "candidate set exceeds in-memory replay budget");
  }

  const uint64_t candidate_wave_bytes = std::max<uint64_t>(
      uint64_t{64} << 20u,
      requested_memory > (uint64_t{256} << 20u)
          ? requested_memory * 2u / 5u
          : requested_memory / 2u);
  const uint64_t candidate_wave_items =
      std::max<uint64_t>(1u, candidate_wave_bytes / sizeof(uint64_t));
  unsigned num_read_order_buckets;
  if (streamed_local_query) {
    const uint64_t rough_waves =
        DivCeiling(candidate_upper_bound, candidate_wave_items);
    const uint64_t desired = std::max<uint64_t>(
        kMinimumReadOrderBuckets,
        rough_waves * static_cast<uint64_t>(std::max(1, options.num_threads)) *
            8u);
    size_t buckets = 1u;
    while (buckets < desired) {
      if (buckets > std::numeric_limits<unsigned>::max() / 2u) {
        throw std::length_error("streamed read-order bucket count overflow");
      }
      buckets *= 2u;
    }
    num_read_order_buckets = static_cast<unsigned>(buckets);
  } else {
    num_read_order_buckets =
        ReadOrderBucketCount(options, metadata, candidate_upper_bound);
  }
  const uint64_t read_order_bucket_words = std::max<uint64_t>(
      1u, DivCeiling(metadata.parsed_words,
                     static_cast<uint64_t>(num_read_order_buckets)));
  if (streamed_local_query) {
    const uint64_t local_anchor_count = local_queries.anchors.size();
    // Directory matching is complete. Exact window verification needs the
    // canonical seed vector, but no longer needs the much larger
    // anchor/offset expansion or directory query buckets.
    std::vector<CompactLocalEndpointQueries::Entry>().swap(
        local_queries.anchors);
    std::vector<std::vector<ActiveAnchorEntry>>().swap(query_buckets);
    const StreamedLocalCandidateStats streamed =
        QueryLocalCandidatesInReadOrderWaves(
            options, metadata, &index_files, groups_by_file, reads.data(),
            local_queries, candidate_upper_bound, requested_memory,
            num_read_order_buckets);
    std::vector<std::vector<MatchedPostingGroup>>().swap(groups_by_file);
    xinfo("Compact index query: {} / {} active anchor keys matched, {} "
          "posting locators, {} upper-bound / {} unique candidate starts, "
          "{} window-filter hits ({} bases); active {.4}, directory {.4}, "
          "candidate count {.4}, gather {.4}, read-order sort {.4}, window "
          "filter {.4} s; directory touched {} sparse blocks / {.3} GiB "
          "sparse + {.3} GiB dense, {} verify workers\n",
          matched_keys, local_anchor_count, matched_postings,
          candidate_upper_bound, streamed.unique_candidates,
          streamed.window_filter_hits, metadata.window_len,
          active_timer.elapsed(), lookup_timer.elapsed(),
          streamed.count_seconds, streamed.gather_seconds,
          streamed.sort_seconds, streamed.verify_seconds,
          sparse_directory_blocks,
          static_cast<double>(sparse_directory_bytes) /
              (uint64_t{1} << 30u),
          static_cast<double>(dense_directory_bytes) /
              (uint64_t{1} << 30u),
          streamed.verify_threads);
    return;
  }
  CandidateBuckets candidate_buckets(num_read_order_buckets);
  const size_t average_reserve = static_cast<size_t>(DivCeiling(
      candidate_upper_bound, static_cast<uint64_t>(num_read_order_buckets)));
  for (auto &bucket : candidate_buckets) {
    bucket.reserve(average_reserve + average_reserve / 16u + 1u);
  }
  const uint64_t position_mask =
      (uint64_t{1} << metadata.position_bits) - 1u;
  const unsigned locator_value_bits =
      metadata.position_bits +
      BitsNeeded(metadata.parsed_words == 0 ? 0 : metadata.parsed_words - 1u);
  const uint64_t locator_value_mask = LowBitMask(locator_value_bits);
  const uint64_t locator_context_mask =
      LowBitMask(metadata.locator_context_bits);
  const uint32_t *read_words = reads.data();
  SimpleTimer gather_timer;
  gather_timer.start();
  std::vector<std::mutex> candidate_mutexes(num_read_order_buckets);
  std::atomic<bool> gather_ok(true);
  const uint64_t requested_stage_items =
      requested_memory / (uint64_t{64} * sizeof(uint64_t));
  const uint64_t workers =
      static_cast<uint64_t>(std::max(1, options.num_threads));
  const size_t candidate_stage_items = static_cast<size_t>(std::max<uint64_t>(
      1024u, std::min<uint64_t>(DivCeiling(candidate_upper_bound, workers),
                                requested_stage_items / workers)));
  xinfo("Read-order replay plan: {} buckets, at most {} staged candidates "
        "per gather worker\n",
        num_read_order_buckets, candidate_stage_items);
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
  for (int file_id = 0; file_id < static_cast<int>(metadata.num_files);
       ++file_id) {
    const unsigned file = static_cast<unsigned>(file_id);
    std::vector<uint8_t> posting_buffer;
    CandidateBuckets local_candidates(num_read_order_buckets);
    size_t local_staged = 0u;
    const auto flush_bucket = [&](unsigned bucket) {
      std::vector<uint64_t> &source = local_candidates[bucket];
      if (source.empty()) return;
      {
        std::lock_guard<std::mutex> lock(candidate_mutexes[bucket]);
        candidate_buckets[bucket].insert(candidate_buckets[bucket].end(),
                                         source.begin(), source.end());
      }
      source.clear();
    };
    const auto flush_all = [&]() {
      for (unsigned bucket = 0; bucket < num_read_order_buckets; ++bucket) {
        flush_bucket(bucket);
      }
      local_staged = 0u;
    };
    for (const MatchedPostingGroup &group : groups_by_file[file]) {
      if (!gather_ok.load(std::memory_order_relaxed)) break;
      const uint64_t posting_file_size = index_files.posting_size(file);
      const uint64_t group_bytes = group.count * metadata.locator_bytes;
      if (group.posting_byte_offset > posting_file_size ||
          group_bytes > posting_file_size - group.posting_byte_offset) {
        gather_ok.store(false, std::memory_order_relaxed);
        break;
      }
      posting_buffer.resize(static_cast<size_t>(group_bytes));
      if (!posting_buffer.empty() &&
          !PreadFully(index_files.posting_fd(file), posting_buffer.data(),
                      posting_buffer.size(), group.posting_byte_offset)) {
        gather_ok.store(false, std::memory_order_relaxed);
        break;
      }
      const uint8_t *group_postings = posting_buffer.data();
      for (uint64_t posting = 0; posting < group.count; ++posting) {
        const uint64_t encoded_locator = LoadLowBytes(
            group_postings + posting * metadata.locator_bytes,
            metadata.locator_bytes);
        const uint64_t locator = encoded_locator & locator_value_mask;
        const uint8_t locator_context = static_cast<uint8_t>(
            (encoded_locator >> locator_value_bits) & locator_context_mask);
        const uint64_t anchor_pos = locator & position_mask;
        const uint64_t word_offset = locator >> metadata.position_bits;
        uint64_t offsets = group.offset_mask;
        while (offsets != 0) {
          const unsigned offset =
              static_cast<unsigned>(__builtin_ctzll(offsets));
          offsets &= offsets - 1u;
          if (!local_candidate_mode &&
              !active_queries.MayMatchContext(
                  group.anchor_key, offset, locator_context)) {
            continue;
          }
          if (anchor_pos < offset || word_offset >= metadata.parsed_words) {
            continue;
          }
          const uint64_t candidate_start = anchor_pos - offset;
          const uint64_t candidate =
              (word_offset << metadata.position_bits) | candidate_start;
          const unsigned read_bucket = std::min(
              static_cast<unsigned>(word_offset / read_order_bucket_words),
              num_read_order_buckets - 1u);
          local_candidates[read_bucket].push_back(candidate);
          ++local_staged;
          if (local_staged >= candidate_stage_items) flush_all();
        }
      }
    }
    flush_all();
  }
  gather_timer.stop();
  if (!gather_ok.load(std::memory_order_relaxed)) {
    throw std::runtime_error("failed while gathering read-index postings");
  }
  std::vector<std::vector<MatchedPostingGroup>>().swap(groups_by_file);

  SimpleTimer sort_timer;
  sort_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(options.num_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(candidate_buckets.size()); ++bucket) {
    auto &candidates = candidate_buckets[static_cast<size_t>(bucket)];
    kmlib::kmsort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
  }
  sort_timer.stop();

  uint64_t unique_candidates = 0;
  for (const auto &bucket : candidate_buckets) {
    unique_candidates += bucket.size();
  }
  const unsigned verify_working_threads =
      ReadWorkingSetThreads(options, metadata, candidate_buckets);
  uint64_t window_filter_hits = 0;
  SimpleTimer verify_timer;
  verify_timer.start();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : window_filter_hits) \
    num_threads(verify_working_threads)
  for (int64_t bucket = 0;
       bucket < static_cast<int64_t>(candidate_buckets.size()); ++bucket) {
    auto &candidates = candidate_buckets[static_cast<size_t>(bucket)];
    size_t output = 0;
    for (uint64_t candidate : candidates) {
      const uint64_t start = candidate & position_mask;
      const uint64_t word_offset = candidate >> metadata.position_bits;
      if (start + metadata.window_len > read_words[word_offset]) {
        continue;
      }
      AnchorWindow window;
      window.InitFromPtr(read_words + word_offset + 1u,
                         static_cast<unsigned>(start), metadata.window_len);
      const bool exact = local_candidate_mode
                             ? local_queries.Contains(window)
                             : active_windows.Contains(window);
      if (exact) {
        candidates[output++] = candidate;
        ++window_filter_hits;
      }
    }
    DiscardCandidateReadPages(candidates, metadata, read_words);
    candidates.resize(output);
    // Replay only needs exact candidates.  Drop complete pages in the unused
    // tail without allocating a second compact vector: capacity remains a
    // virtual reservation, while resident memory and the replay budget track
    // only size().  No subsequent operation appends to these buckets.
    if (candidates.size() < candidates.capacity()) {
      DiscardMemoryPages(candidates.data() + candidates.size(),
                         (candidates.capacity() - candidates.size()) *
                             sizeof(uint64_t));
    }
  }
  verify_timer.stop();
  const unsigned replay_working_threads =
      ReadWorkingSetThreads(options, metadata, candidate_buckets);
  uint64_t aligned_reads = 0;
  if (local_candidate_mode) {
    PublishLocalCandidates(options, metadata, &candidate_buckets);
  } else if (!options.edge_output_prefix.empty()) {
    ReplayCandidateEdges(options, metadata, candidate_buckets, read_words,
                         window_filter_hits, replay_working_threads,
                         prepared_flank_index.get(),
                         &aligned_reads);
  }
  xinfo("Compact index query: {} / {} active anchor keys matched, {} posting "
        "locators, {} upper-bound / {} unique candidate starts, {} "
        "window-filter hits ({} bases); active {.4}, directory {.4}, "
        "gather {.4}, "
        "read-order sort {.4}, window filter {.4} s; directory touched {} "
        "sparse blocks / {.3} GiB sparse + {.3} GiB dense, {} verify / {} "
        "replay workers\n",
        matched_keys,
        local_candidate_mode ? local_queries.anchors.size()
                             : active_queries.key_count(),
        matched_postings,
        candidate_upper_bound, unique_candidates, window_filter_hits,
        metadata.window_len, active_timer.elapsed(), lookup_timer.elapsed(),
        gather_timer.elapsed(), sort_timer.elapsed(), verify_timer.elapsed(),
        sparse_directory_blocks,
        static_cast<double>(sparse_directory_bytes) / (uint64_t{1} << 30u),
        static_cast<double>(dense_directory_bytes) / (uint64_t{1} << 30u),
        verify_working_threads, replay_working_threads);
}

uint64_t FileSize(int fd) {
  struct stat status;
  if (fstat(fd, &status) != 0 || status.st_size < 0) {
    throw std::runtime_error("cannot stat read-index file");
  }
  return static_cast<uint64_t>(status.st_size);
}

void VerifyReadIndex(const ReadIndexOptions &options) {
  ReadIndexMetadata metadata = LoadReadIndexMetadata(options.index_prefix);
  MappedWords reads(options.read_file);
  if (metadata.parsed_words > reads.word_size()) {
    throw std::runtime_error("read index refers past packed-read input");
  }
  ReadIndexShardFiles index_files(options.index_prefix, metadata, false);
  const uint64_t position_mask =
      metadata.position_bits == 64
          ? std::numeric_limits<uint64_t>::max()
          : (uint64_t{1} << metadata.position_bits) - 1u;
  const unsigned locator_value_bits =
      metadata.position_bits +
      BitsNeeded(metadata.parsed_words == 0 ? 0 : metadata.parsed_words - 1u);
  const uint64_t locator_value_mask = LowBitMask(locator_value_bits);
  const uint64_t locator_context_mask =
      LowBitMask(metadata.locator_context_bits);
  const uint32_t *read_words = reads.data();
  uint64_t observed_keys = 0;
  uint64_t observed_postings = 0;
  SimpleTimer timer;
  timer.start();
  for (size_t bucket = 0; bucket < metadata.buckets.size(); ++bucket) {
    const ReadIndexBucketMeta &meta = metadata.buckets[bucket];
    const uint64_t directory_file_size =
        index_files.directory_size(meta.file_id);
    const uint64_t posting_file_size =
        index_files.posting_size(meta.file_id);
    const uint64_t directory_bytes =
        meta.directory_count * metadata.directory_entry_bytes;
    const uint64_t posting_bytes =
        meta.posting_count * metadata.locator_bytes;
    if (meta.directory_offset > directory_file_size ||
        directory_bytes > directory_file_size - meta.directory_offset ||
        meta.posting_offset > posting_file_size ||
        posting_bytes > posting_file_size - meta.posting_offset) {
      throw std::runtime_error("read-index bucket range exceeds file");
    }
    std::vector<uint8_t> directory(static_cast<size_t>(directory_bytes));
    std::vector<uint8_t> postings(static_cast<size_t>(posting_bytes));
    if ((!directory.empty() &&
         !PreadFully(index_files.directory_fd(meta.file_id), directory.data(),
                     directory.size(), meta.directory_offset)) ||
        (!postings.empty() &&
         !PreadFully(index_files.posting_fd(meta.file_id), postings.data(),
                     postings.size(), meta.posting_offset))) {
      throw std::runtime_error("cannot read compact read-index bucket");
    }
    uint64_t expected_begin = 0;
    uint64_t previous_key = 0;
    bool have_previous_key = false;
    for (uint64_t group = 0; group < meta.directory_count; ++group) {
      const uint8_t *entry =
          directory.data() + group * metadata.directory_entry_bytes;
      const uint64_t key = LoadLowBytes(entry, metadata.key_bytes);
      const uint64_t begin = LoadLowBytes(
          entry + metadata.key_bytes, metadata.bucket_begin_bytes);
      const uint64_t count = LoadLowBytes(
          entry + metadata.key_bytes + metadata.bucket_begin_bytes, 4u);
      const bool invalid_key =
          metadata.version == 1
              ? ReadIndexBucket(key, metadata.bucket_bits) != bucket
              : key > LowBitMask(metadata.anchor_len * 2u -
                                 metadata.bucket_bits);
      if (invalid_key || count == 0 ||
          begin != expected_begin ||
          (have_previous_key && key <= previous_key) ||
          begin > meta.posting_count ||
          count > meta.posting_count - begin) {
        throw std::runtime_error("invalid compact read-index directory");
      }
      uint64_t previous_locator = 0;
      for (uint64_t i = 0; i < count; ++i) {
        const uint64_t encoded_locator = LoadLowBytes(
            postings.data() + (begin + i) * metadata.locator_bytes,
            metadata.locator_bytes);
        const uint64_t locator = encoded_locator & locator_value_mask;
        if (metadata.postings_sorted && i != 0 &&
            locator <= previous_locator) {
          throw std::runtime_error("read-index postings are not unique/sorted");
        }
        const uint64_t anchor_pos = locator & position_mask;
        const uint64_t word_offset = locator >> metadata.position_bits;
        if (word_offset >= metadata.parsed_words) {
          throw std::runtime_error("read-index locator exceeds input");
        }
        const uint32_t read_len = read_words[word_offset];
        if (anchor_pos + metadata.anchor_len > read_len) {
          throw std::runtime_error("read-index locator has invalid position");
        }
        const uint32_t *sequence = read_words + word_offset + 1u;
        uint8_t expected_context = 0u;
        if (metadata.locator_context_bits >= 2u &&
            anchor_pos >= metadata.locator_context_distance) {
          expected_context |= PackedBase(
              sequence, static_cast<unsigned>(anchor_pos) -
                            metadata.locator_context_distance);
        }
        if (metadata.locator_context_bits >= 4u &&
            anchor_pos + metadata.anchor_len +
                    metadata.locator_context_distance <=
                read_len) {
          expected_context |= static_cast<uint8_t>(
              PackedBase(sequence,
                         static_cast<unsigned>(anchor_pos) +
                             metadata.anchor_len +
                             metadata.locator_context_distance - 1u)
              << 2u);
        }
        const uint8_t observed_context = static_cast<uint8_t>(
            (encoded_locator >> locator_value_bits) & locator_context_mask);
        if (observed_context != expected_context) {
          throw std::runtime_error("read-index posting/context mismatch");
        }
        uint64_t observed_key = 0;
        for (unsigned base = 0; base < metadata.anchor_len; ++base) {
          observed_key =
              (observed_key << 2u) |
              PackedBase(sequence, static_cast<unsigned>(anchor_pos) + base);
        }
        if (ReadIndexBucketForMetadata(observed_key, metadata) != bucket ||
            ReadIndexDirectoryKey(observed_key, metadata) != key) {
          throw std::runtime_error("read-index posting/key mismatch");
        }
        previous_locator = locator;
      }
      expected_begin += count;
      previous_key = key;
      have_previous_key = true;
    }
    if (expected_begin != meta.posting_count) {
      throw std::runtime_error("read-index bucket lost postings");
    }
    observed_keys += meta.directory_count;
    observed_postings += meta.posting_count;
  }
  timer.stop();
  if (observed_keys != metadata.num_keys ||
      observed_postings != metadata.num_occurrences) {
    throw std::runtime_error("read-index global counts do not match metadata");
  }
  xinfo("Verified exact read index: {} keys / {} postings against packed "
        "reads in {.4} s\n",
        observed_keys, observed_postings, timer.elapsed());
}

unsigned ChooseReadIndexBucketBits(size_t chunk_count, int num_threads,
                                   uint64_t requested_memory) {
  // Physical shards are scheduling units, not an algorithmic k threshold.
  // Keep enough of them to expose parallelism after the input scan, but never
  // create more fine-grained work than the producer chunking can justify.
  // Per-worker HLL and vector metadata is kept within a fixed fraction of the
  // caller's memory budget.  All quantities scale with the requested input
  // and resources; no socket count or particular CPU layout is assumed.
  constexpr uint64_t kTasksPerWorker = 32u;
  constexpr uint64_t kMetadataBudgetFraction = 16u;
  constexpr uint64_t kPerWorkerBucketMetadata =
      kBuildHllRegisters + 2u * sizeof(std::vector<uint8_t>);
  const uint64_t workers = static_cast<uint64_t>(std::max(1, num_threads));
  const uint64_t minimum_buckets = uint64_t{1}
                                   << kMinimumReadIndexBucketBits;
  const uint64_t worker_target =
      workers > std::numeric_limits<uint64_t>::max() / kTasksPerWorker
          ? std::numeric_limits<uint64_t>::max()
          : workers * kTasksPerWorker;
  const uint64_t input_target = std::max<uint64_t>(
      minimum_buckets, static_cast<uint64_t>(chunk_count));
  const uint64_t desired_buckets =
      std::max(minimum_buckets, std::min(worker_target, input_target));
  const uint64_t metadata_bytes_per_bucket =
      workers > std::numeric_limits<uint64_t>::max() /
                    kPerWorkerBucketMetadata
          ? std::numeric_limits<uint64_t>::max()
          : workers * kPerWorkerBucketMetadata;
  const uint64_t metadata_budget =
      std::max<uint64_t>(uint64_t{1} << kMinimumReadIndexBucketBits,
                         requested_memory / kMetadataBudgetFraction);
  const uint64_t maximum_buckets = std::max<uint64_t>(
      minimum_buckets, metadata_budget / metadata_bytes_per_bucket);
  const uint64_t target_buckets =
      std::min(desired_buckets, maximum_buckets);
  unsigned bits = kMinimumReadIndexBucketBits;
  while (bits < kMaximumReadIndexBucketBits &&
         (uint64_t{1} << (bits + 1u)) <= target_buckets) {
    ++bits;
  }
  return bits;
}

void RunIndexBuild(const ReadIndexOptions &options) {
  constexpr size_t kBuildChunkBytes = size_t{1} << 20u;
  MappedWords mapping(options.read_file);
  uint64_t num_reads = 0;
  uint64_t num_bases = 0;
  unsigned max_read_len = 0;
  size_t parsed_words = 0;
  SimpleTimer boundary_timer;
  const double boundary_cpu_begin = ProcessCpuSeconds();
  boundary_timer.start();
  std::vector<BinaryChunk> chunks = ParseChunks(
      mapping,
      options.max_reads > 0 ? static_cast<uint64_t>(options.max_reads) : 0,
      options.num_threads, &num_reads, &num_bases, &max_read_len,
      &parsed_words, kBuildChunkBytes);
  boundary_timer.stop();
  const double boundary_cpu_seconds =
      ProcessCpuSeconds() - boundary_cpu_begin;

  std::unique_ptr<MappedAnchorPositions> stored_anchor_positions;
  if (options.max_reads == 0 &&
      std::getenv("MEGAHIT_DISABLE_READ_ANCHOR_POSITIONS") == nullptr) {
    stored_anchor_positions.reset(new MappedAnchorPositions(
        options.read_file, options.anchor_len, options.window_len,
        mapping.byte_size(), num_reads, num_bases));
    if (!stored_anchor_positions->valid()) {
      stored_anchor_positions.reset();
    } else {
      xinfo("Read-index build reusing buildlib anchor positions: {.3} GiB "
            "compressed side stream\n",
            static_cast<double>(stored_anchor_positions->payload_bytes()) /
                (uint64_t{1} << 30u));
    }
  }

  const unsigned position_bits =
      BitsNeeded(max_read_len == 0 ? 0 : max_read_len - 1u);
  const unsigned word_offset_bits =
      BitsNeeded(parsed_words == 0 ? 0 : parsed_words - 1u);
  const unsigned locator_bits = position_bits + word_offset_bits;
  if (locator_bits > 64u) {
    throw std::length_error("read-index locator exceeds 64 bits");
  }
  const unsigned locator_bytes = DivCeiling(locator_bits, 8u);
  const unsigned original_key_bits = options.anchor_len * 2u;
  const int num_threads = std::max(1, options.num_threads);
  const uint64_t requested_memory =
      options.memory_bytes >=
              static_cast<double>(std::numeric_limits<uint64_t>::max())
          ? std::numeric_limits<uint64_t>::max()
          : std::max<uint64_t>(1u,
                               static_cast<uint64_t>(options.memory_bytes));
  const unsigned bucket_bits = ChooseReadIndexBucketBits(
      chunks.size(), num_threads, requested_memory);
  const unsigned num_index_buckets = 1u << bucket_bits;
  const bool use_permuted_suffix = original_key_bits > bucket_bits;
  const unsigned directory_key_bits =
      use_permuted_suffix ? original_key_bits - bucket_bits
                          : original_key_bits;
  const unsigned key_bytes = DivCeiling(directory_key_bits, 8u);
  const unsigned locator_spare_bits = locator_bytes * 8u - locator_bits;
  const unsigned unaligned_key_bits = directory_key_bits % 8u;
  const unsigned temporary_key_extra_bits =
      use_permuted_suffix && unaligned_key_bits != 0 &&
              locator_spare_bits >= unaligned_key_bits
          ? unaligned_key_bits
          : 0u;
  // The locator is byte-packed, so the high padding bits otherwise go to
  // disk unused. Preserve the bases immediately adjacent to the selected
  // anchor there: an active exact window can reject an incompatible posting
  // before materializing a read candidate, without increasing the index or
  // changing exact verification. The number of usable bits follows solely
  // from the portable on-disk representation, not a dataset/machine cutoff.
  const unsigned available_context_bits =
      locator_spare_bits - temporary_key_extra_bits;
  const unsigned context_span =
      options.window_len > options.anchor_len
          ? options.window_len - options.anchor_len
          : 0u;
  const unsigned locator_context_bits =
      context_span == 0u
          ? 0u
          : (available_context_bits >= 4u
                 ? 4u
                 : (available_context_bits >= 2u ? 2u : 0u));
  // Immediate flanking bases are defined for every occurrence independently
  // of the query window and maximize the number of offsets that can consume
  // both stored bases.
  const unsigned locator_context_distance = 1u;
  const uint64_t locator_context_mask =
      LowBitMask(locator_context_bits);

  // When requested, prepare the first iteration's immutable query state
  // before touching reads.  The index occurrence stream and iterative edge
  // stream can then consume the same decoded read while it is hot, instead
  // of independently traversing a tens-of-gigabytes library.
  std::unique_ptr<ActiveQuerySet> fused_active_queries;
  std::unique_ptr<ActiveWindowSet> fused_active_windows;
  std::unique_ptr<PreparedFlankIndexBase> fused_prepared_flank_index;
  std::unique_ptr<FusedEdgeReplayBase> fused_edge_replay;
  SimpleTimer fused_prepare_timer;
  if (!options.edge_output_prefix.empty()) {
    fused_prepare_timer.start();
    fused_active_queries.reset(new ActiveQuerySet());
    fused_active_queries->ConfigureContext(
        locator_context_bits, options.anchor_len, options.window_len,
        locator_context_distance);
    fused_active_windows.reset(new ActiveWindowSet());
    BuildActiveQueries(options, fused_active_queries.get(),
                       fused_active_windows.get(),
                       &fused_prepared_flank_index);
    fused_active_queries->PrepareRandomLookup();
    fused_edge_replay = MakeFusedEdgeReplay(
        options, fused_prepared_flank_index.get());
    fused_prepare_timer.stop();
    xinfo("Fused first-iteration preparation: {.4} s; exact-window filter "
          "{} bytes\n",
          fused_prepare_timer.elapsed(),
          fused_active_windows->byte_size());
  }

  const unsigned temporary_key_bits =
      directory_key_bits - temporary_key_extra_bits;
  const unsigned temporary_key_bytes = DivCeiling(temporary_key_bits, 8u);
  const unsigned temporary_record_bytes =
      temporary_key_bytes + locator_bytes;
  const uint64_t locator_mask = LowBitMask(locator_bits);
  const uint32_t *words = mapping.data();
  const int num_writer_threads =
      num_threads > 1
          ? std::min(num_threads / 2,
                     std::max(1, static_cast<int>(
                                     std::ceil(std::sqrt(num_threads)))))
          : 0;
  const int num_compute_threads =
      std::max(1, num_threads - num_writer_threads);
  const unsigned num_temporary_files =
      static_cast<unsigned>(std::max(1, num_writer_threads));
  const std::string temporary_path = options.output_prefix + ".ridx.tmp";
  std::vector<std::string> temporary_paths(num_temporary_files);
  std::vector<int> temporary_fds(num_temporary_files, -1);
  for (unsigned file_id = 0; file_id < num_temporary_files; ++file_id) {
    temporary_paths[file_id] =
        file_id == 0 ? temporary_path
                     : temporary_path + "." + std::to_string(file_id);
    temporary_fds[file_id] =
        open(temporary_paths[file_id].c_str(),
             O_CREAT | O_TRUNC | O_RDWR | O_CLOEXEC, 0664);
    if (temporary_fds[file_id] < 0) {
      for (unsigned opened = 0; opened < file_id; ++opened) {
        close(temporary_fds[opened]);
        std::remove(temporary_paths[opened].c_str());
      }
      throw std::runtime_error("cannot create sharded read-index temporary "
                               "files");
    }
  }
  const auto close_temporary_files = [&]() {
    for (int &fd : temporary_fds) {
      if (fd >= 0) {
        close(fd);
        fd = -1;
      }
    }
  };
  const auto remove_temporary_files = [&]() {
    for (const std::string &path : temporary_paths) {
      std::remove(path.c_str());
    }
  };

  struct BucketExtent {
    uint32_t file_id;
    uint64_t byte_offset;
    uint64_t records;
  };
  struct BucketWriteBuffer {
    std::vector<BucketExtent> extents;
    uint64_t records_written{0};
  };
  std::vector<uint64_t> temporary_write_cursors(num_temporary_files, 0);
  using BucketExtentLists = std::vector<std::vector<BucketExtent>>;
  std::vector<BucketExtentLists> thread_bucket_extents(
      static_cast<size_t>(num_compute_threads),
      BucketExtentLists(num_index_buckets));
  std::vector<BucketExtentLists> writer_bucket_extents(
      static_cast<size_t>(num_writer_threads),
      BucketExtentLists(num_index_buckets));

  std::vector<ReadIndexBuildScratch> scratch(
      static_cast<size_t>(num_compute_threads));
  for (ReadIndexBuildScratch &local : scratch) {
    local.bucket_hll.assign(
        static_cast<size_t>(num_index_buckets) * kBuildHllRegisters, 0);
    local.bucket_records.resize(num_index_buckets);
  }

  // One read pass performs anchor discovery, HLL cardinality estimation,
  // physical partitioning, and sequential external output.  The former
  // builder scanned all reads once for counts and again for fill.  Fine input
  // chunks keep producer memory bounded.  Every producer owns its bucket
  // buffers and emits all non-empty buffers in one vectored write.  Only the
  // append offset is shared; bucket extent metadata remains thread-local until
  // the scan has completed.
  std::atomic<bool> fill_ok(true);
  // Producer bucket buffers are a random-scatter working set.  Giving them a
  // large fraction of total memory increases resident pages and TLB/page-fault
  // cost without reducing output traffic: writer-side coalescing already
  // controls the physical I/O size.  Keep this first-level cache deliberately
  // small and scale it only with the caller's memory budget; the bounded batch
  // pool provides backpressure on machines with either few or many workers.
  const uint64_t total_local_buffer_budget =
      std::max<uint64_t>(temporary_record_bytes, requested_memory / 128u);
  const size_t local_flush_bytes = static_cast<size_t>(std::max<uint64_t>(
      temporary_record_bytes,
      total_local_buffer_budget /
          static_cast<uint64_t>(num_compute_threads) /
          temporary_record_bytes * temporary_record_bytes));
  const uint64_t local_flush_records =
      std::max<uint64_t>(1, local_flush_bytes / temporary_record_bytes);
  const uint64_t mean_bucket_records = DivCeiling(
      local_flush_records, static_cast<uint64_t>(num_index_buckets));
  // AnchorBucketHash makes bucket ownership uniform.  Reserve once for the
  // expected occupancy plus a distribution-independent tail margin.  This
  // replaces repeated geometric vector growth (and its copied bytes/page
  // faults) while retaining a normal vector fallback for extreme skew.
  const uint64_t bucket_record_margin = static_cast<uint64_t>(std::ceil(
      6.0 * std::sqrt(static_cast<double>(mean_bucket_records + 1u)))) +
      32u;
  const uint64_t reserved_bucket_records =
      mean_bucket_records + bucket_record_margin;
  const size_t bucket_reserve_bytes = static_cast<size_t>(std::min<uint64_t>(
      std::numeric_limits<size_t>::max(),
      reserved_bucket_records > std::numeric_limits<uint64_t>::max() /
                                    temporary_record_bytes
          ? std::numeric_limits<uint64_t>::max()
          : reserved_bucket_records * temporary_record_bytes));
  const uint64_t total_writer_coalesce_budget = std::max<uint64_t>(
      temporary_record_bytes, requested_memory / 8u);
  const size_t writer_flush_bytes = static_cast<size_t>(std::max<uint64_t>(
      local_flush_bytes,
      num_writer_threads == 0
          ? temporary_record_bytes
          : total_writer_coalesce_budget /
                static_cast<uint64_t>(num_writer_threads) /
                temporary_record_bytes * temporary_record_bytes));
  const uint64_t writer_flush_records =
      std::max<uint64_t>(1, writer_flush_bytes / temporary_record_bytes);
  const uint64_t writer_mean_bucket_records = DivCeiling(
      writer_flush_records, static_cast<uint64_t>(num_index_buckets));
  const uint64_t writer_bucket_margin = static_cast<uint64_t>(std::ceil(
      6.0 * std::sqrt(static_cast<double>(writer_mean_bucket_records + 1u)))) +
      32u;
  const size_t writer_bucket_reserve_bytes = static_cast<size_t>(
      (writer_mean_bucket_records + writer_bucket_margin) *
      temporary_record_bytes);
  struct WriteBatch {
    explicit WriteBatch(unsigned buckets) : bucket_records(buckets) {}
    std::vector<std::vector<uint8_t>> bucket_records;
  };
  std::mutex write_queue_mutex;
  std::condition_variable write_queue_ready;
  std::condition_variable reusable_batch_ready;
  std::condition_variable pending_write_space;
  std::deque<std::unique_ptr<WriteBatch>> pending_writes;
  std::deque<std::unique_ptr<WriteBatch>> reusable_writes;
  const size_t maximum_pending_writes = static_cast<size_t>(
      std::max(1, num_writer_threads * 2));
  const size_t maximum_write_batches =
      maximum_pending_writes + static_cast<size_t>(num_writer_threads);
  size_t allocated_write_batches = 0;
  bool producers_done = false;
  std::atomic<uint64_t> producer_cpu_nanos(0);
  std::atomic<uint64_t> producer_wait_nanos(0);
  std::atomic<uint64_t> writer_cpu_nanos(0);
  std::atomic<uint64_t> writer_copy_nanos(0);
  std::atomic<uint64_t> writer_io_nanos(0);
  std::atomic<uint64_t> writer_batches(0);
  std::atomic<uint64_t> writer_flushes(0);
  std::atomic<uint64_t> writer_bytes(0);
  std::atomic<uint64_t> allocation_nanos(0);
  std::atomic<uint64_t> allocation_calls(0);
  std::vector<uint64_t> allocated_file_bytes(num_temporary_files, 0);
  std::vector<uint8_t> incremental_allocation_supported(
      num_temporary_files, 1);
  const uint64_t allocation_granularity = std::max<uint64_t>(
      uint64_t{16} << 20u,
      std::min<uint64_t>(uint64_t{256} << 20u,
                         total_writer_coalesce_budget / 4u));

  const auto ensure_temporary_space = [&](unsigned file_id,
                                          uint64_t required_bytes) {
    // Each helper writer owns one temporary shard, so allocation is strictly
    // sequential within a file and independent across files.  This avoids a
    // global inode/allocation lock while retaining a normal single-file path
    // when only one worker is requested.
    if (!incremental_allocation_supported[file_id] ||
        required_bytes <= allocated_file_bytes[file_id]) {
      return true;
    }
    const uint64_t target =
        DivCeiling(required_bytes, allocation_granularity) *
        allocation_granularity;
    const uint64_t begin = ClockNanoseconds(CLOCK_MONOTONIC);
    const int error = posix_fallocate(
        temporary_fds[file_id],
        static_cast<off_t>(allocated_file_bytes[file_id]),
        static_cast<off_t>(target - allocated_file_bytes[file_id]));
    allocation_nanos.fetch_add(ClockNanoseconds(CLOCK_MONOTONIC) - begin,
                               std::memory_order_relaxed);
    allocation_calls.fetch_add(1, std::memory_order_relaxed);
    if (error == 0) {
      allocated_file_bytes[file_id] = target;
      return true;
    }
    if (error == EOPNOTSUPP || error == ENOSYS || error == EINVAL) {
      incremental_allocation_supported[file_id] = 0;
      return true;
    }
    errno = error;
    return false;
  };

  const auto write_records = [&](std::vector<std::vector<uint8_t>> *records,
                                 unsigned file_id, uint64_t file_offset) {
    std::vector<struct iovec> vectors;
    vectors.reserve(num_index_buckets);
    for (std::vector<uint8_t> &bucket_records : *records) {
      if (bucket_records.empty()) continue;
      struct iovec vector;
      vector.iov_base = bucket_records.data();
      vector.iov_len = bucket_records.size();
      vectors.push_back(vector);
    }
    return PwritevFully(temporary_fds[file_id], &vectors, file_offset);
  };

  const auto flush_records = [&](std::vector<std::vector<uint8_t>> *records,
                                 BucketExtentLists *extents,
                                 unsigned file_id) {
    size_t total_bytes = 0;
    for (const std::vector<uint8_t> &bucket_records : *records) {
      if (bucket_records.size() >
          std::numeric_limits<size_t>::max() - total_bytes) {
        return false;
      }
      total_bytes += bucket_records.size();
    }
    if (total_bytes == 0) return true;
    const uint64_t file_offset = temporary_write_cursors[file_id];
    if (total_bytes > std::numeric_limits<uint64_t>::max() - file_offset ||
        !ensure_temporary_space(file_id, file_offset + total_bytes)) {
      return false;
    }
    temporary_write_cursors[file_id] += total_bytes;
    size_t bucket_offset = 0;
    for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
      const std::vector<uint8_t> &bucket_records = (*records)[bucket];
      if (bucket_records.empty()) continue;
      (*extents)[bucket].push_back(BucketExtent{
          file_id, file_offset + bucket_offset,
          static_cast<uint64_t>(bucket_records.size() /
                                temporary_record_bytes)});
      bucket_offset += bucket_records.size();
    }
    const uint64_t io_begin = ClockNanoseconds(CLOCK_MONOTONIC);
    const bool ok = write_records(records, file_id, file_offset);
    writer_io_nanos.fetch_add(
        ClockNanoseconds(CLOCK_MONOTONIC) - io_begin,
        std::memory_order_relaxed);
    writer_flushes.fetch_add(1, std::memory_order_relaxed);
    writer_bytes.fetch_add(total_bytes, std::memory_order_relaxed);
    return ok;
  };

  const auto enqueue_local = [&](ReadIndexBuildScratch *local,
                                 BucketExtentLists *extents) {
    size_t total_bytes = 0;
    for (const std::vector<uint8_t> &records : local->bucket_records) {
      if (records.empty()) continue;
      if (records.size() >
          std::numeric_limits<size_t>::max() - total_bytes) {
        return false;
      }
      total_bytes += records.size();
    }
    if (total_bytes == 0) return true;
    if (num_writer_threads == 0) {
      const bool ok = flush_records(&local->bucket_records, extents, 0);
      for (std::vector<uint8_t> &records : local->bucket_records) {
        records.clear();
      }
      return ok;
    }

    std::unique_ptr<WriteBatch> batch;
    {
      std::unique_lock<std::mutex> lock(write_queue_mutex);
      const uint64_t wait_begin = ClockNanoseconds(CLOCK_MONOTONIC);
      reusable_batch_ready.wait(lock, [&] {
        return !reusable_writes.empty() ||
               allocated_write_batches < maximum_write_batches ||
               !fill_ok.load(std::memory_order_relaxed);
      });
      producer_wait_nanos.fetch_add(
          ClockNanoseconds(CLOCK_MONOTONIC) - wait_begin,
          std::memory_order_relaxed);
      if (!fill_ok.load(std::memory_order_relaxed)) return false;
      if (!reusable_writes.empty()) {
        batch = std::move(reusable_writes.front());
        reusable_writes.pop_front();
      } else {
        ++allocated_write_batches;
      }
    }
    if (!batch) batch.reset(new WriteBatch(num_index_buckets));
    batch->bucket_records.swap(local->bucket_records);
    {
      std::unique_lock<std::mutex> lock(write_queue_mutex);
      const uint64_t wait_begin = ClockNanoseconds(CLOCK_MONOTONIC);
      pending_write_space.wait(lock, [&] {
        return pending_writes.size() < maximum_pending_writes ||
               !fill_ok.load(std::memory_order_relaxed);
      });
      producer_wait_nanos.fetch_add(
          ClockNanoseconds(CLOCK_MONOTONIC) - wait_begin,
          std::memory_order_relaxed);
      if (!fill_ok.load(std::memory_order_relaxed)) return false;
      pending_writes.push_back(std::move(batch));
    }
    write_queue_ready.notify_one();
    return true;
  };
  std::vector<std::thread> write_threads;
  write_threads.reserve(static_cast<size_t>(num_writer_threads));
  for (int writer_id = 0; writer_id < num_writer_threads; ++writer_id) {
    write_threads.emplace_back([&, writer_id] {
      // If the parent is an OpenMP-bound master, a std::thread otherwise
      // inherits that single place.  Widen helper writers back to the
      // process's startup cpuset; this remains portable under taskset/cgroups.
      ResetThreadAffinityToStartupMask();
      const uint64_t writer_cpu_begin =
          ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID);
      std::vector<std::vector<uint8_t>> accumulated(num_index_buckets);
      size_t accumulated_bytes = 0;
      while (true) {
        std::unique_ptr<WriteBatch> batch;
        {
          std::unique_lock<std::mutex> lock(write_queue_mutex);
          write_queue_ready.wait(lock, [&] {
            return !pending_writes.empty() || producers_done;
          });
          if (pending_writes.empty()) {
            if (producers_done) break;
            continue;
          }
          batch = std::move(pending_writes.front());
          pending_writes.pop_front();
        }
        pending_write_space.notify_one();
        if (fill_ok.load(std::memory_order_relaxed)) {
          const uint64_t copy_begin = ClockNanoseconds(CLOCK_MONOTONIC);
          for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
            std::vector<uint8_t> &source = batch->bucket_records[bucket];
            if (source.empty()) continue;
            std::vector<uint8_t> &destination = accumulated[bucket];
            if (destination.capacity() == 0) {
              destination.reserve(writer_bucket_reserve_bytes);
            }
            destination.insert(destination.end(), source.begin(),
                               source.end());
            accumulated_bytes += source.size();
          }
          writer_copy_nanos.fetch_add(
              ClockNanoseconds(CLOCK_MONOTONIC) - copy_begin,
              std::memory_order_relaxed);
          writer_batches.fetch_add(1, std::memory_order_relaxed);
        }
        for (std::vector<uint8_t> &records : batch->bucket_records) {
          records.clear();
        }
        {
          std::lock_guard<std::mutex> lock(write_queue_mutex);
          reusable_writes.push_back(std::move(batch));
        }
        reusable_batch_ready.notify_one();
        if (fill_ok.load(std::memory_order_relaxed) &&
            accumulated_bytes >= writer_flush_bytes) {
          if (!flush_records(
                  &accumulated,
                  &writer_bucket_extents[static_cast<size_t>(writer_id)],
                  static_cast<unsigned>(writer_id))) {
            fill_ok.store(false, std::memory_order_relaxed);
            reusable_batch_ready.notify_all();
            pending_write_space.notify_all();
          }
          for (std::vector<uint8_t> &records : accumulated) records.clear();
          accumulated_bytes = 0;
        }
      }
      if (fill_ok.load(std::memory_order_relaxed) &&
          !flush_records(
              &accumulated,
              &writer_bucket_extents[static_cast<size_t>(writer_id)],
              static_cast<unsigned>(writer_id))) {
        fill_ok.store(false, std::memory_order_relaxed);
        reusable_batch_ready.notify_all();
        pending_write_space.notify_all();
      }
      for (std::vector<uint8_t> &records : accumulated) {
        DiscardMemoryPages(records.data(), records.capacity());
      }
      writer_cpu_nanos.fetch_add(
          ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID) - writer_cpu_begin,
          std::memory_order_relaxed);
    });
  }
  SimpleTimer fill_timer;
  const double fill_cpu_begin = ProcessCpuSeconds();
  fill_timer.start();
  size_t replay_wave_chunks = chunks.size();
  uint64_t replay_target_batch_bytes = 0;
  size_t replay_commit_batches = 0;
  double replay_commit_seconds = 0.0;
  if (fused_edge_replay) {
    const uint64_t candidate_bytes =
        fused_edge_replay->candidate_record_bytes();
    // A batch and its radix-partitioned copy coexist briefly.  Keep the
    // candidate stream to one sixteenth of the explicit memory budget; the
    // persistent edge sets, flank queries and index writer consume the rest.
    // The first wave uses a one-edge-per-base upper bound, then adapts from
    // observed candidate density.  This is input/resource driven rather than
    // tied to a core count or a particular host.
    replay_target_batch_bytes =
        std::max<uint64_t>(candidate_bytes, requested_memory / 16u);
    const uint64_t worst_chunk_candidate_bytes = std::max<uint64_t>(
        candidate_bytes,
        static_cast<uint64_t>(kBuildChunkBytes) * 4u * candidate_bytes);
    replay_wave_chunks = static_cast<size_t>(std::max<uint64_t>(
        1u, replay_target_batch_bytes / worst_chunk_candidate_bytes));
    xinfo("Fused edge collection: {}-byte records, {.3} GiB bounded "
          "candidate batches, initial {} read chunks\n",
          candidate_bytes,
          static_cast<double>(replay_target_batch_bytes) /
              (uint64_t{1} << 30u),
          replay_wave_chunks);
  }

  for (size_t wave_begin = 0; wave_begin < chunks.size();) {
    const size_t wave_end = std::min(
        chunks.size(), wave_begin + replay_wave_chunks);
#pragma omp parallel num_threads(num_compute_threads)
    {
      const uint64_t producer_cpu_begin =
          ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID);
      const size_t tid = static_cast<size_t>(omp_get_thread_num());
      ReadIndexBuildScratch &local = scratch[tid];
      BucketExtentLists &local_extents = thread_bucket_extents[tid];
      size_t &buffered_bytes = local.buffered_record_bytes;
      unsigned &reads_until_flush_check = local.reads_until_flush_check;
#pragma omp for schedule(dynamic, 1)
      for (size_t chunk_id = wave_begin; chunk_id < wave_end; ++chunk_id) {
        const BinaryChunk &chunk = chunks[chunk_id];
        const uint8_t *anchor_cursor = nullptr;
        const uint8_t *anchor_end = nullptr;
        if (stored_anchor_positions &&
            !stored_anchor_positions->Range(chunk, &anchor_cursor,
                                            &anchor_end)) {
          fill_ok.store(false, std::memory_order_relaxed);
        }
        size_t word_pos = chunk.word_begin;
        while (word_pos < chunk.word_end &&
               fill_ok.load(std::memory_order_relaxed)) {
          const size_t read_word_offset = word_pos;
          const uint32_t length = words[word_pos++];
          const uint32_t *sequence = words + word_pos;
          std::vector<uint32_t> &replay_positions =
              local.replay_candidate_positions;
          replay_positions.clear();
          const auto emit_anchor = [&](uint64_t key, uint32_t anchor_pos) {
            const uint64_t permuted_key =
                use_permuted_suffix
                    ? PermuteReadIndexKey(key, original_key_bits)
                    : key;
            const unsigned bucket =
                use_permuted_suffix
                    ? static_cast<unsigned>(permuted_key >>
                                            directory_key_bits)
                    : ReadIndexBucket(key, bucket_bits);
            const uint64_t directory_key =
                use_permuted_suffix
                    ? permuted_key & LowBitMask(directory_key_bits)
                    : key;
            UpdateBuildHll(ReadIndexCardinalityHash(key), bucket,
                           &local.bucket_hll);
            const uint64_t locator =
                (static_cast<uint64_t>(read_word_offset) << position_bits) |
                anchor_pos;
            uint8_t locator_context = 0u;
            if (locator_context_bits >= 2u &&
                anchor_pos >= locator_context_distance) {
              locator_context |= PackedBase(
                  sequence, anchor_pos - locator_context_distance);
            }
            if (locator_context_bits >= 4u &&
                anchor_pos + options.anchor_len +
                        locator_context_distance <=
                    length) {
              locator_context |= static_cast<uint8_t>(
                  PackedBase(sequence,
                             anchor_pos + options.anchor_len +
                                 locator_context_distance - 1u)
                  << 2u);
            }
            if (fused_edge_replay) {
              uint64_t offsets = fused_active_queries->QueryOffsets(key);
              while (offsets != 0u) {
                const unsigned offset =
                    static_cast<unsigned>(__builtin_ctzll(offsets));
                offsets &= offsets - 1u;
                if (anchor_pos < offset ||
                    !fused_active_queries->MayMatchContext(
                        key, offset, locator_context)) {
                  continue;
                }
                const unsigned start = anchor_pos - offset;
                if (start + options.window_len > length) continue;
                AnchorWindow window;
                window.InitFromPtr(sequence, start, options.window_len);
                if (fused_active_windows->Contains(window)) {
                  replay_positions.push_back(start);
                }
              }
            }
            const uint64_t temporary_locator =
                locator |
                (((directory_key >> temporary_key_bits) &
                  LowBitMask(temporary_key_extra_bits))
                 << locator_bits) |
                ((static_cast<uint64_t>(locator_context) &
                  locator_context_mask)
                 << (locator_bits + temporary_key_extra_bits));
            std::vector<uint8_t> &records = local.bucket_records[bucket];
            if (records.capacity() == 0) {
              records.reserve(bucket_reserve_bytes);
            }
            const size_t record_offset = records.size();
            records.resize(record_offset + temporary_record_bytes);
            buffered_bytes += temporary_record_bytes;
            StoreLowBytes(records.data() + record_offset, directory_key,
                          temporary_key_bytes);
            StoreLowBytes(records.data() + record_offset + temporary_key_bytes,
                          temporary_locator, locator_bytes);
          };

          if (stored_anchor_positions) {
            uint32_t count = 0;
            bool position_stream_ok = ReadReadAnchorVarint(
                &anchor_cursor, anchor_end, &count);
            uint32_t anchor_pos = 0u;
            for (uint32_t i = 0; position_stream_ok && i < count; ++i) {
              uint32_t delta = 0u;
              position_stream_ok = ReadReadAnchorVarint(
                  &anchor_cursor, anchor_end, &delta);
              if (!position_stream_ok || delta > UINT32_MAX - anchor_pos) {
                position_stream_ok = false;
                break;
              }
              anchor_pos += delta;
              if ((i != 0u && delta == 0u) ||
                  anchor_pos > length ||
                  options.anchor_len > length - anchor_pos) {
                position_stream_ok = false;
                break;
              }
              emit_anchor(ExtractPackedReadKey(
                              sequence, anchor_pos, options.anchor_len),
                          anchor_pos);
            }
            if (!position_stream_ok) {
              fill_ok.store(false, std::memory_order_relaxed);
              break;
            }
          } else {
            ForEachReadAnchor(sequence, length, options, &local.queue_pos,
                              &local.queue_hash, &local.queue_key,
                              emit_anchor);
          }
          if (!replay_positions.empty() &&
              fill_ok.load(std::memory_order_relaxed)) {
            std::sort(replay_positions.begin(), replay_positions.end());
            replay_positions.erase(
                std::unique(replay_positions.begin(), replay_positions.end()),
                replay_positions.end());
            ++local.replay_candidate_reads;
            local.replay_candidate_windows += replay_positions.size();
            uint64_t generated_edges = 0;
            if (fused_edge_replay->ReplayRead(
                    sequence, length, replay_positions, &local.replay_state,
                    &generated_edges)) {
              ++local.replay_aligned_reads;
            }
            local.replay_generated_edges += generated_edges;
          }
          word_pos += DivCeiling(static_cast<size_t>(length),
                                 SeqPackage::kBasesPerWord);
          if (--reads_until_flush_check == 0) {
            reads_until_flush_check = 64;
            if (buffered_bytes >= local_flush_bytes &&
                !enqueue_local(&local, &local_extents)) {
              fill_ok.store(false, std::memory_order_relaxed);
            } else if (buffered_bytes >= local_flush_bytes) {
              buffered_bytes = 0;
            }
          }
        }
        if (stored_anchor_positions && anchor_cursor != anchor_end) {
          fill_ok.store(false, std::memory_order_relaxed);
        }
        if (fill_ok.load(std::memory_order_relaxed)) {
          assert(word_pos == chunk.word_end);

          if (buffered_bytes >= local_flush_bytes &&
              !enqueue_local(&local, &local_extents)) {
            fill_ok.store(false, std::memory_order_relaxed);
          } else if (buffered_bytes >= local_flush_bytes) {
            buffered_bytes = 0;
          }
        }
        if (!fill_ok.load(std::memory_order_relaxed)) {
          for (std::vector<uint8_t> &records : local.bucket_records) {
            records.clear();
          }
          buffered_bytes = 0;
        }
        DiscardMemoryPages(
            const_cast<uint32_t *>(words + chunk.word_begin),
            (chunk.word_end - chunk.word_begin) * sizeof(uint32_t));
      }
      producer_cpu_nanos.fetch_add(
          ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID) - producer_cpu_begin,
          std::memory_order_relaxed);
    }

    if (fused_edge_replay && fill_ok.load(std::memory_order_relaxed)) {
      const double commit_begin = omp_get_wtime();
      const uint64_t batch_records = fused_edge_replay->CommitBatch();
      replay_commit_seconds += omp_get_wtime() - commit_begin;
      ++replay_commit_batches;
      const size_t chunks_in_wave = wave_end - wave_begin;
      if (batch_records != 0u && chunks_in_wave != 0u) {
        const uint64_t batch_bytes =
            batch_records * fused_edge_replay->candidate_record_bytes();
        const uint64_t bytes_per_chunk = std::max<uint64_t>(
            1u, DivCeiling(batch_bytes,
                           static_cast<uint64_t>(chunks_in_wave)));
        const size_t desired = static_cast<size_t>(std::max<uint64_t>(
            1u, replay_target_batch_bytes / bytes_per_chunk));
        const size_t growth_limit =
            replay_wave_chunks > std::numeric_limits<size_t>::max() / 4u
                ? std::numeric_limits<size_t>::max()
                : replay_wave_chunks * 4u;
        replay_wave_chunks = std::min(desired, growth_limit);
      }
    }
    wave_begin = wave_end;
  }

  // Publish each producer's final partial index buffer after the last replay
  // batch.  Keeping these tails across wave boundaries avoids one tiny write
  // batch per producer per replay commit.
#pragma omp parallel num_threads(num_compute_threads)
  {
    const uint64_t producer_cpu_begin =
        ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID);
    const size_t tid = static_cast<size_t>(omp_get_thread_num());
    ReadIndexBuildScratch &local = scratch[tid];
    BucketExtentLists &local_extents = thread_bucket_extents[tid];
    if (fill_ok.load(std::memory_order_relaxed) &&
        !enqueue_local(&local, &local_extents)) {
      fill_ok.store(false, std::memory_order_relaxed);
    }
    producer_cpu_nanos.fetch_add(
        ClockNanoseconds(CLOCK_THREAD_CPUTIME_ID) - producer_cpu_begin,
        std::memory_order_relaxed);
  }
  {
    std::lock_guard<std::mutex> lock(write_queue_mutex);
    producers_done = true;
  }
  write_queue_ready.notify_all();
  for (std::thread &writer : write_threads) writer.join();
  fill_timer.stop();
  const double fill_cpu_seconds = ProcessCpuSeconds() - fill_cpu_begin;
  xinfo("Read-index fill detail: producer CPU {.3} s / queue wait {.3} "
        "thread-s; writer CPU {.3} s / copy {.3} thread-s / I/O {.3} "
        "thread-s; {} batches, {} flushes, {.3} GiB written; allocation "
        "{.3} thread-s / {} calls\n",
        producer_cpu_nanos.load(std::memory_order_relaxed) * 1e-9,
        producer_wait_nanos.load(std::memory_order_relaxed) * 1e-9,
        writer_cpu_nanos.load(std::memory_order_relaxed) * 1e-9,
        writer_copy_nanos.load(std::memory_order_relaxed) * 1e-9,
        writer_io_nanos.load(std::memory_order_relaxed) * 1e-9,
        writer_batches.load(std::memory_order_relaxed),
        writer_flushes.load(std::memory_order_relaxed),
        static_cast<double>(writer_bytes.load(std::memory_order_relaxed)) /
            (uint64_t{1} << 30u),
        allocation_nanos.load(std::memory_order_relaxed) * 1e-9,
        allocation_calls.load(std::memory_order_relaxed));

  if (fused_edge_replay) {
    uint64_t candidate_reads = 0;
    uint64_t candidate_windows = 0;
    uint64_t aligned_reads = 0;
    uint64_t generated_edges = 0;
    for (const ReadIndexBuildScratch &local : scratch) {
      candidate_reads += local.replay_candidate_reads;
      candidate_windows += local.replay_candidate_windows;
      aligned_reads += local.replay_aligned_reads;
      generated_edges += local.replay_generated_edges;
    }
    const size_t unique_edges = fused_edge_replay->unique_edges();
    const uint64_t batch_unique_edges =
        fused_edge_replay->batch_unique_candidates();
    SimpleTimer edge_flush_timer;
    if (fill_ok.load(std::memory_order_relaxed)) {
      edge_flush_timer.start();
      try {
        fused_edge_replay->Flush();
      } catch (...) {
        close_temporary_files();
        remove_temporary_files();
        throw;
      }
      edge_flush_timer.stop();
      xinfo("Fused first iteration: {} exact-window candidates in {} gated "
            "/ {} aligned reads, {} "
            "generated / {} batch-distinct / {} globally unique edges; {} "
            "bounded radix commits in {.4} "
            "s, shared read/index fill {.4}, edge output {.4} s\n",
            candidate_windows, candidate_reads, aligned_reads,
            generated_edges,
            batch_unique_edges, unique_edges, replay_commit_batches,
            replay_commit_seconds,
            fill_timer.elapsed(), edge_flush_timer.elapsed());
    }
    // The compact-index phase has a separate explicit memory budget.  Drop
    // the query tables and edge set before admitting its large bucket tasks.
    fused_edge_replay.reset();
    fused_prepared_flank_index.reset();
    fused_active_windows.reset();
    fused_active_queries.reset();
  }

  SimpleTimer finalize_timer;
  finalize_timer.start();
  // RabbitFX-style bounded reuse: only a fixed number of write blocks exists,
  // and all retained pages are discarded once the producer/consumer phase is
  // complete.  Keeping an unbounded batch cache previously doubled the fill
  // arena and carried several GiB of dead capacity into compact.
  for (ReadIndexBuildScratch &local : scratch) {
    for (std::vector<uint8_t> &records : local.bucket_records) {
      DiscardMemoryPages(records.data(), records.capacity());
    }
  }
  for (std::unique_ptr<WriteBatch> &batch : reusable_writes) {
    for (std::vector<uint8_t> &records : batch->bucket_records) {
      DiscardMemoryPages(records.data(), records.capacity());
    }
  }
  reusable_writes.clear();
  pending_writes.clear();
  if (!fill_ok.load()) {
    close_temporary_files();
    remove_temporary_files();
    throw std::runtime_error("failed while streaming read-index occurrences");
  }

  std::vector<BucketWriteBuffer> bucket_writers(num_index_buckets);
#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int bucket_id = 0;
       bucket_id < static_cast<int>(num_index_buckets); ++bucket_id) {
    const unsigned bucket = static_cast<unsigned>(bucket_id);
    BucketWriteBuffer &writer = bucket_writers[bucket];
    size_t extent_count = 0;
    for (const BucketExtentLists &thread_extents : thread_bucket_extents) {
      extent_count += thread_extents[bucket].size();
    }
    for (const BucketExtentLists &writer_extents : writer_bucket_extents) {
      extent_count += writer_extents[bucket].size();
    }
    writer.extents.reserve(extent_count);
    for (BucketExtentLists &thread_extents : thread_bucket_extents) {
      std::vector<BucketExtent> &extents = thread_extents[bucket];
      for (const BucketExtent &extent : extents) {
        writer.records_written += extent.records;
        writer.extents.push_back(extent);
      }
      std::vector<BucketExtent>().swap(extents);
    }
    for (BucketExtentLists &writer_extents : writer_bucket_extents) {
      std::vector<BucketExtent> &extents = writer_extents[bucket];
      for (const BucketExtent &extent : extents) {
        writer.records_written += extent.records;
        writer.extents.push_back(extent);
      }
      std::vector<BucketExtent>().swap(extents);
    }
    // Lists are appended in file-id order, and the sole owner of each file
    // emits monotonically increasing offsets.  A second per-bucket sort used
    // to reorder millions of entries that are already globally ordered.
    assert(std::is_sorted(
        writer.extents.begin(), writer.extents.end(),
        [](const BucketExtent &lhs, const BucketExtent &rhs) {
          return lhs.file_id < rhs.file_id ||
                 (lhs.file_id == rhs.file_id &&
                  lhs.byte_offset < rhs.byte_offset);
        }));
  }
  std::vector<BucketExtentLists>().swap(thread_bucket_extents);
  std::vector<BucketExtentLists>().swap(writer_bucket_extents);

  std::vector<uint64_t> bucket_offsets(num_index_buckets + 1u, 0);
  std::vector<uint64_t> estimated_unique_keys(num_index_buckets, 0);
  uint64_t max_bucket_count = 0;
  for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
    const uint64_t count = bucket_writers[bucket].records_written;
    bucket_offsets[bucket + 1u] = bucket_offsets[bucket] + count;
    max_bucket_count = std::max(max_bucket_count, count);
  }

#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int bucket_id = 0;
       bucket_id < static_cast<int>(num_index_buckets); ++bucket_id) {
    const unsigned bucket = static_cast<unsigned>(bucket_id);
    const uint64_t count = bucket_writers[bucket].records_written;
    long double harmonic = 0.0L;
    unsigned zero_registers = 0;
    for (unsigned register_id = 0; register_id < kBuildHllRegisters;
         ++register_id) {
      uint8_t value = 0;
      for (const ReadIndexBuildScratch &local : scratch) {
        value = std::max(
            value,
            local.bucket_hll[static_cast<size_t>(bucket) *
                                 kBuildHllRegisters +
                             register_id]);
      }
      zero_registers += value == 0;
      harmonic += std::ldexp(1.0L, -static_cast<int>(value));
    }
    long double estimate =
        0.709L * kBuildHllRegisters * kBuildHllRegisters / harmonic;
    if (zero_registers != 0 &&
        estimate <= 2.5L * kBuildHllRegisters) {
      estimate = kBuildHllRegisters *
                 std::log(static_cast<long double>(kBuildHllRegisters) /
                          zero_registers);
    }
    estimated_unique_keys[bucket] = std::min<uint64_t>(
        count, static_cast<uint64_t>(std::ceil(estimate)));
  }
  const uint64_t total_occurrences = bucket_offsets[num_index_buckets];
  if (total_occurrences >
      std::numeric_limits<uint64_t>::max() / temporary_record_bytes) {
    throw std::length_error("temporary read index is too large");
  }
  const uint64_t temporary_bytes =
      total_occurrences * temporary_record_bytes;
  uint64_t observed_temporary_bytes = 0;
  for (uint64_t file_bytes : temporary_write_cursors) {
    if (file_bytes > std::numeric_limits<uint64_t>::max() -
                         observed_temporary_bytes) {
      observed_temporary_bytes = std::numeric_limits<uint64_t>::max();
      break;
    }
    observed_temporary_bytes += file_bytes;
  }
  if (observed_temporary_bytes != temporary_bytes) {
    close_temporary_files();
    remove_temporary_files();
    throw std::runtime_error("streamed read-index size mismatch");
  }
  for (unsigned file_id = 0; file_id < num_temporary_files; ++file_id) {
    if (ftruncate(temporary_fds[file_id],
                  static_cast<off_t>(temporary_write_cursors[file_id])) !=
        0) {
      close_temporary_files();
      remove_temporary_files();
      throw std::runtime_error("cannot trim sharded read-index allocation");
    }
  }
  uint64_t num_temporary_extents = 0;
  for (const BucketWriteBuffer &writer : bucket_writers) {
    num_temporary_extents += writer.extents.size();
  }
  const size_t num_fill_chunks = chunks.size();
  std::vector<ReadIndexBuildScratch>().swap(scratch);
  std::vector<BinaryChunk>().swap(chunks);

  const unsigned bucket_begin_bits = BitsNeeded(max_bucket_count);
  const unsigned bucket_begin_bytes = DivCeiling(bucket_begin_bits, 8u);
  const unsigned directory_entry_bytes = key_bytes + bucket_begin_bytes + 4u;
  const unsigned directory_fence_stride =
      NativeDirectoryFenceStride(directory_entry_bytes);
  long open_file_limit = sysconf(_SC_OPEN_MAX);
  if (open_file_limit <= 64) open_file_limit = 68;
  const unsigned maximum_output_files = static_cast<unsigned>(
      std::max<long>(1, (open_file_limit - 64) / 4));
  constexpr uint64_t kMinimumOutputShardWork = uint64_t{16} << 20u;
  const unsigned volume_output_files = static_cast<unsigned>(
      std::min<uint64_t>(
          num_index_buckets,
          std::max<uint64_t>(
              1u, DivCeiling(temporary_bytes, kMinimumOutputShardWork))));
  const unsigned num_output_files =
      use_permuted_suffix
          ? std::max<unsigned>(
                1u, std::min<unsigned>(
                        std::min<unsigned>(
                            static_cast<unsigned>(num_threads),
                            volume_output_files),
                        maximum_output_files))
          : 1u;
  std::vector<std::string> directory_paths(num_output_files);
  std::vector<std::string> posting_paths(num_output_files);
  std::vector<int> directory_fds(num_output_files, -1);
  std::vector<int> posting_fds(num_output_files, -1);
  for (unsigned file_id = 0; file_id < num_output_files; ++file_id) {
    directory_paths[file_id] =
        ReadIndexShardPath(options.output_prefix, "keys", file_id);
    posting_paths[file_id] =
        ReadIndexShardPath(options.output_prefix, "postings", file_id);
    directory_fds[file_id] =
        open(directory_paths[file_id].c_str(),
             O_CREAT | O_TRUNC | O_RDWR | O_CLOEXEC, 0664);
    posting_fds[file_id] =
        open(posting_paths[file_id].c_str(),
             O_CREAT | O_TRUNC | O_RDWR | O_CLOEXEC, 0664);
    if (directory_fds[file_id] < 0 || posting_fds[file_id] < 0) {
      for (int fd : directory_fds) {
        if (fd >= 0) close(fd);
      }
      for (int fd : posting_fds) {
        if (fd >= 0) close(fd);
      }
      close_temporary_files();
      remove_temporary_files();
      throw std::runtime_error("cannot create compact read-index shards");
    }
  }
  const auto close_output_files = [&]() {
    for (int &fd : directory_fds) {
      if (fd >= 0) {
        close(fd);
        fd = -1;
      }
    }
    for (int &fd : posting_fds) {
      if (fd >= 0) {
        close(fd);
        fd = -1;
      }
    }
  };

  const auto saturated_add = [](uint64_t lhs, uint64_t rhs) {
    return rhs > std::numeric_limits<uint64_t>::max() - lhs
               ? std::numeric_limits<uint64_t>::max()
               : lhs + rhs;
  };
  const auto saturated_multiply = [](uint64_t lhs, uint64_t rhs) {
    return lhs != 0 && rhs > std::numeric_limits<uint64_t>::max() / lhs
               ? std::numeric_limits<uint64_t>::max()
               : lhs * rhs;
  };
  const uint64_t task_slack_denominator =
      static_cast<uint64_t>(num_threads) * 8u;
  const uint64_t allocator_slack = std::max<uint64_t>(
      uint64_t{1} << 20u, requested_memory / task_slack_denominator);
  std::vector<uint64_t> unique_capacity(num_index_buckets, 0);
  std::vector<uint64_t> compact_task_memory(num_index_buckets, 1);
  std::vector<unsigned> compact_order(num_index_buckets, 0);
  uint64_t maximum_task_memory = 1;
  uint64_t max_estimated_unique = 0;
  for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
    const uint64_t count =
        bucket_offsets[bucket + 1u] - bucket_offsets[bucket];
    // p=6 HLL has ~13% relative error.  A 50% reserve covers estimation
    // variance, phmap growth, and the 12.5% reserve used by each task.
    unique_capacity[bucket] = std::min<uint64_t>(
        count, estimated_unique_keys[bucket] >
                       (std::numeric_limits<uint64_t>::max() - 1u) / 3u * 2u
                   ? count
                   : (estimated_unique_keys[bucket] * 3u + 1u) / 2u);
    max_estimated_unique =
        std::max(max_estimated_unique, estimated_unique_keys[bucket]);
    const bool use_direct_partition =
        directory_key_bits > 16u && directory_key_bits <= 28u &&
        count >= (uint64_t{1} << 16u);
    const uint64_t occurrence_memory = saturated_multiply(
        count, temporary_record_bytes + locator_bytes +
                   (use_direct_partition ? temporary_record_bytes : 0u));
    const uint64_t unique_memory = saturated_multiply(
        unique_capacity[bucket],
        sizeof(uint64_t) + directory_entry_bytes + 32u);
    compact_task_memory[bucket] = saturated_add(
        saturated_add(occurrence_memory, unique_memory), allocator_slack);
    maximum_task_memory =
        std::max(maximum_task_memory, compact_task_memory[bucket]);
    compact_order[bucket] = bucket;
  }
  // Largest-first scheduling prevents a skewed shard becoming a serial tail.
  // A byte-granular limiter below admits as many independent tasks as their
  // actual estimates allow instead of globally capping all work by the single
  // largest shard.  Keep one eighth of the requested memory for allocator,
  // extent, fence, and runtime state not owned by a compact task.
  std::sort(compact_order.begin(), compact_order.end(),
            [&](unsigned lhs, unsigned rhs) {
              return compact_task_memory[lhs] > compact_task_memory[rhs];
            });
  const uint64_t task_memory_budget = std::max<uint64_t>(
      maximum_task_memory, requested_memory - requested_memory / 8u);
  xinfo("Building exact read index: {} reads / {} bases, {} anchors in {} "
        "buckets; {}-byte temp {.3} GiB, {} producer chunks / {} sequential "
        "extents, max bucket {} records / ~{} unique, "
        "up to {} memory-token compact workers / {.3} GiB task budget\n",
        num_reads, num_bases, total_occurrences, num_index_buckets,
        temporary_record_bytes,
        static_cast<double>(temporary_bytes) / (uint64_t{1} << 30u),
        num_fill_chunks, num_temporary_extents, max_bucket_count,
        max_estimated_unique, num_threads,
        static_cast<double>(task_memory_budget) / (uint64_t{1} << 30u));

  std::vector<uint64_t> directory_cursors(num_output_files, 0);
  std::vector<uint64_t> posting_cursors(num_output_files, 0);
  std::vector<std::mutex> output_file_mutexes(num_output_files);
  std::atomic<bool> compact_ok(true);
  std::vector<ReadIndexBucketMeta> bucket_meta(num_index_buckets);
  std::vector<std::vector<uint8_t>> bucket_fence_keys(num_index_buckets);
  std::atomic<uint64_t> total_keys(0);
  std::atomic<uint64_t> compact_read_nanos(0);
  std::atomic<uint64_t> compact_count_nanos(0);
  std::atomic<uint64_t> compact_fill_nanos(0);
  std::atomic<uint64_t> compact_write_nanos(0);
  std::atomic<uint64_t> compact_read_calls(0);
  std::atomic<uint64_t> compact_read_bytes(0);
  CompactMemoryLimiter compact_memory_limiter(task_memory_budget);
  SimpleTimer compact_timer;
  const double compact_cpu_begin = ProcessCpuSeconds();
  finalize_timer.stop();
  compact_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
  for (int task_id = 0; task_id < static_cast<int>(num_index_buckets);
       ++task_id) {
    if (!compact_ok.load(std::memory_order_relaxed)) {
      continue;
    }
    const unsigned bucket = compact_order[static_cast<size_t>(task_id)];
    CompactMemoryLease memory_lease(&compact_memory_limiter,
                                    compact_task_memory[bucket]);
    const uint64_t count =
        bucket_offsets[bucket + 1u] - bucket_offsets[bucket];
    if (count > std::numeric_limits<size_t>::max() /
                    temporary_record_bytes) {
      compact_ok.store(false, std::memory_order_relaxed);
      continue;
    }
    std::vector<uint8_t> packed_records(
        static_cast<size_t>(count) * temporary_record_bytes);
    const uint64_t read_begin = ClockNanoseconds(CLOCK_MONOTONIC);
    uint64_t task_read_calls = 0;
    uint64_t task_read_bytes = 0;
    size_t packed_offset = 0;
    for (const BucketExtent &extent : bucket_writers[bucket].extents) {
      const size_t extent_bytes =
          static_cast<size_t>(extent.records) * temporary_record_bytes;
      ++task_read_calls;
      task_read_bytes += extent_bytes;
      if (extent_bytes > packed_records.size() - packed_offset ||
          extent.file_id >= temporary_fds.size() ||
          !PreadFully(temporary_fds[extent.file_id],
                      packed_records.data() + packed_offset, extent_bytes,
                      extent.byte_offset)) {
        compact_ok.store(false, std::memory_order_relaxed);
        break;
      }
      packed_offset += extent_bytes;
    }
    compact_read_nanos.fetch_add(
        ClockNanoseconds(CLOCK_MONOTONIC) - read_begin,
        std::memory_order_relaxed);
    compact_read_calls.fetch_add(task_read_calls, std::memory_order_relaxed);
    compact_read_bytes.fetch_add(task_read_bytes, std::memory_order_relaxed);
    if (!compact_ok.load(std::memory_order_relaxed) ||
        packed_offset != packed_records.size()) {
      compact_ok.store(false, std::memory_order_relaxed);
      continue;
    }
    const auto load_temporary_locator = [&](const uint8_t *source) {
      return LoadLowBytes(source + temporary_key_bytes, locator_bytes);
    };
    const auto load_temporary_key = [&](const uint8_t *source) {
      const uint64_t low = LoadLowBytes(source, temporary_key_bytes);
      if (temporary_key_extra_bits == 0) return low;
      return low |
             (((load_temporary_locator(source) >> locator_bits) &
               LowBitMask(temporary_key_extra_bits))
              << temporary_key_bits);
    };
    const auto load_final_locator = [&](const uint8_t *source) {
      const uint64_t temporary = load_temporary_locator(source);
      const uint64_t context =
          (temporary >> (locator_bits + temporary_key_extra_bits)) &
          locator_context_mask;
      return (temporary & locator_mask) | (context << locator_bits);
    };
    // RabbitTClust's CSR construction is the right join order here: count
    // unique keys first, sort only that compact key set, then fill one flat
    // postings array.  The old path radix-sorted every occurrence even though
    // most records share a key (7.03B occurrences versus 618M keys on CAMI).
    phmap::flat_hash_map<uint64_t, uint64_t> key_counts;
    std::vector<uint64_t> keys;
    std::vector<uint32_t> group_counts;
    std::vector<size_t> partition_offsets;
    std::vector<size_t> partition_group_offsets;
    std::vector<uint64_t> direct_table;
    const uint64_t count_begin = ClockNanoseconds(CLOCK_MONOTONIC);
    const bool use_direct_partition =
        directory_key_bits > 16u && directory_key_bits <= 28u &&
        count >= (uint64_t{1} << 16u);
    if (use_direct_partition) {
      // One stable high-bit partition turns two global random hash probes per
      // occurrence into direct accesses to a 16-bit table that fits in cache.
      // Within each high partition, sorting only the touched low keys emits
      // globally ordered directory keys without sorting all occurrences.
      constexpr unsigned kDirectLowBits = 16u;
      constexpr uint64_t kDirectLowMask = (uint64_t{1} << kDirectLowBits) - 1u;
      const unsigned partition_bits = directory_key_bits - kDirectLowBits;
      const size_t num_partitions = size_t{1} << partition_bits;
      partition_offsets.assign(num_partitions + 1u, 0);
      for (size_t record = 0;
           record < packed_records.size() / temporary_record_bytes; ++record) {
        const uint8_t *source =
            packed_records.data() + record * temporary_record_bytes;
        const size_t partition = static_cast<size_t>(
            load_temporary_key(source) >> kDirectLowBits);
        ++partition_offsets[partition + 1u];
      }
      for (size_t partition = 0; partition < num_partitions; ++partition) {
        partition_offsets[partition + 1u] += partition_offsets[partition];
      }
      std::vector<size_t> partition_cursors = partition_offsets;
      std::vector<uint8_t> partitioned_records(packed_records.size());
      for (size_t record = 0;
           record < packed_records.size() / temporary_record_bytes; ++record) {
        const uint8_t *source =
            packed_records.data() + record * temporary_record_bytes;
        const size_t partition = static_cast<size_t>(
            load_temporary_key(source) >> kDirectLowBits);
        const size_t destination = partition_cursors[partition]++;
        std::memcpy(partitioned_records.data() +
                        destination * temporary_record_bytes,
                    source, temporary_record_bytes);
      }
      packed_records.swap(partitioned_records);
      std::vector<uint8_t>().swap(partitioned_records);
      std::vector<size_t>().swap(partition_cursors);

      direct_table.assign(size_t{1} << kDirectLowBits, 0);
      std::vector<uint32_t> touched;
      touched.reserve(size_t{1} << kDirectLowBits);
      keys.reserve(static_cast<size_t>(std::min<uint64_t>(
          unique_capacity[bucket], std::numeric_limits<size_t>::max())));
      group_counts.reserve(keys.capacity());
      partition_group_offsets.resize(num_partitions + 1u, 0);
      for (size_t partition = 0; partition < num_partitions; ++partition) {
        partition_group_offsets[partition] = keys.size();
        touched.clear();
        for (size_t record = partition_offsets[partition];
             record < partition_offsets[partition + 1u]; ++record) {
          const uint8_t *source =
              packed_records.data() + record * temporary_record_bytes;
          const uint32_t low = static_cast<uint32_t>(
              load_temporary_key(source) & kDirectLowMask);
          if (direct_table[low] == 0) touched.push_back(low);
          if (direct_table[low] == std::numeric_limits<uint32_t>::max()) {
            compact_ok.store(false, std::memory_order_relaxed);
            break;
          }
          ++direct_table[low];
        }
        if (!compact_ok.load(std::memory_order_relaxed)) break;
        std::sort(touched.begin(), touched.end());
        for (uint32_t low : touched) {
          keys.push_back((static_cast<uint64_t>(partition) <<
                          kDirectLowBits) |
                         low);
          group_counts.push_back(static_cast<uint32_t>(direct_table[low]));
          direct_table[low] = 0;
        }
      }
      partition_group_offsets[num_partitions] = keys.size();
    } else {
      // The previous first-64K extrapolation saw almost exclusively new keys
      // and reserved near the occurrence count.  The order-independent HLL
      // drives reserve for the hash fallback used by genuinely small tasks.
      if (unique_capacity[bucket] <=
          static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        key_counts.reserve(static_cast<size_t>(unique_capacity[bucket]));
      }
      for (size_t record = 0;
           record < packed_records.size() / temporary_record_bytes; ++record) {
        const uint8_t *source =
            packed_records.data() + record * temporary_record_bytes;
        const uint64_t key = load_temporary_key(source);
        auto inserted = key_counts.emplace(key, 0u);
        ++inserted.first->second;
      }
      keys.reserve(key_counts.size());
      for (const auto &entry : key_counts) keys.push_back(entry.first);
      SortReadIndexKeys(&keys, key_bytes);
      group_counts.reserve(keys.size());
      for (uint64_t key : keys) {
        const auto cursor = key_counts.find(key);
        if (cursor == key_counts.end() || cursor->second == 0 ||
            cursor->second > std::numeric_limits<uint32_t>::max()) {
          compact_ok.store(false, std::memory_order_relaxed);
          break;
        }
        group_counts.push_back(static_cast<uint32_t>(cursor->second));
      }
    }
    if (!compact_ok.load(std::memory_order_relaxed) ||
        group_counts.size() != keys.size()) {
      compact_ok.store(false, std::memory_order_relaxed);
      continue;
    }
    compact_count_nanos.fetch_add(
        ClockNanoseconds(CLOCK_MONOTONIC) - count_begin,
        std::memory_order_relaxed);

    const uint64_t posting_fill_begin = ClockNanoseconds(CLOCK_MONOTONIC);
    const uint64_t groups = keys.size();
    std::vector<uint8_t> directory(
        static_cast<size_t>(groups) * directory_entry_bytes);
    std::vector<uint8_t> postings(
        static_cast<size_t>(count) * locator_bytes);
    size_t posting_index = 0;
    for (size_t group_index = 0; group_index < keys.size(); ++group_index) {
      const uint64_t key = keys[group_index];
      const uint64_t group_count = group_counts[group_index];
      if (group_count == 0 || group_count > count - posting_index) {
        compact_ok.store(false, std::memory_order_relaxed);
        break;
      }
      uint8_t *entry =
          directory.data() + group_index * directory_entry_bytes;
      StoreLowBytes(entry, key, key_bytes);
      StoreLowBytes(entry + key_bytes, posting_index, bucket_begin_bytes);
      StoreLowBytes(entry + key_bytes + bucket_begin_bytes, group_count,
                    4u);
      if (!use_direct_partition) {
        auto cursor = key_counts.find(key);
        if (cursor == key_counts.end()) {
          compact_ok.store(false, std::memory_order_relaxed);
          break;
        }
        cursor->second = posting_index;
      }
      posting_index += group_count;
    }
    if (!compact_ok.load(std::memory_order_relaxed) || posting_index != count) {
      compact_ok.store(false, std::memory_order_relaxed);
      continue;
    }

    // Posting order is immaterial to exact replay: candidates are sorted and
    // uniqued in read order before verification.  Avoiding a per-key locator
    // sort saves another full pass over this much larger array.
    if (use_direct_partition) {
      constexpr unsigned kDirectLowBits = 16u;
      constexpr uint64_t kDirectLowMask =
          (uint64_t{1} << kDirectLowBits) - 1u;
      const size_t num_partitions = partition_offsets.size() - 1u;
      for (size_t partition = 0; partition < num_partitions; ++partition) {
        const size_t group_begin = partition_group_offsets[partition];
        const size_t group_end = partition_group_offsets[partition + 1u];
        for (size_t group = group_begin; group < group_end; ++group) {
          const uint32_t low =
              static_cast<uint32_t>(keys[group] & kDirectLowMask);
          const uint8_t *entry =
              directory.data() + group * directory_entry_bytes;
          direct_table[low] =
              LoadLowBytes(entry + key_bytes, bucket_begin_bytes);
        }
        for (size_t packed = partition_offsets[partition];
             packed < partition_offsets[partition + 1u]; ++packed) {
          const uint8_t *source =
              packed_records.data() + packed * temporary_record_bytes;
          const uint32_t low = static_cast<uint32_t>(
              load_temporary_key(source) & kDirectLowMask);
          const uint64_t cursor = direct_table[low]++;
          if (cursor >= count) {
            compact_ok.store(false, std::memory_order_relaxed);
            break;
          }
          const uint64_t locator = load_final_locator(source);
          StoreLowBytes(postings.data() + cursor * locator_bytes, locator,
                        locator_bytes);
        }
        if (!compact_ok.load(std::memory_order_relaxed)) break;
        for (size_t group = group_begin; group < group_end; ++group) {
          const uint32_t low =
              static_cast<uint32_t>(keys[group] & kDirectLowMask);
          const uint8_t *entry =
              directory.data() + group * directory_entry_bytes;
          const uint64_t begin =
              LoadLowBytes(entry + key_bytes, bucket_begin_bytes);
          const uint64_t group_count = LoadLowBytes(
              entry + key_bytes + bucket_begin_bytes, 4u);
          if (direct_table[low] != begin + group_count) {
            compact_ok.store(false, std::memory_order_relaxed);
            break;
          }
          direct_table[low] = 0;
        }
        if (!compact_ok.load(std::memory_order_relaxed)) break;
      }
    } else {
      for (size_t packed = 0;
           packed < packed_records.size() / temporary_record_bytes;
           ++packed) {
        const uint8_t *source =
            packed_records.data() + packed * temporary_record_bytes;
        const uint64_t key = load_temporary_key(source);
        const uint64_t locator = load_final_locator(source);
        auto cursor = key_counts.find(key);
        if (cursor == key_counts.end() || cursor->second >= count) {
          compact_ok.store(false, std::memory_order_relaxed);
          break;
        }
        StoreLowBytes(postings.data() + cursor->second * locator_bytes,
                      locator, locator_bytes);
        ++cursor->second;
      }
      if (compact_ok.load(std::memory_order_relaxed)) {
        posting_index = 0;
        for (size_t group_index = 0; group_index < keys.size();
             ++group_index) {
          posting_index += group_counts[group_index];
          const auto cursor = key_counts.find(keys[group_index]);
          if (cursor == key_counts.end() ||
              cursor->second != posting_index) {
            compact_ok.store(false, std::memory_order_relaxed);
            break;
          }
        }
      }
    }
    if (!compact_ok.load(std::memory_order_relaxed)) continue;
    phmap::flat_hash_map<uint64_t, uint64_t>().swap(key_counts);
    std::vector<uint8_t>().swap(packed_records);
    compact_fill_nanos.fetch_add(
        ClockNanoseconds(CLOCK_MONOTONIC) - posting_fill_begin,
        std::memory_order_relaxed);

    std::vector<uint8_t> &fences = bucket_fence_keys[bucket];
    fences.resize(static_cast<size_t>(DivCeiling(
                      groups, static_cast<uint64_t>(directory_fence_stride))) *
                  key_bytes);
    for (uint64_t group = 0, fence = 0; group < groups;
         group += directory_fence_stride, ++fence) {
      std::memcpy(fences.data() + fence * key_bytes,
                  directory.data() + group * directory_entry_bytes,
                  key_bytes);
    }
    std::vector<uint64_t>().swap(keys);

    const unsigned output_file_id =
        static_cast<unsigned>(omp_get_thread_num()) % num_output_files;
    const uint64_t output_begin = ClockNanoseconds(CLOCK_MONOTONIC);
    uint64_t directory_offset = 0;
    uint64_t posting_offset = 0;
    bool write_ok = true;
    {
      std::lock_guard<std::mutex> lock(
          output_file_mutexes[output_file_id]);
      directory_offset = directory_cursors[output_file_id];
      posting_offset = posting_cursors[output_file_id];
      directory_cursors[output_file_id] += directory.size();
      posting_cursors[output_file_id] += postings.size();
      write_ok =
          (directory.empty() ||
           PwriteFully(directory_fds[output_file_id], directory.data(),
                       directory.size(), directory_offset)) &&
          (postings.empty() ||
           PwriteFully(posting_fds[output_file_id], postings.data(),
                       postings.size(), posting_offset));
    }
    if (!write_ok) {
      compact_ok.store(false, std::memory_order_relaxed);
      continue;
    }
    compact_write_nanos.fetch_add(
        ClockNanoseconds(CLOCK_MONOTONIC) - output_begin,
        std::memory_order_relaxed);
    bucket_meta[bucket] = ReadIndexBucketMeta{
        directory_offset, groups, posting_offset, count, output_file_id};
    total_keys.fetch_add(groups, std::memory_order_relaxed);
  }
  compact_timer.stop();
  SimpleTimer publish_timer;
  publish_timer.start();
  const double compact_cpu_seconds =
      ProcessCpuSeconds() - compact_cpu_begin;
  uint64_t total_directory_bytes = 0;
  uint64_t total_posting_bytes = 0;
  for (unsigned file_id = 0; file_id < num_output_files; ++file_id) {
    total_directory_bytes += directory_cursors[file_id];
    total_posting_bytes += posting_cursors[file_id];
  }
  xinfo("Read-index compact detail: read {.3} thread-s / {} calls / {.3} "
        "GiB; partition+count {.3} thread-s; CSR fill {.3} thread-s; "
        "output {.3} thread-s\n",
        compact_read_nanos.load(std::memory_order_relaxed) * 1e-9,
        compact_read_calls.load(std::memory_order_relaxed),
        static_cast<double>(
            compact_read_bytes.load(std::memory_order_relaxed)) /
            (uint64_t{1} << 30u),
        compact_count_nanos.load(std::memory_order_relaxed) * 1e-9,
        compact_fill_nanos.load(std::memory_order_relaxed) * 1e-9,
        compact_write_nanos.load(std::memory_order_relaxed) * 1e-9);
  close_output_files();
  if (!compact_ok.load()) {
    close_temporary_files();
    remove_temporary_files();
    throw std::runtime_error("failed while compacting read index");
  }

  const std::string info_path = options.output_prefix + ".ridx.info";
  std::ofstream info(info_path.c_str());
  if (!info) {
    remove_temporary_files();
    throw std::runtime_error("cannot write read-index metadata");
  }
  const unsigned output_version =
      locator_context_bits != 0u
          ? 4u
          : (num_output_files > 1 ? 3u
                                  : (use_permuted_suffix ? 2u : 1u));
  info << "version " << output_version << '\n'
       << "num_files " << num_output_files << '\n'
       << "anchor_len " << options.anchor_len << '\n'
       << "window_len " << options.window_len << '\n'
       << "bucket_bits " << bucket_bits << '\n'
       << "position_bits " << position_bits << '\n'
       << "locator_bytes " << locator_bytes << '\n'
       << "locator_context_bits " << locator_context_bits << '\n'
       << "locator_context_distance " << locator_context_distance << '\n'
       << "key_bytes " << key_bytes << '\n'
       << "bucket_begin_bytes " << bucket_begin_bytes << '\n'
       << "directory_entry_bytes " << directory_entry_bytes << '\n'
       << "postings_sorted 0\n"
       << "num_reads " << num_reads << '\n'
       << "num_bases " << num_bases << '\n'
       << "parsed_words " << parsed_words << '\n'
       << "num_keys " << total_keys.load() << '\n'
       << "num_occurrences " << total_occurrences << '\n';
  for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
    const ReadIndexBucketMeta &meta = bucket_meta[bucket];
    info << "bucket " << bucket << ' ';
    if (output_version >= 3) info << meta.file_id << ' ';
    info << meta.directory_offset << ' '
         << meta.directory_count << ' ' << meta.posting_offset << ' '
         << meta.posting_count << '\n';
  }
  info.close();

  DirectoryFenceIndex built_fences;
  built_fences.stride = directory_fence_stride;
  built_fences.bucket_offsets.assign(num_index_buckets + 1u, 0);
  for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
    if (bucket_fence_keys[bucket].size() % key_bytes != 0) {
      throw std::logic_error("invalid constructed directory fences");
    }
    built_fences.bucket_offsets[bucket + 1u] =
        built_fences.bucket_offsets[bucket] +
        bucket_fence_keys[bucket].size() / key_bytes;
  }
  built_fences.keys.resize(static_cast<size_t>(
      built_fences.bucket_offsets.back() * key_bytes));
  for (unsigned bucket = 0; bucket < num_index_buckets; ++bucket) {
    if (!bucket_fence_keys[bucket].empty()) {
      std::memcpy(built_fences.keys.data() +
                      built_fences.bucket_offsets[bucket] * key_bytes,
                  bucket_fence_keys[bucket].data(),
                  bucket_fence_keys[bucket].size());
    }
  }
  PublishDirectoryFenceIndex(
      options.output_prefix + ".ridx.fences", key_bytes, total_keys.load(),
      total_directory_bytes, built_fences);
  std::vector<uint8_t> temporary_cleanup_failed(num_temporary_files, 0);
  const int cleanup_threads = std::max<int>(
      1, std::min<int>(num_threads, num_temporary_files));
#pragma omp parallel for schedule(static) num_threads(cleanup_threads)
  for (int file = 0; file < static_cast<int>(num_temporary_files); ++file) {
    const unsigned file_id = static_cast<unsigned>(file);
    // Unlink while the descriptor is still open, then let independent close
    // operations reclaim each shard's blocks in parallel.  Sequential unlink
    // of a large, closed ext4 file otherwise serializes block-group cleanup.
    const bool removed = std::remove(temporary_paths[file_id].c_str()) == 0;
    const bool closed = temporary_fds[file_id] < 0 ||
                        close(temporary_fds[file_id]) == 0;
    temporary_fds[file_id] = -1;
    temporary_cleanup_failed[file_id] = !removed || !closed;
  }
  for (unsigned file_id = 0; file_id < num_temporary_files; ++file_id) {
    if (temporary_cleanup_failed[file_id] != 0) {
      // A failed unlink may have been transient; retry after the descriptor
      // has closed before leaving a large intermediate behind.
      if (std::remove(temporary_paths[file_id].c_str()) != 0 &&
          errno != ENOENT) {
        xwarn("Could not remove completed read-index temporary file {}\n",
              temporary_paths[file_id].c_str());
      }
    }
  }
  publish_timer.stop();
  xinfo("Exact read index complete: {} keys, directory {.3} GiB, postings "
        "{.3} GiB, fences {.3} MiB; boundary {.4}, one-pass partition/fill "
        "{.4}, unique-key compact {.4} s\n",
        total_keys.load(),
        static_cast<double>(total_directory_bytes) /
            (uint64_t{1} << 30u),
        static_cast<double>(total_posting_bytes) / (uint64_t{1} << 30u),
        static_cast<double>(built_fences.keys.size()) / (uint64_t{1} << 20u),
        boundary_timer.elapsed(), fill_timer.elapsed(), compact_timer.elapsed());
  xinfo("Read-index phase parallelism: boundary {.2} cores, fill {.2} cores "
        "({} compute + {} writers), compact {.2} cores; compact admitted "
        "{} workers / peak charged {.3} GiB\n",
        boundary_timer.elapsed() > 0.0
            ? boundary_cpu_seconds / boundary_timer.elapsed()
            : 0.0,
        fill_timer.elapsed() > 0.0
            ? fill_cpu_seconds / fill_timer.elapsed()
            : 0.0,
        num_compute_threads, num_writer_threads,
        compact_timer.elapsed() > 0.0
            ? compact_cpu_seconds / compact_timer.elapsed()
            : 0.0,
        compact_memory_limiter.peak_workers(),
        static_cast<double>(compact_memory_limiter.peak_used()) /
            (uint64_t{1} << 30u));
  xinfo("Read-index orchestration: extent/HLL/output setup {.4} s; metadata, "
        "close, and temporary cleanup {.4} s\n",
        finalize_timer.elapsed(), publish_timer.elapsed());
}

void RunProfile(const ReadIndexOptions &options) {
  ActiveQuerySet active_queries;
  const ActiveQuerySet *active_query_pointer = nullptr;
  if (!options.contig_file.empty() || !options.bubble_file.empty()) {
    SimpleTimer query_timer;
    query_timer.start();
    BuildActiveQueries(options, &active_queries);
    // The profiling scan probes arbitrary read-side keys.  Production index
    // replay instead enumerates the immutable compact query entries directly,
    // so it deliberately avoids materializing this second hash-table view.
    active_queries.PrepareRandomLookup();
    query_timer.stop();
    xinfo("Active query construction time elapsed: {.6}\n",
          query_timer.elapsed());
    active_query_pointer = &active_queries;
  }
  if (!options.cache_contig_file.empty() ||
      !options.cache_bubble_file.empty()) {
    ReadIndexOptions cache_options = options;
    cache_options.contig_file = options.cache_contig_file;
    cache_options.bubble_file = options.cache_bubble_file;
    cache_options.kmer_k = options.cache_kmer_k;
    cache_options.step = options.cache_step;
    cache_options.cache_contig_file.clear();
    cache_options.cache_bubble_file.clear();
    ActiveQuerySet cached_queries;
    SimpleTimer cache_timer;
    cache_timer.start();
    BuildActiveQueries(cache_options, &cached_queries);
    cache_timer.stop();
    uint64_t covered_keys = 0;
    uint64_t covered_queries = 0;
    active_queries.ForEach([&](uint64_t key,
                               const ActiveQuerySet::QueryInfo &info) {
      if (cached_queries.Contains(key)) {
        ++covered_keys;
        covered_queries += info.query_count;
      }
    });
    xinfo("Cached-anchor coverage: {} / {} active keys ({.6}), {} / {} "
          "oriented queries ({.6}); cache query construction {.6} s\n",
          covered_keys, active_queries.key_count(),
          active_queries.key_count() == 0
              ? 0.0
              : static_cast<double>(covered_keys) /
                    active_queries.key_count(),
          covered_queries, active_queries.query_count(),
          active_queries.query_count() == 0
              ? 0.0
              : static_cast<double>(covered_queries) /
                    active_queries.query_count(),
          cache_timer.elapsed());
  }
  MappedWords mapping(options.read_file);
  SimpleTimer parse_timer;
  parse_timer.start();
  uint64_t num_reads = 0;
  uint64_t num_bases = 0;
  unsigned max_read_len = 0;
  size_t parsed_words = 0;
  const uint64_t max_reads =
      options.max_reads > 0 ? static_cast<uint64_t>(options.max_reads) : 0;
  std::vector<BinaryChunk> chunks =
      ParseChunks(mapping, max_reads, options.num_threads, &num_reads,
                  &num_bases, &max_read_len, &parsed_words);
  parse_timer.stop();
  xinfo("Read-index input: {} bytes, {} parsed reads, {} bases, max length {}, "
        "{} chunks; boundary scan {.6} s\n",
        parsed_words * sizeof(uint32_t), num_reads, num_bases, max_read_len,
        chunks.size(), parse_timer.elapsed());

  const int num_threads = std::max(1, options.num_threads);
  std::vector<ThreadProfile> profiles(static_cast<size_t>(num_threads));
  const uint64_t sample_mask =
      (uint64_t{1} << options.sample_bits) - 1u;
  const uint32_t *words = mapping.data();
  SimpleTimer profile_timer;
  profile_timer.start();
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
  for (size_t chunk_id = 0; chunk_id < chunks.size(); ++chunk_id) {
    ThreadProfile &profile =
        profiles[static_cast<size_t>(omp_get_thread_num())];
    const BinaryChunk &chunk = chunks[chunk_id];
    size_t word_pos = chunk.word_begin;
    while (word_pos < chunk.word_end) {
      const uint32_t length = words[word_pos++];
      const uint32_t *sequence = words + word_pos;
      ProfileRead(sequence, length, options, sample_mask,
                  active_query_pointer, &profile);
      word_pos +=
          DivCeiling(static_cast<size_t>(length), SeqPackage::kBasesPerWord);
    }
    assert(word_pos == chunk.word_end);
    DiscardMemoryPages(
        const_cast<uint32_t *>(words + chunk.word_begin),
        (chunk.word_end - chunk.word_begin) * sizeof(uint32_t));
  }
  profile_timer.stop();

  uint64_t profile_reads = 0;
  uint64_t profile_bases = 0;
  uint64_t total_windows = 0;
  uint64_t total_anchors = 0;
  uint64_t active_postings = 0;
  uint64_t candidate_pairs = 0;
  std::array<uint8_t, kHllRegisters> merged_hll{};
  size_t sampled_key_reserve = 0;
  for (const ThreadProfile &profile : profiles) {
    profile_reads += profile.reads;
    profile_bases += profile.bases;
    total_windows += profile.windows;
    total_anchors += profile.anchors;
    active_postings += profile.active_postings;
    candidate_pairs += profile.candidate_pairs;
    sampled_key_reserve += profile.sampled_postings.size();
    for (size_t i = 0; i < kHllRegisters; ++i) {
      merged_hll[i] = std::max(merged_hll[i], profile.hll[i]);
    }
  }
  if (profile_reads != num_reads || profile_bases != num_bases) {
    throw std::runtime_error("parallel read-index profile lost input records");
  }

  phmap::flat_hash_map<uint64_t, uint64_t> sampled_postings;
  sampled_postings.reserve(sampled_key_reserve);
  for (const ThreadProfile &profile : profiles) {
    for (const auto &entry : profile.sampled_postings) {
      sampled_postings[entry.first] += entry.second;
    }
  }
  std::vector<uint64_t> sampled_counts;
  sampled_counts.reserve(sampled_postings.size());
  uint64_t sampled_occurrences = 0;
  for (const auto &entry : sampled_postings) {
    sampled_counts.push_back(entry.second);
    sampled_occurrences += entry.second;
  }
  std::sort(sampled_counts.begin(), sampled_counts.end());

  double estimated_unique = EstimateHll(merged_hll);
  estimated_unique = std::min(estimated_unique,
                              static_cast<double>(total_anchors));
  const unsigned word_offset_bits =
      BitsNeeded(parsed_words == 0 ? 0 : parsed_words - 1u);
  const unsigned read_position_bits =
      BitsNeeded(max_read_len == 0 ? 0 : max_read_len - 1u);
  const unsigned locator_bits = word_offset_bits + read_position_bits;
  const unsigned locator_bytes = DivCeiling(locator_bits, 8u);
  const long double posting_bytes =
      static_cast<long double>(total_anchors) * locator_bytes;
  const long double directory_bytes =
      static_cast<long double>(estimated_unique) * 16.0L;
  constexpr long double kGiB =
      static_cast<long double>(uint64_t{1} << 30u);

  xinfo("Exact anchor profile: seed {}, window {}, {} windows, {} stored "
        "occurrences, density {.6}; scan {.6} s\n",
        options.anchor_len, options.window_len, total_windows, total_anchors,
        total_windows == 0
            ? 0.0
            : static_cast<double>(total_anchors) / total_windows,
        profile_timer.elapsed());
  xinfo("Estimated unique anchor keys (HLL-{}): {.0}\n", kHllBits,
        estimated_unique);
  xinfo("Posting sample 1/{}, {} keys / {} occurrences: p50 {}, p95 {}, "
        "p99 {}, max {}, mean {.3}\n",
        uint64_t{1} << options.sample_bits, sampled_counts.size(),
        sampled_occurrences, Quantile(sampled_counts, 0.50),
        Quantile(sampled_counts, 0.95), Quantile(sampled_counts, 0.99),
        sampled_counts.empty() ? 0 : sampled_counts.back(),
        sampled_counts.empty()
            ? 0.0
            : static_cast<double>(sampled_occurrences) /
                  sampled_counts.size());
  xinfo("Packed locator estimate: {} word-offset bits + {} position bits = "
        "{} bits / {} bytes; postings {.3} GiB, 16-byte key directory {.3} "
        "GiB, total {.3} GiB\n",
        word_offset_bits, read_position_bits, locator_bits, locator_bytes,
        static_cast<double>(posting_bytes / kGiB),
        static_cast<double>(directory_bytes / kGiB),
        static_cast<double>((posting_bytes + directory_bytes) / kGiB));
  if (active_query_pointer != nullptr) {
    const long double active_bytes =
        static_cast<long double>(active_postings) * locator_bytes;
    xinfo("Active-index traffic estimate: {} postings ({.6} of index, {.3} "
          "GiB packed), {} unique candidate read starts\n",
          active_postings,
          total_anchors == 0
              ? 0.0
              : static_cast<double>(active_postings) / total_anchors,
          static_cast<double>(active_bytes / kGiB), candidate_pairs);
  }
}

}  // namespace

int main_read_index(int argc, char **argv) {
  AutoMaxRssRecorder recorder;
  ReadIndexOptions options;
  OptionsDescription description;
  description.AddOption("read_file", "r", options.read_file,
                        "packed reads.lib.bin file");
  description.AddOption("contig_file", "c", options.contig_file,
                        "optional current-k contig file for active queries");
  description.AddOption("bubble_file", "b", options.bubble_file,
                        "optional current-k bubble file for active queries");
  description.AddOption("cache_contig_file", "", options.cache_contig_file,
                        "optional earlier-k contigs whose anchors are cached");
  description.AddOption("cache_bubble_file", "", options.cache_bubble_file,
                        "optional earlier-k bubbles whose anchors are cached");
  description.AddOption("output_prefix", "o", options.output_prefix,
                        "build a persistent compact exact-anchor index");
  description.AddOption("index_prefix", "", options.index_prefix,
                        "verify an existing compact exact-anchor index");
  description.AddOption("edge_output_prefix", "e", options.edge_output_prefix,
                        "write exact replayed iterative edges");
  description.AddOption("local_candidate_output", "",
                        options.local_candidate_output,
                        "write exact candidate read offsets for local assembly");
  description.AddOption("anchor_len", "a", options.anchor_len,
                        "exact anchor length, at most 31 bases");
  description.AddOption("window_len", "w", options.window_len,
                        "query window covered by each selected anchor");
  description.AddOption("sample_bits", "", options.sample_bits,
                        "sample one key per 2^N for posting statistics");
  description.AddOption("kmer_k", "k", options.kmer_k,
                        "current k when profiling active flanks");
  description.AddOption("step", "s", options.step,
                        "next-k step when profiling active flanks");
  description.AddOption("cache_kmer_k", "", options.cache_kmer_k,
                        "k of the optional cached-anchor source");
  description.AddOption("cache_step", "", options.cache_step,
                        "step of the optional cached-anchor source");
  description.AddOption("local_endpoint_range", "",
                        options.local_endpoint_range,
                        "maximum local-assembly endpoint range");
  description.AddOption("local_seed_len", "", options.local_seed_len,
                        "exact local mapper seed length");
  description.AddOption("local_sparsity", "", options.local_sparsity,
                        "local mapper contig-index sparsity");
  description.AddOption("local_min_contig_len", "",
                        options.local_min_contig_len,
                        "minimum contig length used by local mapping");
  description.AddOption("num_cpu_threads", "t", options.num_threads,
                        "number of profiling workers");
  description.AddOption("max_reads", "", options.max_reads,
                        "profile at most this many reads; 0 means all");
  description.AddOption("memory", "m", options.memory_bytes,
                        "bounded bytes for concurrent index sorting");

  try {
    description.Parse(argc, argv);
    if (options.read_file.empty()) {
      throw std::logic_error("No packed read library specified");
    }
    if (options.anchor_len == 0 || options.anchor_len > 31) {
      throw std::logic_error("Anchor length must be in [1, 31]");
    }
    if (options.window_len < options.anchor_len) {
      throw std::logic_error("Window length must not be shorter than anchor");
    }
    if (options.window_len - options.anchor_len >= 64u) {
      throw std::logic_error(
          "Window/anchor span must fit the 64-bit candidate-offset mask");
    }
    if (options.local_candidate_output.empty() &&
        (!options.contig_file.empty() || !options.bubble_file.empty()) &&
        (options.kmer_k <= 0 || options.step <= 0 ||
         options.window_len > static_cast<unsigned>(options.kmer_k + 1))) {
      throw std::logic_error(
          "Active profiling needs positive k/step and window <= k+1");
    }
    if ((!options.cache_contig_file.empty() ||
         !options.cache_bubble_file.empty()) &&
        (options.cache_kmer_k <= 0 || options.cache_step <= 0 ||
         options.window_len >
             static_cast<unsigned>(options.cache_kmer_k + 1))) {
      throw std::logic_error(
          "Cached-anchor profiling needs positive cache k/step and window "
          "<= cache k+1");
    }
    if (options.sample_bits < 4 || options.sample_bits > 24) {
      throw std::logic_error("Sample bits must be in [4, 24]");
    }
    if (options.max_reads < 0) {
      throw std::logic_error("Maximum read count must not be negative");
    }
    if (options.num_threads == 0) {
      options.num_threads = omp_get_max_threads();
    }
    if (options.num_threads < 1) {
      throw std::logic_error("Number of workers must be positive");
    }
    if (!options.output_prefix.empty() &&
        (!std::isfinite(options.memory_bytes) || options.memory_bytes <= 0)) {
      throw std::logic_error("Index-build memory budget must be positive");
    }
    if (!options.output_prefix.empty() && !options.index_prefix.empty()) {
      throw std::logic_error(
          "Choose either index construction or index verification");
    }
    if (!options.edge_output_prefix.empty() &&
        ((options.index_prefix.empty() && options.output_prefix.empty()) ||
         options.contig_file.empty() || options.bubble_file.empty() ||
         options.kmer_k <= 0 || options.step <= 0 ||
         options.window_len > AnchorWindow::max_size())) {
      throw std::logic_error(
          "Edge replay needs an index build/query, contigs, bubbles, current "
          "k/step and a window of at most 64 bases");
    }
    if (!options.local_candidate_output.empty() &&
        (options.index_prefix.empty() || options.contig_file.empty() ||
         !options.bubble_file.empty() ||
         !options.edge_output_prefix.empty() ||
         options.local_endpoint_range < 0 || options.local_seed_len <= 0 ||
         options.local_seed_len > 32 || options.local_sparsity <= 0 ||
         options.local_min_contig_len < 0)) {
      throw std::logic_error(
          "Local candidates need an index and contigs plus valid local "
          "mapper parameters; bubble/edge replay cannot be combined");
    }
    omp_set_num_threads(options.num_threads);
    if (!options.index_prefix.empty()) {
      if (!options.contig_file.empty() || !options.bubble_file.empty()) {
        ProfileIndexQuery(options);
      } else {
        VerifyReadIndex(options);
      }
    } else if (options.output_prefix.empty()) {
      RunProfile(options);
    } else {
      RunIndexBuild(options);
    }
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    std::cerr << "Usage: read-index -r reads.lib.bin -a 21 -w 40 -t 32\n";
    std::cerr << description << std::endl;
    return 1;
  }
  return 0;
}
