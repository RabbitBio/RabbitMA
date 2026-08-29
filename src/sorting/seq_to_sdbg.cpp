/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include "seq_to_sdbg.h"

#include <omp.h>
#include <array>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fcntl.h>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "sequence/copy_substr.h"
#include "sequence/io/async_sequence_reader.h"
#include "sequence/io/edge/edge_reader.h"
#include "sequence/kmer.h"
#include "kmlib/kmsort.h"
#include "utils/mutex.h"
#include "utils/startup_affinity.h"
#include "utils/utils.h"

const unsigned SeqToSdbg::kMaxLookUpPrefixLength;

namespace {

/**
 * @brief encode seq_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t seq_id, int offset, int strand,
                                const SeqPackage &p) {
  return ((p.GetSeqView(seq_id).full_offset_in_pkg() + offset) << 1) | strand;
}

struct DecodedSeqLocator {
  DecodedSeqLocator(SeqPackage::SeqView seq_view, unsigned offset,
                    unsigned strand)
      : seq_view(seq_view), offset(offset), strand(strand) {}

  SeqPackage::SeqView seq_view;
  unsigned offset;
  unsigned strand;
};

template <bool PackedOffsets>
struct SeqLocatorDecoder;

template <>
struct SeqLocatorDecoder<false> {
  static inline __attribute__((always_inline)) DecodedSeqLocator Decode(
      uint64_t locator, unsigned, uint64_t, const SeqPackage &seq_pkg) {
    auto seq_view = seq_pkg.GetSeqViewByOffset(locator >> 1u);
    return DecodedSeqLocator(
        seq_view,
        static_cast<unsigned>((locator >> 1u) -
                              seq_view.full_offset_in_pkg()),
        locator & 1u);
  }
};

template <>
struct SeqLocatorDecoder<true> {
  static inline __attribute__((always_inline)) DecodedSeqLocator Decode(
      uint64_t locator, unsigned seq_locator_shift, uint64_t seq_offset_mask,
      const SeqPackage &seq_pkg) {
    return DecodedSeqLocator(
        seq_pkg.GetSeqView(locator >> seq_locator_shift),
        static_cast<unsigned>((locator >> 1u) & seq_offset_mask),
        locator & 1u);
  }
};

inline bool IsDiffKMinusOneMer(uint32_t *item1, uint32_t *item2,
                               int64_t spacing, int k) {
  // mask extra bits
  int chars_in_last_word = (k - 1) % kCharsPerEdgeWord;
  int num_full_words = (k - 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    uint32_t w1 = item1[num_full_words * spacing];
    uint32_t w2 = item2[num_full_words * spacing];

    if ((w1 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar) !=
        (w2 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar)) {
      return true;
    }
  }

  for (int i = num_full_words - 1; i >= 0; --i) {
    if (item1[i * spacing] != item2[i * spacing]) {
      return true;
    }
  }

  return false;
}

inline int Extract_a(uint32_t *item, int num_words, int64_t spacing,
                     unsigned k) {
  int non_dollar = (item[(num_words - 1) * spacing] >>
                    (SeqToSdbg::kBWTCharNumBits + kBitsPerMul)) &
                   1;

  if (non_dollar) {
    unsigned which_word = (k - 1) / kCharsPerEdgeWord;
    unsigned word_index = (k - 1) % kCharsPerEdgeWord;
    return (item[which_word * spacing] >>
            (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) &
           kEdgeCharMask;
  } else {
    return SeqToSdbg::kSentinelValue;
  }
}

inline int Extract_b(uint32_t *item, int num_words, int64_t spacing) {
  return (item[(num_words - 1) * spacing] >> kBitsPerMul) &
         ((1 << SeqToSdbg::kBWTCharNumBits) - 1);
}

inline int ExtractCounting(uint32_t *item, int num_words, int64_t spacing) {
  return item[(num_words - 1) * spacing] & kMaxMul;
}

/**
 * @brief build lkt for faster binary search for mercy
 */
struct MercyLookupTable {
  MercyLookupTable(uint64_t num_edges, unsigned query_length) {
    prefix_length = 1;
    const unsigned max_prefix_length =
        std::min(SeqToSdbg::kMaxLookUpPrefixLength, query_length);
    while (prefix_length < max_prefix_length &&
           (uint64_t{1} << (prefix_length * 2u)) < num_edges) {
      ++prefix_length;
    }
    shift = 32 - prefix_length * 2;
    const size_t num_starts =
        (size_t{1} << (prefix_length * 2u)) + 1;
    compact_starts =
        num_edges <= std::numeric_limits<uint32_t>::max() &&
        std::getenv("MEGAHIT_DISABLE_COMPACT_MERCY_STARTS") == nullptr;
    if (compact_starts) {
      starts32.reset(new uint32_t[num_starts]);
    } else {
      starts64.reset(new uint64_t[num_starts]);
    }
    starts_size = num_starts;
  }

  uint32_t Prefix(const uint32_t *kmer) const { return kmer[0] >> shift; }

  inline void Bounds(uint32_t prefix, int64_t *left,
                     int64_t *right) const {
    if (__builtin_expect(compact_starts, true)) {
      *left = starts32[prefix];
      *right = static_cast<int64_t>(starts32[prefix + 1]) - 1;
    } else {
      *left = static_cast<int64_t>(starts64[prefix]);
      *right = static_cast<int64_t>(starts64[prefix + 1]) - 1;
    }
  }

  size_t StartsSize() const {
    return starts_size;
  }

  size_t StartsBytes() const {
    return starts_size * (compact_starts ? sizeof(uint32_t)
                                         : sizeof(uint64_t));
  }

  unsigned prefix_length{};
  unsigned shift{};
  bool compact_starts{};
  size_t starts_size{};
  std::unique_ptr<uint32_t[]> starts32;
  std::unique_ptr<uint64_t[]> starts64;
  unsigned suffix_length{};
  unsigned suffix_words{};
  std::vector<uint32_t> suffixes;
};

// Exact, canonical-order lookup for unordered count output.  For the initial
// graph (the only stage that uses mercy reads), practical k-min values fit in
// at most 40 bases.  The first eight bases are implicit in a 16-bit bucket and
// the remaining bases occupy one uint64_t.  This is both smaller than keeping
// the raw edge+coverage records and independent of their physical file order.
struct CompactMercyLookup {
  enum : unsigned {
    kImplicitPrefixBases = 8u,
    kCompactBuckets = 1u << 16u
  };

  CompactMercyLookup(uint64_t num_edges, unsigned edge_length,
                     unsigned query_length)
      : edge_length(edge_length),
        suffix_bases(edge_length - kImplicitPrefixBases),
        num_edges(num_edges),
        bucket_offsets(kCompactBuckets + 1u, 0) {
    assert(edge_length > kImplicitPrefixBases && edge_length <= 40u);
    // Suffixes are sorted independently inside the implicit 8-base buckets.
    // A directory range shorter than that prefix spans several independently
    // sorted buckets, which is not itself sorted by the stored suffix and
    // therefore cannot be binary-searched.  Keep the lookup prefix at least
    // as long as the physical bucket prefix; larger prefixes may still split
    // a bucket into smaller sorted ranges.
    assert(query_length >= kImplicitPrefixBases);
    prefix_length = kImplicitPrefixBases;
    // Endpoint verification issues hundreds of millions of exact probes on
    // large inputs.  Size this compact-only directory by an algorithmic
    // occupancy target instead of stopping at the legacy 12-base table: no
    // more than about eight suffixes per directory range, capped at 14 bases
    // (1 GiB of uint32 starts at the cap).  The choice depends only on index
    // cardinality and remains identical across machines.
    const unsigned max_prefix_length = std::min(14u, query_length);
    constexpr uint64_t kTargetRangeOccupancy = 8u;
    while (prefix_length < max_prefix_length &&
           (uint64_t{1} << (prefix_length * 2u)) *
                   kTargetRangeOccupancy <
               num_edges) {
      ++prefix_length;
    }
    prefix_shift = 32u - prefix_length * 2u;
    starts_size = (size_t{1} << (prefix_length * 2u)) + 1u;
    compact_starts = num_edges <= std::numeric_limits<uint32_t>::max();
    if (compact_starts) {
      starts32.reset(new uint32_t[starts_size]);
    } else {
      starts64.reset(new uint64_t[starts_size]);
    }
    suffixes.reset(new (std::nothrow) uint64_t[num_edges]);
    if (!suffixes && num_edges != 0) {
      xfatal("Cannot allocate compact mercy edge index\n");
    }
  }

  static uint64_t ExtractSuffix(const uint32_t *packed,
                                unsigned length) {
    assert(length > kImplicitPrefixBases && length <= 40u);
    uint64_t suffix = 0;
    unsigned position = kImplicitPrefixBases;
    while (position < length) {
      const unsigned word = position / kCharsPerEdgeWord;
      const unsigned within = position % kCharsPerEdgeWord;
      const unsigned take = std::min(
          length - position, kCharsPerEdgeWord - within);
      const unsigned shift =
          (kCharsPerEdgeWord - within - take) * kBitsPerEdgeChar;
      const unsigned bits = take * kBitsPerEdgeChar;
      const uint64_t mask =
          bits == 32u ? std::numeric_limits<uint32_t>::max()
                      : (uint64_t{1} << bits) - 1u;
      suffix = (suffix << bits) | ((packed[word] >> shift) & mask);
      position += take;
    }
    return suffix;
  }

  uint32_t Prefix(const uint32_t *packed) const {
    return packed[0] >> prefix_shift;
  }

  void Bounds(uint32_t prefix, int64_t *left, int64_t *right) const {
    if (compact_starts) {
      *left = starts32[prefix];
      *right = static_cast<int64_t>(starts32[prefix + 1u]) - 1;
    } else {
      *left = static_cast<int64_t>(starts64[prefix]);
      *right = static_cast<int64_t>(starts64[prefix + 1u]) - 1;
    }
  }

  template <class KmerType>
  int64_t Search(const KmerType &kmer, unsigned query_length) const {
    assert(query_length <= edge_length && query_length + 1u >= edge_length);
    const uint32_t prefix = Prefix(kmer.data());
    int64_t left;
    int64_t right;
    Bounds(prefix, &left, &right);
    const uint64_t query_suffix =
        ExtractSuffix(kmer.data(), query_length);
    const unsigned trim = (edge_length - query_length) * kBitsPerEdgeChar;
    while (left <= right) {
      const int64_t mid = (left + right) / 2;
      const uint64_t edge_suffix = suffixes[mid] >> trim;
      if (query_suffix > edge_suffix) {
        left = mid + 1;
      } else if (query_suffix < edge_suffix) {
        right = mid - 1;
      } else {
        return mid;
      }
    }
    return -1;
  }

  unsigned LastBase(uint64_t edge_id) const {
    assert(edge_id < num_edges);
    return static_cast<unsigned>(suffixes[edge_id] & kEdgeCharMask);
  }

  size_t StartsBytes() const {
    return starts_size * (compact_starts ? sizeof(uint32_t)
                                         : sizeof(uint64_t));
  }

  unsigned edge_length;
  unsigned suffix_bases;
  uint64_t num_edges;
  unsigned prefix_length{};
  unsigned prefix_shift{};
  bool compact_starts{};
  size_t starts_size{};
  std::unique_ptr<uint32_t[]> starts32;
  std::unique_ptr<uint64_t[]> starts64;
  std::unique_ptr<uint64_t[]> suffixes;
  std::vector<uint64_t> bucket_offsets;
};

template <class OffsetType>
void InitLookupTableSerialFor(MercyLookupTable *lookup_table, SeqPackage &p,
                              OffsetType *starts, size_t starts_size) {

  if (p.seq_count() == 0) {
    std::fill(starts, starts + starts_size, 0);
    return;
  }

  Kmer<1, uint32_t> kmer;
  uint32_t previous_prefix = 0;
  uint32_t next_prefix_to_fill = 0;
  for (int64_t i = 0, num_edges = p.seq_count(); i < num_edges; ++i) {
    auto ptr_and_offset = p.GetSeqView(i).raw_address();
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                     lookup_table->prefix_length);
    uint32_t prefix = lookup_table->Prefix(kmer.data());

    if (i == 0 || prefix != previous_prefix) {
      assert(i == 0 || prefix > previous_prefix);
      std::fill(starts + next_prefix_to_fill, starts + prefix + 1, i);
      next_prefix_to_fill = prefix + 1;
      previous_prefix = prefix;
    }
  }

  std::fill(starts + next_prefix_to_fill, starts + starts_size,
            p.seq_count());
}

template <class OffsetType>
void InitLookupTableParallelFor(MercyLookupTable *lookup_table, SeqPackage &p,
                                OffsetType *starts, size_t starts_size,
                                int num_threads) {
  const OffsetType sentinel = std::numeric_limits<OffsetType>::max();
  const uint64_t num_edges = p.seq_count();

  // Leave the allocation untouched until the worker team is active.  The
  // parallel initialization both removes the serial O(E) setup bottleneck and
  // first-touches the prefix table across the NUMA nodes running the team.
#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int64_t i = 0; i < static_cast<int64_t>(starts_size); ++i) {
    starts[i] = sentinel;
  }

  if (num_edges != 0) {
#pragma omp parallel num_threads(num_threads)
    {
      const uint64_t tid = omp_get_thread_num();
      const uint64_t team_size = omp_get_num_threads();
      const uint64_t edge_begin = num_edges * tid / team_size;
      const uint64_t edge_end = num_edges * (tid + 1) / team_size;
      Kmer<1, uint32_t> kmer;
      uint32_t previous_prefix = std::numeric_limits<uint32_t>::max();

      if (edge_begin != 0) {
        auto raw = p.GetSeqView(edge_begin - 1).raw_address();
        kmer.InitFromPtr(raw.first, raw.second,
                         lookup_table->prefix_length);
        previous_prefix = lookup_table->Prefix(kmer.data());
      }

      for (uint64_t edge_id = edge_begin; edge_id < edge_end; ++edge_id) {
        auto raw = p.GetSeqView(edge_id).raw_address();
        kmer.InitFromPtr(raw.first, raw.second,
                         lookup_table->prefix_length);
        const uint32_t prefix = lookup_table->Prefix(kmer.data());
        if (prefix != previous_prefix) {
          // Prefix runs are contiguous, so exactly one worker observes each
          // global transition (including a transition on a chunk boundary).
          starts[prefix] = static_cast<OffsetType>(edge_id);
          previous_prefix = prefix;
        }
      }
    }
  }

  starts[starts_size - 1] = static_cast<OffsetType>(num_edges);
  OffsetType next_start = static_cast<OffsetType>(num_edges);
  // Only the transition entries need synchronization.  Back-filling 16M
  // compact offsets is a short streaming pass and avoids a parallel prefix
  // algorithm (and its barriers) on small inputs.
  for (size_t i = starts_size - 1; i-- > 0;) {
    if (starts[i] == sentinel) {
      starts[i] = next_start;
    } else {
      next_start = starts[i];
    }
  }
}

void InitLookupTable(MercyLookupTable *lookup_table, SeqPackage &p,
                     int num_threads) {
  const bool parallel =
      num_threads > 1 &&
      std::getenv("MEGAHIT_DISABLE_PARALLEL_MERCY_LOOKUP_INIT") == nullptr;
  if (lookup_table->compact_starts) {
    if (parallel) {
      InitLookupTableParallelFor(lookup_table, p,
                                 lookup_table->starts32.get(),
                                 lookup_table->StartsSize(), num_threads);
    } else {
      InitLookupTableSerialFor(lookup_table, p,
                               lookup_table->starts32.get(),
                               lookup_table->StartsSize());
    }
  } else {
    if (parallel) {
      InitLookupTableParallelFor(lookup_table, p,
                                 lookup_table->starts64.get(),
                                 lookup_table->StartsSize(), num_threads);
    } else {
      InitLookupTableSerialFor(lookup_table, p,
                               lookup_table->starts64.get(),
                               lookup_table->StartsSize());
    }
  }
}

template <class KmerType>
void InitMercySuffixIndex(MercyLookupTable *lookup_table, SeqPackage &p,
                          unsigned edge_length, int num_threads,
                          int64_t host_mem) {
  if (std::getenv("MEGAHIT_DISABLE_MERCY_SUFFIX") != nullptr) {
    return;
  }
  // The compact suffix array improves locality with one worker per physical
  // core, but when SMT siblings share every core its extra memory stream can
  // compete with the already bandwidth-bound mercy scan.  Discover that
  // boundary from the CPUs actually available to the process; do not encode
  // a socket count, NUMA layout, or core count.
  const unsigned physical_cores =
      GetNumaTopology().total_physical_core_count();
  if (physical_cores != 0 &&
      static_cast<unsigned>(num_threads) > physical_cores &&
      std::getenv("MEGAHIT_FORCE_MERCY_SUFFIX") == nullptr) {
    xinfo(
        "Mercy suffix index skipped: {} workers exceed {} available "
        "physical cores\n",
        num_threads, physical_cores);
    return;
  }
  if (lookup_table->prefix_length >= edge_length) {
    return;
  }

  lookup_table->suffix_length = edge_length - lookup_table->prefix_length;
  lookup_table->suffix_words =
      DivCeiling(lookup_table->suffix_length, kCharsPerEdgeWord);
  const uint64_t num_words =
      p.seq_count() * static_cast<uint64_t>(lookup_table->suffix_words);
  const uint64_t num_bytes = num_words * sizeof(uint32_t);

  // The index duplicates only the suffix needed by mercy queries.  Keep a
  // conservative cap for unusually large k/data sets and fall back to the
  // original binary search if the extra sequential array would consume too
  // much of the user's memory budget.
  if (num_bytes > static_cast<uint64_t>(host_mem / 8)) {
    xinfo("Mercy suffix index skipped: {} bytes exceeds budget\n", num_bytes);
    lookup_table->suffix_length = 0;
    lookup_table->suffix_words = 0;
    return;
  }

  lookup_table->suffixes.resize(num_words);
  omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(static)
  for (uint64_t edge_id = 0; edge_id < p.seq_count(); ++edge_id) {
    auto seq_view = p.GetSeqView(edge_id);
    auto raw = seq_view.raw_address(lookup_table->prefix_length);
    KmerType suffix(raw.first, raw.second, lookup_table->suffix_length);
    std::memcpy(lookup_table->suffixes.data() +
                    edge_id * lookup_table->suffix_words,
                suffix.data(), lookup_table->suffix_words * sizeof(uint32_t));
  }
}

/**
 * @brief search mercy kmer
 */
template <class KmerType>
int64_t BinarySearchKmer(KmerType &kmer,
                         const MercyLookupTable &lookup_table, SeqPackage &p,
                         int kmer_size) {
  uint32_t prefix = lookup_table.Prefix(kmer.data());
  int64_t l;
  int64_t r;
  lookup_table.Bounds(prefix, &l, &r);
  KmerType mid_kmer;
  KmerType query_suffix;
  const bool use_suffix =
      !lookup_table.suffixes.empty() &&
      static_cast<unsigned>(kmer_size) > lookup_table.prefix_length;
  const unsigned query_suffix_length =
      use_suffix ? kmer_size - lookup_table.prefix_length : 0;
  if (use_suffix) {
    query_suffix.InitFromPtr(kmer.data(), lookup_table.prefix_length,
                             query_suffix_length);
  }

  while (l <= r) {
    int64_t mid = (l + r) / 2;
    int cmp = 0;
    if (use_suffix) {
      const uint32_t *edge_suffix =
          lookup_table.suffixes.data() + mid * lookup_table.suffix_words;
      const unsigned used_words =
          DivCeiling(query_suffix_length, kCharsPerEdgeWord);
      for (unsigned word = 0; word < used_words; ++word) {
        uint32_t query_value = query_suffix.data()[word];
        uint32_t edge_value = edge_suffix[word];
        if (word + 1 == used_words &&
            query_suffix_length % kCharsPerEdgeWord != 0) {
          const unsigned clean_shift =
              (kCharsPerEdgeWord -
               query_suffix_length % kCharsPerEdgeWord) *
              kBitsPerEdgeChar;
          edge_value = edge_value >> clean_shift << clean_shift;
        }
        if (query_value < edge_value) {
          cmp = -1;
          break;
        }
        if (query_value > edge_value) {
          cmp = 1;
          break;
        }
      }
    } else {
      auto ptr_and_offset = p.GetSeqView(mid).raw_address();
      mid_kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                           kmer_size);
      cmp = kmer.cmp(mid_kmer, kmer_size);
    }

    if (cmp > 0) {
      l = mid + 1;
    } else if (cmp < 0) {
      r = mid - 1;
    } else {
      return mid;
    }
  }

  return -1;
}

}  // namespace

// sorting core functions
SeqToSdbg::~SeqToSdbg() {
  for (int fd : stream_edge_fds_) {
    if (fd >= 0) {
      ::close(fd);
    }
  }
}

bool SeqToSdbg::ConfigureStreamedEdgeInput(
    const EdgeIoMetadata &metadata) {
  if (metadata.num_edges < 0 || metadata.kmer_size != opt_.k ||
      metadata.num_files == 0 ||
      metadata.words_per_edge == 0 || opt_.k < kBucketPrefixLength) {
    return false;
  }
  const uint64_t expected_words = DivCeiling(
      (static_cast<uint64_t>(opt_.k) + 1u) * kBitsPerEdgeChar + kBitsPerMul,
      kBitsPerEdgeWord);
  if (metadata.words_per_edge != expected_words) {
    return false;
  }

  struct FileRange {
    uint64_t begin;
    uint64_t end;
  };
  const uint64_t record_bytes =
      static_cast<uint64_t>(metadata.words_per_edge) * sizeof(uint32_t);
  std::vector<uint64_t> file_records(metadata.num_files, 0);
  uint64_t total_records = 0;
  if (metadata.is_sorted) {
    std::vector<std::vector<FileRange>> ranges(metadata.num_files);
    for (const auto &bucket : metadata.buckets) {
      if (bucket.total_number < 0 || bucket.file_offset < 0 ||
          bucket.file_id < 0 ||
          bucket.file_id >= static_cast<int>(metadata.num_files)) {
        if (bucket.total_number == 0 && bucket.file_id < 0) {
          continue;
        }
        return false;
      }
      if (bucket.total_number == 0) {
        continue;
      }
      const uint64_t begin = static_cast<uint64_t>(bucket.file_offset);
      const uint64_t count = static_cast<uint64_t>(bucket.total_number);
      if (begin > std::numeric_limits<uint64_t>::max() - count ||
          total_records > std::numeric_limits<uint64_t>::max() - count) {
        return false;
      }
      ranges[bucket.file_id].push_back({begin, begin + count});
      total_records += count;
    }
    if (total_records != static_cast<uint64_t>(metadata.num_edges)) {
      return false;
    }
    for (unsigned file_id = 0; file_id < metadata.num_files; ++file_id) {
      auto &file_ranges = ranges[file_id];
      std::sort(file_ranges.begin(), file_ranges.end(),
                [](const FileRange &lhs, const FileRange &rhs) {
                  return lhs.begin < rhs.begin;
                });
      for (const auto &range : file_ranges) {
        // EdgeWriter appends bucket runs to each shard.  Requiring complete,
        // non-overlapping coverage lets the hot loops read every shard as a
        // simple sequential record stream, independent of bucket assignment.
        if (range.begin != file_records[file_id] || range.end < range.begin) {
          return false;
        }
        file_records[file_id] = range.end;
      }
    }
  }

  std::vector<int> fds(metadata.num_files, -1);
  for (unsigned file_id = 0; file_id < metadata.num_files; ++file_id) {
    const std::string path =
        opt_.input_prefix + ".edges." + std::to_string(file_id);
    fds[file_id] = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fds[file_id] < 0) {
      for (int fd : fds) {
        if (fd >= 0) {
          ::close(fd);
        }
      }
      return false;
    }
    (void)::posix_fadvise(fds[file_id], 0, 0, POSIX_FADV_SEQUENTIAL);

    struct stat file_stat {};
    if (::fstat(fds[file_id], &file_stat) != 0 || file_stat.st_size < 0) {
      for (int fd : fds) {
        if (fd >= 0) {
          ::close(fd);
        }
      }
      return false;
    }
    const uint64_t file_bytes = static_cast<uint64_t>(file_stat.st_size);
    if (metadata.is_sorted) {
      if (file_records[file_id] >
              std::numeric_limits<uint64_t>::max() / record_bytes ||
          file_bytes != file_records[file_id] * record_bytes) {
        for (int fd : fds) {
          if (fd >= 0) {
            ::close(fd);
          }
        }
        return false;
      }
    } else {
      if (file_bytes % record_bytes != 0) {
        for (int fd : fds) {
          if (fd >= 0) {
            ::close(fd);
          }
        }
        return false;
      }
      file_records[file_id] = file_bytes / record_bytes;
      if (total_records > std::numeric_limits<uint64_t>::max() -
                              file_records[file_id]) {
        for (int fd : fds) {
          if (fd >= 0) {
            ::close(fd);
          }
        }
        return false;
      }
      total_records += file_records[file_id];
    }
  }
  if (total_records != static_cast<uint64_t>(metadata.num_edges)) {
    for (int fd : fds) {
      if (fd >= 0) {
        ::close(fd);
      }
    }
    return false;
  }

  const uint64_t records_per_chunk = std::max<uint64_t>(
      1u, (uint64_t{2} << 20u) / record_bytes);
  std::vector<StreamEdgeChunk> chunks;
  std::vector<size_t> file_chunk_offsets(metadata.num_files + 1u, 0);
  for (unsigned file_id = 0; file_id < metadata.num_files; ++file_id) {
    file_chunk_offsets[file_id] = chunks.size();
    for (uint64_t first = 0; first < file_records[file_id];) {
      const uint64_t count =
          std::min(records_per_chunk, file_records[file_id] - first);
      chunks.push_back({file_id, first, static_cast<uint32_t>(count)});
      first += count;
    }
  }
  file_chunk_offsets[metadata.num_files] = chunks.size();

  stream_edge_length_ = opt_.k + 1u;
  stream_words_per_edge_ = metadata.words_per_edge;
  stream_num_edges_ = metadata.num_edges;
  stream_edge_chunks_.swap(chunks);
  stream_file_chunk_offsets_.swap(file_chunk_offsets);
  stream_edge_fds_.swap(fds);
  stream_input_edges_ = true;
  stream_input_unordered_ = !metadata.is_sorted;

  EdgeBucketHistogram precomputed;
  if (std::getenv("MEGAHIT_DISABLE_COUNT_SDBG_HISTOGRAM") == nullptr &&
      precomputed.Read(opt_.input_prefix) &&
      precomputed.kmer_size == opt_.k &&
      precomputed.num_files == metadata.num_files &&
      precomputed.num_edges == static_cast<uint64_t>(metadata.num_edges)) {
    bool rows_match_files = true;
    for (unsigned file_id = 0; file_id < metadata.num_files; ++file_id) {
      const uint64_t *row = precomputed.counts.data() +
                            static_cast<size_t>(file_id) * kNumBuckets;
      const uint64_t row_records =
          std::accumulate(row, row + kNumBuckets, uint64_t{0});
      if (row_records != 2u * file_records[file_id]) {
        rows_match_files = false;
        break;
      }
    }
    if (rows_match_files) {
      stream_bucket_histogram_ = std::move(precomputed);
      use_stream_bucket_histogram_ = true;
    }
  }
  xinfo("Streaming {} fixed edge input: {} records in {} sequential chunks; "
        "no resident edge copy\n",
        metadata.is_sorted ? "sorted" : "unordered", stream_num_edges_,
        stream_edge_chunks_.size());
  if (use_stream_bucket_histogram_) {
    xinfo("Loaded count-stage radix histogram; the first edge scan can be "
          "elided when no additional sequences are present\n");
  }
  return true;
}

void SeqToSdbg::RetainMercyEdgesForStreamedInput(
    size_t original_edge_count) {
  if (original_edge_count > seq_pkg_.seq_count() ||
      original_edge_count > multiplicity.size()) {
    xfatal("Invalid original edge count while compacting mercy edges\n");
  }
  const size_t num_mercy_edges = seq_pkg_.seq_count() - original_edge_count;
  const unsigned edge_words =
      DivCeiling(stream_edge_length_, kCharsPerEdgeWord);
  if (num_mercy_edges >
      std::numeric_limits<size_t>::max() / edge_words) {
    xfatal("Mercy edge compaction size overflow\n");
  }

  std::vector<uint32_t> packed(num_mercy_edges * edge_words);
#pragma omp parallel for schedule(static) num_threads(opt_.n_threads)
  for (int64_t i = 0; i < static_cast<int64_t>(num_mercy_edges); ++i) {
    const auto seq_view = seq_pkg_.GetSeqView(original_edge_count + i);
    const auto raw = seq_view.raw_address();
    const unsigned source_words =
        DivCeiling(raw.second + seq_view.length(), kCharsPerEdgeWord);
    CopySubstring(packed.data() + static_cast<size_t>(i) * edge_words,
                  raw.first, raw.second, seq_view.length(), 1, source_words,
                  edge_words);
  }

  std::vector<mul_t> mercy_multiplicity(
      multiplicity.begin() + original_edge_count, multiplicity.end());
  seq_pkg_.ReleaseStorage();
  seq_pkg_.ReserveBases(num_mercy_edges *
                        static_cast<size_t>(stream_edge_length_));
  seq_pkg_.AssignFixedLengthCompactSequences(
      packed.data(), num_mercy_edges, stream_edge_length_, edge_words,
      opt_.n_threads);
  multiplicity.swap(mercy_multiplicity);
  std::vector<mul_t>().swap(mercy_multiplicity);
  std::vector<uint32_t>().swap(packed);
  xinfo("Released resident source edges after mercy lookup; retained {} "
        "new mercy edges\n",
        num_mercy_edges);
}

void SeqToSdbg::ReadStreamedEdgeChunk(
    const StreamEdgeChunk &chunk, std::vector<uint32_t> *records) const {
  assert(stream_input_edges_ && chunk.file_id < stream_edge_fds_.size());
  const uint64_t words =
      static_cast<uint64_t>(chunk.num_records) * stream_words_per_edge_;
  if (words > std::numeric_limits<size_t>::max()) {
    xfatal("Streamed edge chunk is too large\n");
  }
  records->resize(static_cast<size_t>(words));
  char *destination = reinterpret_cast<char *>(records->data());
  size_t bytes_left = static_cast<size_t>(words * sizeof(uint32_t));
  uint64_t byte_offset = chunk.first_record *
                         static_cast<uint64_t>(stream_words_per_edge_) *
                         sizeof(uint32_t);
  while (bytes_left != 0) {
    const ssize_t got = ::pread(stream_edge_fds_[chunk.file_id], destination,
                                bytes_left, static_cast<off_t>(byte_offset));
    if (got < 0 && errno == EINTR) {
      continue;
    }
    if (got <= 0) {
      xfatal("Failed to read streamed edge shard {} at byte {}\n",
             chunk.file_id, byte_offset);
    }
    destination += got;
    byte_offset += static_cast<uint64_t>(got);
    bytes_left -= static_cast<size_t>(got);
  }
}

uint16_t SeqToSdbg::ExtractRawBaseWindow8(const uint32_t *edge,
                                          unsigned base_offset) {
  const unsigned word = base_offset / kCharsPerEdgeWord;
  const unsigned shift =
      (base_offset % kCharsPerEdgeWord) * kBitsPerEdgeChar;
  const uint64_t window =
      (static_cast<uint64_t>(edge[word]) << kBitsPerEdgeWord) |
      edge[word + 1u];
  return static_cast<uint16_t>((window << shift) >> 48u);
}

uint16_t SeqToSdbg::ReverseComplementWindow8(uint16_t bases) {
  uint16_t value = bases;
  value = static_cast<uint16_t>(((value & 0x3333u) << 2u) |
                                ((value & 0xCCCCu) >> 2u));
  value = static_cast<uint16_t>(((value & 0x0F0Fu) << 4u) |
                                ((value & 0xF0F0u) >> 4u));
  value = static_cast<uint16_t>((value << 8u) | (value >> 8u));
  return static_cast<uint16_t>(value ^ 0xFFFFu);
}

unsigned SeqToSdbg::ExtractRawBase(const uint32_t *edge, unsigned offset) {
  return (edge[offset / kCharsPerEdgeWord] >>
          ((kCharsPerEdgeWord - 1u - offset % kCharsPerEdgeWord) *
           kBitsPerEdgeChar)) &
         kEdgeCharMask;
}

void SeqToSdbg::ConfigurePackedSeqOffsets() {
  packed_seq_offsets_ = false;
  seq_locator_shift_ = 0;
  seq_offset_mask_ = 0;

  if (std::getenv("MEGAHIT_DISABLE_PACKED_SEQ_OFFSETS") != nullptr) {
    xinfo("Packed sequence offsets disabled by environment\n");
    return;
  }

  const uint64_t max_seq_offset =
      seq_pkg_.max_length() == 0 ? 0 : seq_pkg_.max_length() - 1;
  unsigned offset_bits = 0;
  for (uint64_t value = max_seq_offset; value != 0; value >>= 1u) {
    ++offset_bits;
  }
  seq_locator_shift_ = offset_bits + 1;  // local offset plus strand
  seq_offset_mask_ =
      offset_bits == 0 ? 0 : (uint64_t{1} << offset_bits) - 1;

  const uint64_t max_seq_id =
      seq_pkg_.seq_count() == 0 ? 0 : seq_pkg_.seq_count() - 1;
  const uint64_t max_locator =
      static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
  if (seq_locator_shift_ >= 63 ||
      max_seq_id > (max_locator >> seq_locator_shift_)) {
    seq_locator_shift_ = 0;
    seq_offset_mask_ = 0;
    xinfo("Packed sequence offsets unavailable: sequence IDs exceed locator "
          "width\n");
    return;
  }

  const uint64_t last_locator =
      (max_seq_id << seq_locator_shift_) | (max_seq_offset << 1u) | 1u;
  if (last_locator > max_locator) {
    seq_locator_shift_ = 0;
    seq_offset_mask_ = 0;
    xinfo("Packed sequence offsets unavailable: local offsets exceed locator "
          "width\n");
    return;
  }

  assert(((max_seq_offset << 1u) | 1u) <
         (uint64_t{1} << seq_locator_shift_));

  // A sequence-id/local-offset locator is cheap to decode, but its fixed
  // power-of-two slot is wasteful when a few long contigs are mixed with
  // millions of much shorter edge records.  The expanded coordinate domain
  // then creates large differential-offset escape tables and can disable the
  // compact cursor even though the packed bases themselves are much smaller.
  // Select from the input layout, not from a data-size or machine-specific
  // threshold: retain the fast form unless it consumes at least two extra
  // locator bits (a fourfold coordinate expansion) over dense base offsets.
  const uint64_t max_signed_locator =
      static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
  const uint64_t total_bases = seq_pkg_.base_count();
  const bool dense_available =
      total_bases == 0 || total_bases - 1u <= (max_signed_locator >> 1u);
  const uint64_t dense_last_locator =
      total_bases == 0 ? 0 : ((total_bases - 1u) << 1u) | 1u;
  if (dense_available && dense_last_locator != 0 &&
      last_locator / dense_last_locator >= 4u) {
    seq_locator_shift_ = 0;
    seq_offset_mask_ = 0;
    xinfo(
        "Using dense global-base offsets: packed locator domain {} is {}x "
        "the dense domain {}\n",
        last_locator, last_locator / dense_last_locator,
        dense_last_locator);
    return;
  }

  packed_seq_offsets_ = true;
  xinfo("Using packed sequence offsets: {} locator bits for sequences of "
        "length <= {}\n",
        seq_locator_shift_, seq_pkg_.max_length());
}

int64_t SeqToSdbg::BoundedTransientWorkspaceLimit(
    uint64_t retained_bytes, bool report) const {
  if (std::getenv("MEGAHIT_FORCE_DIRECT_SEQ_ITEMS") != nullptr ||
      opt_.mem_flag != 1) {
    return std::numeric_limits<int64_t>::max();
  }
  if (opt_.host_mem <= 0) {
    return 0;
  }

  long double target_gib = 20.0L;
  if (const char *value = std::getenv("MEGAHIT_SEQ2SDBG_RSS_TARGET_GIB")) {
    char *end = nullptr;
    const long double parsed = std::strtold(value, &end);
    if (end != value && *end == '\0' && std::isfinite(parsed) &&
        parsed > 0.0L) {
      target_gib = parsed;
    }
  }
  const long double requested =
      target_gib * static_cast<long double>(uint64_t{1} << 30u);
  const uint64_t rss_target = static_cast<uint64_t>(std::min<long double>(
      requested, static_cast<long double>(std::max<double>(0, opt_.host_mem))));

  // The sorting engine keeps three bucket tables per worker while filling
  // either direct records or compact locators.  Account for them here so the
  // same process-wide contract applies to every k, rather than only to the
  // first streamed count input.
  const unsigned __int128 bucket_tables =
      static_cast<uint64_t>(std::max(1, opt_.n_threads)) * kNumBuckets *
      sizeof(int64_t) * 3u;
  const unsigned __int128 tracked =
      static_cast<unsigned __int128>(retained_bytes) + bucket_tables;
  const uint64_t untracked_headroom = std::max<uint64_t>(
      uint64_t{1} << 30u, rss_target / 10u);
  if (tracked >= rss_target ||
      untracked_headroom >= rss_target - static_cast<uint64_t>(tracked)) {
    return 0;
  }
  const uint64_t workspace =
      rss_target - static_cast<uint64_t>(tracked) - untracked_headroom;
  if (report) {
    xinfo("Bounded seq2sdbg workspace: target {.3} GiB, retained {.3} "
          "GiB, safety {.3} GiB, workspace {.3} GiB\n",
          static_cast<double>(rss_target) / (uint64_t{1} << 30u),
          static_cast<double>(tracked) / (uint64_t{1} << 30u),
          static_cast<double>(untracked_headroom) / (uint64_t{1} << 30u),
          static_cast<double>(workspace) / (uint64_t{1} << 30u));
  }
  return workspace >
                 static_cast<uint64_t>(std::numeric_limits<int64_t>::max())
             ? std::numeric_limits<int64_t>::max()
             : static_cast<int64_t>(workspace);
}

uint64_t SeqToSdbg::CurrentRetainedBytes() const {
  const unsigned __int128 retained =
      static_cast<uint64_t>(seq_pkg_.size_in_byte()) +
      static_cast<unsigned __int128>(multiplicity.capacity()) *
          sizeof(mul_t) +
      static_cast<unsigned __int128>(endpoint_tip_items_.capacity()) *
          sizeof(uint32_t) +
      static_cast<unsigned __int128>(
          stream_bucket_histogram_.counts.capacity()) *
          sizeof(uint64_t);
  return retained > std::numeric_limits<uint64_t>::max()
             ? std::numeric_limits<uint64_t>::max()
             : static_cast<uint64_t>(retained);
}

int64_t SeqToSdbg::AutomaticDirectWorkspaceLimit(
    uint64_t retained_bytes) const {
  return BoundedTransientWorkspaceLimit(retained_bytes, false);
}

int64_t SeqToSdbg::EncodeSeqOffset(int64_t seq_id, unsigned offset,
                                   unsigned strand) const {
  if (packed_seq_offsets_) {
    assert(strand <= 1);
    assert(offset <= seq_offset_mask_);
    return static_cast<int64_t>(
        (static_cast<uint64_t>(seq_id) << seq_locator_shift_) |
        (static_cast<uint64_t>(offset) << 1u) | strand);
  }
  return EncodeEdgeOffset(seq_id, offset, strand, seq_pkg_);
}

int64_t SeqToSdbg::Lv0EncodeDiffBase(int64_t seq_id) {
  const int64_t num_stream_chunks =
      static_cast<int64_t>(stream_edge_chunks_.size());
  if (seq_id < num_stream_chunks) {
    // Direct fillers do not differential-code locators, but BaseEngine still
    // requests a monotone partition base.  Chunk IDs are sufficient.
    return seq_id;
  }
  seq_id -= num_stream_chunks;
  assert(seq_id < static_cast<int64_t>(seq_pkg_.seq_count()));
  return EncodeSeqOffset(seq_id, 0, 0);
}

bool SeqToSdbg::Lv1CanUseCompactCursor(int64_t seq_from,
                                       int64_t seq_to) const {
  if (seq_from >= seq_to) {
    return true;
  }
  if (stream_input_edges_) {
    return false;
  }
  const auto last_seq = seq_pkg_.GetSeqView(seq_to - 1);
  const unsigned last_offset =
      last_seq.length() == 0 ? 0 : last_seq.length() - 1;
  const uint64_t begin =
      static_cast<uint64_t>(EncodeSeqOffset(seq_from, 0, 0));
  const uint64_t end = static_cast<uint64_t>(
      EncodeSeqOffset(seq_to - 1, last_offset, 1));
  return end >= begin &&
         end - begin <= std::numeric_limits<uint32_t>::max();
}

bool SeqToSdbg::Lv0BuildBalancedRanges(
    std::vector<std::pair<int64_t, int64_t>> *ranges) const {
  const unsigned num_threads =
      static_cast<unsigned>(std::max(1, opt_.n_threads));
  const int64_t num_stream_chunks =
      static_cast<int64_t>(stream_edge_chunks_.size());
  const int64_t num_sequences = num_stream_chunks +
      static_cast<int64_t>(seq_pkg_.seq_count());
  if (use_stream_bucket_histogram_ && seq_pkg_.seq_count() == 0 &&
      num_threads == stream_bucket_histogram_.num_files &&
      stream_file_chunk_offsets_.size() == num_threads + 1u) {
    ranges->resize(num_threads);
    for (unsigned tid = 0; tid < num_threads; ++tid) {
      (*ranges)[tid] = {
          static_cast<int64_t>(stream_file_chunk_offsets_[tid]),
          static_cast<int64_t>(stream_file_chunk_offsets_[tid + 1u])};
    }
    xinfo("Count-aligned edge partitions: reused {} per-shard radix "
          "histograms\n",
          num_threads);
    return true;
  }
  if (num_sequences == 0 || num_threads <= 1 ||
      (!stream_input_edges_ &&
       seq_pkg_.max_length() <= 4u * (opt_.k + 1u))) {
    return false;
  }

  auto sequence_work = [this, num_stream_chunks](int64_t seq_id) -> uint64_t {
    if (seq_id < num_stream_chunks) {
      return 2u * stream_edge_chunks_[seq_id].num_records;
    }
    const uint64_t length =
        seq_pkg_.GetSeqView(seq_id - num_stream_chunks).length();
    if (length < opt_.k + 1u) {
      return 0u;
    }
    const uint64_t records_per_orientation =
        real_only_records_ ? length - opt_.k : length - opt_.k + 2u;
    return 2u * records_per_orientation;
  };

  uint64_t total_work = 0;
  for (int64_t seq_id = 0; seq_id < num_sequences; ++seq_id) {
    total_work += sequence_work(seq_id);
  }
  if (total_work == 0) {
    return false;
  }

  ranges->assign(num_threads, {num_sequences, num_sequences});
  std::vector<uint64_t> range_work(num_threads, 0);
  int64_t range_begin = 0;
  uint64_t accumulated = 0;
  unsigned range_id = 0;
  for (int64_t seq_id = 0;
       seq_id < num_sequences && range_id + 1 < num_threads; ++seq_id) {
    accumulated += sequence_work(seq_id);
    const uint64_t target = static_cast<uint64_t>(
        (static_cast<unsigned __int128>(total_work) * (range_id + 1u) +
         num_threads - 1u) /
        num_threads);
    if (accumulated >= target) {
      (*ranges)[range_id] = {range_begin, seq_id + 1};
      range_work[range_id] = accumulated -
          static_cast<uint64_t>(
              (static_cast<unsigned __int128>(total_work) * range_id) /
              num_threads);
      range_begin = seq_id + 1;
      ++range_id;
    }
  }
  if (range_id < num_threads) {
    (*ranges)[range_id] = {range_begin, num_sequences};
  }
  for (unsigned i = range_id + 1; i < num_threads; ++i) {
    (*ranges)[i] = {num_sequences, num_sequences};
  }

  uint64_t max_range_work = 0;
  uint64_t prefix_work = 0;
  for (unsigned i = 0; i < num_threads; ++i) {
    uint64_t work = 0;
    for (int64_t seq_id = (*ranges)[i].first;
         seq_id < (*ranges)[i].second; ++seq_id) {
      work += sequence_work(seq_id);
    }
    range_work[i] = work;
    prefix_work += work;
    max_range_work = std::max(max_range_work, work);
  }
  assert(prefix_work == total_work);
  xinfo("Work-balanced sequence partitions: {} records, max/avg: {.3}\n",
        total_work,
        static_cast<double>(max_range_work) * num_threads / total_work);
  return true;
}

template <unsigned NWords>
void SeqToSdbg::GenMercyEdgesForK() {
  using MercyKmer = Kmer<NWords, uint32_t>;
  if (opt_.k + 1 > MercyKmer::max_size()) {
    xfatal("Invalid mercy edge length: {}\n", opt_.k + 1);
  }

  // The candidate reader is asynchronous and sleeps while the compute team
  // owns a batch.  Let the full requested team search that batch; retaining
  // the legacy one-worker reservation is useful for architecture A/B tests.
  int num_threads = std::max(
      1, opt_.n_threads -
             (std::getenv("MEGAHIT_RESERVE_MERCY_IO_THREAD") != nullptr));
  omp_set_num_threads(num_threads);

  SimpleTimer phase_timer;
  phase_timer.start();
  std::unique_ptr<MercyLookupTable> edge_lookup;
  std::unique_ptr<CompactMercyLookup> compact_lookup;
  if (stream_compact_mercy_index_) {
    compact_lookup.reset(new CompactMercyLookup(
        static_cast<uint64_t>(stream_num_edges_), stream_edge_length_,
        opt_.k));
    const size_t matrix_items =
        static_cast<size_t>(num_threads) * kNumBuckets;
    std::vector<uint64_t> thread_bucket_counts(matrix_items, 0);
    DiscardMemoryPages(thread_bucket_counts.data(),
                       thread_bucket_counts.size() * sizeof(uint64_t));
    AdviseHugePages(thread_bucket_counts.data(),
                    thread_bucket_counts.size() * sizeof(uint64_t));

#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      uint64_t *local_counts =
          thread_bucket_counts.data() +
          static_cast<size_t>(tid) * kNumBuckets;
      std::fill(local_counts, local_counts + kNumBuckets, uint64_t{0});
      std::vector<uint32_t> records;
#pragma omp for schedule(static)
      for (int64_t chunk_id = 0;
           chunk_id < static_cast<int64_t>(stream_edge_chunks_.size());
           ++chunk_id) {
        const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
        ReadStreamedEdgeChunk(chunk, &records);
        for (uint32_t i = 0; i < chunk.num_records; ++i) {
          const uint32_t *edge =
              records.data() + static_cast<size_t>(i) *
                                   stream_words_per_edge_;
          ++local_counts[edge[0] >> 16u];
        }
      }
    }

    for (unsigned prefix = 0; prefix < kNumBuckets; ++prefix) {
      uint64_t count = 0;
      for (int tid = 0; tid < num_threads; ++tid) {
        count += thread_bucket_counts[
            static_cast<size_t>(tid) * kNumBuckets + prefix];
      }
      compact_lookup->bucket_offsets[prefix + 1u] =
          compact_lookup->bucket_offsets[prefix] + count;
    }
    if (compact_lookup->bucket_offsets.back() !=
        static_cast<uint64_t>(stream_num_edges_)) {
      xfatal("Compact mercy index edge count mismatch\n");
    }
    const size_t suffix_bytes =
        static_cast<size_t>(stream_num_edges_) * sizeof(uint64_t);
    InterleaveMemoryPages(compact_lookup->suffixes.get(), suffix_bytes);
    AdviseHugePages(compact_lookup->suffixes.get(), suffix_bytes);

    for (unsigned prefix = 0; prefix < kNumBuckets; ++prefix) {
      uint64_t cursor = compact_lookup->bucket_offsets[prefix];
      for (int tid = 0; tid < num_threads; ++tid) {
        uint64_t &entry = thread_bucket_counts[
            static_cast<size_t>(tid) * kNumBuckets + prefix];
        const uint64_t count = entry;
        entry = cursor;
        cursor += count;
      }
      assert(cursor == compact_lookup->bucket_offsets[prefix + 1u]);
    }

#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      uint64_t *local_cursors =
          thread_bucket_counts.data() +
          static_cast<size_t>(tid) * kNumBuckets;
      std::vector<uint32_t> records;
#pragma omp for schedule(static)
      for (int64_t chunk_id = 0;
           chunk_id < static_cast<int64_t>(stream_edge_chunks_.size());
           ++chunk_id) {
        const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
        ReadStreamedEdgeChunk(chunk, &records);
        for (uint32_t i = 0; i < chunk.num_records; ++i) {
          const uint32_t *edge =
              records.data() + static_cast<size_t>(i) *
                                   stream_words_per_edge_;
          const unsigned prefix = edge[0] >> 16u;
          compact_lookup->suffixes[local_cursors[prefix]++] =
              CompactMercyLookup::ExtractSuffix(edge,
                                                 stream_edge_length_);
        }
      }
    }
    std::vector<uint64_t>().swap(thread_bucket_counts);

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int prefix = 0; prefix < static_cast<int>(kNumBuckets); ++prefix) {
      uint64_t *begin = compact_lookup->suffixes.get() +
                        compact_lookup->bucket_offsets[prefix];
      uint64_t *end = compact_lookup->suffixes.get() +
                      compact_lookup->bucket_offsets[prefix + 1u];
      if (end - begin > 1) kmlib::kmsort(begin, end);
    }

    if (compact_lookup->prefix_length <= kBucketPrefixLength) {
      const unsigned prefix_bits = compact_lookup->prefix_length * 2u;
      const size_t logical_prefixes = size_t{1} << prefix_bits;
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int64_t prefix = 0;
           prefix <= static_cast<int64_t>(logical_prefixes); ++prefix) {
        const unsigned bucket = static_cast<unsigned>(prefix) <<
                                (16u - prefix_bits);
        const uint64_t offset = compact_lookup->bucket_offsets[bucket];
        if (compact_lookup->compact_starts) {
          compact_lookup->starts32[prefix] =
              static_cast<uint32_t>(offset);
        } else {
          compact_lookup->starts64[prefix] = offset;
        }
      }
    } else {
      if (compact_lookup->compact_starts) {
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t i = 0;
             i < static_cast<int64_t>(compact_lookup->starts_size); ++i) {
          compact_lookup->starts32[i] =
              std::numeric_limits<uint32_t>::max();
        }
      } else {
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t i = 0;
             i < static_cast<int64_t>(compact_lookup->starts_size); ++i) {
          compact_lookup->starts64[i] =
              std::numeric_limits<uint64_t>::max();
        }
      }
      const unsigned extra_bases =
          compact_lookup->prefix_length - kBucketPrefixLength;
      const unsigned extra_bits = extra_bases * kBitsPerEdgeChar;
      const unsigned suffix_shift =
          (compact_lookup->suffix_bases - extra_bases) *
          kBitsPerEdgeChar;
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int prefix = 0; prefix < static_cast<int>(kNumBuckets); ++prefix) {
        const uint64_t begin = compact_lookup->bucket_offsets[prefix];
        const uint64_t end = compact_lookup->bucket_offsets[prefix + 1u];
        unsigned previous = std::numeric_limits<unsigned>::max();
        for (uint64_t i = begin; i < end; ++i) {
          const unsigned extra = static_cast<unsigned>(
              compact_lookup->suffixes[i] >> suffix_shift);
          if (extra == previous) continue;
          const size_t logical_prefix =
              (static_cast<size_t>(prefix) << extra_bits) | extra;
          if (compact_lookup->compact_starts) {
            compact_lookup->starts32[logical_prefix] =
                static_cast<uint32_t>(i);
          } else {
            compact_lookup->starts64[logical_prefix] = i;
          }
          previous = extra;
        }
      }
      if (compact_lookup->compact_starts) {
        compact_lookup->starts32[compact_lookup->starts_size - 1u] =
            static_cast<uint32_t>(stream_num_edges_);
        uint32_t next = static_cast<uint32_t>(stream_num_edges_);
        for (size_t i = compact_lookup->starts_size - 1u; i-- > 0;) {
          if (compact_lookup->starts32[i] ==
              std::numeric_limits<uint32_t>::max()) {
            compact_lookup->starts32[i] = next;
          } else {
            next = compact_lookup->starts32[i];
          }
        }
      } else {
        compact_lookup->starts64[compact_lookup->starts_size - 1u] =
            static_cast<uint64_t>(stream_num_edges_);
        uint64_t next = static_cast<uint64_t>(stream_num_edges_);
        for (size_t i = compact_lookup->starts_size - 1u; i-- > 0;) {
          if (compact_lookup->starts64[i] ==
              std::numeric_limits<uint64_t>::max()) {
            compact_lookup->starts64[i] = next;
          } else {
            next = compact_lookup->starts64[i];
          }
        }
      }
    }
    phase_timer.stop();
    xinfo("Compact mercy index: {} exact suffixes, {} bytes/suffix, "
          "{}-base lookup prefix, {} start bytes; build {.4} sec\n",
          stream_num_edges_, sizeof(uint64_t),
          compact_lookup->prefix_length, compact_lookup->StartsBytes(),
          phase_timer.elapsed());
    if (std::getenv("MEGAHIT_VALIDATE_COMPACT_MERCY") != nullptr) {
      uint64_t checked = 0;
      uint64_t missing_full = 0;
      uint64_t missing_prefix = 0;
      uint64_t wrong_last = 0;
#pragma omp parallel for schedule(static) num_threads(num_threads) \
    reduction(+ : checked, missing_full, missing_prefix, wrong_last)
      for (int64_t chunk_id = 0;
           chunk_id < std::min<int64_t>(
                          stream_edge_chunks_.size(), num_threads);
           ++chunk_id) {
        std::vector<uint32_t> records;
        const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
        ReadStreamedEdgeChunk(chunk, &records);
        for (uint32_t i = 0; i < chunk.num_records; ++i) {
          const uint32_t *edge =
              records.data() + static_cast<size_t>(i) *
                                   stream_words_per_edge_;
          MercyKmer edge_kmer;
          edge_kmer.InitFromPtr(edge, 0, stream_edge_length_);
          const int64_t full =
              compact_lookup->Search(edge_kmer, stream_edge_length_);
          const int64_t prefix = compact_lookup->Search(edge_kmer, opt_.k);
          ++checked;
          missing_full += full < 0;
          missing_prefix += prefix < 0;
          if (prefix >= 0) {
            wrong_last += compact_lookup->LastBase(prefix) !=
                          ExtractRawBase(edge, opt_.k);
          }
        }
      }
      xinfo("Compact mercy self-check: {} edges, missing full {}, missing "
            "prefix {}, prefix last-base differences {}\n",
            checked, missing_full, missing_prefix, wrong_last);
    }
  } else {
    edge_lookup.reset(new MercyLookupTable(seq_pkg_.seq_count(), opt_.k));
    InitLookupTable(edge_lookup.get(), seq_pkg_, num_threads);
    phase_timer.stop();
    const double lookup_init_seconds = phase_timer.elapsed();

    phase_timer.reset();
    phase_timer.start();
    InitMercySuffixIndex<MercyKmer>(edge_lookup.get(), seq_pkg_, opt_.k + 1,
                                    num_threads, opt_.host_mem);
    phase_timer.stop();
    xinfo(
        "Mercy lookup prefix: {}, memory: {} bytes, suffix index: {} bytes\n",
        edge_lookup->prefix_length, edge_lookup->StartsBytes(),
        edge_lookup->suffixes.size() * sizeof(edge_lookup->suffixes[0]));
    xinfo("Mercy lookup init: {.4} sec, suffix init: {.4} sec, workers: {}\n",
          lookup_init_seconds, phase_timer.elapsed(), num_threads);
  }

  const auto search_edge = [&](MercyKmer &kmer,
                               unsigned query_length) -> int64_t {
    if (compact_lookup) {
      return compact_lookup->Search(kmer, query_length);
    }
    return BinarySearchKmer(kmer, *edge_lookup, seq_pkg_, query_length);
  };
  const auto edge_last_base = [&](uint64_t edge_id) -> unsigned {
    if (compact_lookup) return compact_lookup->LastBase(edge_id);
    return seq_pkg_.GetSeqView(edge_id).base_at(opt_.k);
  };

  // Count already has the complete occurrence-context histogram for each
  // solid edge.  A missing context is a necessary condition for a graph
  // endpoint, so its compact .endcand shards form a strict candidate
  // superset.  Materialize both endpoint orientations from that small stream;
  // the existing per-group solid-side masks perform the final exact filter
  // after mercy edges are present.  This removes the only semantic reason
  // that every real record had to remain resident at once.
  if (stream_input_edges_ && stream_words_per_edge_ == NWords &&
      std::getenv("MEGAHIT_DISABLE_PRECOMPUTED_ENDPOINTS") == nullptr) {
    struct CandidateTip {
      uint16_t bucket;
      uint16_t reserved;
      std::array<uint32_t, NWords> data;
    };

    const unsigned item_words = static_cast<unsigned>(words_per_substr_);
    const unsigned label_words = DivCeiling(opt_.k, kCharsPerEdgeWord);
    const size_t candidate_record_bytes =
        static_cast<size_t>(stream_words_per_edge_) * sizeof(uint32_t) + 1u;
    std::vector<std::vector<CandidateTip>> thread_tips(num_threads);
    std::atomic<bool> invalid_candidate_stream{false};
    uint64_t candidate_records = 0;
    uint64_t exact_source_nodes = 0;
    uint64_t exact_sink_nodes = 0;
    int candidate_files_found = 0;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads) \
    reduction(+ : candidate_records, exact_source_nodes, exact_sink_nodes, \
                  candidate_files_found)
    for (int file_id = 0;
         file_id < static_cast<int>(stream_edge_fds_.size()); ++file_id) {
      const std::string path =
          opt_.input_prefix + ".endcand." + std::to_string(file_id);
      std::ifstream input(path, std::ios::binary | std::ios::ate);
      if (!input) continue;
      ++candidate_files_found;
      const std::streampos end = input.tellg();
      if (end < std::streampos{0} ||
          static_cast<uint64_t>(end) % candidate_record_bytes != 0) {
        invalid_candidate_stream.store(true, std::memory_order_relaxed);
        continue;
      }
      const uint64_t bytes = static_cast<uint64_t>(end);
      const uint64_t num_records = bytes / candidate_record_bytes;
      candidate_records += num_records;

      auto &tips = thread_tips[omp_get_thread_num()];
      const auto mask_label = [&](uint32_t *label) {
        const unsigned used = opt_.k % kCharsPerEdgeWord;
        if (used != 0) {
          label[label_words - 1u] &=
              std::numeric_limits<uint32_t>::max()
              << ((kCharsPerEdgeWord - used) * kBitsPerEdgeChar);
        }
      };
      const auto append_tip = [&](const MercyKmer &node, bool source) {
        CandidateTip tip{};
        std::array<uint32_t, NWords> bases{};
        std::copy(node.data(), node.data() + label_words, bases.begin());
        uint32_t metadata;
        if (source) {
          tip.bucket = static_cast<uint16_t>(bases[0] >> 16u);
          metadata =
              (uint32_t{1} << (kBWTCharNumBits + kBitsPerMul)) |
              (uint32_t{kSentinelValue} << kBitsPerMul) | kMaxMul;
        } else {
          const uint32_t first_base = bases[0] >> 30u;
          for (unsigned word = 0; word < label_words; ++word) {
            bases[word] <<= kBitsPerEdgeChar;
            if (word + 1u < label_words) {
              bases[word] |=
                  node.data()[word + 1u] >>
                  (kBitsPerEdgeWord - kBitsPerEdgeChar);
            }
          }
          mask_label(bases.data());
          tip.bucket = static_cast<uint16_t>(bases[0] >> 16u);
          metadata = (first_base << kBitsPerMul) | kMaxMul;
        }

        if (bucket_packed_records_) {
          constexpr unsigned kPrefixBits =
              kBucketPrefixLength * kBitsPerEdgeChar;
          for (unsigned word = 0; word < item_words; ++word) {
            tip.data[word] = bases[word] << kPrefixBits;
            if (word + 1u < label_words) {
              tip.data[word] |=
                  bases[word + 1u] >> (kBitsPerEdgeWord - kPrefixBits);
            }
          }
          AppendPackedMetadata(tip.data.data(), metadata);
        } else {
          std::copy(bases.begin(), bases.begin() + label_words,
                    tip.data.begin());
          tip.data[item_words - 1u] |= metadata;
        }
        tips.push_back(tip);
      };
      const auto process_record = [&](const uint8_t *record) {
        std::array<uint32_t, NWords> aligned_edge{};
        std::memcpy(aligned_edge.data(), record,
                    stream_words_per_edge_ * sizeof(uint32_t));
        const uint32_t *edge = aligned_edge.data();
        const uint8_t flags = record[stream_words_per_edge_ *
                                     sizeof(uint32_t)];
        // Context flags are expressed in the canonical count orientation.
        // Bits 0/1 prove a missing side and therefore emit a conservative
        // pair of physical endpoints.  Bits 2/3 mark the only ambiguous case:
        // support below twice the solid threshold can have been collected
        // from two physical orientations.  Resolve those sides against the
        // exact solid-edge index before emitting a sentinel.  The regular
        // per-(k-1)-mer masks below still make the final decision after mercy
        // edges have been merged.
        if ((flags & 3u) != 0) {
          if ((flags & 1u) != 0) {
            MercyKmer prefix(edge, 0, opt_.k);
            append_tip(prefix, true);
            prefix.ReverseComplement(opt_.k);
            append_tip(prefix, false);
            ++exact_source_nodes;
            ++exact_sink_nodes;
          }
          if ((flags & 2u) != 0) {
            MercyKmer suffix(edge, 1u, opt_.k);
            append_tip(suffix, false);
            suffix.ReverseComplement(opt_.k);
            append_tip(suffix, true);
            ++exact_source_nodes;
            ++exact_sink_nodes;
          }
        } else {
          if ((flags & 4u) != 0) {
            MercyKmer prefix(edge, 0, opt_.k);
            const unsigned base = (flags >> 4u) & 3u;
            MercyKmer query(prefix);
            query.ShiftPreappend(static_cast<uint8_t>(base), opt_.k + 1u);
            query = query.unique_format(opt_.k + 1u);
            const bool has_exact_in =
                search_edge(query, opt_.k + 1u) >= 0;
            if (!has_exact_in) {
              append_tip(prefix, true);
              prefix.ReverseComplement(opt_.k);
              append_tip(prefix, false);
              ++exact_source_nodes;
              ++exact_sink_nodes;
            }
          }
          if ((flags & 8u) != 0) {
            MercyKmer suffix(edge, 1u, opt_.k);
            const unsigned base = (flags >> 6u) & 3u;
            MercyKmer query(suffix);
            // `suffix` already occupies positions [0, k).  Growing it to a
            // (k+1)-mer only sets the newly exposed final base.
            query.SetBase(opt_.k, static_cast<uint8_t>(base));
            query = query.unique_format(opt_.k + 1u);
            const bool has_exact_out =
                search_edge(query, opt_.k + 1u) >= 0;
            if (!has_exact_out) {
              append_tip(suffix, false);
              suffix.ReverseComplement(opt_.k);
              append_tip(suffix, true);
              ++exact_source_nodes;
              ++exact_sink_nodes;
            }
          }
        }
      };

      // Candidate shards can be many GiB on a full metagenome.  Reading a
      // complete shard per worker made the temporary input buffers alone
      // scale with the data set and could exceed the intended RSS envelope.
      // Stream fixed-size record blocks instead; only the exact sentinel
      // candidates survive in thread_tips.
      constexpr uint64_t kRecordsPerBlock = uint64_t{1} << 16u;
      std::vector<uint8_t> records(
          static_cast<size_t>(std::min<uint64_t>(num_records,
                                                 kRecordsPerBlock)) *
          candidate_record_bytes);
      input.seekg(0);
      for (uint64_t first = 0; first < num_records;) {
        const uint64_t count =
            std::min<uint64_t>(kRecordsPerBlock, num_records - first);
        const size_t block_bytes =
            static_cast<size_t>(count) * candidate_record_bytes;
        input.read(reinterpret_cast<char *>(records.data()),
                   static_cast<std::streamsize>(block_bytes));
        if (!input) {
          invalid_candidate_stream.store(true, std::memory_order_relaxed);
          break;
        }
        for (uint64_t i = 0; i < count; ++i) {
          process_record(records.data() +
                         static_cast<size_t>(i) * candidate_record_bytes);
        }
        first += count;
      }
    }

    if (invalid_candidate_stream.load(std::memory_order_relaxed)) {
      xfatal("Invalid endpoint-candidate stream\n");
    }
    if (candidate_files_found ==
        static_cast<int>(stream_edge_fds_.size())) {
      const size_t matrix_items =
          static_cast<size_t>(num_threads) * kNumBuckets;
      std::vector<uint64_t> thread_bucket_counts(matrix_items, 0u);
      uint64_t total_tips = 0;
      for (unsigned tid = 0; tid < static_cast<unsigned>(num_threads);
           ++tid) {
        uint64_t *counts = thread_bucket_counts.data() +
                           static_cast<size_t>(tid) * kNumBuckets;
        for (const CandidateTip &tip : thread_tips[tid]) {
          ++counts[tip.bucket];
          ++total_tips;
        }
      }
      std::vector<uint64_t> all_bucket_begin(kNumBuckets + 1u, 0u);
      for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
        uint64_t count = 0;
        for (int tid = 0; tid < num_threads; ++tid) {
          count += thread_bucket_counts[
              static_cast<size_t>(tid) * kNumBuckets + bucket];
        }
        all_bucket_begin[bucket + 1u] = all_bucket_begin[bucket] + count;
      }
      std::vector<uint64_t> thread_bucket_cursor(matrix_items, 0u);
      for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
        uint64_t cursor = all_bucket_begin[bucket];
        for (int tid = 0; tid < num_threads; ++tid) {
          const size_t index =
              static_cast<size_t>(tid) * kNumBuckets + bucket;
          thread_bucket_cursor[index] = cursor;
          cursor += thread_bucket_counts[index];
        }
      }
      std::vector<uint32_t> all_tip_items(
          static_cast<size_t>(total_tips) * item_words);
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int tid = 0; tid < num_threads; ++tid) {
        uint64_t *cursor = thread_bucket_cursor.data() +
                           static_cast<size_t>(tid) * kNumBuckets;
        for (const CandidateTip &tip : thread_tips[tid]) {
          const uint64_t index = cursor[tip.bucket]++;
          std::copy(tip.data.begin(), tip.data.begin() + item_words,
                    all_tip_items.data() + index * item_words);
        }
        std::vector<CandidateTip>().swap(thread_tips[tid]);
      }
      std::vector<std::vector<CandidateTip>>().swap(thread_tips);
      std::vector<uint64_t>().swap(thread_bucket_counts);
      std::vector<uint64_t>().swap(thread_bucket_cursor);

      auto tip_sort = SelectSortingFunc(
          item_words, 0, Lv2SortIgnoredLowBytes(),
          Lv2SortIgnoredHighBytes());
#pragma omp parallel for schedule(dynamic, 16) num_threads(num_threads)
      for (int bucket_int = 0;
           bucket_int < static_cast<int>(kNumBuckets); ++bucket_int) {
        const unsigned bucket = static_cast<unsigned>(bucket_int);
        const uint64_t begin = all_bucket_begin[bucket];
        const uint64_t count = all_bucket_begin[bucket + 1u] - begin;
        if (count > 1) {
          tip_sort(all_tip_items.data() + begin * item_words, count);
        }
      }

      endpoint_tip_bucket_begin_.assign(kNumBuckets + 1u, 0u);
      std::vector<uint64_t> unique_counts(kNumBuckets, 0u);
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int bucket_int = 0;
           bucket_int < static_cast<int>(kNumBuckets); ++bucket_int) {
        const unsigned bucket = static_cast<unsigned>(bucket_int);
        const uint64_t begin = all_bucket_begin[bucket];
        const uint64_t end = all_bucket_begin[bucket + 1u];
        uint64_t count = 0;
        const uint32_t *previous = nullptr;
        for (uint64_t i = begin; i < end; ++i) {
          const uint32_t *item =
              all_tip_items.data() + i * item_words;
          if (previous == nullptr ||
              std::memcmp(previous, item,
                          item_words * sizeof(uint32_t)) != 0) {
            ++count;
            previous = item;
          }
        }
        unique_counts[bucket] = count;
      }
      for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
        endpoint_tip_bucket_begin_[bucket + 1u] =
            endpoint_tip_bucket_begin_[bucket] + unique_counts[bucket];
      }
      endpoint_tip_items_.resize(
          static_cast<size_t>(endpoint_tip_bucket_begin_.back()) *
          item_words);
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int bucket_int = 0;
           bucket_int < static_cast<int>(kNumBuckets); ++bucket_int) {
        const unsigned bucket = static_cast<unsigned>(bucket_int);
        const uint64_t begin = all_bucket_begin[bucket];
        const uint64_t end = all_bucket_begin[bucket + 1u];
        uint64_t output = endpoint_tip_bucket_begin_[bucket];
        const uint32_t *previous = nullptr;
        for (uint64_t i = begin; i < end; ++i) {
          const uint32_t *item =
              all_tip_items.data() + i * item_words;
          if (previous == nullptr ||
              std::memcmp(previous, item,
                          item_words * sizeof(uint32_t)) != 0) {
            std::copy(item, item + item_words,
                      endpoint_tip_items_.data() + output * item_words);
            ++output;
            previous = item;
          }
        }
        assert(output == endpoint_tip_bucket_begin_[bucket + 1u]);
      }
      std::vector<uint32_t>().swap(all_tip_items);
      std::vector<uint64_t>().swap(all_bucket_begin);
      std::vector<uint64_t>().swap(unique_counts);
      precomputed_endpoint_tips_ = true;
      xinfo("Endpoint candidates: {} records -> {} unique sentinel "
            "records (raw source/sink candidates {}/{}); all buckets are "
            "independently processable\n",
            candidate_records, endpoint_tip_bucket_begin_.back(),
            exact_source_nodes, exact_sink_nodes);
    } else if (candidate_files_found != 0) {
      xwarn("Incomplete endpoint-candidate stream ({}/{} shards); retaining "
            "global endpoint difference\n",
            candidate_files_found, stream_edge_fds_.size());
    }
  }

  BinaryReader binary_reader(opt_.input_prefix + ".cand");
  AsyncSequenceReader reader(&binary_reader, false, 1 << 16, 1 << 23);

  std::vector<std::vector<uint64_t>> thread_mercy_offsets(num_threads);
  std::vector<size_t> mercy_begins(num_threads + 1);
  std::vector<uint32_t> packed_mercy_edges;

  int64_t num_mercy_edges = 0;
  int64_t num_mercy_reads = 0;
  double query_seconds = 0;
  double append_seconds = 0;
  phase_timer.reset();
  phase_timer.start();

  while (true) {
    const auto &batch_reads = reader.Next();
    if (batch_reads.seq_count() == 0) {
      break;
    }
    xinfo("Read {} reads to search for mercy k-mers\n",
          batch_reads.seq_count());

    num_mercy_reads += batch_reads.seq_count();
    for (auto &offsets : thread_mercy_offsets) {
      offsets.clear();
    }

    SimpleTimer batch_timer;
    batch_timer.start();
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      auto &mercy_offsets = thread_mercy_offsets[tid];
      std::vector<uint8_t> state;

#pragma omp for schedule(static)
      for (unsigned read_id = 0; read_id < batch_reads.seq_count(); ++read_id) {
        auto seq_view = batch_reads.GetSeqView(read_id);
        unsigned read_len = seq_view.length();

        if (read_len < opt_.k + 2) {
          continue;
        }

        // bit 0: has incoming; bit 1: has outgoing.  Reuse one byte array per
        // worker instead of allocating two vector<bool> objects per read.
        state.assign(read_len, 0);
        MercyKmer kmer, rev_kmer;

        auto ptr_and_offset = seq_view.raw_address();
        kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, opt_.k);
        rev_kmer = kmer;
        rev_kmer.ReverseComplement(opt_.k);

        for (int i = 0; i + opt_.k <= read_len; ++i) {
          if (!(state[i] & 1u)) {
            if (search_edge(rev_kmer, opt_.k) != -1) {
              state[i] |= 1u;
            } else {
              rev_kmer.SetBase(opt_.k, 3);
              kmer.ShiftPreappend(0, opt_.k + 1);

              for (int c = 0; c < 4; ++c) {
                kmer.SetBase(0, c);
                if (kmer.cmp(rev_kmer, opt_.k + 1) > 0) {
                  break;
                }
                if (search_edge(kmer, opt_.k + 1) != -1) {
                  state[i] |= 1u;
                  break;
                }
              }

              rev_kmer.SetBase(opt_.k, 0);
              kmer.ShiftAppend(0, opt_.k + 1);
            }
          }

          int64_t edge_id = search_edge(kmer, opt_.k);

          if (edge_id != -1) {
            state[i] |= 2u;
            if (i + opt_.k < read_len &&
                edge_last_base(edge_id) ==
                    seq_view.base_at(i + opt_.k)) {
              state[i + 1] |= 1u;
            }
          } else {
            kmer.SetBase(opt_.k, 3);
            int next_char =
                i + opt_.k < read_len ? 3 - seq_view.base_at(i + opt_.k) : 0;
            rev_kmer.ShiftPreappend(next_char, opt_.k + 1);

            if (rev_kmer.cmp(kmer, opt_.k + 1) <= 0 &&
                search_edge(rev_kmer, opt_.k + 1) != -1) {
              state[i] |= 2u;
              state[i + 1] |= 1u;
            } else {
              for (int c = 0; c < 4; ++c) {
                if (c == next_char) {
                  continue;
                }
                rev_kmer.SetBase(0, c);
                if (rev_kmer.cmp(kmer, opt_.k + 1) > 0) {
                  break;
                }
                if (search_edge(rev_kmer, opt_.k + 1) != -1) {
                  state[i] |= 2u;
                  break;
                }
              }
            }

            kmer.SetBase(opt_.k, 0);
            rev_kmer.ShiftAppend(0, opt_.k + 1);
          }

          if (i + opt_.k < read_len) {
            int next_char = seq_view.base_at(i + opt_.k);
            kmer.ShiftAppend(next_char, opt_.k);
            rev_kmer.ShiftPreappend(3 - next_char, opt_.k);
          }
        }

        int last_no_out = -1;
        for (int i = 0; i + opt_.k <= read_len; ++i) {
          switch (state[i] & 3u) {
            case 1:
              last_no_out = i;
              break;
            case 2:
              if (last_no_out >= 0) {
                for (int j = last_no_out; j < i; ++j) {
                  mercy_offsets.emplace_back(
                      (uint64_t{read_id} << 32u) | static_cast<uint32_t>(j));
                }
              }
              last_no_out = -1;
              break;
            case 3:
              last_no_out = -1;
              break;
            default:
              break;
          }
        }
      }
    }
    batch_timer.stop();
    query_seconds += batch_timer.elapsed();

    batch_timer.reset();
    batch_timer.start();
    size_t batch_mercy_edges = 0;
    if (std::getenv("MEGAHIT_DISABLE_PARALLEL_MERCY_APPEND") != nullptr) {
      for (const auto &offsets : thread_mercy_offsets) {
        batch_mercy_edges += offsets.size();
        for (uint64_t encoded : offsets) {
          const unsigned read_id = encoded >> 32u;
          const unsigned offset = static_cast<uint32_t>(encoded);
          auto seq_view = batch_reads.GetSeqView(read_id);
          auto raw_address = seq_view.raw_address();
          MercyKmer mercy_edge(raw_address.first,
                               raw_address.second + offset, opt_.k + 1);
          seq_pkg_.AppendCompactSequence(mercy_edge.data(), opt_.k + 1);
        }
      }
    } else {
      mercy_begins[0] = 0;
      for (int tid = 0; tid < num_threads; ++tid) {
        mercy_begins[tid + 1] =
            mercy_begins[tid] + thread_mercy_offsets[tid].size();
      }
      batch_mercy_edges = mercy_begins.back();
      if (batch_mercy_edges >
          std::numeric_limits<size_t>::max() / NWords) {
        xfatal("Mercy-edge batch size overflow\n");
      }
      packed_mercy_edges.resize(batch_mercy_edges * NWords);

#pragma omp parallel for schedule(static)
      for (int tid = 0; tid < num_threads; ++tid) {
        size_t out_id = mercy_begins[tid];
        for (uint64_t encoded : thread_mercy_offsets[tid]) {
          const unsigned read_id = encoded >> 32u;
          const unsigned offset = static_cast<uint32_t>(encoded);
          auto seq_view = batch_reads.GetSeqView(read_id);
          auto raw_address = seq_view.raw_address();
          MercyKmer mercy_edge(raw_address.first,
                               raw_address.second + offset, opt_.k + 1);
          std::memcpy(packed_mercy_edges.data() + out_id * NWords,
                      mercy_edge.data(), sizeof(uint32_t) * NWords);
          ++out_id;
        }
      }
      seq_pkg_.AppendFixedLengthCompactSequences32(
          packed_mercy_edges.data(), batch_mercy_edges, opt_.k + 1, NWords,
          num_threads);
    }
    num_mercy_edges += batch_mercy_edges;
    batch_timer.stop();
    append_seconds += batch_timer.elapsed();
  }

  multiplicity.insert(multiplicity.end(), num_mercy_edges, 1);
  phase_timer.stop();
  xinfo("Mercy candidate scan and append: {.4} sec\n",
        phase_timer.elapsed());
  xinfo("Mercy query: {.4} sec, append: {.4} sec\n", query_seconds,
        append_seconds);
  xinfo("Number of reads: {}, Number of mercy edges: {}\n", num_mercy_reads,
        num_mercy_edges);
}

void SeqToSdbg::GenMercyEdges() {
  const unsigned words = DivCeiling(opt_.k + 1, kCharsPerEdgeWord);
#define DISPATCH_MERCY_WORDS(N) \
  case N:                        \
    return GenMercyEdgesForK<N>()
  switch (words) {
    DISPATCH_MERCY_WORDS(1);
    DISPATCH_MERCY_WORDS(2);
    DISPATCH_MERCY_WORDS(3);
    DISPATCH_MERCY_WORDS(4);
    DISPATCH_MERCY_WORDS(5);
    DISPATCH_MERCY_WORDS(6);
    DISPATCH_MERCY_WORDS(7);
    DISPATCH_MERCY_WORDS(8);
    DISPATCH_MERCY_WORDS(9);
    DISPATCH_MERCY_WORDS(10);
    DISPATCH_MERCY_WORDS(11);
    DISPATCH_MERCY_WORDS(12);
    DISPATCH_MERCY_WORDS(13);
    DISPATCH_MERCY_WORDS(14);
    DISPATCH_MERCY_WORDS(15);
    DISPATCH_MERCY_WORDS(16);
    default:
      xfatal("Invalid number of mercy k-mer words: {}\n", words);
  }
#undef DISPATCH_MERCY_WORDS
}

SeqToSdbg::MemoryStat SeqToSdbg::Initialize() {
  full_words_per_substr_ =
      DivCeiling(opt_.k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMul,
                 kBitsPerEdgeWord);
  const unsigned bucket_prefix_bits =
      kBucketPrefixLength * kBitsPerEdgeChar;
  bucket_packed_base_bits_ =
      opt_.k > kBucketPrefixLength
          ? opt_.k * kBitsPerEdgeChar - bucket_prefix_bits
          : 0;
  bucket_packed_total_bits_ =
      bucket_packed_base_bits_ + kBWTCharNumBits + 1 + kBitsPerMul;
  const int64_t bucket_packed_words =
      DivCeiling(bucket_packed_total_bits_, kBitsPerEdgeWord);
  bucket_packed_records_ =
      opt_.k > kBucketPrefixLength &&
      bucket_packed_words < full_words_per_substr_ &&
      std::getenv("MEGAHIT_DISABLE_BUCKET_PACKED_SEQ_RECORDS") == nullptr;
  words_per_substr_ = bucket_packed_records_ ? bucket_packed_words
                                              : full_words_per_substr_;

  const bool real_only_requested =
      opt_.mem_flag != 0 &&
      (std::getenv("MEGAHIT_ENABLE_REAL_ONLY_SEQ_RECORDS") != nullptr ||
       std::getenv("MEGAHIT_DISABLE_REAL_ONLY_SEQ_RECORDS") == nullptr);
  real_only_records_ = real_only_requested;

  EdgeIoMetadata edge_metadata;
  bool have_edge_metadata = false;
  bool mercy_candidates_empty = false;
  bool stream_requires_mercy_index = false;
  int64_t edge_reserve_count = 0;
  int64_t edge_bases_to_reserve = 0;
  int64_t edge_length = 0;

  // reserve space
  {
    int64_t bases_to_reserve = 0;
    int64_t num_contigs_to_reserve = 0;
    int64_t num_multiplicities_to_reserve = 0;
    auto checked_reserve_add = [](int64_t *total, int64_t value,
                                  const char *description) {
      if (value < 0 ||
          *total > std::numeric_limits<int64_t>::max() - value) {
        xfatal("Invalid or overflowing {s} reserve\n", description);
      }
      *total += value;
    };
    auto checked_reserve_multiply = [](int64_t lhs, int64_t rhs,
                                       const char *description) -> int64_t {
      if (lhs < 0 || rhs < 0 ||
          (lhs != 0 && rhs > std::numeric_limits<int64_t>::max() / lhs)) {
        xfatal("Invalid or overflowing {s} reserve\n", description);
      }
      return lhs * rhs;
    };

    if (!opt_.input_prefix.empty()) {
      std::ifstream edge_info(opt_.input_prefix + ".edges.info");
      edge_metadata.Deserialize(edge_info);
      have_edge_metadata = true;
      int64_t num_edges = edge_metadata.num_edges;
      xinfo("Number edges: {}\n", num_edges);
      if (num_edges < 0) {
        xfatal("Invalid negative edge count in metadata: {}\n", num_edges);
      }

      if (opt_.need_mercy) {
        const char *mercy_factor =
            std::getenv("MEGAHIT_NUM_MERCY_FACTOR");
        bool custom_factor_applied = false;
        if (mercy_factor != nullptr) {
          errno = 0;
          char *end = nullptr;
          const double factor = std::strtod(mercy_factor, &end);
          if (end != mercy_factor && *end == '\0' && errno == 0 &&
              std::isfinite(factor) && factor >= 0) {
            const long double scaled_edges =
                static_cast<long double>(num_edges) *
                (1.0L + static_cast<long double>(factor));
            // 2^63 is exactly representable even when long double aliases
            // double; every non-negative int64_t value is strictly below it.
            const long double int64_limit = std::ldexp(1.0L, 63);
            if (scaled_edges < int64_limit) {
              num_edges = static_cast<int64_t>(scaled_edges);
              custom_factor_applied = true;
            }
          }
          if (!custom_factor_applied) {
            xwarn(
                "Invalid MEGAHIT_NUM_MERCY_FACTOR '{s}'; using default "
                "25% reserve\n",
                mercy_factor);
          }
        }
        if (!custom_factor_applied) {
          const int64_t extra_edges = num_edges >> 2;
          if (extra_edges > std::numeric_limits<int64_t>::max() - num_edges) {
            xfatal("Mercy-edge reserve overflows int64_t\n");
          }
          num_edges += extra_edges;  // it is rare that # mercy > 25%
        }
      }

      edge_length = static_cast<int64_t>(edge_metadata.kmer_size) + 1;
      edge_reserve_count = num_edges;
      edge_bases_to_reserve =
          checked_reserve_multiply(num_edges, edge_length, "edge bases");

      if (opt_.need_mercy) {
        std::ifstream candidate(opt_.input_prefix + ".cand",
                                std::ios::binary | std::ios::ate);
        mercy_candidates_empty =
            candidate.is_open() && candidate.tellg() == std::streampos{0};
      }
    }

    if (!opt_.contig.empty()) {
      auto sizes = ContigReader(opt_.contig).GetNumContigsAndBases();
      checked_reserve_add(&bases_to_reserve, sizes.second, "total bases");
      checked_reserve_add(&num_contigs_to_reserve, sizes.first, "contig");
      checked_reserve_add(&num_multiplicities_to_reserve, sizes.first,
                          "multiplicity");
    }

    if (!opt_.bubble_seq.empty()) {
      auto sizes = ContigReader(opt_.bubble_seq).GetNumContigsAndBases();
      checked_reserve_add(&bases_to_reserve, sizes.second, "total bases");
      checked_reserve_add(&num_contigs_to_reserve, sizes.first, "contig");
      checked_reserve_add(&num_multiplicities_to_reserve, sizes.first,
                          "multiplicity");
    }

    if (!opt_.addi_contig.empty()) {
      auto sizes = ContigReader(opt_.addi_contig).GetNumContigsAndBases();
      checked_reserve_add(&bases_to_reserve, sizes.second, "total bases");
      checked_reserve_add(&num_contigs_to_reserve, sizes.first, "contig");
      checked_reserve_add(&num_multiplicities_to_reserve, sizes.first,
                          "multiplicity");
    }

    if (!opt_.local_contig.empty()) {
      auto sizes = ContigReader(opt_.local_contig).GetNumContigsAndBases();
      checked_reserve_add(&bases_to_reserve, sizes.second, "total bases");
      checked_reserve_add(&num_contigs_to_reserve, sizes.first, "contig");
      checked_reserve_add(&num_multiplicities_to_reserve, sizes.first,
                          "multiplicity");
    }

    const auto estimate_data_bytes = [&](int64_t retained_edges) {
      unsigned __int128 bases = static_cast<uint64_t>(bases_to_reserve);
      unsigned __int128 multiplicities =
          static_cast<uint64_t>(num_multiplicities_to_reserve);
      if (retained_edges > 0) {
        bases += static_cast<uint64_t>(retained_edges) *
                 static_cast<uint64_t>(edge_length);
        multiplicities += static_cast<uint64_t>(retained_edges);
      }
      // CompactVector stores four DNA bases per byte.  Include sequence
      // boundary capacity and multiplicities as conservative resident data.
      return (bases + 3u) / 4u +
             multiplicities * sizeof(mul_t) +
             static_cast<uint64_t>(num_contigs_to_reserve + 1u) *
                 sizeof(uint64_t);
    };
    const auto estimate_real_items = [&](int64_t input_edges) {
      return static_cast<unsigned __int128>(2u) *
             (static_cast<uint64_t>(input_edges) +
              static_cast<uint64_t>(bases_to_reserve));
    };
    const auto direct_fits = [&](unsigned __int128 items,
                                 unsigned __int128 data_bytes) {
      if (!real_only_requested || opt_.host_mem <= 0) {
        return false;
      }
      const unsigned __int128 direct_bytes =
          items * static_cast<uint64_t>(words_per_substr_) * sizeof(uint32_t);
      if (direct_bytes > static_cast<unsigned __int128>(
                             AutomaticDirectWorkspaceLimit(static_cast<uint64_t>(
                                 std::min<unsigned __int128>(
                                     data_bytes,
                                     std::numeric_limits<uint64_t>::max()))))) {
        return false;
      }
      const uint64_t fixed_headroom = std::max<uint64_t>(
          uint64_t{256} << 20u, static_cast<uint64_t>(opt_.host_mem) / 100u);
      const uint64_t per_thread_headroom =
          static_cast<uint64_t>(std::max(1, opt_.n_threads)) * kNumBuckets *
          sizeof(int64_t) * 3u;
      return data_bytes + direct_bytes + fixed_headroom +
                 per_thread_headroom <=
             static_cast<uint64_t>(opt_.host_mem);
    };

    const bool stream_semantically_available =
        have_edge_metadata && real_only_requested &&
        std::getenv("MEGAHIT_DISABLE_STREAMED_EDGE_INPUT") == nullptr;
    // Streamed fixed edges live in the sortable records, not in SeqPackage.
    // Only supplemental contigs remain resident alongside that workspace.
    // Charging the packed edge copy here would reject a layout that never
    // allocates that copy and unnecessarily fall back to the slower resident
    // representation.
    const bool streamed_direct_fits = direct_fits(
        estimate_real_items(edge_reserve_count),
        estimate_data_bytes(0));

    // The generic compact-locator path can reread resident SeqPackage
    // sequences, whereas streamed fixed edges currently have no stable
    // random-access locator once the sequential chunk buffer is released.
    // A pure edge stream is handled by the bounded partial-direct path and a
    // mixed stream is handled when the complete direct workspace fits.  If a
    // mixed input cannot use direct records, retain the exact resident-edge
    // representation instead of selecting an unsupported join of streamed
    // edges and compact sequence locators.  This is a representation
    // capability check, not a workload- or machine-specific threshold.
    const bool has_supplemental_sequences =
        bases_to_reserve != 0 || num_contigs_to_reserve != 0;
    const bool stream_layout_supported =
        !has_supplemental_sequences || streamed_direct_fits;
    if (stream_semantically_available && stream_layout_supported) {
      stream_input_edges_ = ConfigureStreamedEdgeInput(edge_metadata);
      stream_requires_mercy_index =
          stream_input_edges_ && opt_.need_mercy && !mercy_candidates_empty;
      stream_compact_mercy_index_ =
          stream_requires_mercy_index && stream_input_unordered_ &&
          edge_length > static_cast<int64_t>(kBucketPrefixLength) &&
          edge_length <= 40;
      if (stream_compact_mercy_index_) {
        xinfo("Unordered mercy input will use an exact canonical-order "
              "compact edge index\n");
      }
    } else if (stream_semantically_available &&
               has_supplemental_sequences && !streamed_direct_fits) {
      xinfo("Streaming fixed edges skipped: mixed input requires the exact "
            "resident-locator path when the complete direct workspace does "
            "not fit\n");
    }
    if (!stream_input_edges_ && !streamed_direct_fits) {
      real_only_records_ = false;
      if (real_only_requested) {
        xinfo("Real-edge-only sorting skipped: complete exact workspace does "
              "not fit the requested memory budget\n");
      }
    }
    if (!stream_input_edges_ ||
        (stream_requires_mercy_index && !stream_compact_mercy_index_)) {
      checked_reserve_add(&bases_to_reserve, edge_bases_to_reserve,
                          "total bases");
      checked_reserve_add(&num_multiplicities_to_reserve,
                          edge_reserve_count, "multiplicity");
    } else if (opt_.need_mercy) {
      xinfo("Mercy candidate stream is empty; skipping edge lookup and "
            "retaining the raw edge stream\n");
    }

    xinfo("Bases to reserve: {}, number contigs: {}, number multiplicity: {}\n",
          bases_to_reserve, num_contigs_to_reserve,
          num_multiplicities_to_reserve);
    seq_pkg_.ReserveSequences(num_contigs_to_reserve);
    seq_pkg_.ReserveBases(bases_to_reserve);
    multiplicity.reserve(num_multiplicities_to_reserve);
  }

  xinfo("Before reading, sizeof seq_package: {}, multiplicity vector: {}\n",
        seq_pkg_.size_in_byte(), multiplicity.capacity());

  if (!opt_.input_prefix.empty() &&
      (!stream_input_edges_ ||
       (stream_requires_mercy_index && !stream_compact_mercy_index_))) {
    EdgeReader reader(opt_.input_prefix);
    reader.SetMultiplicityVec(&multiplicity);
    auto n_read = reader.ReadSortedBulk(&seq_pkg_, &multiplicity,
                                        opt_.n_threads, opt_.host_mem);
    if (n_read < 0) {
      n_read = reader.ReadAll(&seq_pkg_, false);
    }
    xinfo("Read {} edges.\n", n_read);
    xinfo(
        "After reading, sizeof seq_package: {}/{}/{}, multiplicity vector: "
        "{}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (opt_.need_mercy &&
      (!stream_input_edges_ || stream_requires_mercy_index)) {
    SimpleTimer timer;
    timer.reset();
    timer.start();
    xinfo("Adding mercy edges...\n");

    GenMercyEdges();
    timer.stop();
    xinfo("Done. Time elapsed: {.4}\n", timer.elapsed());
    xinfo(
        "After adding mercy, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
    if (stream_requires_mercy_index && !stream_compact_mercy_index_) {
      RetainMercyEdgesForStreamedInput(
          static_cast<size_t>(edge_metadata.num_edges));
    }
  }

  if (!opt_.contig.empty()) {
    ContigReader reader(opt_.contig);
    reader.SetExtendLoop(opt_.k_from, opt_.k)->SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());

    // read bubble
    ContigReader bubble_reader(opt_.bubble_seq);
    bubble_reader.SetMinLen(opt_.k + 1);
    n_read = bubble_reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                   contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.bubble_seq.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (!opt_.addi_contig.empty()) {
    ContigReader reader(opt_.addi_contig);
    reader.SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.addi_contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  if (!opt_.local_contig.empty()) {
    ContigReader reader(opt_.local_contig);
    reader.SetMinLen(opt_.k + 1);
    bool contig_reverse = true;
    auto n_read = reader.ReadAllWithMultiplicity(&seq_pkg_, &multiplicity,
                                                 contig_reverse);
    xinfo("Read {} contigs from {}.\n", n_read, opt_.local_contig.c_str());
    xinfo(
        "After reading contigs, sizeof seq_package: {}/{}/{}, multiplicity "
        "vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());
  }

  xinfo("Finally, sizeof seq_package: {}/{}/{}, multiplicity vector: {}/{}\n",
        seq_pkg_.size_in_byte(), seq_pkg_.seq_count(), seq_pkg_.base_count(),
        multiplicity.size(), multiplicity.capacity());

  ConfigurePackedSeqOffsets();
  // Recovering a sequence ID from a global base offset requires pos_to_id.
  // The packed fast path carries the ID directly and can omit this index.
  if (!packed_seq_offsets_) {
    seq_pkg_.BuildIndex();
  }
  if (bucket_packed_records_) {
    xinfo(
        "Bucket-packed sequence records: {} -> {} words/item ({} shared "
        "prefix bits omitted)\n",
        full_words_per_substr_, words_per_substr_, bucket_prefix_bits);
  }
  if (real_only_records_) {
    xinfo("Real-edge-only sequence sorting enabled; necessary sentinels "
          "will be synthesized by exact endpoint set difference\n");
  }

  // For a pure streamed edge input the direct-record cardinality is exact.
  // Fault its pages in one linear parallel walk before the unordered producer
  // starts, instead of paying millions of first-touch faults on scattered
  // bucket destinations in the hot Lv1 loop.  AdjustMemory revalidates the
  // allocation against the actual histogram and memory budget before adopting
  // it, so mixed/mercy inputs retain the regular sizing path.
  if (stream_input_edges_ && seq_pkg_.seq_count() == 0 &&
      real_only_records_ && stream_num_edges_ > 0 &&
      std::getenv("MEGAHIT_DISABLE_STREAM_DIRECT_PREFAULT") == nullptr &&
      stream_num_edges_ <= std::numeric_limits<int64_t>::max() / 2) {
    PreallocateDirectItems(
        stream_num_edges_ * 2,
        static_cast<int64_t>(seq_pkg_.size_in_byte() +
                             multiplicity.capacity() * sizeof(mul_t)));
  }

  sdbg_writer_.set_num_threads(opt_.n_threads);
  sdbg_writer_.set_kmer_size(opt_.k);
  sdbg_writer_.set_num_buckets(kNumBuckets);
  sdbg_writer_.set_file_prefix(opt_.output_prefix);
  sdbg_writer_.InitFiles();

  return {
      static_cast<int64_t>(stream_edge_chunks_.size() +
                           seq_pkg_.seq_count()),
      static_cast<int64_t>(seq_pkg_.size_in_byte() +
                           multiplicity.capacity() * sizeof(mul_t)),
      words_per_substr_,
      0,
  };
}

void SeqToSdbg::Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                  std::array<int64_t, kNumBuckets> *out) {
  auto &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);

  if (use_stream_bucket_histogram_ && seq_pkg_.seq_count() == 0) {
    for (unsigned file_id = 0;
         file_id + 1u < stream_file_chunk_offsets_.size(); ++file_id) {
      if (seq_from !=
              static_cast<int64_t>(stream_file_chunk_offsets_[file_id]) ||
          seq_to != static_cast<int64_t>(
                        stream_file_chunk_offsets_[file_id + 1u])) {
        continue;
      }
      const uint64_t *row = stream_bucket_histogram_.counts.data() +
                            static_cast<size_t>(file_id) * kNumBuckets;
      for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
        if (row[bucket] >
            static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
          xfatal("Precomputed bucket {} in edge shard {} overflows int64\n",
                 bucket, file_id);
        }
        bucket_sizes[bucket] = static_cast<int64_t>(row[bucket]);
      }
      return;
    }
  }

  const int64_t num_stream_chunks =
      static_cast<int64_t>(stream_edge_chunks_.size());
  if (seq_from < num_stream_chunks) {
    std::vector<uint32_t> records;
    const int64_t chunk_end = std::min(seq_to, num_stream_chunks);
    const unsigned reverse_window_offset =
        stream_edge_length_ - 1u - kBucketPrefixLength;
    for (int64_t chunk_id = seq_from; chunk_id < chunk_end; ++chunk_id) {
      const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
      ReadStreamedEdgeChunk(chunk, &records);
      const uint32_t *edge = records.data();
      for (uint32_t i = 0; i < chunk.num_records; ++i,
                    edge += stream_words_per_edge_) {
        const uint16_t key = ExtractRawBaseWindow8(edge, 1u);
        const uint16_t reverse_key = ReverseComplementWindow8(
            ExtractRawBaseWindow8(edge, reverse_window_offset));
        ++bucket_sizes[key];
        ++bucket_sizes[reverse_key];
      }
    }
  }

  const int64_t loaded_from =
      std::max(seq_from, num_stream_chunks) - num_stream_chunks;
  const int64_t loaded_to =
      std::max(std::min(seq_to,
                        num_stream_chunks +
                            static_cast<int64_t>(seq_pkg_.seq_count())),
               num_stream_chunks) -
      num_stream_chunks;
  for (int64_t seq_id = loaded_from; seq_id < loaded_to; ++seq_id) {
    auto seq_view = seq_pkg_.GetSeqView(seq_id);
    unsigned seq_len = seq_view.length();

    if (seq_len < opt_.k + 1) {
      continue;
    }

    uint32_t key = 0;  // $$$$$$$$

    // build initial partial key
    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase + seq_view.base_at(i);
    }

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    unsigned local_offset = 0;
    const unsigned final_offset = seq_len - opt_.k + 1u;
    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = (key * kBucketBase + seq_view.base_at(i)) % kNumBuckets;
      if (!real_only_records_ ||
          (local_offset > 0 && local_offset < final_offset)) {
        bucket_sizes[key]++;
      }
      ++local_offset;
    }

    // reverse complement
    key = 0;

    for (int i = 0; i < static_cast<int>(kBucketPrefixLength) - 1; ++i) {
      key = key * kBucketBase +
            (3 - seq_view.base_at(seq_len - 1 - i));  // complement
    }

    local_offset = 0;
    for (int i = kBucketPrefixLength - 1;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = key * kBucketBase + (3 - seq_view.base_at(seq_len - 1 - i));
      key %= kNumBuckets;
      if (!real_only_records_ ||
          (local_offset > 0 && local_offset < final_offset)) {
        bucket_sizes[key]++;
      }
      ++local_offset;
    }
  }
}

template <bool PackedOffsets>
void SeqToSdbg::Lv1FillOffsetsFor(OffsetFiller &filler, int64_t seq_from,
                                  int64_t seq_to) {
  const int64_t num_stream_chunks =
      static_cast<int64_t>(stream_edge_chunks_.size());
  if (seq_from < num_stream_chunks) {
    xfatal("Streamed edge input unexpectedly entered compact-offset mode\n");
  }
  seq_from -= num_stream_chunks;
  seq_to -= num_stream_chunks;
// ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(key, strand)                     \
  do {                                                         \
    if (filler.IsHandling(key)) {                              \
      filler.WriteNextOffset(key, encoded_offset | (strand)); \
    }                                                          \
  } while (0)
  // =========== end macro ==========================

  for (int64_t seq_id = seq_from; seq_id < seq_to; ++seq_id) {
    auto seq_view = seq_pkg_.GetSeqView(seq_id);
    unsigned seq_len = seq_view.length();
    if (seq_len < opt_.k + 1) {
      continue;
    }

    // build initial partial key
    Kmer<1, uint32_t> kmer, rev_kmer;
    auto ptr_and_offset = seq_view.raw_address();
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                     kBucketPrefixLength);
    auto rev_ptr_and_offset =
        seq_view.raw_address(seq_len - kBucketPrefixLength);
    rev_kmer.InitFromPtr(rev_ptr_and_offset.first, rev_ptr_and_offset.second,
                         kBucketPrefixLength);
    rev_kmer.ReverseComplement(kBucketPrefixLength);

    int key = kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    int rev_key = rev_kmer.data()[0] >> (32 - kBucketPrefixLength * 2);
    int64_t encoded_offset =
        PackedOffsets
            ? static_cast<int64_t>(static_cast<uint64_t>(seq_id)
                                   << seq_locator_shift_)
            : static_cast<int64_t>(seq_view.full_offset_in_pkg() << 1u);
    unsigned local_offset = 0;
    const unsigned final_offset = seq_len - opt_.k + 1u;
    if (!real_only_records_) {
      CHECK_AND_SAVE_OFFSET(key, 0);
      CHECK_AND_SAVE_OFFSET(rev_key, 1);
    }

    // sequence = xxxxxxxxx
    // edges = $xxxx, xxxxx, ..., xxxx$
    for (int i = kBucketPrefixLength;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = (key * kBucketBase + seq_view.base_at(i)) % kNumBuckets;
      rev_key = rev_key * kBucketBase + (3 - seq_view.base_at(seq_len - 1 - i));
      rev_key %= kNumBuckets;
      encoded_offset += 2;
      ++local_offset;
      if (!real_only_records_ || local_offset < final_offset) {
        CHECK_AND_SAVE_OFFSET(key, 0);
        CHECK_AND_SAVE_OFFSET(rev_key, 1);
      }
    }
  }
#undef CHECK_AND_SAVE_OFFSET
}

void SeqToSdbg::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                               int64_t seq_to) {
  if (UsingDirectLv1Items()) {
#define DISPATCH_DIRECT_WORDS(N)                                      \
  case N:                                                             \
    if (bucket_packed_records_) {                                     \
      return Lv1FillDirectItemsFor<N, true>(filler, seq_from, seq_to); \
    }                                                                 \
    return Lv1FillDirectItemsFor<N, false>(filler, seq_from, seq_to)
    switch (full_words_per_substr_) {
      case 1:
        return Lv1FillDirectItemsFor<1, false>(filler, seq_from, seq_to);
      DISPATCH_DIRECT_WORDS(2);
      DISPATCH_DIRECT_WORDS(3);
      DISPATCH_DIRECT_WORDS(4);
      DISPATCH_DIRECT_WORDS(5);
      DISPATCH_DIRECT_WORDS(6);
      DISPATCH_DIRECT_WORDS(7);
      DISPATCH_DIRECT_WORDS(8);
      DISPATCH_DIRECT_WORDS(9);
      DISPATCH_DIRECT_WORDS(10);
      DISPATCH_DIRECT_WORDS(11);
      DISPATCH_DIRECT_WORDS(12);
      DISPATCH_DIRECT_WORDS(13);
      DISPATCH_DIRECT_WORDS(14);
      DISPATCH_DIRECT_WORDS(15);
      DISPATCH_DIRECT_WORDS(16);
      DISPATCH_DIRECT_WORDS(17);
      default:
        xfatal("Invalid number of direct substring words: {}\n",
               words_per_substr_);
    }
#undef DISPATCH_DIRECT_WORDS
  }

  if (packed_seq_offsets_) {
    return Lv1FillOffsetsFor<true>(filler, seq_from, seq_to);
  }
  return Lv1FillOffsetsFor<false>(filler, seq_from, seq_to);
}

bool SeqToSdbg::Lv1SupportsDirectItems() const {
  return std::getenv("MEGAHIT_DISABLE_DIRECT_SEQ_ITEMS") == nullptr;
}

int64_t SeqToSdbg::Lv1DirectWordsPerItem() const {
  return words_per_substr_;
}

int64_t SeqToSdbg::Lv1DirectAuxWordsPerItem() const {
  return 0;
}

int64_t SeqToSdbg::Lv1DirectMemoryLimit() const {
  return BoundedTransientWorkspaceLimit(CurrentRetainedBytes(), true);
}

bool SeqToSdbg::Lv1AllowsPartialDirectItems() const {
  return stream_input_edges_ && precomputed_endpoint_tips_ &&
         std::getenv("MEGAHIT_VALIDATE_PRECOMPUTED_ENDPOINTS") == nullptr;
}

int64_t SeqToSdbg::Lv1AutoWorkspaceLimit() const {
  if (opt_.mem_flag != 1) return 0;
  return BoundedTransientWorkspaceLimit(CurrentRetainedBytes(), false);
}

bool SeqToSdbg::Lv1UseWriteCombine() const {
  // Unordered fixed-edge streams are processed in cache/TLB-sized macro
  // partitions below.  Their active destination set is already bounded, so
  // allocating one 192-byte staging block for every bucket and worker would
  // add ~1.5 GiB without improving locality.
  return !stream_input_unordered_ &&
         BaseSequenceSortingEngine::Lv1UseWriteCombine();
}

int SeqToSdbg::Lv2SortIgnoredLowBytes() const {
  if (!bucket_packed_records_) {
    return 0;
  }
  return static_cast<int>(
      (words_per_substr_ * kBitsPerEdgeWord - bucket_packed_total_bits_) / 8u);
}

int SeqToSdbg::Lv2SortIgnoredHighBytes() const {
  return bucket_packed_records_ ? 0 : 2;
}

bool SeqToSdbg::Lv2NeedsEmptyBucketPostprocess(unsigned bucket) const {
  return precomputed_endpoint_tips_ &&
         bucket + 1u < endpoint_tip_bucket_begin_.size() &&
         endpoint_tip_bucket_begin_[bucket] !=
             endpoint_tip_bucket_begin_[bucket + 1u];
}

void SeqToSdbg::AppendPackedMetadata(uint32_t *item,
                                     uint32_t metadata) const {
  assert(bucket_packed_records_);
  constexpr unsigned kMetadataBits =
      kBWTCharNumBits + 1u + kBitsPerMul;
  constexpr uint32_t kMetadataMask = (uint32_t{1} << kMetadataBits) - 1u;
  const unsigned bit_offset = bucket_packed_base_bits_;
  const unsigned word = bit_offset / kBitsPerEdgeWord;
  const unsigned in_word = bit_offset % kBitsPerEdgeWord;
  assert(word < static_cast<unsigned>(words_per_substr_));
  assert(in_word + kMetadataBits <= 2u * kBitsPerEdgeWord);
  const uint64_t field = static_cast<uint64_t>(metadata & kMetadataMask)
                         << (2u * kBitsPerEdgeWord - in_word - kMetadataBits);
  item[word] |= static_cast<uint32_t>(field >> kBitsPerEdgeWord);
  if (in_word + kMetadataBits > kBitsPerEdgeWord) {
    assert(word + 1u < static_cast<unsigned>(words_per_substr_));
    item[word + 1u] |= static_cast<uint32_t>(field);
  }
}

uint32_t SeqToSdbg::ExtractPackedField(const uint32_t *item,
                                       unsigned bit_offset,
                                       unsigned bit_width) const {
  assert(bit_width > 0 && bit_width < 32);
  const unsigned word = bit_offset / kBitsPerEdgeWord;
  const unsigned in_word = bit_offset % kBitsPerEdgeWord;
  assert(word < static_cast<unsigned>(words_per_substr_));
  uint64_t window = static_cast<uint64_t>(item[word]) << kBitsPerEdgeWord;
  if (word + 1u < static_cast<unsigned>(words_per_substr_)) {
    window |= item[word + 1u];
  }
  const unsigned shift =
      2u * kBitsPerEdgeWord - in_word - bit_width;
  return static_cast<uint32_t>(window >> shift) &
         ((uint32_t{1} << bit_width) - 1u);
}

uint32_t SeqToSdbg::ExtractPackedMetadata(const uint32_t *item) const {
  constexpr unsigned kMetadataBits =
      kBWTCharNumBits + 1u + kBitsPerMul;
  return ExtractPackedField(item, bucket_packed_base_bits_, kMetadataBits);
}

bool SeqToSdbg::IsDiffPackedKMinusOneMer(const uint32_t *lhs,
                                         const uint32_t *rhs) const {
  // The shared bucket prefix is equal by construction.  Compare the rest of
  // the (k-1)-mer, stopping before the final 'a' base and metadata.
  const unsigned group_bits =
      (opt_.k - 1u - kBucketPrefixLength) * kBitsPerEdgeChar;
  const unsigned full_words = group_bits / kBitsPerEdgeWord;
  const unsigned remaining_bits = group_bits % kBitsPerEdgeWord;
  for (unsigned i = 0; i < full_words; ++i) {
    if (lhs[i] != rhs[i]) {
      return true;
    }
  }
  if (remaining_bits != 0) {
    const uint32_t mask = std::numeric_limits<uint32_t>::max()
                          << (kBitsPerEdgeWord - remaining_bits);
    if ((lhs[full_words] & mask) != (rhs[full_words] & mask)) {
      return true;
    }
  }
  return false;
}

int SeqToSdbg::ExtractPackedA(const uint32_t *item) const {
  const uint32_t metadata = ExtractPackedMetadata(item);
  const bool non_dollar =
      (metadata >> (kBWTCharNumBits + kBitsPerMul)) & 1u;
  if (!non_dollar) {
    return kSentinelValue;
  }
  const unsigned last_base_offset =
      (opt_.k - 1u - kBucketPrefixLength) * kBitsPerEdgeChar;
  return static_cast<int>(ExtractPackedField(
      item, last_base_offset, kBitsPerEdgeChar));
}

void SeqToSdbg::BuildPackedTipLabel(const uint32_t *item,
                                    unsigned bucket_id,
                                    uint32_t *label) const {
  constexpr unsigned kPrefixBits =
      kBucketPrefixLength * kBitsPerEdgeChar;
  const unsigned label_words = DivCeiling(opt_.k, kCharsPerEdgeWord);
  for (unsigned i = 0; i < label_words; ++i) {
    const uint32_t left =
        i == 0 ? static_cast<uint32_t>(bucket_id << kPrefixBits)
               : item[i - 1u] << kPrefixBits;
    const uint32_t right =
        i < static_cast<unsigned>(words_per_substr_)
            ? item[i] >> (kBitsPerEdgeWord - kPrefixBits)
            : 0u;
    label[i] = left | right;
  }
  const unsigned used_in_last = opt_.k % kCharsPerEdgeWord;
  if (used_in_last != 0) {
    label[label_words - 1u] &=
        std::numeric_limits<uint32_t>::max()
        << ((kCharsPerEdgeWord - used_in_last) * kBitsPerEdgeChar);
  }
}

template <unsigned NWords, bool BucketPacked>
inline __attribute__((always_inline)) void SeqToSdbg::MaterializeItem(
    const SeqPackage::SeqView &seq_view, int offset, unsigned strand,
    uint32_t *item) const {
  const unsigned seq_len = seq_view.length();
  const unsigned num_chars_to_copy =
      opt_.k - (offset + static_cast<int>(opt_.k) >
                static_cast<int>(seq_len));
  constexpr unsigned omitted_prefix_chars =
      BucketPacked ? kBucketPrefixLength : 0u;
  constexpr unsigned stored_words = BucketPacked ? NWords - 1u : NWords;
  static_assert(!BucketPacked || NWords > 1,
                "a packed record must retain at least one word");
  assert(num_chars_to_copy >= omitted_prefix_chars);
  int counting = 0;
  if (offset > 0 && offset + static_cast<int>(opt_.k) <=
                        static_cast<int>(seq_len)) {
    counting = multiplicity[seq_view.id()];
  }

  auto ptr_and_offset = seq_view.raw_address();
  const unsigned start_offset = ptr_and_offset.second;
  const unsigned words_this_seq = DivCeiling(start_offset + seq_len, 16);
  const uint32_t *edge_p = ptr_and_offset.first;
  unsigned prev_char;

  if (strand == 0) {
    prev_char =
        offset == 0 ? kSentinelValue : seq_view.base_at(offset - 1);
    CopySubstring(item, edge_p,
                  offset + start_offset + omitted_prefix_chars,
                  num_chars_to_copy - omitted_prefix_chars, 1,
                  words_this_seq, stored_words);
  } else {
    prev_char = offset == 0
                    ? kSentinelValue
                    : 3 - seq_view.base_at(seq_len - offset);
    int copy_offset = static_cast<int>(seq_len) - 1 - offset -
                      (static_cast<int>(opt_.k) - 1);
    if (copy_offset < 0) {
      assert(num_chars_to_copy == opt_.k - 1);
      copy_offset = 0;
    }
    // Skipping the prefix of a reverse-complemented result removes the same
    // number of bases from the high end of its source interval.  The source
    // start therefore stays unchanged while the copied length shrinks.
    CopySubstringRC(item, edge_p, copy_offset + start_offset,
                    num_chars_to_copy - omitted_prefix_chars, 1,
                    words_this_seq, stored_words);
  }

  const uint32_t metadata =
      (unsigned(num_chars_to_copy == opt_.k)
       << (kBWTCharNumBits + kBitsPerMul)) |
      (prev_char << kBitsPerMul) |
      static_cast<unsigned>(std::max(0, kMaxMul - counting));
  if (BucketPacked) {
    AppendPackedMetadata(item, metadata);
  } else {
    item[NWords - 1u] |= metadata;
  }
}

template <unsigned NWords, bool BucketPacked>
inline __attribute__((always_inline)) void SeqToSdbg::MaterializeStreamedEdge(
    const uint32_t *edge, unsigned strand, uint32_t *item) const {
  constexpr unsigned omitted_prefix_chars =
      BucketPacked ? kBucketPrefixLength : 0u;
  constexpr unsigned stored_words = BucketPacked ? NWords - 1u : NWords;
  static_assert(!BucketPacked || NWords > 1,
                "a packed record must retain at least one word");
  const unsigned copied_chars = opt_.k - omitted_prefix_chars;
  unsigned previous_base;
  if (strand == 0) {
    previous_base = ExtractRawBase(edge, 0u);
    CopySubstring(item, edge, 1u + omitted_prefix_chars, copied_chars, 1,
                  stream_words_per_edge_, stored_words);
  } else {
    previous_base = 3u - ExtractRawBase(edge, opt_.k);
    // The real reverse-complement record begins at offset one.  After its
    // shared prefix is omitted, its source interval is the low end of the raw
    // edge and remains anchored at zero.
    CopySubstringRC(item, edge, 0u, copied_chars, 1,
                    stream_words_per_edge_, stored_words);
  }
  const unsigned counting = edge[stream_words_per_edge_ - 1u] & kMaxMul;
  const uint32_t metadata =
      (uint32_t{1} << (kBWTCharNumBits + kBitsPerMul)) |
      (previous_base << kBitsPerMul) | (kMaxMul - counting);
  if (BucketPacked) {
    AppendPackedMetadata(item, metadata);
  } else {
    item[NWords - 1u] |= metadata;
  }
}

template <unsigned NWords, bool BucketPacked>
void SeqToSdbg::Lv1FillStreamedEdgesFor(OffsetFiller &filler,
                                        int64_t chunk_from,
                                        int64_t chunk_to) {
  std::vector<uint32_t> records;
  const unsigned reverse_window_offset =
      stream_edge_length_ - 1u - kBucketPrefixLength;
  // The macro-descriptor path deliberately postpones materialization until
  // after a whole chunk has been partitioned.  It is ideal for the one-pass
  // case where every descriptor is live.  A bounded arena handles an
  // arbitrary rank-selected subset of biological buckets; scan that subset
  // directly so a descriptor can never outlive or be confused with a bucket
  // outside the current reusable-arena pass.
  bool handles_all_buckets = stream_input_unordered_;
  if (handles_all_buckets) {
    for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
      if (!filler.IsHandling(bucket)) {
        handles_all_buckets = false;
        break;
      }
    }
  }
  if (stream_input_unordered_ && handles_all_buckets) {
    // A full 16-bit radix partition emits mostly tiny runs for a 2 MiB input
    // chunk (about five records per bucket here), which turns into millions of
    // small copies.  Partition only by the high byte instead.  During the
    // final pass at most 256 bucket streams are active, keeping their current
    // destination pages and cache lines within the per-core TLB/cache budget,
    // while every record is materialized exactly once into the final array.
    constexpr unsigned kMacroBits = 8u;
    constexpr unsigned kBucketBits = 2u * kBucketPrefixLength;
    static_assert(kBucketBits >= kMacroBits, "invalid macro radix");
    constexpr unsigned kMacroShift = kBucketBits - kMacroBits;
    constexpr unsigned kNumMacros = 1u << kMacroBits;
    std::array<uint32_t, kNumMacros> macro_counts{};
    std::array<uint32_t, kNumMacros + 1u> macro_begin{};
    std::array<uint32_t, kNumMacros> macro_cursor{};
    std::vector<uint32_t> descriptors;

    for (int64_t chunk_id = chunk_from; chunk_id < chunk_to; ++chunk_id) {
      const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
      ReadStreamedEdgeChunk(chunk, &records);
      macro_counts.fill(0);

      const uint32_t *edge = records.data();
      for (uint32_t i = 0; i < chunk.num_records; ++i,
                    edge += stream_words_per_edge_) {
        const uint16_t key = ExtractRawBaseWindow8(edge, 1u);
        const uint16_t reverse_key = ReverseComplementWindow8(
            ExtractRawBaseWindow8(edge, reverse_window_offset));
        if (filler.IsHandling(key)) {
          ++macro_counts[key >> kMacroShift];
        }
        if (filler.IsHandling(reverse_key)) {
          ++macro_counts[reverse_key >> kMacroShift];
        }
      }

      macro_begin[0] = 0;
      for (unsigned macro = 0; macro < kNumMacros; ++macro) {
        macro_begin[macro + 1u] =
            macro_begin[macro] + macro_counts[macro];
        macro_cursor[macro] = macro_begin[macro];
      }
      descriptors.resize(macro_begin[kNumMacros]);

      edge = records.data();
      for (uint32_t i = 0; i < chunk.num_records; ++i,
                    edge += stream_words_per_edge_) {
        const uint16_t key = ExtractRawBaseWindow8(edge, 1u);
        const uint16_t reverse_key = ReverseComplementWindow8(
            ExtractRawBaseWindow8(edge, reverse_window_offset));
        if (filler.IsHandling(key)) {
          descriptors[macro_cursor[key >> kMacroShift]++] = i << 1u;
        }
        if (filler.IsHandling(reverse_key)) {
          descriptors[macro_cursor[reverse_key >> kMacroShift]++] =
              (i << 1u) | 1u;
        }
      }

      for (unsigned macro = 0; macro < kNumMacros; ++macro) {
        for (uint32_t j = macro_begin[macro];
             j < macro_begin[macro + 1u]; ++j) {
          const uint32_t descriptor = descriptors[j];
          const unsigned strand = descriptor & 1u;
          const uint32_t *source =
              records.data() +
              static_cast<size_t>(descriptor >> 1u) * stream_words_per_edge_;
          const uint16_t key =
              strand == 0u
                  ? ExtractRawBaseWindow8(source, 1u)
                  : ReverseComplementWindow8(ExtractRawBaseWindow8(
                        source, reverse_window_offset));
          assert((key >> kMacroShift) == macro);
          MaterializeStreamedEdge<NWords, BucketPacked>(
              source, strand, filler.ReserveNextItem(key));
        }
      }
    }
    return;
  }

  for (int64_t chunk_id = chunk_from; chunk_id < chunk_to; ++chunk_id) {
    const StreamEdgeChunk &chunk = stream_edge_chunks_[chunk_id];
    ReadStreamedEdgeChunk(chunk, &records);
    const uint32_t *edge = records.data();
    for (uint32_t i = 0; i < chunk.num_records; ++i,
                  edge += stream_words_per_edge_) {
      const uint16_t key = ExtractRawBaseWindow8(edge, 1u);
      const uint16_t reverse_key = ReverseComplementWindow8(
          ExtractRawBaseWindow8(edge, reverse_window_offset));
      if (filler.IsHandling(key)) {
        MaterializeStreamedEdge<NWords, BucketPacked>(
            edge, 0u, filler.ReserveNextItem(key));
      }
      if (filler.IsHandling(reverse_key)) {
        MaterializeStreamedEdge<NWords, BucketPacked>(
            edge, 1u, filler.ReserveNextItem(reverse_key));
      }
    }
  }
}

template <unsigned NWords, bool BucketPacked>
void SeqToSdbg::Lv1FillDirectItemsFor(OffsetFiller &filler, int64_t seq_from,
                                      int64_t seq_to) {
  const int64_t num_stream_chunks =
      static_cast<int64_t>(stream_edge_chunks_.size());
  if (seq_from < num_stream_chunks) {
    Lv1FillStreamedEdgesFor<NWords, BucketPacked>(
        filler, seq_from, std::min(seq_to, num_stream_chunks));
  }
  const int64_t loaded_from =
      std::max(seq_from, num_stream_chunks) - num_stream_chunks;
  const int64_t loaded_to =
      std::max(std::min(seq_to,
                        num_stream_chunks +
                            static_cast<int64_t>(seq_pkg_.seq_count())),
               num_stream_chunks) -
      num_stream_chunks;
  for (int64_t seq_id = loaded_from; seq_id < loaded_to; ++seq_id) {
    auto seq_view = seq_pkg_.GetSeqView(seq_id);
    const unsigned seq_len = seq_view.length();
    if (seq_len < opt_.k + 1) {
      continue;
    }

    Kmer<1, uint32_t> kmer, rev_kmer;
    auto ptr_and_offset = seq_view.raw_address();
    kmer.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second,
                     kBucketPrefixLength);
    auto rev_ptr_and_offset =
        seq_view.raw_address(seq_len - kBucketPrefixLength);
    rev_kmer.InitFromPtr(rev_ptr_and_offset.first, rev_ptr_and_offset.second,
                         kBucketPrefixLength);
    rev_kmer.ReverseComplement(kBucketPrefixLength);

    unsigned key = kmer.data()[0] >>
                   (kBitsPerEdgeWord -
                    kBucketPrefixLength * kBitsPerEdgeChar);
    unsigned rev_key = rev_kmer.data()[0] >>
                       (kBitsPerEdgeWord -
                        kBucketPrefixLength * kBitsPerEdgeChar);
    unsigned local_offset = 0;
    const unsigned final_offset = seq_len - opt_.k + 1u;
    if (!real_only_records_) {
      if (filler.IsHandling(key)) {
        MaterializeItem<NWords, BucketPacked>(seq_view, local_offset, 0,
                                              filler.ReserveNextItem(key));
      }
      if (filler.IsHandling(rev_key)) {
        MaterializeItem<NWords, BucketPacked>(
            seq_view, local_offset, 1, filler.ReserveNextItem(rev_key));
      }
    }

    for (int i = kBucketPrefixLength;
         i - (static_cast<int>(kBucketPrefixLength) - 1) + opt_.k - 1 <=
         seq_len;
         ++i) {
      key = (key * kBucketBase + seq_view.base_at(i)) % kNumBuckets;
      rev_key =
          (rev_key * kBucketBase + (3 - seq_view.base_at(seq_len - 1 - i))) %
          kNumBuckets;
      ++local_offset;
      if (!real_only_records_ || local_offset < final_offset) {
        if (filler.IsHandling(key)) {
          MaterializeItem<NWords, BucketPacked>(
              seq_view, local_offset, 0, filler.ReserveNextItem(key));
        }
        if (filler.IsHandling(rev_key)) {
          MaterializeItem<NWords, BucketPacked>(
              seq_view, local_offset, 1,
              filler.ReserveNextItem(rev_key));
        }
      }
    }
  }
}

template <unsigned NWords, bool PackedOffsets, bool BucketPacked>
void SeqToSdbg::Lv2ExtractSubStringFor(OffsetFetcher &fetcher,
                                       SubstrPtr substr) {
  // Bucket order makes locators monotone within each producer partition, but
  // consecutive records can still land on different packed-sequence cache
  // lines.  Keep an independent read-only cursor ahead of the consumer so the
  // sequence words and multiplicity arrive while earlier records are decoded
  // and copied.  This is layout-driven and applies to every packed-locator
  // input; the legacy global-offset path remains unchanged.
  static const unsigned kPrefetchDistance = []() {
    const char *value = std::getenv("MEGAHIT_SEQ_READ_PREFETCH_DISTANCE");
    if (value == nullptr) {
      return 24u;
    }
    char *end = nullptr;
    const unsigned long parsed = std::strtoul(value, &end, 10);
    return end != value && *end == '\0' && parsed <= 128
               ? static_cast<unsigned>(parsed)
               : 24u;
  }();
  static const bool prefetch_packed_sequences =
      PackedOffsets && kPrefetchDistance != 0 &&
      std::getenv("MEGAHIT_DISABLE_SEQ_READ_PREFETCH") == nullptr;
  OffsetFetcher prefetcher = fetcher;
  auto prefetch_sequence = [&](uint64_t locator) {
    const auto decoded = SeqLocatorDecoder<PackedOffsets>::Decode(
        locator, seq_locator_shift_, seq_offset_mask_, seq_pkg_);
    const auto &seq_view = decoded.seq_view;
    const auto raw = seq_view.raw_address();
    const unsigned first_word = raw.second / kCharsPerEdgeWord;
    const unsigned last_word =
        (raw.second + seq_view.length() - 1u) / kCharsPerEdgeWord;
    __builtin_prefetch(raw.first + first_word, 0, 1);
    __builtin_prefetch(raw.first + last_word, 0, 1);
    __builtin_prefetch(multiplicity.data() + seq_view.id(), 0, 1);
  };

  if (prefetch_packed_sequences) {
    for (unsigned i = 0;
         i < kPrefetchDistance && prefetcher.HasNext(); ++i) {
      prefetch_sequence(prefetcher.Next());
    }
  }

  while (fetcher.HasNext()) {
    if (prefetch_packed_sequences && prefetcher.HasNext()) {
      prefetch_sequence(prefetcher.Next());
    }
    const uint64_t locator = fetcher.Next();
    uint32_t *const output_item = &*substr;
    const auto decoded =
        SeqLocatorDecoder<PackedOffsets>::Decode(locator, seq_locator_shift_,
                                                  seq_offset_mask_, seq_pkg_);
    MaterializeItem<NWords, BucketPacked>(
        decoded.seq_view, decoded.offset, decoded.strand, output_item);
    substr += words_per_substr_;
  }
}

void SeqToSdbg::Lv2ExtractSubString(OffsetFetcher &fetcher, SubstrPtr substr) {
#define DISPATCH_LV2_WORDS(N)                                         \
  case N:                                                             \
    if (packed_seq_offsets_) {                                        \
      if (bucket_packed_records_) {                                   \
        return Lv2ExtractSubStringFor<N, true, true>(fetcher, substr); \
      }                                                               \
      return Lv2ExtractSubStringFor<N, true, false>(fetcher, substr);  \
    }                                                                 \
    if (bucket_packed_records_) {                                     \
      return Lv2ExtractSubStringFor<N, false, true>(fetcher, substr);  \
    }                                                                 \
    return Lv2ExtractSubStringFor<N, false, false>(fetcher, substr)

  switch (full_words_per_substr_) {
    case 1:
      if (packed_seq_offsets_) {
        return Lv2ExtractSubStringFor<1, true, false>(fetcher, substr);
      }
      return Lv2ExtractSubStringFor<1, false, false>(fetcher, substr);
    DISPATCH_LV2_WORDS(2);
    DISPATCH_LV2_WORDS(3);
    DISPATCH_LV2_WORDS(4);
    DISPATCH_LV2_WORDS(5);
    DISPATCH_LV2_WORDS(6);
    DISPATCH_LV2_WORDS(7);
    DISPATCH_LV2_WORDS(8);
    DISPATCH_LV2_WORDS(9);
    DISPATCH_LV2_WORDS(10);
    DISPATCH_LV2_WORDS(11);
    DISPATCH_LV2_WORDS(12);
    DISPATCH_LV2_WORDS(13);
    DISPATCH_LV2_WORDS(14);
    DISPATCH_LV2_WORDS(15);
    DISPATCH_LV2_WORDS(16);
    DISPATCH_LV2_WORDS(17);
    default:
      xfatal("Invalid number of full substring words: {}\n",
             full_words_per_substr_);
  }
#undef DISPATCH_LV2_WORDS
}

bool SeqToSdbg::Lv2DeferPostprocess() const {
  return real_only_records_ &&
         (!precomputed_endpoint_tips_ ||
          std::getenv("MEGAHIT_VALIDATE_PRECOMPUTED_ENDPOINTS") != nullptr);
}

void SeqToSdbg::Lv2PostprocessDeferred() {
#define DISPATCH_REAL_ONLY(N)                                  \
  case N:                                                      \
    if (bucket_packed_records_) {                              \
      return Lv2PostprocessRealOnly<N, true>();                \
    }                                                          \
    return Lv2PostprocessRealOnly<N, false>()
  switch (full_words_per_substr_) {
    case 1:
      return Lv2PostprocessRealOnly<1, false>();
    DISPATCH_REAL_ONLY(2);
    DISPATCH_REAL_ONLY(3);
    DISPATCH_REAL_ONLY(4);
    DISPATCH_REAL_ONLY(5);
    DISPATCH_REAL_ONLY(6);
    DISPATCH_REAL_ONLY(7);
    DISPATCH_REAL_ONLY(8);
    DISPATCH_REAL_ONLY(9);
    DISPATCH_REAL_ONLY(10);
    DISPATCH_REAL_ONLY(11);
    DISPATCH_REAL_ONLY(12);
    DISPATCH_REAL_ONLY(13);
    DISPATCH_REAL_ONLY(14);
    DISPATCH_REAL_ONLY(15);
    DISPATCH_REAL_ONLY(16);
    DISPATCH_REAL_ONLY(17);
    default:
      xfatal("Invalid number of real-only substring words: {}\n",
             full_words_per_substr_);
  }
#undef DISPATCH_REAL_ONLY
}

template <unsigned NWords, bool BucketPacked>
void SeqToSdbg::Lv2PostprocessRealOnly() {
  static constexpr unsigned kStoredWords =
      BucketPacked ? NWords - 1u : NWords;
  static_assert(kStoredWords > 0, "real-only records cannot be empty");

  struct TipRecord {
    uint16_t bucket;
    uint16_t reserved;
    uint32_t data[kStoredWords];
  };
  struct PairEntry {
    int a;
    int b;
    int counting;
    const uint32_t *label;
  };

  const unsigned label_words = DivCeiling(opt_.k, kCharsPerEdgeWord);
  const auto mask_label = [&](uint32_t *label) {
    const unsigned used = opt_.k % kCharsPerEdgeWord;
    if (used != 0) {
      label[label_words - 1u] &=
          std::numeric_limits<uint32_t>::max()
          << ((kCharsPerEdgeWord - used) * kBitsPerEdgeChar);
    }
  };
  const auto build_full_label = [&](const uint32_t *item,
                                    unsigned bucket,
                                    uint32_t *label) {
    std::fill(label, label + label_words, 0u);
    if (BucketPacked) {
      BuildPackedTipLabel(item, bucket, label);
    } else {
      std::copy(item, item + label_words, label);
      mask_label(label);
    }
  };
  const auto group_diff = [&](const uint32_t *lhs, const uint32_t *rhs) {
    return BucketPacked
               ? IsDiffPackedKMinusOneMer(lhs, rhs)
               : IsDiffKMinusOneMer(const_cast<uint32_t *>(lhs),
                                    const_cast<uint32_t *>(rhs), 1, opt_.k);
  };
  const auto extract_a = [&](const uint32_t *item) {
    return BucketPacked
               ? ExtractPackedA(item)
               : Extract_a(const_cast<uint32_t *>(item), kStoredWords, 1,
                           opt_.k);
  };
  const auto extract_b = [&](const uint32_t *item) {
    return BucketPacked
               ? static_cast<int>((ExtractPackedMetadata(item) >>
                                   kBitsPerMul) &
                                  ((1u << kBWTCharNumBits) - 1u))
               : Extract_b(const_cast<uint32_t *>(item), kStoredWords, 1);
  };
  const auto extract_counting = [&](const uint32_t *item) {
    return BucketPacked
               ? static_cast<int>(ExtractPackedMetadata(item) & kMaxMul)
               : ExtractCounting(const_cast<uint32_t *>(item),
                                 kStoredWords, 1);
  };
  const auto node_compare = [&](const uint32_t *lhs, const uint32_t *rhs) {
    for (unsigned i = 0; i < label_words; ++i) {
      if (lhs[i] < rhs[i]) {
        return -1;
      }
      if (lhs[i] > rhs[i]) {
        return 1;
      }
    }
    return 0;
  };
  const auto packed_dest_prefix_compare =
      [&](const uint32_t *dest, const uint32_t *source,
          unsigned source_bucket) {
        assert(BucketPacked && bucket_packed_base_bits_ >= 2u);
        const unsigned words =
            DivCeiling(bucket_packed_base_bits_, kBitsPerEdgeWord);
        const unsigned tail_bits =
            bucket_packed_base_bits_ % kBitsPerEdgeWord;
        for (unsigned i = 0; i < words; ++i) {
          uint32_t prefix_word =
              i == 0
                  ? ((source_bucket & 3u) << 30u) | (source[0] >> 2u)
                  : (source[i - 1u] << 30u) | (source[i] >> 2u);
          uint32_t dest_word = dest[i];
          if (i + 1u == words && tail_bits != 0) {
            const uint32_t mask = std::numeric_limits<uint32_t>::max()
                                  << (kBitsPerEdgeWord - tail_bits);
            prefix_word &= mask;
            dest_word &= mask;
          }
          if (dest_word < prefix_word) {
            return -1;
          }
          if (dest_word > prefix_word) {
            return 1;
          }
        }
        return 0;
      };
  const auto record_less = [&](const uint32_t *lhs, const uint32_t *rhs) {
    for (unsigned i = 0; i < kStoredWords; ++i) {
      if (lhs[i] < rhs[i]) {
        return true;
      }
      if (lhs[i] > rhs[i]) {
        return false;
      }
    }
    return false;
  };

  const unsigned n_threads = std::max(1, opt_.n_threads);
  std::vector<std::vector<TipRecord>> thread_tips(n_threads);
  std::vector<uint32_t> thread_bucket_counts(
      static_cast<size_t>(n_threads) * kNumBuckets, 0u);
  std::vector<uint64_t> thread_source_tips(n_threads, 0u);
  std::vector<uint64_t> thread_sink_tips(n_threads, 0u);

  const auto append_tip = [&](const uint32_t *node, bool source,
                              unsigned tid) {
    TipRecord tip{};
    uint32_t bases[32]{};
    uint32_t metadata;
    if (source) {
      std::copy(node, node + label_words, bases);
      tip.bucket = static_cast<uint16_t>(bases[0] >> 16u);
      metadata = (uint32_t{1} << (kBWTCharNumBits + kBitsPerMul)) |
                 (uint32_t{kSentinelValue} << kBitsPerMul) | kMaxMul;
      ++thread_source_tips[tid];
    } else {
      const uint32_t first_base = node[0] >> 30u;
      for (unsigned i = 0; i < label_words; ++i) {
        bases[i] = node[i] << kBitsPerEdgeChar;
        if (i + 1u < label_words) {
          bases[i] |= node[i + 1u] >>
                      (kBitsPerEdgeWord - kBitsPerEdgeChar);
        }
      }
      mask_label(bases);
      tip.bucket = static_cast<uint16_t>(bases[0] >> 16u);
      metadata = (first_base << kBitsPerMul) | kMaxMul;
      ++thread_sink_tips[tid];
    }

    if (BucketPacked) {
      constexpr unsigned kPrefixBits =
          kBucketPrefixLength * kBitsPerEdgeChar;
      for (unsigned i = 0; i < kStoredWords; ++i) {
        tip.data[i] = bases[i] << kPrefixBits;
        if (i + 1u < label_words) {
          tip.data[i] |=
              bases[i + 1u] >> (kBitsPerEdgeWord - kPrefixBits);
        }
      }
      AppendPackedMetadata(tip.data, metadata);
    } else {
      std::copy(bases, bases + label_words, tip.data);
      tip.data[kStoredWords - 1u] |= metadata;
    }

    auto &count = thread_bucket_counts[
        static_cast<size_t>(tid) * kNumBuckets + tip.bucket];
    if (count == std::numeric_limits<uint32_t>::max()) {
      xfatal("Too many deferred sentinel records in bucket {}\n",
             tip.bucket);
    }
    ++count;
    thread_tips[tid].push_back(tip);
  };

  xinfo("Computing exact prefix/suffix k-mer set difference...\n");
  SimpleTimer endpoint_timer;
  endpoint_timer.start();
  omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(dynamic, 1)
  for (int output_suffix_int = 0; output_suffix_int < (1 << 14);
       ++output_suffix_int) {
    const unsigned output_suffix =
        static_cast<unsigned>(output_suffix_int);
    const unsigned source_prefix = output_suffix << 2u;
    const unsigned tid = omp_get_thread_num();

    struct DestState {
      unsigned bucket;
      const uint32_t *data;
      int64_t size;
      int64_t index;
      bool valid;
      const uint32_t *item;
      uint32_t node[32];
    };
    DestState destinations[4]{};
    const auto advance_dest = [&](unsigned b) {
      DestState &dest = destinations[b];
      if (dest.index >= dest.size) {
        dest.valid = false;
        return;
      }
      const uint32_t *item = dest.data + dest.index * kStoredWords;
      dest.item = item;
      const int a = extract_a(item);
      if (!BucketPacked) {
        build_full_label(item, dest.bucket, dest.node);
      }
      ++dest.index;
      while (dest.index < dest.size) {
        const uint32_t *next = dest.data + dest.index * kStoredWords;
        if (group_diff(item, next) || extract_a(next) != a) {
          break;
        }
        ++dest.index;
      }
      dest.valid = true;
    };
    for (unsigned b = 0; b < 4u; ++b) {
      DestState &dest = destinations[b];
      dest.bucket = (b << 14u) | output_suffix;
      dest.data = SortedDirectBucketData(dest.bucket);
      dest.size = SortedDirectBucketSize(dest.bucket);
      advance_dest(b);
    }

    uint32_t scratch_label[32]{};
    uint32_t prefix_node[32]{};
    const auto append_sink = [&](DestState &dest) {
      if (BucketPacked) {
        build_full_label(dest.item, dest.bucket, dest.node);
      }
      append_tip(dest.node, false, tid);
    };
    for (unsigned source_suffix = 0; source_suffix < 4u;
         ++source_suffix) {
      const unsigned source_bucket = source_prefix | source_suffix;
      const uint32_t *source_data = SortedDirectBucketData(source_bucket);
      const int64_t source_size = SortedDirectBucketSize(source_bucket);
      int64_t source_index = 0;
      while (source_index < source_size) {
        const uint32_t *group_item =
            source_data + source_index * kStoredWords;
        uint32_t b_mask = 0;
        int64_t group_end = source_index;
        do {
          const int b =
              extract_b(source_data + group_end * kStoredWords);
          assert(b >= 0 && b < 4);
          b_mask |= uint32_t{1} << b;
          ++group_end;
        } while (group_end < source_size &&
                 !group_diff(group_item,
                             source_data + group_end * kStoredWords));
        source_index = group_end;

        bool source_label_ready = false;
        const auto build_prefix_node = [&](unsigned b) {
          if (!source_label_ready) {
            build_full_label(group_item, source_bucket, scratch_label);
            source_label_ready = true;
          }
          prefix_node[0] = (b << 30u) | (scratch_label[0] >> 2u);
          for (unsigned i = 1; i < label_words; ++i) {
            prefix_node[i] = (scratch_label[i - 1u] << 30u) |
                             (scratch_label[i] >> 2u);
          }
          mask_label(prefix_node);
        };
        for (unsigned b = 0; b < 4u; ++b) {
          if ((b_mask & (uint32_t{1} << b)) == 0) {
            continue;
          }
          if (!BucketPacked) {
            build_prefix_node(b);
          }
          DestState &dest = destinations[b];
          int comparison =
              dest.valid
                  ? (BucketPacked
                         ? packed_dest_prefix_compare(dest.item, group_item,
                                                      source_bucket)
                         : node_compare(dest.node, prefix_node))
                  : 1;
          while (dest.valid && comparison < 0) {
            append_sink(dest);
            advance_dest(b);
            if (dest.valid) {
              comparison =
                  BucketPacked
                      ? packed_dest_prefix_compare(dest.item, group_item,
                                                   source_bucket)
                      : node_compare(dest.node, prefix_node);
            }
          }
          if (dest.valid && comparison == 0) {
            advance_dest(b);
          } else {
            if (BucketPacked) {
              build_prefix_node(b);
            }
            assert((prefix_node[0] >> 16u) == dest.bucket);
            append_tip(prefix_node, true, tid);
          }
        }
      }
    }
    for (unsigned b = 0; b < 4u; ++b) {
      DestState &dest = destinations[b];
      while (dest.valid) {
        append_sink(dest);
        advance_dest(b);
      }
    }
  }
  endpoint_timer.stop();

  uint64_t source_tips = 0;
  uint64_t sink_tips = 0;
  for (unsigned tid = 0; tid < n_threads; ++tid) {
    source_tips += thread_source_tips[tid];
    sink_tips += thread_sink_tips[tid];
  }
  const uint64_t total_tips = source_tips + sink_tips;
  xinfo("Exact endpoint difference: {} source, {} sink, time {.4}\n",
        source_tips, sink_tips, endpoint_timer.elapsed());

  std::vector<uint64_t> tip_bucket_begin(kNumBuckets + 1u, 0u);
  for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
    uint64_t count = 0;
    for (unsigned tid = 0; tid < n_threads; ++tid) {
      count += thread_bucket_counts[
          static_cast<size_t>(tid) * kNumBuckets + bucket];
    }
    tip_bucket_begin[bucket + 1u] = tip_bucket_begin[bucket] + count;
  }
  assert(tip_bucket_begin.back() == total_tips);
  if (total_tips >
      std::numeric_limits<size_t>::max() / kStoredWords) {
    xfatal("Deferred sentinel array is too large\n");
  }
  std::vector<uint32_t> tip_items(
      static_cast<size_t>(total_tips) * kStoredWords);
  std::vector<uint64_t> thread_bucket_cursor(
      static_cast<size_t>(n_threads) * kNumBuckets, 0u);
  for (unsigned bucket = 0; bucket < kNumBuckets; ++bucket) {
    uint64_t cursor = tip_bucket_begin[bucket];
    for (unsigned tid = 0; tid < n_threads; ++tid) {
      const size_t index = static_cast<size_t>(tid) * kNumBuckets + bucket;
      thread_bucket_cursor[index] = cursor;
      cursor += thread_bucket_counts[index];
    }
  }

#pragma omp parallel for schedule(static)
  for (int tid_int = 0; tid_int < static_cast<int>(n_threads); ++tid_int) {
    const unsigned tid = static_cast<unsigned>(tid_int);
    uint64_t *cursor = thread_bucket_cursor.data() +
                       static_cast<size_t>(tid) * kNumBuckets;
    for (const TipRecord &tip : thread_tips[tid]) {
      const uint64_t index = cursor[tip.bucket]++;
      std::copy(tip.data, tip.data + kStoredWords,
                tip_items.data() + index * kStoredWords);
    }
    std::vector<TipRecord>().swap(thread_tips[tid]);
  }
  thread_tips.clear();
  thread_tips.shrink_to_fit();
  thread_bucket_counts.clear();
  thread_bucket_counts.shrink_to_fit();
  thread_bucket_cursor.clear();
  thread_bucket_cursor.shrink_to_fit();

  auto tip_sort = SelectSortingFunc(kStoredWords, 0,
                                    Lv2SortIgnoredLowBytes(),
                                    Lv2SortIgnoredHighBytes());
#pragma omp parallel for schedule(dynamic, 16)
  for (int bucket_int = 0; bucket_int < static_cast<int>(kNumBuckets);
       ++bucket_int) {
    const unsigned bucket = static_cast<unsigned>(bucket_int);
    const uint64_t begin = tip_bucket_begin[bucket];
    const uint64_t size = tip_bucket_begin[bucket + 1u] - begin;
    if (size > 1) {
      tip_sort(tip_items.data() + begin * kStoredWords, size);
    }
  }

  if (precomputed_endpoint_tips_ &&
      std::getenv("MEGAHIT_VALIDATE_PRECOMPUTED_ENDPOINTS") != nullptr) {
    uint64_t missing_precomputed = 0;
#pragma omp parallel for schedule(static) reduction(+ : missing_precomputed)
    for (int bucket_int = 0;
         bucket_int < static_cast<int>(kNumBuckets); ++bucket_int) {
      const unsigned bucket = static_cast<unsigned>(bucket_int);
      uint64_t expected = tip_bucket_begin[bucket];
      const uint64_t expected_end = tip_bucket_begin[bucket + 1u];
      uint64_t candidate = endpoint_tip_bucket_begin_[bucket];
      const uint64_t candidate_end =
          endpoint_tip_bucket_begin_[bucket + 1u];
      while (expected < expected_end) {
        const uint32_t *expected_item =
            tip_items.data() + expected * kStoredWords;
        while (candidate < candidate_end) {
          const uint32_t *candidate_item =
              endpoint_tip_items_.data() + candidate * kStoredWords;
          if (!record_less(candidate_item, expected_item)) break;
          ++candidate;
        }
        if (candidate >= candidate_end ||
            record_less(expected_item,
                        endpoint_tip_items_.data() +
                            candidate * kStoredWords)) {
          ++missing_precomputed;
        }
        ++expected;
      }
    }
    xinfo("Precomputed endpoint containment check: {} of {} exact global "
          "sentinels missing\n",
          missing_precomputed, total_tips);
    if (missing_precomputed != 0) {
      uint64_t printed = 0;
      for (unsigned bucket = 0;
           bucket < kNumBuckets && printed < 64u; ++bucket) {
        uint64_t expected = tip_bucket_begin[bucket];
        const uint64_t expected_end = tip_bucket_begin[bucket + 1u];
        uint64_t candidate = endpoint_tip_bucket_begin_[bucket];
        const uint64_t candidate_end =
            endpoint_tip_bucket_begin_[bucket + 1u];
        while (expected < expected_end && printed < 64u) {
          const uint32_t *expected_item =
              tip_items.data() + expected * kStoredWords;
          while (candidate < candidate_end &&
                 record_less(endpoint_tip_items_.data() +
                                 candidate * kStoredWords,
                             expected_item)) {
            ++candidate;
          }
          if (candidate >= candidate_end ||
              record_less(expected_item,
                          endpoint_tip_items_.data() +
                              candidate * kStoredWords)) {
            xinfo("Missing endpoint sentinel bucket {} item", bucket);
            for (unsigned word = 0; word < kStoredWords; ++word) {
              xinfoc(" {}", expected_item[word]);
            }
            xinfoc("{s}", "\n");
            ++printed;
          }
          ++expected;
        }
      }
    }
  }

  xinfo("Merging real records with necessary sentinels...\n");
  SimpleTimer merge_timer;
  merge_timer.start();
#pragma omp parallel for schedule(dynamic, 1)
  for (int bucket_int = 0; bucket_int < static_cast<int>(kNumBuckets);
       ++bucket_int) {
    const unsigned bucket = static_cast<unsigned>(bucket_int);
    const unsigned tid = omp_get_thread_num();
    const uint32_t *real = SortedDirectBucketData(bucket);
    const int64_t n_real = SortedDirectBucketSize(bucket);
    const uint64_t tip_begin = tip_bucket_begin[bucket];
    const uint32_t *tips = tip_items.data() + tip_begin * kStoredWords;
    const int64_t n_tips = static_cast<int64_t>(
        tip_bucket_begin[bucket + 1u] - tip_begin);
    if (n_real == 0 && n_tips == 0) {
      continue;
    }

    int64_t real_index = 0;
    int64_t tip_index = 0;
    const auto peek = [&](bool *is_tip) -> const uint32_t * {
      if (real_index >= n_real) {
        *is_tip = true;
        return tips + tip_index * kStoredWords;
      }
      if (tip_index >= n_tips) {
        *is_tip = false;
        return real + real_index * kStoredWords;
      }
      const uint32_t *real_item = real + real_index * kStoredWords;
      const uint32_t *tip_item = tips + tip_index * kStoredWords;
      *is_tip = record_less(tip_item, real_item);
      return *is_tip ? tip_item : real_item;
    };
    const auto pop = [&](bool *is_tip) -> const uint32_t * {
      const uint32_t *item = peek(is_tip);
      if (*is_tip) {
        ++tip_index;
      } else {
        ++real_index;
      }
      return item;
    };

    SdbgWriter::Snapshot snapshot;
    uint32_t tip_label[32]{};
    while (real_index < n_real || tip_index < n_tips) {
      bool from_tip = false;
      const uint32_t *group_item = pop(&from_tip);
      PairEntry pairs[25];
      unsigned n_pairs = 0;

      const auto add_pair = [&](const uint32_t *item) {
        const int a = extract_a(item);
        const int b = extract_b(item);
        const int counting = extract_counting(item);
        if (n_pairs != 0 && pairs[n_pairs - 1u].a == a &&
            pairs[n_pairs - 1u].b == b) {
          if (counting < pairs[n_pairs - 1u].counting) {
            pairs[n_pairs - 1u].counting = counting;
            pairs[n_pairs - 1u].label = item;
          }
          return;
        }
        assert(n_pairs < 25u);
        pairs[n_pairs++] = {a, b, counting, item};
      };
      add_pair(group_item);

      while (real_index < n_real || tip_index < n_tips) {
        bool next_is_tip = false;
        const uint32_t *next = peek(&next_is_tip);
        if (group_diff(group_item, next)) {
          break;
        }
        pop(&next_is_tip);
        add_pair(next);
      }

      uint32_t outputed_b = 0;
      for (unsigned i = 0; i < n_pairs; ++i) {
        const PairEntry &pair = pairs[i];
        const bool is_dollar = pair.a == kSentinelValue;
        const int w = pair.b == kSentinelValue
                          ? 0
                          : ((outputed_b & (1u << pair.b))
                                 ? pair.b + 5
                                 : pair.b + 1);
        const int last =
            pair.a == kSentinelValue
                ? 0
                : (i + 1u == n_pairs || pairs[i + 1u].a != pair.a);
        outputed_b |= 1u << pair.b;
        uint32_t *label = const_cast<uint32_t *>(pair.label);
        if (is_dollar && BucketPacked) {
          BuildPackedTipLabel(pair.label, bucket, tip_label);
          label = tip_label;
        }
        sdbg_writer_.Write(tid, bucket, w, last, is_dollar,
                           kMaxMul - pair.counting, label, &snapshot);
      }
    }
    sdbg_writer_.SaveSnapshot(snapshot);
  }
  merge_timer.stop();
  xinfo("Real/sentinel merge time: {.4}\n", merge_timer.elapsed());
}

void SeqToSdbg::Lv2Postprocess(int64_t from, int64_t to, int tid,
                               uint32_t *substr, unsigned bucket_id) {
  std::vector<uint32_t> merged_items;
  if (precomputed_endpoint_tips_) {
    assert(bucket_id + 1u < endpoint_tip_bucket_begin_.size());
    const uint64_t tip_begin = endpoint_tip_bucket_begin_[bucket_id];
    const uint64_t tip_end = endpoint_tip_bucket_begin_[bucket_id + 1u];
    const int64_t num_real = to - from;
    const int64_t num_tips = static_cast<int64_t>(tip_end - tip_begin);
    if (num_tips != 0) {
      merged_items.resize(
          static_cast<size_t>(num_real + num_tips) * words_per_substr_);
      const uint32_t *real =
          num_real == 0 ? nullptr : substr + from * words_per_substr_;
      const uint32_t *tips = endpoint_tip_items_.data() +
                             tip_begin * words_per_substr_;
      const auto less = [&](const uint32_t *lhs, const uint32_t *rhs) {
        for (int64_t word = 0; word < words_per_substr_; ++word) {
          if (lhs[word] < rhs[word]) return true;
          if (lhs[word] > rhs[word]) return false;
        }
        return false;
      };
      int64_t real_index = 0;
      int64_t tip_index = 0;
      int64_t output = 0;
      while (real_index < num_real || tip_index < num_tips) {
        const bool take_tip =
            real_index >= num_real ||
            (tip_index < num_tips &&
             less(tips + tip_index * words_per_substr_,
                  real + real_index * words_per_substr_));
        const uint32_t *item =
            take_tip ? tips + tip_index++ * words_per_substr_
                     : real + real_index++ * words_per_substr_;
        std::copy(item, item + words_per_substr_,
                  merged_items.data() + output * words_per_substr_);
        ++output;
      }
      substr = merged_items.data();
      from = 0;
      to = num_real + num_tips;
    }
  }

  int64_t start_idx, end_idx;
  int has_solid_a = 0;  // has solid (k+1)-mer aSb
  int has_solid_b = 0;  // has solid aSb
  int64_t last_a[4], outputed_b;
  uint32_t tip_label[32];
  SdbgWriter::Snapshot snapshot;

  const auto is_diff_k_minus_one = [&](const uint32_t *lhs,
                                       const uint32_t *rhs) {
    return bucket_packed_records_
               ? IsDiffPackedKMinusOneMer(lhs, rhs)
               : IsDiffKMinusOneMer(const_cast<uint32_t *>(lhs),
                                    const_cast<uint32_t *>(rhs), 1, opt_.k);
  };
  const auto extract_a = [&](const uint32_t *item) {
    return bucket_packed_records_
               ? ExtractPackedA(item)
               : Extract_a(const_cast<uint32_t *>(item), words_per_substr_, 1,
                           opt_.k);
  };
  const auto extract_b = [&](const uint32_t *item) {
    return bucket_packed_records_
               ? static_cast<int>((ExtractPackedMetadata(item) >> kBitsPerMul) &
                                  ((1u << kBWTCharNumBits) - 1u))
               : Extract_b(const_cast<uint32_t *>(item), words_per_substr_, 1);
  };
  const auto extract_counting = [&](const uint32_t *item) {
    return bucket_packed_records_
               ? static_cast<int>(ExtractPackedMetadata(item) & kMaxMul)
               : ExtractCounting(const_cast<uint32_t *>(item),
                                 words_per_substr_, 1);
  };

  for (start_idx = from; start_idx < to; start_idx = end_idx) {
    end_idx = start_idx + 1;
    uint32_t *item = substr + start_idx * words_per_substr_;

    while (end_idx < to &&
           !is_diff_k_minus_one(item,
                                substr + end_idx * words_per_substr_)) {
      ++end_idx;
    }

    // clean marking
    has_solid_a = has_solid_b = 0;
    outputed_b = 0;

    for (int64_t i = start_idx; i < end_idx; ++i) {
      uint32_t *cur_item = substr + i * words_per_substr_;
      int a = extract_a(cur_item);
      int b = extract_b(cur_item);

      if (a != kSentinelValue && b != kSentinelValue) {
        has_solid_a |= 1 << a;
        has_solid_b |= 1 << b;
      }

      if (a != kSentinelValue &&
          (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
        last_a[a] = i;
      }
    }

    for (int64_t i = start_idx, j; i < end_idx; i = j) {
      uint32_t *cur_item = substr + i * words_per_substr_;
      int a = extract_a(cur_item);
      int b = extract_b(cur_item);

      j = i + 1;

      while (j < end_idx) {
        uint32_t *next_item = substr + j * words_per_substr_;

        if (extract_a(next_item) != a || extract_b(next_item) != b) {
          break;
        } else {
          ++j;
        }
      }

      int w, last, is_dollar = 0;

      if (a == kSentinelValue) {
        assert(b != kSentinelValue);

        if (has_solid_b & (1 << b)) {
          continue;
        }

        is_dollar = 1;
      }

      if (b == kSentinelValue) {
        assert(a != kSentinelValue);

        if (has_solid_a & (1 << a)) {
          continue;
        }
      }

      w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
      last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);
      outputed_b |= 1 << b;

      uint32_t *label = cur_item;
      if (is_dollar && bucket_packed_records_) {
        BuildPackedTipLabel(cur_item, bucket_id, tip_label);
        label = tip_label;
      }
      sdbg_writer_.Write(tid, bucket_id, w, last, is_dollar,
                         kMaxMul - extract_counting(cur_item), label,
                         &snapshot);
    }
  }
  sdbg_writer_.SaveSnapshot(snapshot);
}

void SeqToSdbg::Lv0Postprocess() {
  sdbg_writer_.Finalize();
  xinfo("Number of $ A C G T A- C- G- T-:\n");
  xinfo("");
  for (int i = 0; i < 9; ++i) {
    xinfoc("{} ", sdbg_writer_.final_meta().w_count(i));
  }

  xinfoc("{s}", "\n");
  xinfo("Total number of edges: {}\n", sdbg_writer_.final_meta().item_count());
  xinfo("Total number of ONEs: {}\n", sdbg_writer_.final_meta().ones_in_last());
  xinfo("Total number of $v edges: {}\n",
        sdbg_writer_.final_meta().tip_count());

  assert(sdbg_writer_.final_meta().w_count(0) ==
         sdbg_writer_.final_meta().tip_count());
}
