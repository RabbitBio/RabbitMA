//
// Created by vout on 11/23/18.
//

#ifndef MEGAHIT_JUNCION_INDEX_H
#define MEGAHIT_JUNCION_INDEX_H

#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <sdbg/sdbg_def.h>
#include <utils/utils.h>
#include "parallel_hashmap/phmap.h"
#include "sequence/kmer_plus.h"
#include "sequence/sequence_package.h"
#include "utils/mutex.h"

template <class KmerType>
class ContigFlankIndex {
 public:
  struct FlankInfo {
    uint64_t ext_seq : 58;
    unsigned ext_len : 6;
    float mul;
  } __attribute__((packed));
  struct StoredFlankInfo {
    uint64_t ext_seq : 56;
    unsigned ext_len : 5;
    unsigned orientation : 1;
    unsigned has_secondary : 1;
    unsigned reserved : 1;
    float mul;
  } __attribute__((packed));
  using Flank = KmerPlus<KmerType, FlankInfo>;
  using HashMap =
      phmap::flat_hash_map<KmerType, StoredFlankInfo, KmerHash>;
  using SecondaryMap =
      phmap::flat_hash_map<KmerType, FlankInfo, KmerHash>;
  static_assert(sizeof(StoredFlankInfo) == sizeof(FlankInfo),
                "canonical flank metadata must remain compact");

  struct PendingFlank {
    KmerType canonical;
    FlankInfo info;
    unsigned orientation;
  };

  struct alignas(64) HashShard {
    HashMap primary;
    SecondaryMap secondary;
  };

  // Bloom membership does not need to share phmap's comparatively expensive
  // XXH3 hash.  A deterministic one-multiply word combiner is sufficient for
  // a probabilistic prefilter: collisions can only cause extra table probes,
  // never false negatives, because construction and lookup use the same key.
  static size_t FilterHash(const KmerType &kmer) {
    const typename KmerType::word_type *words = kmer.data();
    uint64_t hash = UINT64_C(0x9e3779b97f4a7c15);
    for (unsigned i = 0; i < KmerType::kNumWords; ++i) {
      hash ^= static_cast<uint64_t>(words[i]) +
              UINT64_C(0x9e3779b97f4a7c15) + (hash << 6u) + (hash >> 2u);
    }
    hash ^= hash >> 29u;
    hash *= UINT64_C(0x165667919e3779f9);
    hash ^= hash >> 32u;
    return static_cast<size_t>(hash);
  }

  // Most read windows are absent from the flank index.  A one-word blocked
  // Bloom filter rejects those windows with one compact, cache-friendly load
  // before touching the much larger hash table.  Four bits in the selected
  // word give a low false-positive rate while retaining the essential Bloom
  // property: an indexed k-mer is never rejected.
  class BlockedBloomFilter {
   public:
    void Build(const std::vector<HashShard> &shards, size_t index_size) {
      if (index_size == 0) {
        words_.clear();
        word_mask_ = 0;
        return;
      }

      // At least eight bits per indexed k-mer.  Power-of-two sizing makes the
      // hot query path independent of integer division and naturally scales
      // with the input rather than with any particular machine.
      const size_t target_words = DivCeiling(index_size, size_t{8});
      size_t num_words = 1;
      while (num_words < target_words) {
        if (num_words > std::numeric_limits<size_t>::max() / 2) {
          throw std::length_error("flank Bloom filter is too large");
        }
        num_words *= 2;
      }
      words_.assign(num_words, 0);
      word_mask_ = num_words - 1;

      for (const auto &shard : shards) {
        for (const auto &entry : shard.primary) {
          const size_t hash = FilterHash(entry.first);
          words_[hash & word_mask_] |= BitMask(hash);
        }
      }
    }

    bool MayContain(size_t hash) const {
      if (words_.empty()) {
        return false;
      }
      const uint64_t mask = BitMask(hash);
      return (words_[hash & word_mask_] & mask) == mask;
    }

    size_t byte_size() const { return words_.size() * sizeof(uint64_t); }

   private:
    static uint64_t BitMask(uint64_t hash) {
      // XXH3 already provides well-distributed bits.  Keep the low bits for
      // the word index and draw four independent six-bit positions from the
      // upper half of the hash.
      return (uint64_t{1} << ((hash >> 32u) & 63u)) |
             (uint64_t{1} << ((hash >> 40u) & 63u)) |
             (uint64_t{1} << ((hash >> 48u) & 63u)) |
             (uint64_t{1} << ((hash >> 56u) & 63u));
    }

    std::vector<uint64_t> words_;
    size_t word_mask_{0};
 };

  // A short oriented prefix is a necessary (but deliberately not sufficient)
  // condition for a read-side (k+1)-mer to match a stored flank.  Checking it
  // before canonicalizing the full multiword k-mer avoids the k-dependent
  // comparison and full-key Bloom hash at the overwhelming majority of read
  // positions.  Bloom false positives merely fall through to the original
  // exact lookup; all stored oriented prefixes are inserted, so false
  // negatives and result changes are impossible.
  class OrientedPrefixFilter {
   public:
    void Build(const std::vector<HashShard> &shards, size_t index_size,
               unsigned kmer_length) {
      prefix_len_ = std::min(19u, kmer_length);
      if (index_size == 0 || prefix_len_ == 0 ||
          std::getenv("MEGAHIT_DISABLE_PREFIX_FLANK_FILTER") != nullptr) {
        words_.clear();
        word_mask_ = 0;
        return;
      }

      // There are at most two oriented prefixes per canonical flank.  Use the
      // same eight-bit/key target as the full-key filter and power-of-two
      // blocked layout, derived solely from the live index cardinality.
      const size_t oriented_keys =
          index_size > std::numeric_limits<size_t>::max() / 2u
              ? std::numeric_limits<size_t>::max()
              : index_size * 2u;
      const size_t target_words = DivCeiling(oriented_keys, size_t{8});
      size_t num_words = 1;
      while (num_words < target_words) {
        if (num_words > std::numeric_limits<size_t>::max() / 2u) {
          throw std::length_error("oriented prefix filter is too large");
        }
        num_words *= 2u;
      }
      words_.assign(num_words, 0);
      word_mask_ = num_words - 1u;

      for (const HashShard &shard : shards) {
        for (const auto &entry : shard.primary) {
          Add(PrefixOf(entry.first, prefix_len_));
          KmerType reverse = entry.first;
          reverse.ReverseComplement(kmer_length);
          Add(PrefixOf(reverse, prefix_len_));
        }
      }
    }

    bool MayContain(uint64_t prefix) const {
      if (words_.empty()) {
        return true;
      }
      const uint64_t hash = Hash(prefix);
      const uint64_t mask = BitMask(hash);
      return (words_[hash & word_mask_] & mask) == mask;
    }

    unsigned prefix_len() const { return prefix_len_; }
    size_t byte_size() const { return words_.size() * sizeof(uint64_t); }

   private:
    static uint64_t PrefixOf(const KmerType &kmer, unsigned length) {
      uint64_t prefix = 0;
      for (unsigned i = 0; i < length; ++i) {
        prefix = (prefix << 2u) | kmer.GetBase(i);
      }
      return prefix;
    }

    static uint64_t Hash(uint64_t value) {
      value ^= value >> 30u;
      value *= UINT64_C(0xbf58476d1ce4e5b9);
      value ^= value >> 27u;
      value *= UINT64_C(0x94d049bb133111eb);
      return value ^ (value >> 31u);
    }

    static uint64_t BitMask(uint64_t hash) {
      return (uint64_t{1} << ((hash >> 32u) & 63u)) |
             (uint64_t{1} << ((hash >> 40u) & 63u)) |
             (uint64_t{1} << ((hash >> 48u) & 63u)) |
             (uint64_t{1} << ((hash >> 56u) & 63u));
    }

    void Add(uint64_t prefix) {
      const uint64_t hash = Hash(prefix);
      words_[hash & word_mask_] |= BitMask(hash);
    }

    std::vector<uint64_t> words_;
    size_t word_mask_{0};
    unsigned prefix_len_{0};
  };

 public:
  ContigFlankIndex(unsigned k, unsigned step) : k_(k), step_(step) {
    const size_t workers =
        static_cast<size_t>(std::max(1, omp_get_max_threads()));
    num_shards_ = 1;
    while (num_shards_ < workers) {
      if (num_shards_ > std::numeric_limits<size_t>::max() / 2) {
        throw std::length_error("too many flank-index workers");
      }
      num_shards_ *= 2;
    }
    shard_mask_ = num_shards_ - 1;
    hash_shards_.resize(num_shards_);
  }
  size_t size() const { return index_size_; }

  template <class Visitor>
  void ForEachCanonicalKmer(const Visitor &visitor) const {
    for (const HashShard &shard : hash_shards_) {
      for (const auto &entry : shard.primary) {
        visitor(entry.first);
      }
    }
  }

  template <class Visitor>
  void ForEachCanonicalKmerParallel(const Visitor &visitor) const {
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t shard = 0; shard < hash_shards_.size(); ++shard) {
      for (const auto &entry : hash_shards_[shard].primary) {
        visitor(entry.first);
      }
    }
  }

  void Finalize() {
    flank_filter_.Build(hash_shards_, index_size_);
    oriented_prefix_filter_.Build(hash_shards_, index_size_, k_ + 1u);
    xinfo("Flank membership filter: {} bytes for {} canonical k-mers\n",
          flank_filter_.byte_size(), index_size_);
    xinfo("Oriented {}-mer prefix filter: {} bytes\n",
          oriented_prefix_filter_.prefix_len(),
          oriented_prefix_filter_.byte_size());
  }

  void FeedBatchContigs(SeqPackage &seq_pkg, const std::vector<float> &mul) {
    const int num_threads = omp_get_max_threads();
    const size_t num_outboxes = static_cast<size_t>(num_threads) * num_shards_;
    std::vector<std::vector<PendingFlank>> thread_flanks(num_outboxes);
    const size_t reserve_per_outbox = DivCeiling(
        seq_pkg.seq_count() * 2,
        static_cast<size_t>(num_threads) * num_shards_);
    for (auto &flanks : thread_flanks) {
      flanks.reserve(reserve_per_outbox);
    }

#pragma omp parallel for
    for (size_t i = 0; i < seq_pkg.seq_count(); ++i) {
      auto seq_view = seq_pkg.GetSeqView(i);
      size_t seq_len = seq_view.length();
      if (seq_len < k_ + 1) {
        continue;
      }
      for (int strand = 0; strand < 2; ++strand) {
        auto get_jth_char = [&seq_view, strand,
                             seq_len](unsigned j) -> uint8_t {
          uint8_t c = seq_view.base_at(strand == 0 ? j : (seq_len - 1 - j));
          return strand == 0 ? c : 3u ^ c;
        };

        KmerType kmer;
        for (unsigned j = 0; j < k_ + 1; ++j) {
          kmer.ShiftAppend(get_jth_char(j), k_ + 1);
        }
        if (kmer.IsPalindrome(k_ + 1)) {
          continue;
        }

        unsigned ext_len =
            std::min(static_cast<size_t>(step_ - 1), seq_len - (k_ + 1));
        uint64_t ext_seq = 0;
        for (unsigned j = 0; j < ext_len; ++j) {
          ext_seq |= uint64_t(get_jth_char(k_ + 1 + j)) << j * 2;
        }

        KmerType canonical = kmer;
        KmerType reverse = canonical;
        reverse.ReverseComplement(k_ + 1);
        unsigned orientation = 0;
        if (reverse < canonical) {
          canonical = reverse;
          orientation = 1;
        }
        const size_t hash = MixedHash(canonical);
        auto &flanks =
            thread_flanks[static_cast<size_t>(omp_get_thread_num()) *
                              num_shards_ +
                          (hash & shard_mask_)];
        flanks.push_back(
            PendingFlank{canonical, FlankInfo{ext_seq, ext_len}, orientation});
        if (seq_len == k_ + 1) {
          break;
        }
      }
    }

    // Hash partitioning gives every shard a single writer.  All shard merges
    // therefore run in parallel without locks, while each key still follows
    // the exact max-extension winner rule used by the serial map.
    size_t total_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : total_size)
    for (size_t shard_id = 0; shard_id < num_shards_; ++shard_id) {
      HashShard &shard = hash_shards_[shard_id];
      for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        const auto &flanks =
            thread_flanks[static_cast<size_t>(thread_id) * num_shards_ +
                          shard_id];
        for (const auto &flank : flanks) {
          auto res = shard.primary.emplace(
              flank.canonical,
              StoredFlankInfo{flank.info.ext_seq, flank.info.ext_len,
                              flank.orientation, 0, 0, flank.info.mul});
          if (res.second) {
            continue;
          }

          StoredFlankInfo &stored = res.first->second;
          if (stored.orientation == flank.orientation) {
            if (stored.ext_len < flank.info.ext_len ||
                (stored.ext_len == flank.info.ext_len &&
                 stored.ext_seq < flank.info.ext_seq)) {
              stored.ext_seq = flank.info.ext_seq;
              stored.ext_len = flank.info.ext_len;
              stored.mul = flank.info.mul;
            }
          } else {
            auto secondary =
                shard.secondary.emplace(flank.canonical, flank.info);
            if (!secondary.second) {
              FlankInfo &old = secondary.first->second;
              if (old.ext_len < flank.info.ext_len ||
                  (old.ext_len == flank.info.ext_len &&
                   old.ext_seq < flank.info.ext_seq)) {
                old = flank.info;
              }
            }
            stored.has_secondary = 1;
          }
        }
      }
      total_size += shard.primary.size();
    }
    index_size_ = total_size;
  }

  // Run the original exact state machine for one arbitrary packed-read view.
  // Keeping this primitive independent of SeqPackage lets the persistent
  // index builder fuse its read scan with first-round iterative edge
  // generation without copying or decoding each read into a second owner.
  template <class ReadViewType, class CollectorType>
  bool FindNextKmersFromRead(
      const ReadViewType &seq_view, CollectorType *out,
      std::vector<uint32_t> *state_scratch = nullptr,
      uint64_t *num_generated_edges_out = nullptr,
      const std::vector<uint32_t> *candidate_positions = nullptr) const {
    const size_t length = seq_view.length();
    if (length < k_ + step_ + 1) return false;

    std::vector<uint32_t> local_state;
    std::vector<uint32_t> &kmer_state =
        state_scratch == nullptr ? local_state : *state_scratch;
    kmer_state.assign(length, 0u);
    uint64_t num_generated_edges = 0;
    bool success = false;

    Flank flank, rflank;
    auto &kmer = flank.kmer;
    auto &rkmer = rflank.kmer;

    const unsigned prefix_len = oriented_prefix_filter_.prefix_len();
    const uint64_t prefix_mask =
        prefix_len == 32u
            ? std::numeric_limits<uint64_t>::max()
            : (uint64_t{1} << (2u * prefix_len)) - 1u;
    uint64_t oriented_prefix = 0;
    if (candidate_positions == nullptr) {
      for (unsigned j = 0; j < prefix_len; ++j) {
        oriented_prefix =
            (oriented_prefix << 2u) | seq_view.base_at(j);
      }
    }

    for (unsigned j = 0; j < k_ + 1; ++j) {
      kmer.ShiftAppend(seq_view.base_at(j), k_ + 1);
    }
    rkmer = kmer;
    rkmer.ReverseComplement(k_ + 1);

    unsigned cur_pos = 0;
    size_t candidate_index = 0;
    while (cur_pos + k_ + 1 <= length) {
      unsigned next_pos = cur_pos + 1;
      bool may_match = false;
      if (candidate_positions == nullptr) {
        may_match = oriented_prefix_filter_.MayContain(oriented_prefix);
      } else {
        while (candidate_index < candidate_positions->size() &&
               (*candidate_positions)[candidate_index] < cur_pos) {
          ++candidate_index;
        }
        may_match = candidate_index < candidate_positions->size() &&
                    (*candidate_positions)[candidate_index] == cur_pos;
      }

      if (kmer_state[cur_pos] == 0 && may_match) {
        const bool query_is_reverse = rkmer < kmer;
        const KmerType &canonical = query_is_reverse ? rkmer : kmer;
        const unsigned forward_orientation = query_is_reverse ? 1u : 0u;
        FlankInfo forward_info{};
        FlankInfo reverse_info{};
        bool has_forward = false;
        bool has_reverse = false;
        if (flank_filter_.MayContain(FilterHash(canonical))) {
          // Only Bloom positives need phmap's exact pre-mixed XXH3 value.
          const size_t canonical_hash = MixedHash(canonical);
          const HashShard &shard =
              hash_shards_[canonical_hash & shard_mask_];
          auto iter = shard.primary.find(canonical, canonical_hash);
          if (iter != shard.primary.end()) {
            const StoredFlankInfo &stored = iter->second;
            FlankInfo primary{stored.ext_seq, stored.ext_len, stored.mul};
            if (stored.orientation == forward_orientation) {
              forward_info = primary;
              has_forward = true;
            } else {
              reverse_info = primary;
              has_reverse = true;
            }
            if (stored.has_secondary) {
              auto secondary =
                  shard.secondary.find(canonical, canonical_hash);
              assert(secondary != shard.secondary.end());
              if (stored.orientation == forward_orientation) {
                reverse_info = secondary->second;
                has_reverse = true;
              } else {
                forward_info = secondary->second;
                has_forward = true;
              }
            }
          }
        }
        if (has_forward) {
          const FlankInfo &info = forward_info;
          const uint64_t ext_seq = info.ext_seq;
          const unsigned ext_len = info.ext_len;
          const float mul = info.mul;
          kmer_state[cur_pos] = EncodeMultiplicity(mul);

          for (unsigned j = 0;
               j < ext_len && cur_pos + k_ + 1 + j < length;
               ++j, ++next_pos) {
            if (seq_view.base_at(cur_pos + k_ + 1 + j) ==
                ((ext_seq >> j * 2u) & 3u)) {
              kmer_state[cur_pos + j + 1] = EncodeMultiplicity(mul);
            } else {
              break;
            }
          }
        }
        if (has_reverse) {
          const FlankInfo &info = reverse_info;
          const uint64_t ext_seq = info.ext_seq;
          const unsigned ext_len = info.ext_len;
          const float mul = info.mul;
          kmer_state[cur_pos] =
              kmer_state[cur_pos] != 0
                  ? EncodeMultiplicity(
                        (DecodeMultiplicity(kmer_state[cur_pos]) + mul) / 2)
                  : EncodeMultiplicity(mul);

          for (unsigned j = 0; j < ext_len && cur_pos >= j + 1; ++j) {
            if ((3u ^ seq_view.base_at(cur_pos - 1 - j)) ==
                ((ext_seq >> j * 2u) & 3u)) {
              uint32_t &state = kmer_state[cur_pos - 1 - j];
              state = state != 0
                          ? EncodeMultiplicity(
                                (DecodeMultiplicity(state) + mul) / 2)
                          : EncodeMultiplicity(mul);
            } else {
              break;
            }
          }
        }
      }

      if (next_pos + k_ + 1 <= length) {
        while (cur_pos < next_pos) {
          ++cur_pos;
          const uint8_t c = seq_view.base_at(cur_pos + k_);
          kmer.ShiftAppend(c, k_ + 1);
          rkmer.ShiftPreappend(3u ^ c, k_ + 1);
          if (candidate_positions == nullptr) {
            oriented_prefix =
                ((oriented_prefix << 2u) |
                 seq_view.base_at(cur_pos + prefix_len - 1u)) &
                prefix_mask;
          }
        }
      } else {
        break;
      }
    }

    typename CollectorType::kmer_type new_kmer, new_rkmer;
    float prefix_mul = 0;
    for (unsigned accumulated_len = 0, j = 0, end_pos = 0;
         j + k_ < length; ++j) {
      const bool kmer_exists = kmer_state[j] != 0;
      if (kmer_exists) prefix_mul += DecodeMultiplicity(kmer_state[j]);
      // Earlier positions are no longer queried for membership, so reuse the
      // same compact array for the prefix sums consumed by later windows.
      kmer_state[j] = FloatBits(prefix_mul);
      accumulated_len = kmer_exists ? accumulated_len + 1 : 0;
      if (accumulated_len < step_ + 1) continue;

      const unsigned target_end = j + k_ + 1;
      static const bool direct_edge_windows =
          std::getenv("MEGAHIT_EXPERIMENTAL_DIRECT_ITERATE_WINDOWS") != nullptr;
      if (direct_edge_windows) {
        const unsigned next_length = k_ + step_ + 1u;
        if (end_pos == 0u || target_end - end_pos > 8u) {
          InitReadKmer(seq_view, target_end - next_length, next_length,
                       &new_kmer);
          new_rkmer = new_kmer;
          new_rkmer.ReverseComplement(next_length);
          end_pos = target_end;
        } else {
          while (end_pos < target_end) {
            const auto base = seq_view.base_at(end_pos++);
            new_kmer.ShiftAppend(base, next_length);
            new_rkmer.ShiftPreappend(3u ^ base, next_length);
          }
        }
      } else if (end_pos + 8 < target_end) {
        while (end_pos < target_end) {
          const auto c = seq_view.base_at(end_pos++);
          new_kmer.ShiftAppend(c, k_ + step_ + 1);
          new_rkmer.ShiftPreappend(3u ^ c, k_ + step_ + 1);
        }
      } else {
        if (end_pos + k_ + step_ + 1 < target_end) {
          end_pos = target_end - (k_ + step_ + 1);
        }
        while (end_pos < target_end) {
          new_kmer.ShiftAppend(seq_view.base_at(end_pos++),
                               k_ + step_ + 1);
        }
        new_rkmer = new_kmer;
        new_rkmer.ReverseComplement(k_ + step_ + 1);
      }
      const float previous_prefix =
          j >= step_ + 1 ? BitsFloat(kmer_state[j - (step_ + 1)]) : 0;
      const float mul = (prefix_mul - previous_prefix) / (step_ + 1);
      assert(mul <= kMaxMul + 1);
      out->Insert(new_kmer < new_rkmer ? new_kmer : new_rkmer,
                  static_cast<mul_t>(
                      std::min(kMaxMul, static_cast<int>(mul + 0.5))));
      ++num_generated_edges;
      success = true;
    }
    if (num_generated_edges_out != nullptr) {
      *num_generated_edges_out += num_generated_edges;
    }
    return success;
  }

  template <class CollectorType>
  size_t FindNextKmersFromReads(
      const SeqPackage &seq_pkg, CollectorType *out,
      uint64_t *num_generated_edges_out = nullptr) const {
    std::vector<uint32_t> kmer_state;
    size_t num_aligned_reads = 0;
    uint64_t num_generated_edges = 0;

#pragma omp parallel for reduction(+ : num_aligned_reads, num_generated_edges) \
    private(kmer_state)
    for (unsigned seq_id = 0; seq_id < seq_pkg.seq_count(); ++seq_id) {
      uint64_t read_generated_edges = 0;
      num_aligned_reads += FindNextKmersFromRead(
          seq_pkg.GetSeqView(seq_id), out, &kmer_state,
          &read_generated_edges);
      num_generated_edges += read_generated_edges;
    }
    if (num_generated_edges_out != nullptr) {
      *num_generated_edges_out += num_generated_edges;
    }
    return num_aligned_reads;
  }

  // Replay the original per-read state machine from an exact, complete set of
  // candidate starts supplied in increasing order.  Candidate generation may
  // use a persistent short-seed occurrence index, but every candidate still
  // performs the original full (k+1)-mer lookup and all multiplicity/
  // extension updates occur in the same left-to-right order.  Consequently a
  // false-positive candidate only costs work; it cannot alter the output.
  template <class ReadViewType, class CollectorType>
  bool FindNextKmersFromReadCandidates(
      const ReadViewType &seq_view,
      const std::vector<uint32_t> &candidate_positions,
      CollectorType *out,
      std::vector<uint32_t> *state_scratch = nullptr) const {
    const size_t length = seq_view.length();
    if (length < k_ + step_ + 1) {
      return false;
    }
    std::vector<uint32_t> local_state;
    std::vector<uint32_t> &kmer_state =
        state_scratch == nullptr ? local_state : *state_scratch;
    kmer_state.assign(length, 0u);
    uint32_t previous_position = std::numeric_limits<uint32_t>::max();

    for (uint32_t cur_pos : candidate_positions) {
      if (cur_pos == previous_position) continue;
      previous_position = cur_pos;
      if (cur_pos + k_ + 1 > length || kmer_state[cur_pos] != 0) continue;

      KmerType kmer;
      InitCandidateKmer(seq_view, cur_pos, &kmer);
      KmerType rkmer = kmer;
      rkmer.ReverseComplement(k_ + 1);
      const bool query_is_reverse = rkmer < kmer;
      const KmerType &canonical = query_is_reverse ? rkmer : kmer;
      const unsigned forward_orientation = query_is_reverse ? 1u : 0u;
      FlankInfo forward_info{};
      FlankInfo reverse_info{};
      bool has_forward = false;
      bool has_reverse = false;

      if (flank_filter_.MayContain(FilterHash(canonical))) {
        const size_t canonical_hash = MixedHash(canonical);
        const HashShard &shard =
            hash_shards_[canonical_hash & shard_mask_];
        auto iterator = shard.primary.find(canonical, canonical_hash);
        if (iterator != shard.primary.end()) {
          const StoredFlankInfo &stored = iterator->second;
          const FlankInfo primary{stored.ext_seq, stored.ext_len, stored.mul};
          if (stored.orientation == forward_orientation) {
            forward_info = primary;
            has_forward = true;
          } else {
            reverse_info = primary;
            has_reverse = true;
          }
          if (stored.has_secondary) {
            auto secondary = shard.secondary.find(canonical, canonical_hash);
            assert(secondary != shard.secondary.end());
            if (stored.orientation == forward_orientation) {
              reverse_info = secondary->second;
              has_reverse = true;
            } else {
              forward_info = secondary->second;
              has_forward = true;
            }
          }
        }
      }

      if (has_forward) {
        const FlankInfo &info = forward_info;
        kmer_state[cur_pos] = EncodeMultiplicity(info.mul);
        for (unsigned j = 0;
             j < info.ext_len && cur_pos + k_ + 1 + j < length; ++j) {
          if (seq_view.base_at(cur_pos + k_ + 1 + j) !=
              ((info.ext_seq >> (j * 2u)) & 3u)) {
            break;
          }
          kmer_state[cur_pos + j + 1u] = EncodeMultiplicity(info.mul);
        }
      }
      if (has_reverse) {
        const FlankInfo &info = reverse_info;
        kmer_state[cur_pos] =
            kmer_state[cur_pos] != 0
                ? EncodeMultiplicity(
                      (DecodeMultiplicity(kmer_state[cur_pos]) + info.mul) / 2)
                : EncodeMultiplicity(info.mul);
        for (unsigned j = 0; j < info.ext_len && cur_pos >= j + 1u; ++j) {
          if ((3u ^ seq_view.base_at(cur_pos - 1u - j)) !=
              ((info.ext_seq >> (j * 2u)) & 3u)) {
            break;
          }
          uint32_t &state = kmer_state[cur_pos - 1u - j];
          state = state != 0
                      ? EncodeMultiplicity(
                            (DecodeMultiplicity(state) + info.mul) / 2)
                      : EncodeMultiplicity(info.mul);
        }
      }
    }

    typename CollectorType::kmer_type new_kmer, new_rkmer;
    float prefix_mul = 0;
    bool success = false;
    for (unsigned accumulated_len = 0, j = 0, end_pos = 0;
         j + k_ < length; ++j) {
      const bool kmer_exists = kmer_state[j] != 0;
      if (kmer_exists) prefix_mul += DecodeMultiplicity(kmer_state[j]);
      kmer_state[j] = FloatBits(prefix_mul);
      accumulated_len = kmer_exists ? accumulated_len + 1u : 0u;
      if (accumulated_len < step_ + 1u) continue;

      const unsigned target_end = j + k_ + 1u;
      static const bool direct_edge_windows =
          std::getenv("MEGAHIT_EXPERIMENTAL_DIRECT_ITERATE_WINDOWS") != nullptr;
      if (direct_edge_windows) {
        const unsigned next_length = k_ + step_ + 1u;
        if (end_pos == 0u || target_end - end_pos > 8u) {
          InitReadKmer(seq_view, target_end - next_length, next_length,
                       &new_kmer);
          new_rkmer = new_kmer;
          new_rkmer.ReverseComplement(next_length);
          end_pos = target_end;
        } else {
          while (end_pos < target_end) {
            const auto base = seq_view.base_at(end_pos++);
            new_kmer.ShiftAppend(base, next_length);
            new_rkmer.ShiftPreappend(3u ^ base, next_length);
          }
        }
      } else if (end_pos + 8u < target_end) {
        while (end_pos < target_end) {
          const auto c = seq_view.base_at(end_pos++);
          new_kmer.ShiftAppend(c, k_ + step_ + 1u);
          new_rkmer.ShiftPreappend(3u ^ c, k_ + step_ + 1u);
        }
      } else {
        if (end_pos + k_ + step_ + 1u < target_end) {
          end_pos = target_end - (k_ + step_ + 1u);
        }
        while (end_pos < target_end) {
          new_kmer.ShiftAppend(seq_view.base_at(end_pos++),
                               k_ + step_ + 1u);
        }
        new_rkmer = new_kmer;
        new_rkmer.ReverseComplement(k_ + step_ + 1u);
      }
      const float previous_prefix =
          j >= step_ + 1u ? BitsFloat(kmer_state[j - (step_ + 1u)]) : 0;
      const float mul = (prefix_mul - previous_prefix) / (step_ + 1u);
      out->Insert(new_kmer < new_rkmer ? new_kmer : new_rkmer,
                  static_cast<mul_t>(
                      std::min(kMaxMul, static_cast<int>(mul + 0.5))));
      success = true;
    }
    return success;
  }

 private:
  template <class ReadViewType>
  void InitCandidateKmer(const ReadViewType &seq_view, unsigned position,
                         KmerType *kmer) const {
    InitReadKmer(seq_view, position, k_ + 1u, kmer);
  }

  template <class ReadViewType, class TargetKmer>
  void InitReadKmer(const ReadViewType &seq_view, unsigned position,
                    unsigned length, TargetKmer *kmer) const {
    InitReadKmerImpl(seq_view, position, length, kmer,
        std::integral_constant<bool,
            std::is_same<typename TargetKmer::word_type,
                         typename ReadViewType::word_type>::value>());
  }

  template <class ReadViewType, class TargetKmer>
  void InitReadKmerImpl(const ReadViewType &seq_view, unsigned position,
                        unsigned length, TargetKmer *kmer,
                        std::true_type) const {
    const auto raw = seq_view.raw_address();
    kmer->InitFromPtr(raw.first, raw.second + position, length);
  }

  template <class ReadViewType, class TargetKmer>
  void InitReadKmerImpl(const ReadViewType &seq_view, unsigned position,
                        unsigned length, TargetKmer *kmer,
                        std::false_type) const {
    InitReadKmerCrossWords(seq_view, position, length, kmer,
        std::integral_constant<bool,
            sizeof(typename ReadViewType::word_type) == 4u &&
            sizeof(typename TargetKmer::word_type) == 8u>());
  }

  template <class ReadViewType, class TargetKmer>
  void InitReadKmerCrossWords(const ReadViewType &seq_view,
                              unsigned position, unsigned length,
                              TargetKmer *kmer, std::true_type) const {
    // The uint32 loader checks whether the final source word is needed.
    // Join its exact, zero-padded output into uint64 lanes without reading
    // beyond a packed read at the end of the mapped file.
    const auto raw = seq_view.raw_address();
    Kmer<TargetKmer::kNumWords * 2u, uint32_t> narrow;
    narrow.InitFromPtr(raw.first, raw.second + position, length);
    uint64_t packed[TargetKmer::kNumWords];
    for (unsigned word = 0; word < TargetKmer::kNumWords; ++word) {
      packed[word] = (uint64_t(narrow.data()[word * 2u]) << 32u) |
                     narrow.data()[word * 2u + 1u];
    }
    kmer->InitFromPtr(packed, 0u, length);
  }

  template <class ReadViewType, class TargetKmer>
  void InitReadKmerCrossWords(const ReadViewType &seq_view,
                              unsigned position, unsigned length,
                              TargetKmer *kmer, std::false_type) const {
    for (unsigned base = 0; base < length; ++base) {
      kmer->ShiftAppend(seq_view.base_at(position + base), length);
    }
  }

  static size_t MixedHash(const KmerType &kmer) {
    return phmap::phmap_mix<sizeof(size_t)>()(KmerHash{}(kmer));
  }

  static uint32_t FloatBits(float value) {
    uint32_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
  }

  static float BitsFloat(uint32_t bits) {
    float value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
  }

  static uint32_t EncodeMultiplicity(float value) {
    // Multiplicity is non-negative.  Set the sign bit to distinguish a real
    // zero multiplicity from an absent k-mer without changing its payload.
    return FloatBits(value) | uint32_t{0x80000000u};
  }

  static float DecodeMultiplicity(uint32_t state) {
    return BitsFloat(state & uint32_t{0x7fffffffu});
  }

  std::vector<HashShard> hash_shards_;
  BlockedBloomFilter flank_filter_;
  OrientedPrefixFilter oriented_prefix_filter_;
  size_t num_shards_{0};
  size_t shard_mask_{0};
  size_t index_size_{0};
  unsigned k_{};
  unsigned step_{};
};

#endif  // MEGAHIT_JUNCION_INDEX_H
