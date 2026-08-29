//
// Created by vout on 11/23/18.
//

#ifndef MEGAHIT_JUNCION_INDEX_H
#define MEGAHIT_JUNCION_INDEX_H

#include <algorithm>
#include <cstring>
#include <limits>
#include <omp.h>
#include <stdexcept>
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

  void Finalize() {
    flank_filter_.Build(hash_shards_, index_size_);
    xinfo("Flank membership filter: {} bytes for {} canonical k-mers\n",
          flank_filter_.byte_size(), index_size_);
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

  template <class CollectorType>
  size_t FindNextKmersFromReads(const SeqPackage &seq_pkg,
                                CollectorType *out) const {
    std::vector<uint32_t> kmer_state;
    size_t num_aligned_reads = 0;

#pragma omp parallel for reduction(+ : num_aligned_reads) private(kmer_state)
    for (unsigned seq_id = 0; seq_id < seq_pkg.seq_count(); ++seq_id) {
      auto seq_view = seq_pkg.GetSeqView(seq_id);
      size_t length = seq_view.length();
      if (length < k_ + step_ + 1) {
        continue;
      }

      bool success = false;
      kmer_state.clear();
      kmer_state.resize(length, 0);

      Flank flank, rflank;
      auto &kmer = flank.kmer;
      auto &rkmer = rflank.kmer;

      for (unsigned j = 0; j < k_ + 1; ++j) {
        kmer.ShiftAppend(seq_view.base_at(j), k_ + 1);
      }
      rkmer = kmer;
      rkmer.ReverseComplement(k_ + 1);

      unsigned cur_pos = 0;
      while (cur_pos + k_ + 1 <= length) {
        unsigned next_pos = cur_pos + 1;

        if (kmer_state[cur_pos] == 0) {
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
            uint64_t ext_seq = info.ext_seq;
            unsigned ext_len = info.ext_len;
            float mul = info.mul;
            kmer_state[cur_pos] = EncodeMultiplicity(mul);

            for (unsigned j = 0; j < ext_len && cur_pos + k_ + 1 + j < length;
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
            uint64_t ext_seq = info.ext_seq;
            unsigned ext_len = info.ext_len;
            float mul = info.mul;
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
            uint8_t c = seq_view.base_at(cur_pos + k_);
            kmer.ShiftAppend(c, k_ + 1);
            rkmer.ShiftPreappend(3u ^ c, k_ + 1);
          }
        } else {
          break;
        }
      }

      typename CollectorType::kmer_type new_kmer, new_rkmer;

      float prefix_mul = 0;
      for (unsigned accumulated_len = 0, j = 0, end_pos = 0; j + k_ < length;
           ++j) {
        const bool kmer_exists = kmer_state[j] != 0;
        if (kmer_exists) {
          prefix_mul += DecodeMultiplicity(kmer_state[j]);
        }
        // Earlier positions are no longer queried for membership, so reuse the
        // same compact array for the prefix sums consumed by later windows.
        kmer_state[j] = FloatBits(prefix_mul);
        accumulated_len = kmer_exists ? accumulated_len + 1 : 0;
        if (accumulated_len >= step_ + 1) {
          unsigned target_end = j + k_ + 1;
          if (end_pos + 8 < target_end) {
            while (end_pos < target_end) {
              auto c = seq_view.base_at(end_pos);
              new_kmer.ShiftAppend(c, k_ + step_ + 1);
              new_rkmer.ShiftPreappend(3u ^ c, k_ + step_ + 1);
              end_pos++;
            }
          } else {
            if (end_pos + k_ + step_ + 1 < target_end) {
              end_pos = target_end - (k_ + step_ + 1);
            }
            while (end_pos < target_end) {
              auto c = seq_view.base_at(end_pos);
              new_kmer.ShiftAppend(c, k_ + step_ + 1);
              end_pos++;
            }
            new_rkmer = new_kmer;
            new_rkmer.ReverseComplement(k_ + step_ + 1);
          }
          const float previous_prefix =
              j >= step_ + 1 ? BitsFloat(kmer_state[j - (step_ + 1)]) : 0;
          float mul = (prefix_mul - previous_prefix) / (step_ + 1);
          assert(mul <= kMaxMul + 1);
          out->Insert(new_kmer < new_rkmer ? new_kmer : new_rkmer,
                      static_cast<mul_t>(
                          std::min(kMaxMul, static_cast<int>(mul + 0.5))));
          success = true;
        }
      }
      num_aligned_reads += success;
    }
    return num_aligned_reads;
  }

 private:
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
  size_t num_shards_{0};
  size_t shard_mask_{0};
  size_t index_size_{0};
  unsigned k_{};
  unsigned step_{};
};

#endif  // MEGAHIT_JUNCION_INDEX_H
