/**
 * @file kmer.h
 * @brief IdbaKmer Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __BASIC_KMER_H_

#define __BASIC_KMER_H_

#include <stdint.h>

#include <algorithm>
#include <cstring>
#ifdef USE_BMI2
#include <immintrin.h>
#endif
#include "xxhash/xxh3.h"

#include "bit_operation.h"
#include "definitions.h"

/**
 * @brief It represents a k-mer. The value of k is limited by the number of
 * uint64 words used. The maximum value can be calculated by max_size().
 */
class IdbaKmer {
 public:
  IdbaKmer() { std::memset(data_, 0, sizeof(uint64_t) * kNumUint64); }

  IdbaKmer(const IdbaKmer &kmer) {
    std::memcpy(data_, kmer.data_, sizeof(uint64_t) * kNumUint64);
  }

  explicit IdbaKmer(uint32_t size) {
    std::memset(data_, 0, sizeof(uint64_t) * kNumUint64);
    resize(size);
  }

  IdbaKmer(const uint8_t *bases, uint32_t size) {
    AssignBases(bases, size);
  }

  ~IdbaKmer() {}

  const IdbaKmer &operator=(const IdbaKmer &kmer) {
    std::memcpy(data_, kmer.data_, sizeof(uint64_t) * kNumUint64);
    return *this;
  }

  bool operator<(const IdbaKmer &kmer) const {
    const uint32_t this_size = size();
    const uint32_t other_size = kmer.size();
    if (this_size != other_size) return this_size < other_size;
    for (int i = static_cast<int>(NumBaseWords(this_size)) - 1; i >= 0; --i) {
      if (data_[i] != kmer.data_[i]) return data_[i] < kmer.data_[i];
    }
    return false;
  }

  bool operator>(const IdbaKmer &kmer) const {
    const uint32_t this_size = size();
    const uint32_t other_size = kmer.size();
    if (this_size != other_size) return this_size > other_size;
    for (int i = static_cast<int>(NumBaseWords(this_size)) - 1; i >= 0; --i) {
      if (data_[i] != kmer.data_[i]) return data_[i] > kmer.data_[i];
    }
    return false;
  }

  bool operator==(const IdbaKmer &kmer) const {
    const uint32_t this_size = size();
    if (this_size != kmer.size()) return false;
    const uint32_t used_words = NumBaseWords(this_size);
    for (uint32_t i = 0; i < used_words; ++i) {
      if (data_[i] != kmer.data_[i]) return false;
    }
    return true;
  }

  bool operator!=(const IdbaKmer &kmer) const {
    const uint32_t this_size = size();
    if (this_size != kmer.size()) return true;
    const uint32_t used_words = NumBaseWords(this_size);
    for (uint32_t i = 0; i < used_words; ++i) {
      if (data_[i] != kmer.data_[i]) return true;
    }
    return false;
  }

  const IdbaKmer &ReverseComplement() {
    uint32_t kmer_size = size();
    uint32_t used_words = (kmer_size + 31) >> 5;

    resize(0);

    for (unsigned i = 0; i < used_words; ++i)
      bit_operation::ReverseComplement(data_[i]);

    for (unsigned i = 0; i < (used_words >> 1); ++i)
      std::swap(data_[i], data_[used_words - 1 - i]);

    if ((kmer_size & 31) != 0) {
      unsigned offset = (32 - (kmer_size & 31)) << 1;
      for (unsigned i = 0; i + 1 < used_words; ++i)
        data_[i] = (data_[i] >> offset) | data_[i + 1] << (64 - offset);
      data_[used_words - 1] >>= offset;
    }

    resize(kmer_size);

    return *this;
  }

  void ShiftAppend(uint8_t ch) {
    ch &= 3;
    uint32_t kmer_size = size();
    uint32_t used_words = (kmer_size + 31) >> 5;

    resize(0);

    for (unsigned i = 0; i + 1 < used_words; ++i)
      data_[i] = (data_[i] >> 2) | (data_[i + 1] << 62);
    data_[used_words - 1] = (data_[used_words - 1] >> 2) |
                            (uint64_t(ch) << (((kmer_size - 1) & 31) << 1));

    resize(kmer_size);
  }

  void ShiftPreappend(uint8_t ch) {
    ch &= 3;
    uint32_t kmer_size = size();
    uint32_t used_words = (kmer_size + 31) >> 5;

    resize(0);

    for (int i = used_words - 1; i > 0; --i)
      data_[i] = (data_[i] << 2) | (data_[i - 1] >> 62);
    data_[0] = (data_[0] << 2) | ch;

    if ((kmer_size & 31) != 0)
      data_[used_words - 1] &= (1ULL << ((kmer_size & 31) << 1)) - 1;

    resize(kmer_size);
  }

  bool IsPalindrome() const {
    // No DNA base is its own complement, so an odd-length sequence cannot be
    // equal to its reverse complement.  Local assembly commonly uses only
    // odd k values; avoid constructing and reversing a full IdbaKmer for this
    // mathematically impossible case.
    if ((size() & 1u) != 0) return false;
    IdbaKmer kmer(*this);
    return kmer.ReverseComplement() == *this;
  }

  uint64_t hash() const {
    // Keep the historical full-capacity XXH3 value bit-for-bit: IDBA graph
    // traversal is sensitive to bucket order.  The configured 255-base Idba
    // key occupies 72 bytes, however, and most local-assembly keys leave the
    // upper words permanently zero.  XXH3's 72-byte path is six independent
    // 16-byte mixes.  Fold mixes whose two input words are known constants
    // into a per-k accumulator, while evaluating every mix which contains an
    // active base word normally.  This is algebraically the same XXH3 call,
    // not a replacement hash.
    if (kNumUint64 != 9u || sizeof(data_) != 72u || size() > 255u) {
      return XXH3_64bits(data_, sizeof(data_));
    }

    const uint32_t used = NumBaseWords(size());
    const char *bytes = reinterpret_cast<const char *>(data_);
    const char *key = reinterpret_cast<const char *>(kKey);
    uint64_t acc = FixedHashAccumulators()[size()];
    if (used > 4u) acc += XXH3_mix16B(bytes + 32, key + 64, 0);
    if (used > 3u) acc += XXH3_mix16B(bytes + 24, key + 80, 0);
    if (used > 2u) acc += XXH3_mix16B(bytes + 16, key + 32, 0);
    if (used > 5u) acc += XXH3_mix16B(bytes + 40, key + 48, 0);
    acc += XXH3_mix16B(bytes, key, 0);
    if (used > 7u) acc += XXH3_mix16B(bytes + 56, key + 16, 0);
    return XXH3_avalanche(acc);
  }

  // Hash an active-word rolling key without first materializing and clearing
  // the 72-byte IdbaKmer. `words` provides one zero look-ahead word. The
  // result is exactly hash() for an IdbaKmer containing these active words.
  template <unsigned Words>
  static uint64_t HashPackedWords(const uint64_t *words,
                                  uint32_t kmer_size) {
    return HashPackedWords<Words>(words, kmer_size,
                                  PackedHashFixedAccumulator(kmer_size));
  }

  template <unsigned Words>
  static uint64_t HashPackedWords(const uint64_t *words,
                                  uint32_t kmer_size,
                                  uint64_t fixed_accumulator) {
    const char *bytes = reinterpret_cast<const char *>(words);
    const char *key = reinterpret_cast<const char *>(kKey);
    uint64_t acc = fixed_accumulator;
    if (Words > 4u) acc += XXH3_mix16B(bytes + 32, key + 64, 0);
    if (Words > 3u) acc += XXH3_mix16B(bytes + 24, key + 80, 0);
    if (Words > 2u) acc += XXH3_mix16B(bytes + 16, key + 32, 0);
    if (Words > 5u) acc += XXH3_mix16B(bytes + 40, key + 48, 0);
    acc += XXH3_mix16B(bytes, key, 0);
    if (Words > 7u) {
      const uint64_t final_words[2] = {
          words[7], uint64_t(kmer_size) << (64u - kBitsForSize)};
      acc += XXH3_mix16B(final_words, key + 16, 0);
    }
    return XXH3_avalanche(acc);
  }

  static uint64_t PackedHashFixedAccumulator(uint32_t kmer_size) {
    return FixedHashAccumulators()[kmer_size];
  }

  template <unsigned Words>
  bool EqualsPackedWords(const uint64_t *words, uint32_t kmer_size) const {
    if (size() != kmer_size) return false;
    for (unsigned word = 0; word < Words; ++word) {
      if (data_[word] != words[word]) return false;
    }
    return true;
  }

  IdbaKmer unique_format() const {
    IdbaKmer rev_comp = *this;
    rev_comp.ReverseComplement();
    return (*this < rev_comp ? *this : rev_comp);
  }

  uint8_t operator[](uint32_t index) const {
    return (data_[index >> 5] >> ((index & 31) << 1)) & 3;
  }

  uint8_t get_base(uint32_t index) const {
    return (data_[index >> 5] >> ((index & 31) << 1)) & 3;
  }

  uint64_t word(uint32_t index) const { return data_[index]; }
  const uint64_t *words() const { return data_; }

  void AssignWords(const uint64_t *words, uint32_t new_size) {
    std::memset(data_, 0, sizeof(data_));
    std::memcpy(data_, words,
                NumBaseWords(new_size) * sizeof(uint64_t));
    resize(new_size);
  }

  void set_base(uint32_t index, uint8_t ch) {
    ch &= 3;
    unsigned offset = (index & 31) << 1;
    data_[index >> 5] =
        (data_[index >> 5] & ~(3ULL << offset)) | (uint64_t(ch) << offset);
  }

  void AssignBases(const uint8_t *bases, uint32_t new_size) {
    std::memset(data_, 0, sizeof(data_));
    const uint32_t used_words = NumBaseWords(new_size);
    for (uint32_t word_id = 0; word_id < used_words; ++word_id) {
      const uint32_t begin = word_id << 5u;
      const uint32_t count = std::min<uint32_t>(32u, new_size - begin);
      uint64_t packed = 0;
      uint32_t i = 0;
#ifdef USE_BMI2
      // Sequence stores one two-bit base in each byte.  PEXT packs eight
      // bytes into sixteen contiguous bits, replacing eight dependent
      // read/modify/write operations in the legacy set_base loop.
      constexpr uint64_t kTwoBitsPerByte = UINT64_C(0x0303030303030303);
      for (; i + 8u <= count; i += 8u) {
        uint64_t bytes;
        std::memcpy(&bytes, bases + begin + i, sizeof(bytes));
        packed |= _pext_u64(bytes, kTwoBitsPerByte) << (i << 1u);
      }
#endif
      for (; i < count; ++i) {
        packed |= uint64_t(bases[begin + i] & 3u) << (i << 1u);
      }
      data_[word_id] = packed;
    }
    resize(new_size);
  }

  // Load a k-mer directly from a little-endian 2-bit base stream (base 0 in
  // bits 1:0).  The source must provide one zero-padded look-ahead word.  A
  // local read can therefore be packed once and reused by every inner-k
  // round, instead of replaying k dependent ShiftAppend/ShiftPreappend steps
  // before the first candidate of every round.
  void AssignPackedBases(const uint64_t *words, uint32_t base_offset,
                         uint32_t new_size) {
    std::memset(data_, 0, sizeof(data_));
    const uint32_t used_words = NumBaseWords(new_size);
    const uint32_t source_word = base_offset >> 5u;
    const uint32_t shift = (base_offset & 31u) << 1u;
    for (uint32_t i = 0; i < used_words; ++i) {
      uint64_t value = words[source_word + i] >> shift;
      if (shift != 0) {
        value |= words[source_word + i + 1u] << (64u - shift);
      }
      data_[i] = value;
    }
    const uint32_t tail = new_size & 31u;
    if (tail != 0) {
      data_[used_words - 1u] &= (uint64_t{1} << (tail << 1u)) - 1u;
    }
    resize(new_size);
  }

  void swap(IdbaKmer &kmer) {
    if (this != &kmer) {
      for (unsigned i = 0; i < kNumUint64; ++i)
        std::swap(data_[i], kmer.data_[i]);
    }
  }

  uint32_t size() const { return data_[kNumUint64 - 1] >> (64 - kBitsForSize); }
  void resize(uint32_t new_size) {
    data_[kNumUint64 - 1] =
        ((data_[kNumUint64 - 1] << kBitsForSize) >> kBitsForSize) |
        (uint64_t(new_size) << (64 - kBitsForSize));
  }

  void clear() {
    uint32_t kmer_size = size();
    memset(data_, 0, sizeof(uint64_t) * kNumUint64);
    resize(kmer_size);
  }

  static uint32_t max_size() { return kMaxSize; }

  static const uint32_t kNumUint64 = kUint64PerIdbaKmerMaxK;
  static const uint32_t kBitsForSize =
      ((kNumUint64 <= 2) ? 6 : ((kNumUint64 <= 8) ? 8 : 16));
  static const uint32_t kBitsForIdbaKmer = (kNumUint64 * 64 - kBitsForSize);
  static const uint32_t kMaxSize = kBitsForIdbaKmer / 2;

 private:
  static const uint64_t *FixedHashAccumulators() {
    struct FixedAccumulators {
      FixedAccumulators() {
        const char *key = reinterpret_cast<const char *>(kKey);
        for (uint32_t k = 0; k <= 255u; ++k) {
          uint64_t words[9] = {};
          words[8] = uint64_t(k) << (64u - kBitsForSize);
          const uint32_t used = NumBaseWords(k);
          uint64_t value = uint64_t(sizeof(words)) * PRIME64_1;
          if (used <= 4u)
            value += XXH3_mix16B(
                reinterpret_cast<const char *>(words) + 32, key + 64, 0);
          if (used <= 3u)
            value += XXH3_mix16B(
                reinterpret_cast<const char *>(words) + 24, key + 80, 0);
          if (used <= 2u)
            value += XXH3_mix16B(
                reinterpret_cast<const char *>(words) + 16, key + 32, 0);
          if (used <= 5u)
            value += XXH3_mix16B(
                reinterpret_cast<const char *>(words) + 40, key + 48, 0);
          if (used <= 7u)
            value += XXH3_mix16B(
                reinterpret_cast<const char *>(words) + 56, key + 16, 0);
          values[k] = value;
        }
      }
      uint64_t values[256];
    };
    static const FixedAccumulators fixed;
    return fixed.values;
  }

  static uint32_t NumBaseWords(uint32_t kmer_size) {
    return (kmer_size + 31u) >> 5u;
  }

  uint64_t data_[kNumUint64];
};

namespace std {
template <>
inline void swap(IdbaKmer &kmer1, IdbaKmer &kmer2) {
  kmer1.swap(kmer2);
}
}  // namespace std

#endif
