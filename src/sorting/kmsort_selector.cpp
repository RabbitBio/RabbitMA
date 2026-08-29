//
// Created by vout on 12/8/18.
//

#include "kmsort_selector.h"

#include <definitions.h>
#include <utils/utils.h>
#include <cassert>
#include <limits>
#include "kmlib/kmsort.h"

namespace {
template <int NWords, int NExtraWord, int IgnoredLowBytes,
          int IgnoredHighBytes>
struct Substr {
  uint32_t data[NWords + NExtraWord];
  static const int n_bytes =
      static_cast<int>(sizeof(data[0])) * NWords - IgnoredHighBytes -
      IgnoredLowBytes;
  bool operator<(const Substr &rhs) const {
    for (int i = 0; i < NWords; ++i) {
      uint32_t lhs = data[i];
      uint32_t other = rhs.data[i];
      if (i == 0 && IgnoredHighBytes != 0) {
        const uint32_t high_mask =
            std::numeric_limits<uint32_t>::max() >> (IgnoredHighBytes * 8);
        lhs &= high_mask;
        other &= high_mask;
      }
      if (i == NWords - 1 && IgnoredLowBytes != 0) {
        lhs >>= IgnoredLowBytes * 8;
        other >>= IgnoredLowBytes * 8;
      }
      if (lhs < other) {
        return true;
      } else if (lhs > other) {
        return false;
      }
    }
    return false;
  }
  int kth_byte(int k) const {
    k += IgnoredLowBytes;
    return data[NWords - 1 - k / sizeof(uint32_t)] >>
               (k % sizeof(uint32_t) * 8) &
           0xFF;
  }
} __attribute((packed));

constexpr int kMaxWords =
    (kMaxK * kBitsPerEdgeChar + 3 + 1 + kBitsPerMul + kBitsPerEdgeWord - 1) /
    kBitsPerEdgeWord;

template <int NWords, int NExtraWords, int IgnoredHighBytes>
struct LowByteSelector {
  static std::function<void(uint32_t *, int64_t)> Select(
      int ignored_low_bytes) {
#define RETURN_SORTER(IGNORED)                                                \
  return [](uint32_t *substr_ptr, int64_t n) {                                \
    auto ptr = reinterpret_cast<                                               \
        Substr<NWords, NExtraWords, IGNORED, IgnoredHighBytes> *>(substr_ptr); \
    kmlib::kmsort(ptr, ptr + n);                                               \
  }
  switch (ignored_low_bytes) {
    case 1:
      RETURN_SORTER(1);
    default:
      assert(ignored_low_bytes == 0);
      RETURN_SORTER(0);
  }
#undef RETURN_SORTER
  }
};

// All existing non-packed records have the fixed 16-bit bucket prefix.  Their
// only reachable low-byte layouts are zero (the usual case) and two (packed
// count metadata).  Keeping this specialization narrow avoids instantiating
// every radix sorter for combinations the program can never request.
template <int NWords, int NExtraWords>
struct LowByteSelector<NWords, NExtraWords, 2> {
  static std::function<void(uint32_t *, int64_t)> Select(
      int ignored_low_bytes) {
#define RETURN_SORTER(IGNORED)                                         \
  return [](uint32_t *substr_ptr, int64_t n) {                         \
    auto ptr = reinterpret_cast<Substr<NWords, NExtraWords, IGNORED,  \
                                        2> *>(substr_ptr);             \
    kmlib::kmsort(ptr, ptr + n);                                       \
  }
    switch (ignored_low_bytes) {
      case 2:
        RETURN_SORTER(2);
      default:
        assert(ignored_low_bytes == 0);
        RETURN_SORTER(0);
    }
#undef RETURN_SORTER
  }
};

template <int NWords, int NExtraWords, int IgnoredHighBytes>
std::function<void(uint32_t *, int64_t)> SelectSortingFuncHelper(
    int words_per_substr, int extra_words, int ignored_low_bytes) {
  assert(words_per_substr > 0 && words_per_substr <= NWords);
  assert(extra_words >= 0 && extra_words <= NExtraWords);
  if (words_per_substr < NWords) {
    return SelectSortingFuncHelper<(NWords - 1 > 1 ? NWords - 1 : 1),
                                   NExtraWords, IgnoredHighBytes>(
        words_per_substr, extra_words, ignored_low_bytes);
  }
  if (extra_words < NExtraWords) {
    return SelectSortingFuncHelper<
        NWords, (NExtraWords - 1 > 0 ? NExtraWords - 1 : 0),
        IgnoredHighBytes>(words_per_substr, extra_words, ignored_low_bytes);
  }
  return LowByteSelector<NWords, NExtraWords, IgnoredHighBytes>::Select(
      ignored_low_bytes);
}

template <int IgnoredHighBytes>
std::function<void(uint32_t *, int64_t)> SelectSortingFuncForHighBytes(
    int words_per_substr, int extra_words, int ignored_low_bytes) {
  return SelectSortingFuncHelper<kMaxWords, 2, IgnoredHighBytes>(
      words_per_substr, extra_words, ignored_low_bytes);
}

}  // namespace

std::function<void(uint32_t *, int64_t)> SelectSortingFunc(int words_per_substr,
                                                           int extra_words,
                                                           int ignored_low_bytes,
                                                           int ignored_high_bytes) {
  assert((ignored_high_bytes == 0 &&
          (ignored_low_bytes == 0 || ignored_low_bytes == 1)) ||
         (ignored_high_bytes == 2 &&
          (ignored_low_bytes == 0 || ignored_low_bytes == 2)));
  switch (ignored_high_bytes) {
    case 0:
      return SelectSortingFuncForHighBytes<0>(
          words_per_substr, extra_words, ignored_low_bytes);
    default:
      return SelectSortingFuncForHighBytes<2>(
          words_per_substr, extra_words, ignored_low_bytes);
  }
}
