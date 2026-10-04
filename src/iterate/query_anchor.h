#ifndef MEGAHIT_ITERATE_QUERY_ANCHOR_H_
#define MEGAHIT_ITERATE_QUERY_ANCHOR_H_

#include <array>
#include <cassert>
#include <cstdint>

#include "sequence/io/read_anchor_positions.h"
#include "sequence/kmer.h"

namespace query_anchor {

struct SelectedAnchor {
  uint64_t key;
  unsigned offset;
  uint8_t context;
};

// Kmer stores bases most-significant first in either 32- or 64-bit words.
// An anchor (1..31 bases) fits in one uint64_t; no base-by-base seed build is
// needed, even when it straddles two uint32_t words.
template <class KmerType>
uint64_t PrefixKey(const KmerType &kmer, unsigned length) {
  assert(length > 0u && length <= 31u && length <= KmerType::max_size());
  uint64_t prefix = kmer.data()[0];
  if (sizeof(typename KmerType::word_type) == sizeof(uint32_t)) {
    prefix <<= 32u;
    if (length > 16u) prefix |= kmer.data()[1];
  }
  return prefix >> (64u - 2u * length);
}

template <class KmerType>
uint8_t Context(const KmerType &kmer, unsigned offset, unsigned anchor_len,
                unsigned window_len, unsigned bits, unsigned distance) {
  uint8_t context = 0u;
  if (bits >= 2u && offset >= distance) {
    context |= kmer.GetBase(offset - distance);
  }
  if (bits >= 4u && offset + anchor_len + distance <= window_len) {
    context |= static_cast<uint8_t>(
        kmer.GetBase(offset + anchor_len + distance - 1u) << 2u);
  }
  return context;
}

template <class KmerType>
SelectedAnchor Select(const KmerType &kmer, unsigned anchor_len,
                      unsigned window_len, unsigned context_bits,
                      unsigned context_distance) {
  assert(anchor_len <= window_len && window_len <= KmerType::max_size());
  const uint64_t mask = (uint64_t{1} << (2u * anchor_len)) - 1u;
  uint64_t key = PrefixKey(kmer, anchor_len);
  uint64_t best_hash = ReadAnchorRankHash(key);
  SelectedAnchor selected{key, 0u, 0u};
  for (unsigned i = anchor_len; i < window_len; ++i) {
    key = ((key << 2u) | kmer.GetBase(i)) & mask;
    const uint64_t hash = ReadAnchorRankHash(key);
    // The read-side minimizer chooses the right-most equal minimum.
    if (hash <= best_hash) {
      best_hash = hash;
      selected.key = key;
      selected.offset = i + 1u - anchor_len;
    }
  }
  selected.context = Context(kmer, selected.offset, anchor_len, window_len,
                             context_bits, context_distance);
  return selected;
}

// The reverse query is a prefix of the *full reverse-complement flank*, not
// necessarily the reverse complement of the forward query window. Only fuse
// the candidate streams when both windows cover the entire flank.
template <class KmerType>
std::array<SelectedAnchor, 2> SelectPair(
    const KmerType &forward, const KmerType &reverse, unsigned flank_len,
    unsigned anchor_len, unsigned window_len, unsigned context_bits,
    unsigned context_distance) {
  assert(window_len <= flank_len && flank_len <= KmerType::max_size());
  assert(anchor_len > 0u && anchor_len <= window_len && anchor_len <= 31u);
  if (window_len != flank_len) {
    return {{Select(forward, anchor_len, window_len, context_bits,
                    context_distance),
             Select(reverse, anchor_len, window_len, context_bits,
                    context_distance)}};
  }
  const uint64_t mask = (uint64_t{1} << (2u * anchor_len)) - 1u;
  const unsigned reverse_shift = 2u * (anchor_len - 1u);
  const unsigned last_offset = window_len - anchor_len;
  uint64_t key = PrefixKey(forward, anchor_len);
  uint64_t reverse_key = kmlib::bit::ReverseComplement<2>(key) >>
                         (64u - 2u * anchor_len);
  uint64_t best_hash = ReadAnchorRankHash(key);
  uint64_t best_reverse_hash = ReadAnchorRankHash(reverse_key);
  std::array<SelectedAnchor, 2> selected{{
      {key, 0u, 0u}, {reverse_key, last_offset, 0u}}};
  for (unsigned offset = 1u; offset <= last_offset; ++offset) {
    const uint64_t base = forward.GetBase(offset + anchor_len - 1u);
    key = ((key << 2u) | base) & mask;
    reverse_key = (reverse_key >> 2u) | ((base ^ 3u) << reverse_shift);
    const uint64_t hash = ReadAnchorRankHash(key);
    const uint64_t reverse_hash = ReadAnchorRankHash(reverse_key);
    if (hash <= best_hash) {
      best_hash = hash;
      selected[0].key = key;
      selected[0].offset = offset;
    }
    // Reverse offsets run from right to left. Keep the first equal minimum
    // so this direction also retains the historical right-most tie rule.
    if (reverse_hash < best_reverse_hash) {
      best_reverse_hash = reverse_hash;
      selected[1].key = reverse_key;
      selected[1].offset = last_offset - offset;
    }
  }
  selected[0].context = Context(forward, selected[0].offset, anchor_len,
                                window_len, context_bits, context_distance);
  selected[1].context = Context(reverse, selected[1].offset, anchor_len,
                                window_len, context_bits, context_distance);
  return selected;
}

// Preserve the exact bytes hashed by ActiveWindowSet, including zero padding.
// Copy packed words instead of shifting the whole window once per base.
template <class KmerType>
Kmer<4, uint32_t> PackedWindow(const KmerType &kmer, unsigned length) {
  assert(length > 0u && length <= 64u && length <= KmerType::max_size());
  uint32_t words[4]{};
  const unsigned used_words = (length + 15u) / 16u;
  for (unsigned i = 0; i < used_words; ++i) {
    if (sizeof(typename KmerType::word_type) == sizeof(uint32_t)) {
      words[i] = kmer.data()[i];
    } else {
      words[i] = static_cast<uint32_t>(
          static_cast<uint64_t>(kmer.data()[i / 2u]) >>
          ((i & 1u) == 0u ? 32u : 0u));
    }
  }
  return Kmer<4, uint32_t>(words, 0u, length);
}

}  // namespace query_anchor

#endif  // MEGAHIT_ITERATE_QUERY_ANCHOR_H_
