#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "iterate/query_anchor.h"

namespace {

// Independent reference: reconstruct every candidate from unpacked bases.
query_anchor::SelectedAnchor Reference(const std::vector<uint8_t> &bases,
                                      unsigned a, unsigned w, unsigned bits,
                                      unsigned distance) {
  query_anchor::SelectedAnchor selected{};
  uint64_t best_hash = std::numeric_limits<uint64_t>::max();
  for (unsigned offset = 0; offset + a <= w; ++offset) {
    uint64_t key = 0;
    for (unsigned i = 0; i < a; ++i) key = (key << 2u) | bases[offset + i];
    const uint64_t hash = ReadAnchorRankHash(key);
    if (hash <= best_hash) {
      best_hash = hash;
      selected.key = key;
      selected.offset = offset;
    }
  }
  if (bits >= 2 && selected.offset >= distance) {
    selected.context |= bases[selected.offset - distance];
  }
  if (bits >= 4 && selected.offset + a + distance <= w) {
    selected.context |= bases[selected.offset + a + distance - 1] << 2u;
  }
  return selected;
}

void Require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

void Compare(const query_anchor::SelectedAnchor &expected,
             const query_anchor::SelectedAnchor &actual) {
  Require(expected.key == actual.key, "anchor key differs");
  Require(expected.offset == actual.offset, "anchor offset/tie differs");
  Require(expected.context == actual.context, "anchor context differs");
}

template <class KmerType>
uint64_t Check(unsigned flank_len) {
  uint64_t checked = 0;
  std::mt19937_64 random(1291 + flank_len);
  for (unsigned sample = 0; sample < 80; ++sample) {
    std::vector<uint8_t> bases(flank_len);
    for (unsigned i = 0; i < flank_len; ++i) {
      // Homopolymers and periodic sequences force equal-ranked candidates,
      // including ties at each end of the reverse window.
      bases[i] = sample < 4 ? sample :
                 sample < 8 ? (i + sample) % 4 : random() % 4;
    }
    if (sample >= 8 && sample < 16) {
      for (unsigned i = 0; i < flank_len / 2; ++i) {
        bases[flank_len - i - 1] = bases[i] ^ 3u;
      }
    }
    std::vector<uint8_t> reverse_bases(bases.rbegin(), bases.rend());
    for (auto &base : reverse_bases) base ^= 3u;
    KmerType forward;
    for (uint8_t base : bases) forward.ShiftAppend(base, flank_len);
    KmerType reverse = forward;
    reverse.ReverseComplement(flank_len);
    for (unsigned w : {1u, 15u, 16u, 17u, 19u, 21u, 31u, 32u, 33u,
                       40u, 48u, 63u, 64u, 65u, 80u, 94u, flank_len}) {
      if (w > flank_len) continue;
      if (w <= 64u) {
        for (unsigned strand = 0; strand < 2; ++strand) {
          const auto &oriented = strand == 0 ? forward : reverse;
          Kmer<4, uint32_t> expected;
          const auto &unpacked = strand == 0 ? bases : reverse_bases;
          for (unsigned i = 0; i < w; ++i) expected.ShiftAppend(unpacked[i], w);
          Require(expected == query_anchor::PackedWindow(oriented, w),
                  "window filter bytes or zero padding differ");
        }
      }
      for (unsigned a = 1; a <= std::min(31u, w); ++a) {
        // The public index accepts at most 64 candidate offsets.
        if (w - a >= 64u) continue;
        for (unsigned bits : {0u, 2u, 4u}) {
          for (unsigned distance : {1u, 2u, 7u, w}) {
            const auto expected_forward = Reference(bases, a, w, bits, distance);
            const auto expected_reverse = Reference(reverse_bases, a, w, bits, distance);
            const auto actual = query_anchor::SelectPair(
                forward, reverse, flank_len, a, w, bits, distance);
            Compare(expected_forward, actual[0]);
            Compare(expected_reverse, actual[1]);
            Compare(expected_forward, query_anchor::Select(forward, a, w, bits, distance));
            Compare(expected_reverse, query_anchor::Select(reverse, a, w, bits, distance));
            if (sample < 4) {
              Require(actual[0].offset == w - a && actual[1].offset == w - a,
                      "homopolymer did not keep right-most tie on both strands");
            }
            ++checked;
          }
        }
      }
    }
  }
  return checked;
}
}  // namespace

int main() {
  uint64_t checked = 0;
  checked += Check<Kmer<1, uint32_t>>(16);
  checked += Check<Kmer<1, uint64_t>>(32);
  checked += Check<Kmer<3, uint32_t>>(40);
  checked += Check<Kmer<2, uint64_t>>(40);
  checked += Check<Kmer<4, uint32_t>>(64);
  checked += Check<Kmer<2, uint64_t>>(64);
  checked += Check<Kmer<5, uint32_t>>(80);
  checked += Check<Kmer<3, uint64_t>>(94);
  checked += Check<Kmer<4, uint64_t>>(120);
  checked += Check<Kmer<9, uint32_t>>(140);
  std::cout << checked << " anchor pairs match unpacked reference, including "
               "ties, context, both word widths and distinct endpoint windows\n";
}
