#include "assembly/legacy_order_sort.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

static void Require(bool value, const char *message) {
  if (!value) { std::fprintf(stderr, "%s\n", message); std::exit(2); }
}

int main() {
  std::mt19937_64 random(47292);
  for (size_t size : {size_t{0}, size_t{1}, size_t{2}, size_t{63},
                     size_t{64}, size_t{65}, size_t{4097}, size_t{100000}}) {
    for (uint64_t base : {uint64_t{0}, (uint64_t{1} << 32u) - 129u,
                          uint64_t{1} << 40u}) {
      std::vector<uint64_t> keys(size), handles(size);
      std::iota(handles.begin(), handles.end(), uint64_t{0});
      for (size_t i = 0; i < size; ++i) keys[i] = base + i * 7u;
      std::shuffle(keys.begin(), keys.end(), random);
      std::shuffle(handles.begin(), handles.end(), random);
      auto reference = handles;
      std::sort(reference.begin(), reference.end(), [&](uint64_t a, uint64_t b) {
        return keys[a] < keys[b];
      });
      const uint64_t max_key = keys.empty() ? 0 : *std::max_element(keys.begin(), keys.end());
      Require(legacy_order::TrySortPackedStarts(&handles, max_key,
                  size == 0 ? 0 : size - 1u, [&](uint64_t h) { return keys[h]; }),
              "representable keys rejected");
      Require(handles == reference, "sorted handle order differs");
    }
  }
  // Handles can themselves exceed 32 bits; no truncation of oriented IDs.
  std::vector<uint64_t> handles{(uint64_t{1} << 32u) + 7u,
                                (uint64_t{1} << 32u) + 3u};
  Require(legacy_order::TrySortPackedStarts(&handles, 7u, handles[0],
              [](uint64_t h) { return h & 7u; }), "wide handles rejected");
  Require(handles[0] == (uint64_t{1} << 32u) + 3u, "wide handle truncated");

  // Use the full 64-bit packed word, including its top bit.
  handles = {(uint64_t{1} << 30u) + 1u, 0u};
  Require(legacy_order::TrySortPackedStarts(&handles,
              (uint64_t{1} << 33u) - 1u, handles[0], [](uint64_t h) {
                return h == 0 ? 0u : (uint64_t{1} << 33u) - 1u;
              }), "exact 64-bit fit rejected");
  Require(handles[0] == 0 && handles[1] == (uint64_t{1} << 30u) + 1u,
          "packed top bit changed unsigned ordering");

  handles = {3u, 0u, 2u, 1u};
  const auto original = handles;
  Require(!legacy_order::TrySortPackedStarts(&handles, uint64_t{1} << 63u,
              3u, [](uint64_t h) { return h; }), "overflow not rejected");
  Require(handles == original, "width fallback mutated input");
  Require(!legacy_order::TrySortPackedStarts(&handles, 1u, 3u,
              [](uint64_t h) { return h & 1u; }), "duplicate keys not rejected");
  std::puts("Packed legacy starts: exact order, 32/40-bit boundaries, wide handles, overflow and duplicate fallback passed");
}
