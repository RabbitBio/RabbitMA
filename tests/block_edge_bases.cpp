#include "sdbg/block_edge_bases.h"
#include <cstdio>
#include <cstdlib>
#include <vector>

static void Require(bool condition) {
  if (!condition) std::abort();
}

int main() {
  BlockEdgeBases bases;
  const uint64_t u32 = std::numeric_limits<uint32_t>::max();
  const uint64_t nil = BlockEdgeBases::kNull;
  bases.Allocate(3, false);
  for (size_t i = 0; i < 3; ++i) bases.Store(i, i == 2 ? nil : i * (u32 - 1));
  Require(bases.bytes() == 12 && !bases.wide());
  Require(bases.Load(0) == 0 && bases.Load(1) == u32 - 1 && bases.Load(2) == nil);
  const std::vector<uint64_t> values{0, u32 - 1, u32, u32 + 1,
                                    uint64_t{1} << 40, nil - 253, nil};
  bases.Allocate(values.size(), true);
  for (size_t i = 0; i < values.size(); ++i) bases.Store(i, values[i]);
  for (size_t i = 0; i < values.size(); ++i) Require(bases.Load(i) == values[i]);
  Require(bases.bytes() == values.size() * 8 && bases.wide());
  // Exercise both crossing UINT32_MAX and exact overflow fallback boundaries.
  for (uint64_t base : {uint64_t{0}, u32 - 100, u32, u32 + 1, uint64_t{1} << 40}) {
    for (uint64_t delta = 0; delta < 1024; ++delta) {
      const uint8_t code = EncodeSimplePathDelta(base, base + delta);
      if (delta < 254) {
        Require(code != 0 && code != 255 && base + code - 1 == base + delta);
      } else Require(code == 255);
    }
    if (base) Require(EncodeSimplePathDelta(base, base - 1) == 255);
    Require(EncodeSimplePathDelta(base, nil) == 0);
  }
  Require(EncodeSimplePathDelta(nil, u32) == 255);
  Require(EncodeSimplePathDelta(nil - 253, nil - 1) == 253);
  Require(EncodeSimplePathDelta(0, nil - 1) == 255);
  BlockEdgeBases moved = std::move(bases);
  Require(moved.Load(2) == u32);
  moved.reset();
  Require(!moved && moved.bytes() == 0);
  std::puts("Block edge bases: boundary, null, overflow, and move checks passed");
}
