#ifndef MEGAHIT_ASSEMBLY_PACKED_ENDPOINT_KEY_H
#define MEGAHIT_ASSEMBLY_PACKED_ENDPOINT_KEY_H

#include <cstdint>
#include <cstddef>
#include <utility>
#include "parallel_hashmap/phmap.h"

namespace endpoint_map {

// Retain every address bit while requiring only four-byte alignment. A
// (uint64_t, uint32_t) hash slot otherwise pads its twelve payload bytes to
// sixteen. Two naturally aligned words avoid packed/unaligned C++ objects.
struct PackedKey {
  PackedKey(uint64_t value = 0)
      : low(static_cast<uint32_t>(value)),
        high(static_cast<uint32_t>(value >> 32u)) {}
  uint64_t value() const {
    return (uint64_t{high} << 32u) | low;
  }
  bool operator==(const PackedKey &rhs) const {
    return value() == rhs.value();
  }
  uint32_t low;
  uint32_t high;
};

struct PackedKeyHash {
  using is_transparent = void;
  size_t operator()(const PackedKey &key) const {
    // Preserve the old uint64 key hash exactly, including the high bits.
    return phmap::Hash<uint64_t>{}(key.value());
  }
  size_t operator()(uint64_t key) const {
    return phmap::Hash<uint64_t>{}(key);
  }
};

struct PackedKeyEqual {
  using is_transparent = void;
  bool operator()(const PackedKey &lhs, const PackedKey &rhs) const {
    return lhs == rhs;
  }
  bool operator()(const PackedKey &lhs, uint64_t rhs) const {
    return lhs.value() == rhs;
  }
  bool operator()(uint64_t lhs, const PackedKey &rhs) const {
    return lhs == rhs.value();
  }
};

using PackedMap = phmap::flat_hash_map<PackedKey, uint32_t, PackedKeyHash,
                                       PackedKeyEqual>;
static_assert(sizeof(PackedKey) == 8, "Endpoint keys must retain 64 bits");
static_assert(alignof(PackedKey) == alignof(uint32_t), "Four-byte key alignment");
static_assert(sizeof(PackedMap::value_type) == 12, "No tail padding in endpoint slots");

}  // namespace endpoint_map

#endif  // MEGAHIT_ASSEMBLY_PACKED_ENDPOINT_KEY_H
