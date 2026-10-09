#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "assembly/packed_endpoint_key.h"

int main() {
  phmap::flat_hash_map<uint64_t, uint32_t> reference;
  endpoint_map::PackedMap packed;
  std::vector<uint64_t> keys{0, 1, UINT32_MAX, uint64_t{1} << 32u,
      (uint64_t{1} << 32u) + 1u, uint64_t{1} << 48u, UINT64_MAX};
  // Equal low words with different high words must never alias.
  for (uint64_t high = 0; high < 256; ++high) keys.push_back((high << 32u) | 123u);
  std::mt19937_64 rng(20261005);
  for (unsigned i = 0; i < 200000; ++i) keys.push_back(rng());
  reference.reserve(keys.size());
  packed.reserve(keys.size());
  for (size_t i = 0; i < keys.size(); ++i) {
    const uint64_t key = keys[i];
    assert(endpoint_map::PackedKey(key).value() == key);
    assert(endpoint_map::PackedKeyHash{}(endpoint_map::PackedKey(key)) ==
           phmap::Hash<uint64_t>{}(key));
    assert(endpoint_map::PackedKeyHash{}(endpoint_map::PackedKey(key)) ==
           endpoint_map::PackedKeyHash{}(key));
    reference[key] = static_cast<uint32_t>(i);
    packed[key] = static_cast<uint32_t>(i);
  }
  for (const uint64_t key : keys) assert(packed.at(key) == reference.at(key));
  for (size_t i = 0; i < keys.size(); i += 3) {
    assert(packed.erase(keys[i]) == reference.erase(keys[i]));
    assert(packed.count(keys[i]) == 0);
  }
  for (size_t i = 0; i < keys.size(); i += 3) {
    reference[keys[i]] = static_cast<uint32_t>(i ^ 999u);
    packed[keys[i]] = static_cast<uint32_t>(i ^ 999u);
  }
  // Existing-entry mutation is the path used by endpoint remapping.
  for (const uint64_t key : keys) {
    reference.find(key)->second ^= 42u;
    packed.find(key)->second ^= 42u;
  }
  assert(reference.size() == packed.size());
  for (const auto &entry : reference) assert(packed.at(entry.first) == entry.second);
  packed.clear();
  assert(packed.empty());
  packed[UINT64_MAX] = UINT32_MAX;
  assert(packed.at(UINT64_MAX) == UINT32_MAX);
  std::cout << "PASS: 64-bit identity/hash, lookup, erase/reinsert, remap, clear; "
            << sizeof(decltype(reference)::value_type) << " -> "
            << sizeof(endpoint_map::PackedMap::value_type) << " bytes per slot\n";
}
