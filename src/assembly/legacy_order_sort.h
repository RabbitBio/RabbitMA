#ifndef MEGAHIT_ASSEMBLY_LEGACY_ORDER_SORT_H
#define MEGAHIT_ASSEMBLY_LEGACY_ORDER_SORT_H

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "kmlib/kmsort.h"

namespace legacy_order {

struct SortProfile {
  double pack_seconds{0};
  double radix_seconds{0};
  double finish_seconds{0};
  const char *fallback_reason{"none"};
};

inline unsigned SignificantBits(uint64_t value) {
  return value == 0 ? 0u : 64u - __builtin_clzll(value);
}

// Sort by the exact SDBG start address, keeping the existing eight bytes per
// oriented handle. The address is not restricted to 32 bits: packing is used
// only when both structural bounds fit together in a uint64_t.
//
// A false return requires the caller to regenerate the ORIGINAL handle order
// and apply its historical comparator. In particular, duplicate start keys
// must not acquire a new handle-ID tie-break from the packed representation.
// This function may have changed the array before detecting that fallback.
template <class StartAddress>
bool TrySortPackedStarts(std::vector<uint64_t> *handles,
                         uint64_t max_start_address, uint64_t max_handle,
                         StartAddress start_address,
                         SortProfile *profile = nullptr) {
  using Clock = std::chrono::steady_clock;
  auto begin = profile ? Clock::now() : Clock::time_point{};
  if (profile) *profile = SortProfile{};
  if (handles->empty()) return true;
  const unsigned handle_bits = std::max(1u, SignificantBits(max_handle));
  if (handle_bits >= 64u ||
      max_start_address >
          (std::numeric_limits<uint64_t>::max() >> handle_bits)) {
    if (profile) {
      profile->pack_seconds = std::chrono::duration<double>(Clock::now() - begin).count();
      profile->fallback_reason = "bit-width";
    }
    return false;
  }
  const uint64_t handle_mask = (uint64_t{1} << handle_bits) - 1u;
  for (uint64_t &handle : *handles) {
    const uint64_t key = start_address(handle);
    if (handle > max_handle || key > max_start_address) {
      if (profile) {
        profile->pack_seconds = std::chrono::duration<double>(Clock::now() - begin).count();
        profile->fallback_reason = "bounds";
      }
      return false;
    }
    handle = (key << handle_bits) | handle;
  }
  if (profile) {
    const auto now = Clock::now();
    profile->pack_seconds = std::chrono::duration<double>(now - begin).count();
    begin = now;
  }

  kmlib::kmsort(handles->begin(), handles->end());
  if (profile) {
    const auto now = Clock::now();
    profile->radix_seconds = std::chrono::duration<double>(now - begin).count();
    begin = now;
  }
  for (size_t i = 1; i < handles->size(); ++i) {
    if (((*handles)[i - 1u] >> handle_bits) ==
        ((*handles)[i] >> handle_bits)) {
      if (profile) {
        profile->finish_seconds = std::chrono::duration<double>(Clock::now() - begin).count();
        profile->fallback_reason = "duplicate-start";
      }
      return false;
    }
  }
  for (uint64_t &handle : *handles) handle &= handle_mask;
  if (profile) {
    profile->finish_seconds = std::chrono::duration<double>(Clock::now() - begin).count();
  }
  return true;
}

}  // namespace legacy_order

#endif  // MEGAHIT_ASSEMBLY_LEGACY_ORDER_SORT_H
