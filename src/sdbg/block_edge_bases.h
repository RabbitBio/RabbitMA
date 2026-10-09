#ifndef MEGAHIT_BLOCK_EDGE_BASES_H
#define MEGAHIT_BLOCK_EDGE_BASES_H

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>

// One absolute base per symbol per block, not one wide ID per graph edge.
// Small graphs retain their original four-byte directory entries. Null and
// edge UINT32_MAX must remain distinct in the wide representation.
class BlockEdgeBases {
 public:
  static constexpr uint64_t kNull = std::numeric_limits<uint64_t>::max();

  void reset() {
    narrow_.reset();
    wide_.reset();
    bytes_ = 0;
  }
  void Allocate(size_t entries, bool wide) {
    reset();
    const size_t width = wide ? sizeof(uint64_t) : sizeof(uint32_t);
    if (entries > std::numeric_limits<size_t>::max() / width) {
      throw std::bad_alloc();
    }
    if (wide) wide_.reset(new uint64_t[entries]);
    else narrow_.reset(new uint32_t[entries]);
    bytes_ = entries * width;
  }
  explicit operator bool() const { return narrow_ || wide_; }
  size_t bytes() const { return bytes_; }
  bool wide() const { return wide_ != nullptr; }

  void Store(size_t i, uint64_t value) {
    if (wide_) wide_[i] = value;
    else {
      assert(value == kNull || value < std::numeric_limits<uint32_t>::max());
      narrow_[i] = static_cast<uint32_t>(value);
    }
  }
  uint64_t Load(size_t i) const {
    if (wide_) return wide_[i];
    const uint32_t value = narrow_[i];
    return value == std::numeric_limits<uint32_t>::max() ? kNull : value;
  }

 private:
  std::unique_ptr<uint32_t[]> narrow_;
  std::unique_ptr<uint64_t[]> wide_;
  size_t bytes_{0};
};

inline uint8_t EncodeSimplePathDelta(uint64_t base, uint64_t next) {
  if (next == BlockEdgeBases::kNull) return 0u;
  if (base != BlockEdgeBases::kNull && next >= base && next - base < 254u) {
    return static_cast<uint8_t>(next - base + 1u);
  }
  return 255u;  // Exact live-topology fallback, never a truncated address.
}

#endif
