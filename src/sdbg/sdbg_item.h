//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_ITEM_H
#define MEGAHIT_SDBG_ITEM_H

#include <cstdint>
#include <type_traits>
#include "sdbg_def.h"

/**
 * Pack w, last, tip and multiplicity (small) into one struct of 16 bits
 */
struct SdbgItem {
  SdbgItem() = default;
  SdbgItem(uint8_t w, uint8_t last, uint8_t tip, small_mul_t mul)
      : w(w), last(last), tip(tip), padding(0), mul(mul) {}
  uint8_t w : kBitsPerWChar;
  uint8_t last : 1;
  uint8_t tip : 1;
  // This object is serialized verbatim.  Keep the otherwise-unused bits
  // deterministic instead of leaking indeterminate stack bits into the file.
  uint8_t padding : 6 - kBitsPerWChar;
  small_mul_t mul : sizeof(small_mul_t) * 8;
  static_assert(sizeof(small_mul_t) <= 2, "");
};

static_assert(sizeof(SdbgItem) == sizeof(uint16_t),
              "SDBG on-disk item layout must remain 16 bits");
static_assert(std::is_trivially_copyable<SdbgItem>::value,
              "SDBG items must remain safe to serialize as object bytes");

#endif  // MEGAHIT_SDBG_ITEM_H
