//
// Created by vout on 11/5/18.
//

#ifndef MEGAHIT_SDBG_RAW_CONTENT_H
#define MEGAHIT_SDBG_RAW_CONTENT_H

#include "kmlib/kmcompactvector.h"
#include "sdbg_def.h"
#include "sdbg_meta.h"

/**
 * The raw (non-indexed) data of a SDBG
 */
struct SdbgRawContent {
  // Eight packed words per directory entry keep the rank directory compact.
  // Large multiplicities are rare, so this gives the best measured balance
  // between directory traffic and the occasional local popcount scan.
  static constexpr size_t kLargeMulRankBlockEdges = 512u;
  SdbgRawContent() = default;
  ~SdbgRawContent() = default;
  SdbgMeta meta;
  kmlib::CompactVector<kAlphabetSize, uint64_t> w;
  kmlib::CompactVector<1, uint64_t> last, tip;
  std::vector<small_mul_t> small_mul;
  std::vector<label_word_t> tip_lables;
  // Sparse multiplicities are represented in edge order.  The old hash map
  // required a serial insertion for every exceptional edge after parallel
  // loading and retained substantial empty-table slack.  A bit directory plus
  // block ranks gives exact O(1) lookup with compact contiguous storage.
  std::vector<uint64_t> large_mul_bits;
  std::vector<uint64_t> large_mul_rank;
  std::vector<mul_t> large_mul_values;
  std::vector<mul_t> full_mul;
};

void LoadSdbgRawContent(SdbgRawContent *raw_content,
                        const std::string &file_prefix);

#endif  // MEGAHIT_SDBG_RAW_CONTENT_H
