/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#ifndef MEGAHIT_KMER_COUNTER_H
#define MEGAHIT_KMER_COUNTER_H

#include <stdint.h>
#include <mutex>
#include <string>
#include <vector>
#include "base_engine.h"
#include "definitions.h"
#include "edge_counter.h"
#include "sequence/io/edge/edge_writer.h"
#include "sequence/sequence_package.h"
#include "utils/atomic_wrapper.h"

class MappedCountReads;

struct KmerCounterOption {
  unsigned k{21};
  int solid_threshold{2};
  double host_mem{0};
  int n_threads{0};
  std::string read_lib_file{};
  std::string output_prefix{"out"};
  int mem_flag{1};
};

class KmerCounter : public BaseSequenceSortingEngine {
 public:
  static const unsigned kSentinelValue = 4;
  static const uint32_t kSentinelOffset = 4294967295U;
  explicit KmerCounter(const KmerCounterOption &opt);
  ~KmerCounter() final;

 public:
  MemoryStat Initialize() override;

 protected:
  int64_t Lv0EncodeDiffBase(int64_t) override;
  void Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                         std::array<int64_t, kNumBuckets> *out) override;
  void Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                      int64_t seq_to) override;
  void Lv2ExtractSubString(OffsetFetcher &fetcher,
                           SubstrPtr substr_ptr) override;
  void Lv2Postprocess(int64_t start_index, int64_t end_index, int thread_id,
                      uint32_t *substr_ptr, unsigned bucket_id) override;
  void Lv0Postprocess() override;
  bool RunSpecializedMainLoop() override;
  bool Lv1SupportsDirectItems() const override;
  int64_t Lv1DirectWordsPerItem() const override;
  int64_t Lv1DirectAuxWordsPerItem() const override;
  int64_t Lv1DirectMemoryLimit() const override;
  unsigned Lv1BytesPerOffset() const override;
  int64_t Lv1AutoWorkspaceLimit() const override;
  int Lv2SortIgnoredLowBytes() const override;
  bool Lv1CanUseCompactCursor(int64_t seq_from,
                              int64_t seq_to) const override;

 private:
  // Short-read minimizer runs need only a read ID plus two byte-sized
  // coordinates.  Keeping this descriptor at six bytes is material at scale:
  // it is the only object retained between the single producer scan and the
  // per-shard exact sorts.
  struct CompactCountSegment {
    uint32_t read_id;
    uint8_t start;
    uint8_t length;
  } __attribute__((packed));
  static_assert(sizeof(CompactCountSegment) == 6,
                "compact count segment must not contain padding");

  void PackEdge(uint32_t *dest, uint32_t *item, int64_t counting);
  void ConfigurePackedReadOffsets(uint64_t max_read_len, uint64_t num_reads);
  void ConfigureRecordLayout(uint64_t max_read_len, uint64_t num_reads);
  int64_t EncodeReadOffset(int64_t read_id, unsigned offset,
                           unsigned strand) const;

  template <unsigned NWords>
  void Lv0CalcBucketSizeFor(int64_t seq_from, int64_t seq_to,
                            std::array<int64_t, kNumBuckets> *out);
  template <unsigned NWords, bool PackedOffsets>
  void Lv1FillOffsetsFor(OffsetFiller &filler, int64_t seq_from,
                         int64_t seq_to);
  template <unsigned NWords, bool PackedOffsets, bool CompactItems>
  void Lv1FillDirectItemsFor(OffsetFiller &filler, int64_t seq_from,
                             int64_t seq_to);
  template <unsigned NWords, bool PackedOffsets>
  void Lv1FillWideCompactItemsFor(OffsetFiller &filler, int64_t seq_from,
                                  int64_t seq_to);
  template <bool PackedOffsets, bool CompactItems, bool WideCompactLocator>
  void Lv2ExtractSubStringFor(OffsetFetcher &fetcher,
                              SubstrPtr substr_ptr);
  template <bool PackedOffsets, bool CompactItems, bool WideCompactLocator>
  void Lv2PostprocessFor(int64_t start_index, int64_t end_index,
                         int thread_id, uint32_t *substr_ptr);
  template <bool WideCompactLocator>
  uint64_t DecodeCompactLocator(const uint32_t *item) const;
  template <unsigned NWords>
  void MaterializeSegment(const CompactCountSegment &segment,
                          uint32_t *dest) const;
  template <unsigned NWords>
  bool RunSegmentedCountFor();
  template <unsigned NWords>
  bool RunExternalSegmentedCountFor();
  void UpdateLast0In(uint32_t read_id, uint32_t offset);
  void UpdateFirst0Out(uint32_t read_id, uint32_t offset);
  uint32_t LoadFirst0Out(uint32_t read_id) const;
  uint32_t LoadLast0In(uint32_t read_id) const;
  void WriteEndpointCandidate(const uint32_t *packed_edge,
                              bool missing_in, bool missing_out,
                              uint8_t verify_in_mask,
                              uint8_t verify_out_mask,
                              int thread_id);

 private:
  KmerCounterOption opt_;

  int words_per_edge_{};  // number of (32-bit) words needed to represent a
                          // (k+1)-mer
  int64_t words_per_substr_{};  // substrings to be sorted by GPU
  SeqPackage seq_pkg_;
  std::unique_ptr<MappedCountReads> mapped_count_reads_;
  bool packed_read_offsets_{};
  bool compact_items_{};
  uint32_t compact_key_mask_{};
  uint32_t compact_locator_high_mask_{};
  unsigned read_locator_shift_{};
  uint64_t read_offset_mask_{};
  bool segmented_count_enabled_{};
  std::vector<uint64_t> seq_bucket_histograms_;
  // Uninitialized on purpose: tens of millions of entries are 0xFF-filled in
  // parallel (with NUMA-local first touch) instead of serially constructed.
  struct ByteFreeDeleter {
    void operator()(AtomicWrapper<uint8_t> *p) const { std::free(p); }
  };
  struct WordFreeDeleter {
    void operator()(AtomicWrapper<uint32_t> *p) const { std::free(p); }
  };
  // Short-read offsets fit in one byte (segmented_count_enabled_ requires a
  // maximum read length <= 255).  The two arrays are updated only with
  // monotone min/max CAS operations, so byte atomics preserve the exact
  // result while removing six bytes of per-read resident state.  Generic
  // long-read inputs retain the legacy 32-bit path through separate storage.
  bool compact_read_flags_{};
  std::unique_ptr<AtomicWrapper<uint8_t>[], ByteFreeDeleter> first_0_out8_;
  std::unique_ptr<AtomicWrapper<uint8_t>[], ByteFreeDeleter> last_0_in8_;
  std::unique_ptr<AtomicWrapper<uint32_t>[], WordFreeDeleter> first_0_out32_;
  std::unique_ptr<AtomicWrapper<uint32_t>[], WordFreeDeleter> last_0_in32_;
  // stat
  EdgeMultiplicityRecorder edge_counter_;
  // output
  EdgeWriter edge_writer_;
  // A missing (k+2)-context is a necessary, but not sufficient, condition
  // for a solid (k+1)-edge endpoint.  Persist that strict superset while the
  // count group and its context histogram are hot.  SeqToSdbg later filters
  // it against the exact solid-edge index.  This replaces the all-record
  // residency formerly needed solely for the endpoint set difference.
  std::vector<std::unique_ptr<std::ofstream>> endpoint_candidate_files_;
};

#endif  // MEGAHIT_KMER_COUNTER_H
