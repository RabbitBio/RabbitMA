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

#ifndef MEGAHIT_SEQ_TO_SDBG_H
#define MEGAHIT_SEQ_TO_SDBG_H

#include <sequence/sequence_package.h>
#include <stdint.h>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>
#include "base_engine.h"
#include "sdbg/sdbg_writer.h"
#include "sequence/io/edge/edge_bucket_histogram.h"

struct EdgeIoMetadata;

struct Seq2SdbgOption {
  double host_mem{0};
  int n_threads{0};
  unsigned k{0};
  unsigned k_from{0};
  std::string contig;
  std::string bubble_seq;
  std::string addi_contig;
  std::string local_contig;
  std::string input_prefix;
  std::string output_prefix;
  int mem_flag{1};
  bool need_mercy{false};
};

class SeqToSdbg : public BaseSequenceSortingEngine {
 public:
  // binary search look up table
  static const unsigned kMaxLookUpPrefixLength = 12;
  static const unsigned kSentinelValue = 4;
  static const unsigned kBWTCharNumBits = 3;

  explicit SeqToSdbg(const Seq2SdbgOption &opt)
      : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
        opt_(opt) {}
  ~SeqToSdbg() override;

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
  bool Lv1SupportsDirectItems() const override;
  int64_t Lv1DirectWordsPerItem() const override;
  int64_t Lv1DirectAuxWordsPerItem() const override;
  int64_t Lv1DirectMemoryLimit() const override;
  bool Lv1AllowsPartialDirectItems() const override;
  int64_t Lv1AutoWorkspaceLimit() const override;
  bool Lv1UseWriteCombine() const override;
  int Lv2SortIgnoredLowBytes() const override;
  int Lv2SortIgnoredHighBytes() const override;
  bool Lv2NeedsEmptyBucketPostprocess(unsigned bucket) const override;
  bool Lv0BuildBalancedRanges(
      std::vector<std::pair<int64_t, int64_t>> *ranges) const override;
  bool Lv1CanUseCompactCursor(int64_t seq_from,
                              int64_t seq_to) const override;
  bool Lv2DeferPostprocess() const override;
  void Lv2PostprocessDeferred() override;

 private:
  // input options
  Seq2SdbgOption opt_;
  int64_t full_words_per_substr_{};
  int64_t words_per_substr_{};
  bool bucket_packed_records_{false};
  bool real_only_records_{false};
  unsigned bucket_packed_base_bits_{0};
  unsigned bucket_packed_total_bits_{0};

  // big arrays
  SeqPackage seq_pkg_;
  std::vector<mul_t> multiplicity;

  // Fixed-length edges can be consumed twice as sequential raw streams:
  // once for exact bucket counts and once to materialize the final sortable
  // records.  This avoids retaining a second packed copy of the complete edge
  // set alongside the direct-item array.
  struct StreamEdgeChunk {
    uint32_t file_id;
    uint64_t first_record;
    uint32_t num_records;
  };
  bool stream_input_edges_{false};
  bool stream_input_unordered_{false};
  bool stream_compact_mercy_index_{false};
  uint32_t stream_edge_length_{0};
  uint32_t stream_words_per_edge_{0};
  int64_t stream_num_edges_{0};
  std::vector<StreamEdgeChunk> stream_edge_chunks_;
  std::vector<size_t> stream_file_chunk_offsets_;
  std::vector<int> stream_edge_fds_;
  EdgeBucketHistogram stream_bucket_histogram_;
  bool use_stream_bucket_histogram_{false};

  // Exact source/sink sentinel records, grouped by the same 16-bit bucket as
  // the real records.  Count emits only a strict endpoint-candidate superset;
  // GenMercyEdges filters it through the exact solid-edge index and dedupes
  // the result.  Once this is available every real bucket is independent and
  // the direct workspace can safely be reused in bounded batches.
  bool precomputed_endpoint_tips_{false};
  std::vector<uint64_t> endpoint_tip_bucket_begin_;
  std::vector<uint32_t> endpoint_tip_items_;

  // The sorting engine only requires locators to be monotonically increasing
  // inside each input partition.  Keeping the sequence ID and local offset in
  // the locator avoids recovering an ID from a global base offset for every
  // Lv2 item.  Very large inputs transparently retain the legacy encoding.
  bool packed_seq_offsets_{false};
  unsigned seq_locator_shift_{0};
  uint64_t seq_offset_mask_{0};

  // output
  SdbgWriter sdbg_writer_;
  int64_t AutomaticDirectWorkspaceLimit(uint64_t retained_bytes = 0) const;
  int64_t BoundedTransientWorkspaceLimit(uint64_t retained_bytes,
                                         bool report) const;
  uint64_t CurrentRetainedBytes() const;
  void ConfigurePackedSeqOffsets();
  bool ConfigureStreamedEdgeInput(const EdgeIoMetadata &metadata);
  void RetainMercyEdgesForStreamedInput(size_t original_edge_count);
  void ReadStreamedEdgeChunk(const StreamEdgeChunk &chunk,
                             std::vector<uint32_t> *records) const;
  static uint16_t ExtractRawBaseWindow8(const uint32_t *edge,
                                        unsigned base_offset);
  static uint16_t ReverseComplementWindow8(uint16_t bases);
  static unsigned ExtractRawBase(const uint32_t *edge, unsigned offset);
  template <unsigned NWords, bool BucketPacked>
  void MaterializeStreamedEdge(const uint32_t *edge, unsigned strand,
                               uint32_t *item) const;
  template <unsigned NWords, bool BucketPacked>
  void Lv1FillStreamedEdgesFor(OffsetFiller &filler, int64_t chunk_from,
                               int64_t chunk_to);
  int64_t EncodeSeqOffset(int64_t seq_id, unsigned offset,
                          unsigned strand) const;
  template <bool PackedOffsets>
  void Lv1FillOffsetsFor(OffsetFiller &filler, int64_t seq_from,
                         int64_t seq_to);
  template <unsigned NWords, bool BucketPacked>
  void Lv1FillDirectItemsFor(OffsetFiller &filler, int64_t seq_from,
                             int64_t seq_to);
  template <unsigned NWords, bool BucketPacked>
  void MaterializeItem(const SeqPackage::SeqView &seq_view, int offset,
                       unsigned strand, uint32_t *item) const;
  void AppendPackedMetadata(uint32_t *item, uint32_t metadata) const;
  uint32_t ExtractPackedField(const uint32_t *item, unsigned bit_offset,
                              unsigned bit_width) const;
  uint32_t ExtractPackedMetadata(const uint32_t *item) const;
  bool IsDiffPackedKMinusOneMer(const uint32_t *lhs,
                                const uint32_t *rhs) const;
  int ExtractPackedA(const uint32_t *item) const;
  void BuildPackedTipLabel(const uint32_t *item, unsigned bucket_id,
                           uint32_t *label) const;
  template <unsigned NWords, bool PackedOffsets, bool BucketPacked>
  void Lv2ExtractSubStringFor(OffsetFetcher &fetcher, SubstrPtr substr_ptr);
  template <unsigned NWords, bool BucketPacked>
  void Lv2PostprocessRealOnly();
  void GenMercyEdges();
  template <unsigned NWords>
  void GenMercyEdgesForK();
};
#endif  // MEGAHIT_SEQ_TO_SDBG_H
