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

#ifndef MEGAHIT_SEQUENCE_PACKAGE_H
#define MEGAHIT_SEQUENCE_PACKAGE_H

#include <cassert>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>
#include "kmlib/kmbit.h"
#include "kmlib/kmcompactvector.h"
#include "utils/mutex.h"
#include "utils/startup_affinity.h"
#include "utils/utils.h"

/**
 * @brief hold a set of sequences
 */
template <class WordType = unsigned long>
class SequencePackage {
 public:
  /**
   * The sequence view of a sequence
   */
  class SeqView {
   public:
    using word_type = WordType;

    SeqView(const SequencePackage *pkg, size_t seq_id)
        : package_(pkg), seq_id_(seq_id) {
      assert(seq_id < pkg->seq_count());
      pkg->GetSeqBounds(seq_id, &start_pos_, &length_);
    }

    unsigned length() const { return length_; }

    std::pair<const WordType *, unsigned> raw_address(
        unsigned offset = 0) const {
      assert(offset < length());
      size_t index = start_pos_ + offset;
      return {package_->sequences_.data() + index / kBasesPerWord,
              static_cast<unsigned>(index % kBasesPerWord)};
    }

    size_t base_at(unsigned index) const {
      assert(index < length());
      return package_->sequences_[start_pos_ + index];
    }

    size_t id() const { return seq_id_; }

    size_t full_offset_in_pkg() const { return start_pos_; }

   private:
    const SequencePackage *package_;
    size_t seq_id_;
    uint64_t start_pos_{};
    unsigned length_{};
  };

  using TWord = WordType;
  using TVector = kmlib::CompactVector<2, TWord, kmlib::kBigEndian>;
  using TAddress = std::pair<const TWord *, unsigned>;
  const static unsigned kBasesPerWord = TVector::kBasesPerWord;

 public:
  SequencePackage() {
    Clear();
    for (int i = 0; i < 10; ++i) {
      dna_map_[static_cast<int>("ACGTNacgtn"[i])] = "0123201232"[i] - '0';
    }
  }

  void Clear() {
    sequences_.clear();
    start_pos_.clear();
    start_pos_.push_back(0);
    sparse_gap_block_begin_.clear();
    sparse_gap_before_block_.clear();
    sparse_read_gaps_.clear();
    read_gap16_.clear();
    read_gap32_.clear();
    compact_seq_count_ = 0;
    compact_base_count_ = 0;
    pos_to_id_.clear();
    max_len_ = 0;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
  }

  // Clear both logical content and reserved backing storage.  Large staged
  // pipelines need this stronger operation when a source index has served its
  // purpose and the next phase must not inherit its multi-GiB capacity.
  void ReleaseStorage() {
    sequences_.release();
    std::vector<uint64_t>().swap(start_pos_);
    std::vector<uint32_t>().swap(sparse_gap_block_begin_);
    std::vector<uint32_t>().swap(sparse_gap_before_block_);
    std::vector<SparseReadGap>().swap(sparse_read_gaps_);
    std::vector<uint16_t>().swap(read_gap16_);
    std::vector<uint32_t>().swap(read_gap32_);
    std::vector<uint64_t>().swap(pos_to_id_);
    start_pos_.push_back(0);
    compact_seq_count_ = 0;
    compact_base_count_ = 0;
    max_len_ = 0;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
  }

  void ReserveBases(size_t num_bases) { sequences_.reserve(num_bases); }

  void ReserveSequences(size_t num_seq) { start_pos_.reserve(num_seq + 1); }

  /**
   * Replace an empty package with fixed-length sequences stored in strided
   * compact records.  Each source record begins with the sequence bases in
   * the same big-endian 2-bit representation used by TVector; words after
   * `len` bases (for example edge multiplicity) are ignored.
   *
   * Packing by destination-word chunks gives every worker exclusive output
   * words, so sequence boundaries need neither atomics nor a padded in-memory
   * representation.  This is intentionally an assignment operation: the
   * sorted edge loader is the first producer of seq_pkg_, while subsequent
   * mercy/contig inputs continue to use the normal append path.
   */
  void AssignFixedLengthCompactSequences(const TWord *records,
                                         size_t num_records, unsigned len,
                                         unsigned record_stride,
                                         int num_threads) {
    assert(seq_count() == 0 && base_count() == 0);
    assert(len > 0);
    assert(record_stride >= DivCeiling(len, kBasesPerWord));
    if (num_records == 0) {
      return;
    }

    const size_t total_bases = num_records * static_cast<size_t>(len);
    sequences_.resize(total_bases);
    const size_t sequence_bytes = sequences_.word_count() * sizeof(TWord);
    if (num_threads > 1 &&
        InterleaveMemoryPages(sequences_.data(), sequence_bytes)) {
      // resize() supplied the zero state expected by packed assignments.
      // Discarding after mbind retains that state lazily and avoids pinning
      // every page to the serial allocator thread.
      DiscardMemoryPages(sequences_.data(), sequence_bytes);
    }
    AdviseHugePages(sequences_.data(), sequence_bytes);
    fixed_len_ = len;
    num_fixed_len_ = num_records;
    start_pos_.back() = total_bases;
    max_len_ = len;

    const size_t num_words = TVector::size_to_word_count(total_bases);
    // Amortize the global-base -> (record, offset) division over a sizeable
    // sequential chunk while retaining enough chunks for large OpenMP teams.
    const size_t kWordsPerChunk = 4096;
    const size_t num_chunks = DivCeiling(num_words, kWordsPerChunk);
    num_threads = std::max(1, num_threads);

#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int64_t chunk = 0; chunk < static_cast<int64_t>(num_chunks);
         ++chunk) {
      const size_t word_begin = chunk * kWordsPerChunk;
      const size_t word_end =
          std::min(num_words, word_begin + kWordsPerChunk);
      const size_t first_base = word_begin * kBasesPerWord;
      size_t record_id = first_base / len;
      unsigned record_offset = first_base - record_id * len;

      for (size_t out_word = word_begin; out_word < word_end; ++out_word) {
        const size_t out_base = out_word * kBasesPerWord;
        unsigned remaining = static_cast<unsigned>(
            std::min<size_t>(kBasesPerWord, total_bases - out_base));
        unsigned out_offset = 0;
        TWord packed = 0;

        while (remaining > 0) {
          const unsigned source_word_offset =
              record_offset % kBasesPerWord;
          const unsigned take =
              std::min(remaining,
                       std::min(len - record_offset,
                                kBasesPerWord - source_word_offset));
          const TWord source =
              records[record_id * record_stride +
                      record_offset / kBasesPerWord];
          packed |= TVector::sub_word(source, source_word_offset, take)
                    << TVector::bit_shift(out_offset, take);

          remaining -= take;
          out_offset += take;
          record_offset += take;
          if (record_offset == len) {
            ++record_id;
            record_offset = 0;
          }
        }
        sequences_.data()[out_word] = packed;
      }
    }
  }

  /**
   * Replace an empty package with variable-length records taken from a
   * packed word stream (each record `r` occupies
   * `src[src_word_begin[r] .. +DivCeiling(lens[r],16))`, zero-padded in its
   * final word).  Optionally base-reverses every record, matching
   * AppendReversedCompactSequence.
   *
   * Unlike the fixed-length assign, destination words here may straddle two
   * records processed by different threads, so each record's first and last
   * output words are merged with relaxed atomic OR into the pre-zeroed
   * destination; all interior words are exclusively owned.  Zero-length
   * records are faked as a single 'A' exactly like the append path.
   */
  void AssignCompactRecords(const TWord *src, const uint64_t *src_word_begin,
                            const uint32_t *lens, size_t num_records,
                            bool reverse, int num_threads) {
    assert(seq_count() == 0 && base_count() == 0);
    if (num_records == 0) {
      return;
    }

    // Replicate UpdateLength bookkeeping: a maximal prefix of equal lengths
    // is stored implicitly; the remaining reads get explicit start positions.
    const auto eff_len = [lens](size_t r) -> uint32_t {
      return lens[r] == 0 ? 1u : lens[r];
    };
    size_t fixed_prefix = 1;
    while (fixed_prefix < num_records &&
           eff_len(fixed_prefix) == eff_len(0)) {
      ++fixed_prefix;
    }
    fixed_len_ = eff_len(0);
    num_fixed_len_ = fixed_prefix;
    max_len_ = fixed_len_;

    std::vector<uint64_t> dest_begin(num_records + 1);
    dest_begin[0] = 0;
    for (size_t r = 0; r < num_records; ++r) {
      dest_begin[r + 1] = dest_begin[r] + eff_len(r);
    }
    const uint64_t total_bases = dest_begin[num_records];

    start_pos_.clear();
    start_pos_.reserve(num_records - fixed_prefix + 1);
    start_pos_.push_back(static_cast<uint64_t>(fixed_len_) * fixed_prefix);
    for (size_t r = fixed_prefix; r < num_records; ++r) {
      start_pos_.push_back(dest_begin[r + 1]);
      max_len_ = std::max(max_len_, eff_len(r));
    }

    sequences_.resize(total_bases);  // zero-filled by the underlying vector
    TWord *data = sequences_.data();
    const size_t sequence_bytes = sequences_.word_count() * sizeof(TWord);
    if (num_threads > 1 && InterleaveMemoryPages(data, sequence_bytes)) {
      DiscardMemoryPages(data, sequence_bytes);
    }
    AdviseHugePages(data, sequence_bytes);
    num_threads = std::max(1, num_threads);

#pragma omp parallel for schedule(dynamic, 16384) num_threads(num_threads)
    for (int64_t r = 0; r < static_cast<int64_t>(num_records); ++r) {
      const uint32_t len = lens[r];
      if (len == 0) {
        continue;  // faked 'A' == 0 bits; destination already zero
      }
      const TWord *rec = src + src_word_begin[r];
      const uint64_t base = dest_begin[r];
      const size_t first_word = base / kBasesPerWord;
      const size_t last_word = (base + len - 1) / kBasesPerWord;
      const unsigned off = static_cast<unsigned>(base % kBasesPerWord);
      const size_t n_groups = DivCeiling(len, kBasesPerWord);
      const unsigned tail = len % kBasesPerWord;

      const auto emit = [&](size_t w, TWord bits) {
        if (bits == 0) {
          return;
        }
        if (w == first_word || w == last_word) {
          __atomic_fetch_or(&data[w], bits, __ATOMIC_RELAXED);
        } else {
          data[w] |= bits;
        }
      };

      // Aligned 16-base groups of the (possibly reversed) record.
      const unsigned shift_head = tail ? 2 * (kBasesPerWord - tail) : 0;
      const unsigned shift_tail = 2 * tail;
      TWord rev_cur =
          reverse && tail ? kmlib::bit::Reverse<2>(rec[n_groups - 1]) : 0;
      for (size_t i = 0; i < n_groups; ++i) {
        TWord g;
        if (!reverse) {
          g = rec[i];
          if (i + 1 == n_groups && tail) {
            g &= ~((TWord{1} << (2 * (kBasesPerWord - tail))) - 1);
          }
        } else if (tail == 0) {
          g = kmlib::bit::Reverse<2>(rec[n_groups - 1 - i]);
        } else {
          const TWord next =
              i + 1 < n_groups
                  ? kmlib::bit::Reverse<2>(rec[n_groups - 2 - i])
                  : TWord{0};
          g = (rev_cur << shift_head) | (next >> shift_tail);
          rev_cur = next;
        }
        const size_t w = first_word + i;
        if (off == 0) {
          emit(w, g);
        } else {
          emit(w, g >> (2 * off));
          emit(w + 1, g << (2 * (kBasesPerWord - off)));
        }
      }
    }
  }

  /**
   * Replace an empty package from the on-disk binary-library stream
   * `[uint32 length][packed words]...` without materializing the stream or a
   * per-read source-offset array in anonymous memory.
   *
   * The input is first walked once to validate record boundaries, build the
   * read-only length index, and split it into byte-bounded chunks.  Chunks are
   * then repacked in parallel directly into their disjoint destination base
   * ranges.  Only the two destination words straddling a chunk boundary need
   * atomic ORs; words inside a chunk have a single writer.  This removes the
   * two atomic updates per record used by AssignCompactRecords and bounds the
   * resident source mapping independently of library size.
   *
   * `expected_max_len` comes from the binary library's own metadata.  A
   * mismatch makes the caller fall back to the portable streaming reader.
   */
  bool AssignMappedBinaryRecords(const TWord *records, size_t total_words,
                                 size_t num_records,
                                 uint64_t expected_bases,
                                 unsigned expected_max_len, bool reverse,
                                 int num_threads) {
    static_assert(sizeof(TWord) == sizeof(uint32_t),
                  "binary libraries use uint32 packed words");
    if (records == nullptr || total_words == 0 || num_records == 0 ||
        expected_max_len == 0) {
      return false;
    }
    if (num_records >
        std::numeric_limits<uint64_t>::max() / expected_max_len) {
      return false;
    }

    struct BinaryChunk {
      size_t word_begin;
      size_t word_end;
      size_t record_begin;
      size_t record_end;
      uint64_t base_begin;
      uint64_t base_end;
    };

    const uint64_t nominal_bases =
        static_cast<uint64_t>(num_records) * expected_max_len;
    if (nominal_bases < expected_bases) {
      return false;
    }
    const uint64_t expected_gap = nominal_bases - expected_bases;

    enum LengthIndexKind { kSparseGap, kGap16, kGap32, kAbsolute };
    LengthIndexKind index_kind = kAbsolute;
    if (expected_max_len <= std::numeric_limits<uint16_t>::max() &&
        expected_gap <= std::numeric_limits<uint32_t>::max() &&
        expected_gap <= num_records / 64u) {
      index_kind = kSparseGap;
    } else if (expected_gap <= std::numeric_limits<uint16_t>::max()) {
      index_kind = kGap16;
    } else if (expected_gap <= std::numeric_limits<uint32_t>::max()) {
      index_kind = kGap32;
    }

    constexpr unsigned kMappedSparseGapBlockShift = 10;
    const size_t num_sparse_blocks =
        DivCeiling(num_records,
                   size_t{1} << kMappedSparseGapBlockShift);
    std::vector<uint32_t> sparse_block_begin;
    std::vector<uint32_t> sparse_before_block;
    std::vector<SparseReadGap> sparse_read_gaps;
    std::vector<uint16_t> gap16;
    std::vector<uint32_t> gap32;
    std::vector<uint64_t> absolute_starts;
    if (index_kind == kSparseGap) {
      sparse_block_begin.resize(num_sparse_blocks + 1);
      sparse_before_block.resize(num_sparse_blocks);
      sparse_read_gaps.reserve(static_cast<size_t>(expected_gap));
    } else if (index_kind == kGap16) {
      gap16.resize(num_records + 1);
    } else if (index_kind == kGap32) {
      gap32.resize(num_records + 1);
    } else {
      absolute_starts.resize(num_records + 1);
    }

    // Chunks are bounded by input bytes rather than a machine-specific read
    // count.  This keeps mapped-source residency bounded for both short and
    // unusually long records while leaving thousands of independent chunks
    // on large libraries.
    constexpr size_t kTargetSourceChunkBytes = size_t{16} << 20u;
    constexpr size_t kTargetSourceChunkWords =
        kTargetSourceChunkBytes / sizeof(TWord);
    std::vector<BinaryChunk> chunks;
    chunks.reserve(DivCeiling(total_words, kTargetSourceChunkWords));

    size_t word_pos = 0;
    size_t chunk_word_begin = 0;
    size_t chunk_record_begin = 0;
    uint64_t base_pos = 0;
    uint64_t chunk_base_begin = 0;
    uint64_t raw_base_sum = 0;
    uint64_t cumulative_gap = 0;
    unsigned actual_max_len = 0;
    size_t current_sparse_block = std::numeric_limits<size_t>::max();

    for (size_t record_id = 0; record_id < num_records; ++record_id) {
      if (word_pos >= total_words) {
        return false;
      }
      const uint32_t len = records[word_pos];
      // The normal reader represents an empty input record as one A.  Current
      // library writers already store that effective length in both the file
      // and metadata; reject legacy zero records here so the exact fallback
      // path can apply its historical rule.
      if (len == 0 || len > expected_max_len) {
        return false;
      }
      const size_t record_words =
          DivCeiling(static_cast<size_t>(len), kBasesPerWord);
      if (record_words > total_words - word_pos - 1) {
        return false;
      }

      actual_max_len = std::max(actual_max_len, static_cast<unsigned>(len));
      raw_base_sum += len;
      base_pos += len;
      cumulative_gap += expected_max_len - len;

      if (index_kind == kSparseGap) {
        const size_t block = record_id >> kMappedSparseGapBlockShift;
        if (block != current_sparse_block) {
          current_sparse_block = block;
          sparse_block_begin[block] = sparse_read_gaps.size();
          sparse_before_block[block] =
              static_cast<uint32_t>(cumulative_gap -
                                    (expected_max_len - len));
        }
        if (len != expected_max_len) {
          sparse_read_gaps.push_back(
              {static_cast<uint16_t>(
                   record_id &
                   ((size_t{1} << kMappedSparseGapBlockShift) - 1u)),
               static_cast<uint16_t>(expected_max_len - len)});
        }
      } else if (index_kind == kGap16) {
        gap16[record_id + 1] = static_cast<uint16_t>(cumulative_gap);
      } else if (index_kind == kGap32) {
        gap32[record_id + 1] = static_cast<uint32_t>(cumulative_gap);
      } else {
        absolute_starts[record_id + 1] = base_pos;
      }

      word_pos += 1 + record_words;
      if (word_pos - chunk_word_begin >= kTargetSourceChunkWords ||
          record_id + 1 == num_records) {
        chunks.push_back({chunk_word_begin, word_pos, chunk_record_begin,
                          record_id + 1, chunk_base_begin, base_pos});
        // The page-interior helper is also valid for a private file mapping:
        // it removes completed PTEs without changing file contents.  Boundary
        // pages are deliberately retained because the next record may share
        // one of them.
        DiscardMemoryPages(
            const_cast<TWord *>(records + chunk_word_begin),
            (word_pos - chunk_word_begin) * sizeof(TWord));
        chunk_word_begin = word_pos;
        chunk_record_begin = record_id + 1;
        chunk_base_begin = base_pos;
      }
    }

    if (word_pos != total_words || raw_base_sum != expected_bases ||
        base_pos != expected_bases || actual_max_len != expected_max_len ||
        cumulative_gap != expected_gap) {
      return false;
    }
    if (index_kind == kSparseGap) {
      sparse_block_begin[num_sparse_blocks] = sparse_read_gaps.size();
    }

    Clear();
    compact_seq_count_ = num_records;
    compact_base_count_ = expected_bases;
    max_len_ = expected_max_len;
    fixed_len_ = 0;
    num_fixed_len_ = 0;
    if (index_kind == kSparseGap) {
      sparse_gap_block_begin_.swap(sparse_block_begin);
      sparse_gap_before_block_.swap(sparse_before_block);
      sparse_read_gaps_.swap(sparse_read_gaps);
      std::vector<uint64_t>().swap(start_pos_);
    } else if (index_kind == kGap16) {
      read_gap16_.swap(gap16);
      std::vector<uint64_t>().swap(start_pos_);
    } else if (index_kind == kGap32) {
      read_gap32_.swap(gap32);
      std::vector<uint64_t>().swap(start_pos_);
    } else {
      start_pos_.swap(absolute_starts);
      compact_seq_count_ = 0;
      compact_base_count_ = 0;
    }

    sequences_.resize(expected_bases);
    TWord *destination = sequences_.data();
    const size_t destination_bytes =
        sequences_.word_count() * sizeof(TWord);
    num_threads = std::max(1, num_threads);
    if (num_threads > 1 &&
        InterleaveMemoryPages(destination, destination_bytes)) {
      DiscardMemoryPages(destination, destination_bytes);
    }
    AdviseHugePages(destination, destination_bytes);

#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
      const BinaryChunk &chunk = chunks[chunk_id];
      const size_t first_destination_word =
          chunk.base_begin / kBasesPerWord;
      const size_t last_destination_word =
          (chunk.base_end - 1) / kBasesPerWord;
      uint64_t destination_base = chunk.base_begin;
      size_t source_word = chunk.word_begin;

      const auto emit_group = [&](TWord value, unsigned source_offset,
                                  unsigned count) {
        while (count != 0) {
          const size_t out_word = destination_base / kBasesPerWord;
          const unsigned out_offset =
              static_cast<unsigned>(destination_base % kBasesPerWord);
          const unsigned take =
              std::min(count, kBasesPerWord - out_offset);
          const TWord bits =
              TVector::sub_word(value, source_offset, take)
              << TVector::bit_shift(out_offset, take);
          if (out_word == first_destination_word ||
              out_word == last_destination_word) {
            __atomic_fetch_or(destination + out_word, bits,
                              __ATOMIC_RELAXED);
          } else {
            destination[out_word] |= bits;
          }
          source_offset += take;
          destination_base += take;
          count -= take;
        }
      };

      for (size_t record_id = chunk.record_begin;
           record_id < chunk.record_end; ++record_id) {
        const uint32_t len = records[source_word++];
        const size_t num_groups =
            DivCeiling(static_cast<size_t>(len), kBasesPerWord);
        const unsigned tail = len % kBasesPerWord;
        if (!reverse) {
          const size_t full_groups = len / kBasesPerWord;
          for (size_t group = 0; group < full_groups; ++group) {
            emit_group(records[source_word + group], 0, kBasesPerWord);
          }
          if (tail != 0) {
            emit_group(records[source_word + full_groups], 0, tail);
          }
        } else {
          size_t full_groups = num_groups;
          if (tail != 0) {
            const TWord value =
                kmlib::bit::Reverse<2>(records[source_word + num_groups - 1]);
            emit_group(value, kBasesPerWord - tail, tail);
            --full_groups;
          }
          while (full_groups != 0) {
            --full_groups;
            emit_group(kmlib::bit::Reverse<2>(
                           records[source_word + full_groups]),
                       0, kBasesPerWord);
          }
        }
        source_word += num_groups;
      }
      assert(source_word == chunk.word_end);
      assert(destination_base == chunk.base_end);
      DiscardMemoryPages(const_cast<TWord *>(records + chunk.word_begin),
                         (chunk.word_end - chunk.word_begin) * sizeof(TWord));
    }
    return true;
  }

  /**
   * Build a dense read-only package from an ordered subset of an indexed
   * binary library.  Input chunks are scanned independently, and each worker
   * writes one contiguous destination base interval.  Thus source residency
   * is bounded by active chunks and the destination remains tightly packed;
   * only destination words shared by adjacent chunks require atomic ORs.
   *
   * `Chunk` is intentionally structural (word_begin/end, read_begin/end), so
   * this core container does not depend on an I/O-layer index header.
   */
  template <class Chunk>
  bool AssignSelectedMappedBinaryRecords(
      const TWord *records, size_t total_words, size_t num_source_records,
      const std::vector<uint64_t> &selected_ids,
      const std::vector<Chunk> &chunks, int num_threads) {
    static_assert(sizeof(TWord) == sizeof(uint32_t),
                  "binary libraries use uint32 packed words");
    assert(seq_count() == 0 && base_count() == 0);
    if (selected_ids.empty()) {
      return true;
    }
    if (records == nullptr || chunks.empty() ||
        selected_ids.back() >= num_source_records ||
        chunks.back().word_end != total_words ||
        chunks.back().read_end != num_source_records) {
      return false;
    }

    const size_t num_chunks = chunks.size();
    std::vector<size_t> selected_begin(num_chunks);
    std::vector<size_t> selected_end(num_chunks);
    std::vector<uint64_t> chunk_bases(num_chunks + 1u, 0);
    std::vector<uint32_t> lengths(selected_ids.size());
    std::vector<uint8_t> valid(num_chunks, 1u);
    num_threads = std::max(1, num_threads);

    // First pass obtains exact lengths and per-chunk destination sizes.  It
    // reads only record headers, then drops completed mapped pages.
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(num_chunks); ++chunk_id) {
      const Chunk &chunk = chunks[chunk_id];
      const auto first = std::lower_bound(selected_ids.begin(),
                                          selected_ids.end(),
                                          chunk.read_begin);
      const auto last = std::lower_bound(first, selected_ids.end(),
                                         chunk.read_end);
      const size_t first_id = static_cast<size_t>(first - selected_ids.begin());
      const size_t last_id = static_cast<size_t>(last - selected_ids.begin());
      selected_begin[chunk_id] = first_id;
      selected_end[chunk_id] = last_id;

      const TWord *cursor = records + chunk.word_begin;
      size_t selected_id = first_id;
      uint64_t base_sum = 0;
      for (uint64_t read_id = chunk.read_begin;
           read_id < chunk.read_end; ++read_id) {
        if (cursor >= records + chunk.word_end) {
          valid[chunk_id] = 0;
          break;
        }
        const uint32_t len = *cursor;
        const size_t words = DivCeiling(static_cast<size_t>(len),
                                        kBasesPerWord);
        if (len == 0 || words > static_cast<size_t>(records + chunk.word_end -
                                                    cursor - 1)) {
          valid[chunk_id] = 0;
          break;
        }
        if (selected_id < last_id && selected_ids[selected_id] == read_id) {
          lengths[selected_id++] = len;
          base_sum += len;
        }
        cursor += 1u + words;
      }
      if (cursor != records + chunk.word_end || selected_id != last_id) {
        valid[chunk_id] = 0;
      }
      chunk_bases[chunk_id + 1u] = base_sum;
      DiscardMemoryPages(
          const_cast<TWord *>(records + chunk.word_begin),
          (chunk.word_end - chunk.word_begin) * sizeof(TWord));
    }
    for (uint8_t chunk_valid : valid) {
      if (!chunk_valid) return false;
    }
    for (size_t chunk = 0; chunk < num_chunks; ++chunk) {
      if (chunk_bases[chunk + 1u] >
          std::numeric_limits<uint64_t>::max() - chunk_bases[chunk]) {
        return false;
      }
      chunk_bases[chunk + 1u] += chunk_bases[chunk];
    }
    const uint64_t total_bases = chunk_bases.back();

    // Build the same fixed-prefix/absolute-tail length representation used by
    // normal appends.  The caller may compact rare length gaps afterwards.
    Clear();
    fixed_len_ = lengths[0];
    size_t fixed_prefix = 1;
    while (fixed_prefix < lengths.size() &&
           lengths[fixed_prefix] == fixed_len_) {
      ++fixed_prefix;
    }
    num_fixed_len_ = fixed_prefix;
    max_len_ = fixed_len_;
    start_pos_.clear();
    uint64_t base = static_cast<uint64_t>(fixed_prefix) * fixed_len_;
    start_pos_.push_back(base);
    for (size_t read = fixed_prefix; read < lengths.size(); ++read) {
      base += lengths[read];
      start_pos_.push_back(base);
      max_len_ = std::max(max_len_, lengths[read]);
    }
    if (base != total_bases) {
      return false;
    }

    sequences_.resize(total_bases);
    TWord *destination = sequences_.data();
    const size_t destination_bytes =
        sequences_.word_count() * sizeof(TWord);
    if (num_threads > 1 &&
        InterleaveMemoryPages(destination, destination_bytes)) {
      DiscardMemoryPages(destination, destination_bytes);
    }
    AdviseHugePages(destination, destination_bytes);

    // Second pass copies only selected records while their input chunk is
    // resident, then immediately releases the source pages.
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(num_chunks); ++chunk_id) {
      const Chunk &chunk = chunks[chunk_id];
      const size_t first_selected = selected_begin[chunk_id];
      const size_t last_selected = selected_end[chunk_id];
      if (first_selected == last_selected) {
        DiscardMemoryPages(
            const_cast<TWord *>(records + chunk.word_begin),
            (chunk.word_end - chunk.word_begin) * sizeof(TWord));
        continue;
      }

      uint64_t destination_base = chunk_bases[chunk_id];
      const uint64_t destination_end = chunk_bases[chunk_id + 1u];
      const size_t first_destination_word =
          destination_base / kBasesPerWord;
      const size_t last_destination_word =
          (destination_end - 1u) / kBasesPerWord;
      const TWord *cursor = records + chunk.word_begin;
      size_t selected_id = first_selected;

      const auto emit_group = [&](TWord value, unsigned count) {
        unsigned source_offset = 0;
        while (count != 0) {
          const size_t out_word = destination_base / kBasesPerWord;
          const unsigned out_offset =
              static_cast<unsigned>(destination_base % kBasesPerWord);
          const unsigned take =
              std::min(count, kBasesPerWord - out_offset);
          const TWord bits =
              TVector::sub_word(value, source_offset, take)
              << TVector::bit_shift(out_offset, take);
          if (bits != 0) {
            if (out_word == first_destination_word ||
                out_word == last_destination_word) {
              __atomic_fetch_or(destination + out_word, bits,
                                __ATOMIC_RELAXED);
            } else {
              destination[out_word] |= bits;
            }
          }
          source_offset += take;
          destination_base += take;
          count -= take;
        }
      };

      for (uint64_t read_id = chunk.read_begin;
           read_id < chunk.read_end; ++read_id) {
        const uint32_t len = *cursor++;
        const size_t num_words =
            DivCeiling(static_cast<size_t>(len), kBasesPerWord);
        if (selected_id < last_selected &&
            selected_ids[selected_id] == read_id) {
          const size_t full_words = len / kBasesPerWord;
          for (size_t word = 0; word < full_words; ++word) {
            emit_group(cursor[word], kBasesPerWord);
          }
          const unsigned tail = len % kBasesPerWord;
          if (tail != 0) emit_group(cursor[full_words], tail);
          ++selected_id;
        }
        cursor += num_words;
      }
      assert(selected_id == last_selected);
      assert(destination_base == destination_end);
      DiscardMemoryPages(
          const_cast<TWord *>(records + chunk.word_begin),
          (chunk.word_end - chunk.word_begin) * sizeof(TWord));
    }
    return true;
  }

  /**
   * Append fixed-length records stored as big-endian 2-bit uint32_t words.
   *
   * Mercy-edge discovery naturally produces independent fixed-size k-mers.
   * Appending them one base/record at a time serialized an otherwise parallel
   * stage.  Packing by destination-word chunks keeps every output word owned
   * by exactly one worker, including the word shared with the old tail.
   */
  void AppendFixedLengthCompactSequences32(const uint32_t *records,
                                           size_t num_records, unsigned len,
                                           unsigned record_stride,
                                           int num_threads) {
    assert(len > 0);
    assert(record_stride >= DivCeiling(len, unsigned{16}));
    assert(IsFixedLength());
    assert(num_fixed_len_ == 0 || fixed_len_ == len);
    if (num_records == 0) {
      return;
    }
    if (num_records >
        (std::numeric_limits<size_t>::max() - base_count()) / len) {
      xfatal("Fixed-length sequence append size overflow\n");
    }

    const size_t old_total_bases = base_count();
    const size_t added_bases = num_records * static_cast<size_t>(len);
    const size_t new_total_bases = old_total_bases + added_bases;
    const size_t first_word = old_total_bases / kBasesPerWord;
    const size_t end_word = TVector::size_to_word_count(new_total_bases);

    sequences_.resize(new_total_bases);
    if (num_fixed_len_ == 0) {
      fixed_len_ = len;
    }
    num_fixed_len_ += num_records;
    start_pos_.back() = new_total_bases;
    max_len_ = std::max(max_len_, len);

    const size_t num_words = end_word - first_word;
    const size_t kWordsPerChunk = 4096;
    const size_t num_chunks = DivCeiling(num_words, kWordsPerChunk);
    num_threads = std::max(1, num_threads);

#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int64_t chunk = 0; chunk < static_cast<int64_t>(num_chunks);
         ++chunk) {
      const size_t word_begin =
          first_word + static_cast<size_t>(chunk) * kWordsPerChunk;
      const size_t word_end =
          std::min(end_word, word_begin + kWordsPerChunk);

      for (size_t out_word = word_begin; out_word < word_end; ++out_word) {
        const size_t out_base = out_word * kBasesPerWord;
        const size_t write_begin = std::max(out_base, old_total_bases);
        const size_t write_end =
            std::min(out_base + kBasesPerWord, new_total_bases);
        size_t relative_base = write_begin - old_total_bases;
        size_t record_id = relative_base / len;
        unsigned record_offset = relative_base - record_id * len;
        unsigned out_offset = write_begin - out_base;
        unsigned remaining = write_end - write_begin;
        TWord packed = out_word == first_word ? sequences_.data()[out_word] : 0;

        while (remaining > 0) {
          const unsigned source_word_offset = record_offset & 15u;
          const unsigned take =
              std::min(remaining,
                       std::min(len - record_offset,
                                16u - source_word_offset));
          const uint32_t source =
              records[record_id * record_stride + record_offset / 16u];
          const unsigned source_shift =
              32u - (source_word_offset + take) * 2u;
          uint64_t value = source >> source_shift;
          if (take != 16u) {
            value &= (uint64_t{1} << (take * 2u)) - 1u;
          }
          packed |= static_cast<TWord>(value)
                    << TVector::bit_shift(out_offset, take);

          remaining -= take;
          out_offset += take;
          record_offset += take;
          if (record_offset == len) {
            ++record_id;
            record_offset = 0;
          }
        }
        sequences_.data()[out_word] = packed;
      }
    }
  }

  size_t seq_count() const {
    if (!sparse_gap_block_begin_.empty() || !read_gap16_.empty() ||
        !read_gap32_.empty()) {
      return compact_seq_count_;
    }
    return num_fixed_len_ + start_pos_.size() - 1;
  }

  size_t base_count() const {
    if (!sparse_gap_block_begin_.empty() || !read_gap16_.empty() ||
        !read_gap32_.empty()) {
      return compact_base_count_;
    }
    return start_pos_.back();
  }

  size_t size_in_byte() const {
    return sizeof(TWord) * sequences_.word_capacity() +
           sizeof(uint64_t) * start_pos_.capacity() +
           sizeof(uint32_t) * sparse_gap_block_begin_.capacity() +
           sizeof(uint32_t) * sparse_gap_before_block_.capacity() +
           sizeof(SparseReadGap) * sparse_read_gaps_.capacity() +
           sizeof(uint16_t) * read_gap16_.capacity() +
           sizeof(uint32_t) * read_gap32_.capacity() +
           sizeof(uint64_t) * pos_to_id_.capacity();
  }

  unsigned max_length() const { return max_len_; }

  /**
   * Replace the 64-bit absolute start table with a compact cumulative-gap
   * table for a read-only package.  Read libraries are commonly almost fixed
   * length: after one trimmed read the old representation nevertheless stores
   * an absolute 64-bit start for every remaining read.  For a sequence i,
   *
   *   start(i) = i * max_length - cumulative_gap(i)
   *
   * and adjacent gaps also recover its exact length.  Very sparse trimming is
   * represented as a hot block-prefix table plus only the exceptional reads;
   * otherwise use the narrowest safe dense gap type.  Appending after this
   * operation is not supported; callers intentionally use it only after
   * loading count input.
   */
  unsigned CompactReadOnlyLengthIndex() {
    if (!sparse_gap_block_begin_.empty() || !read_gap16_.empty() ||
        !read_gap32_.empty() ||
        start_pos_.size() <= 1 || max_len_ == 0) {
      return 0;
    }

    const size_t n = seq_count();
    const uint64_t total_bases = base_count();
    if (n > std::numeric_limits<uint64_t>::max() / max_len_) {
      return 0;
    }
    const uint64_t nominal_bases = static_cast<uint64_t>(n) * max_len_;
    if (nominal_bases < total_bases) {
      return 0;
    }
    const uint64_t total_gap = nominal_bases - total_bases;
    if (total_gap > std::numeric_limits<uint32_t>::max()) {
      return 0;
    }

    auto actual_start = [&](size_t boundary) -> uint64_t {
      if (boundary < num_fixed_len_) {
        return static_cast<uint64_t>(boundary) * fixed_len_;
      }
      return start_pos_[boundary - num_fixed_len_];
    };

    // With at most one exceptional read per 64 reads, a two-level sparse
    // index is both smaller and faster than touching a random dense array.
    // Each 1024-read block stores its cumulative gap and exception range; the
    // common no-exception lookup touches only the small block tables.
    if (max_len_ <= std::numeric_limits<uint16_t>::max() &&
        total_gap <= n / 64u) {
      const size_t num_blocks =
          DivCeiling(n, size_t{1} << kSparseGapBlockShift);
      sparse_gap_block_begin_.resize(num_blocks + 1);
      sparse_gap_before_block_.resize(num_blocks);
      sparse_read_gaps_.reserve(static_cast<size_t>(total_gap));
      uint64_t cumulative_gap = 0;
      for (size_t block = 0; block < num_blocks; ++block) {
        const size_t seq_begin = block << kSparseGapBlockShift;
        const size_t seq_end =
            std::min(n, seq_begin + (size_t{1} << kSparseGapBlockShift));
        sparse_gap_block_begin_[block] = sparse_read_gaps_.size();
        sparse_gap_before_block_[block] =
            static_cast<uint32_t>(cumulative_gap);
        uint64_t start = actual_start(seq_begin);
        for (size_t seq_id = seq_begin; seq_id < seq_end; ++seq_id) {
          const uint64_t end = actual_start(seq_id + 1);
          const unsigned len = static_cast<unsigned>(end - start);
          assert(len <= max_len_);
          if (len != max_len_) {
            const unsigned gap = max_len_ - len;
            sparse_read_gaps_.push_back(
                {static_cast<uint16_t>(seq_id - seq_begin),
                 static_cast<uint16_t>(gap)});
            cumulative_gap += gap;
          }
          start = end;
        }
      }
      sparse_gap_block_begin_[num_blocks] = sparse_read_gaps_.size();
      assert(cumulative_gap == total_gap);
    } else if (total_gap <= std::numeric_limits<uint16_t>::max()) {
      read_gap16_.resize(n + 1);
#pragma omp parallel for schedule(static)
      for (int64_t i = 0; i <= static_cast<int64_t>(n); ++i) {
        read_gap16_[i] = static_cast<uint16_t>(
            static_cast<uint64_t>(i) * max_len_ - actual_start(i));
      }
    } else {
      read_gap32_.resize(n + 1);
#pragma omp parallel for schedule(static)
      for (int64_t i = 0; i <= static_cast<int64_t>(n); ++i) {
        read_gap32_[i] = static_cast<uint32_t>(
            static_cast<uint64_t>(i) * max_len_ - actual_start(i));
      }
    }

    compact_seq_count_ = n;
    compact_base_count_ = total_bases;
    std::vector<uint64_t>().swap(start_pos_);
    if (!sparse_gap_block_begin_.empty()) {
      return 1u;
    }
    return read_gap16_.empty() ? 32u : 16u;
  }

  size_t compact_length_exceptions() const {
    return sparse_read_gaps_.size();
  }

  SeqView GetSeqView(size_t seq_id) const { return SeqView(this, seq_id); }

  SeqView GetSeqViewByOffset(size_t offset) const {
    return SeqView(this, GetSeqID(offset));
  }

 private:
  void GetSeqBounds(size_t seq_id, uint64_t *start, unsigned *len) const {
    if (!sparse_gap_block_begin_.empty()) {
      const size_t block = seq_id >> kSparseGapBlockShift;
      const unsigned in_block =
          seq_id & ((size_t{1} << kSparseGapBlockShift) - 1u);
      uint32_t gap = sparse_gap_before_block_[block];
      *len = max_len_;
      const uint32_t exception_end = sparse_gap_block_begin_[block + 1];
      for (uint32_t i = sparse_gap_block_begin_[block]; i < exception_end;
           ++i) {
        const SparseReadGap exception = sparse_read_gaps_[i];
        if (exception.offset < in_block) {
          gap += exception.gap;
        } else {
          if (exception.offset == in_block) {
            *len -= exception.gap;
          }
          break;
        }
      }
      *start = static_cast<uint64_t>(seq_id) * max_len_ - gap;
      return;
    }
    if (!read_gap16_.empty()) {
      const uint16_t gap0 = read_gap16_[seq_id];
      const uint16_t gap1 = read_gap16_[seq_id + 1];
      *start = static_cast<uint64_t>(seq_id) * max_len_ - gap0;
      *len = max_len_ - static_cast<unsigned>(gap1 - gap0);
      return;
    }
    if (!read_gap32_.empty()) {
      const uint32_t gap0 = read_gap32_[seq_id];
      const uint32_t gap1 = read_gap32_[seq_id + 1];
      *start = static_cast<uint64_t>(seq_id) * max_len_ - gap0;
      *len = max_len_ - static_cast<unsigned>(gap1 - gap0);
      return;
    }
    if (seq_id < num_fixed_len_) {
      *start = static_cast<uint64_t>(seq_id) * fixed_len_;
      *len = fixed_len_;
    } else {
      const size_t variable_id = seq_id - num_fixed_len_;
      *start = start_pos_[variable_id];
      *len = start_pos_[variable_id + 1] - start_pos_[variable_id];
    }
  }

  unsigned GetSeqLength(size_t seq_id) const {
    if (!sparse_gap_block_begin_.empty()) {
      uint64_t ignored_start;
      unsigned len;
      GetSeqBounds(seq_id, &ignored_start, &len);
      return len;
    }
    if (!read_gap16_.empty()) {
      return max_len_ -
             static_cast<unsigned>(read_gap16_[seq_id + 1] -
                                   read_gap16_[seq_id]);
    }
    if (!read_gap32_.empty()) {
      return max_len_ -
             static_cast<unsigned>(read_gap32_[seq_id + 1] -
                                   read_gap32_[seq_id]);
    }
    if (seq_id < num_fixed_len_) {
      return fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_ + 1] -
             start_pos_[seq_id - num_fixed_len_];
    }
  }

  uint8_t GetBase(size_t seq_id, unsigned offset) const {
    return sequences_[StartPos(seq_id) + offset];
  }

  uint64_t StartPos(size_t seq_id) const {
    if (!sparse_gap_block_begin_.empty()) {
      if (seq_id == compact_seq_count_) {
        return compact_base_count_;
      }
      uint64_t start;
      unsigned ignored_len;
      GetSeqBounds(seq_id, &start, &ignored_len);
      return start;
    }
    if (!read_gap16_.empty()) {
      return static_cast<uint64_t>(seq_id) * max_len_ - read_gap16_[seq_id];
    }
    if (!read_gap32_.empty()) {
      return static_cast<uint64_t>(seq_id) * max_len_ - read_gap32_[seq_id];
    }
    if (seq_id < num_fixed_len_) {
      return seq_id * fixed_len_;
    } else {
      return start_pos_[seq_id - num_fixed_len_];
    }
  }

  TAddress GetRawAddress(size_t seq_id, unsigned offset = 0) const {
    size_t index = StartPos(seq_id) + offset;
    return {sequences_.data() + index / kBasesPerWord, index % kBasesPerWord};
  }

  uint64_t GetSeqID(size_t full_offset) const {
    assert(sparse_gap_block_begin_.empty() && read_gap16_.empty() &&
           read_gap32_.empty());
    if (full_offset < num_fixed_len_ * fixed_len_) {
      return full_offset / fixed_len_;
    } else {
      size_t look_up_entry =
          (full_offset - num_fixed_len_ * fixed_len_) / kLookupStep;
      size_t l = pos_to_id_[look_up_entry], r = pos_to_id_[look_up_entry + 1];

      while (l < r) {
        size_t mid = (l + r) / 2;
        if (start_pos_[mid - num_fixed_len_] > full_offset) {
          r = mid - 1;
        } else if (start_pos_[mid - num_fixed_len_ + 1] <= full_offset) {
          l = mid + 1;
        } else {
          return mid;
        }
      }
      return l;
    }
  }

 public:
  void AppendStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s, s + len, len);
  }

  void AppendReversedStringSequence(const char *s, unsigned len) {
    AppendStringSequence(s + len - 1, s - 1, len);
  }

  void AppendCompactSequence(const TWord *s, unsigned len) {
    AppendCompactSequence(s, len, false);
  }

  /** Append a possibly word-unaligned view without decoding its bases. */
  void AppendSequenceView(const SeqView &view) {
    unsigned len = view.length();
    if (len == 0) {
      TWord fake_sequence = 0;
      return AppendCompactSequence(&fake_sequence, 1, false);
    }
    UpdateLength(len);
    auto address = view.raw_address();
    const TWord *ptr = address.first;
    unsigned offset = address.second;
    if (offset != 0) {
      const unsigned take = std::min(len, kBasesPerWord - offset);
      sequences_.push_word(*ptr, offset, take);
      len -= take;
      ++ptr;
    }
    while (len >= kBasesPerWord) {
      sequences_.push_word(*ptr++);
      len -= kBasesPerWord;
    }
    if (len != 0) {
      sequences_.push_word(*ptr, 0, len);
    }
  }

  void AppendReversedCompactSequence(const TWord *s, unsigned len) {
    AppendCompactSequence(s, len, true);
  }

  void FetchSequence(size_t seq_id, std::vector<TWord> *s) const {
    TVector cvec(s);
    auto ptr_and_offset = GetRawAddress(seq_id);
    auto ptr = ptr_and_offset.first;
    auto offset = ptr_and_offset.second;
    auto len = GetSeqLength(seq_id);
    if (offset != 0) {
      unsigned remaining_len = std::min(len, kBasesPerWord - offset);
      cvec.push_word(*ptr, offset, remaining_len);
      len -= remaining_len;
      ++ptr;
    }
    unsigned n_full = len / kBasesPerWord;
    for (unsigned i = 0; i < n_full; ++i) {
      cvec.push_word(ptr[i]);
    }
    if (len % kBasesPerWord > 0) {
      cvec.push_word(ptr[n_full], 0, len % kBasesPerWord);
    }
  }

  void BuildIndex() {
    assert(sparse_gap_block_begin_.empty() && read_gap16_.empty() &&
           read_gap32_.empty());
    pos_to_id_.clear();
    pos_to_id_.reserve(start_pos_.back() / kLookupStep + 4);
    size_t abs_offset = num_fixed_len_ * fixed_len_;
    size_t cur_id = num_fixed_len_;

    while (abs_offset <= start_pos_.back()) {
      while (cur_id < seq_count() &&
             start_pos_[cur_id - num_fixed_len_ + 1] <= abs_offset) {
        ++cur_id;
      }

      pos_to_id_.push_back(cur_id);
      abs_offset += kLookupStep;
    }

    pos_to_id_.push_back(seq_count());
    pos_to_id_.push_back(seq_count());
  }

  void WriteSequences(std::ostream &os, int64_t from = 0,
                      int64_t to = -1) const {
    if (to == -1) {
      to = seq_count() - 1;
    }

    uint32_t len;
    std::vector<uint32_t> s;

    for (int64_t i = from; i <= to; ++i) {
      len = GetSeqView(i).length();
      FetchSequence(i, &s);
      os.write(reinterpret_cast<const char *>(&len), sizeof(uint32_t));
      os.write(reinterpret_cast<const char *>(s.data()),
               sizeof(uint32_t) * s.size());
    }
  }

 private:
  bool IsFixedLength() const { return start_pos_.size() == 1; }

  void UpdateLength(unsigned len) {
    assert(sparse_gap_block_begin_.empty() && read_gap16_.empty() &&
           read_gap32_.empty());
    if (num_fixed_len_ == 0) {
      num_fixed_len_ = 1;
      fixed_len_ = len;
      start_pos_.back() += len;
    } else if (IsFixedLength() && len == fixed_len_) {
      num_fixed_len_++;
      start_pos_.back() += len;
    } else {
      start_pos_.push_back(start_pos_.back() + len);
    }
    if (len > max_len_) {
      max_len_ = len;
    }
  }

  void AppendStringSequence(const char *from, const char *to, unsigned len) {
    if (len == 0) {
      // Fake a sequence whose length is 1, as we need all sequences' length > 0
      // to make `GetSeqID` working
      auto fake_sequence = "A";
      return AppendStringSequence(fake_sequence, fake_sequence + 1, 1);
    }
    UpdateLength(len);
    std::ptrdiff_t step = from < to ? 1 : -1;
    for (auto ptr = from; ptr != to; ptr += step) {
      sequences_.push_back(dna_map_[static_cast<int>(*ptr)]);
    }
  }

  void AppendCompactSequence(const TWord *ptr, unsigned len, bool rev) {
    if (len == 0) {
      // Fake a sequence whose length is 1, as we need all sequences' length > 0
      // to make `GetSeqID` working
      TWord fake_sequence = 0;
      return AppendCompactSequence(&fake_sequence, 1, false);
    }
    UpdateLength(len);

    if (rev) {
      auto rptr = ptr + DivCeiling(len, kBasesPerWord) - 1;
      unsigned bases_in_last_word = len % kBasesPerWord;
      if (bases_in_last_word > 0) {
        auto val = kmlib::bit::Reverse<2>(*rptr);
        sequences_.push_word(val, kBasesPerWord - bases_in_last_word,
                             bases_in_last_word);
        --rptr;
      }
      for (auto p = rptr; p >= ptr; --p) {
        sequences_.push_word(kmlib::bit::Reverse<2>(*p));
      }
    } else {
      while (len >= kBasesPerWord) {
        sequences_.push_word(*ptr);
        len -= kBasesPerWord;
        ++ptr;
      }
      if (len > 0) {
        sequences_.push_word(*ptr, 0, len);
      }
    }
  }

 private:
  TVector sequences_;
  unsigned fixed_len_{0};
  size_t num_fixed_len_{0};
  std::vector<uint64_t>
      start_pos_;  // the index of the starting position of a sequence
  struct SparseReadGap {
    uint16_t offset;
    uint16_t gap;
  };
  static constexpr unsigned kSparseGapBlockShift = 10;
  std::vector<uint32_t> sparse_gap_block_begin_;
  std::vector<uint32_t> sparse_gap_before_block_;
  std::vector<SparseReadGap> sparse_read_gaps_;
  std::vector<uint16_t> read_gap16_;
  std::vector<uint32_t> read_gap32_;
  size_t compact_seq_count_{0};
  uint64_t compact_base_count_{0};
  char dna_map_[256]{};
  unsigned max_len_{0};

  // for looking up the seq_id of a full offset
  std::vector<uint64_t> pos_to_id_;
  const static unsigned kLookupStep = 1024;
};

using SeqPackage = SequencePackage<uint32_t>;

#endif
