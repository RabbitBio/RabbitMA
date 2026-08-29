#ifndef MEGAHIT_BASE_ENGINE_H
#define MEGAHIT_BASE_ENGINE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#if defined(__SSE2__)
#include <emmintrin.h>
#endif

#include "kmsort_selector.h"
#include "utils/utils.h"

/**
 * The base class of sequence sorting engine
 */
class BaseSequenceSortingEngine {
 public:
  static const unsigned kBucketBase = 4;
  static const unsigned kBucketPrefixLength = 8;
  static const unsigned kNumBuckets = 65536;
  static const int64_t kDefaultLv1ScanTime = 8;
  static const int64_t kMaxLv1ScanTime = 128;

  static const unsigned kDefaultLv1BytesPerItem =
      4;  // traditional 32-bit differential offset

  struct MemoryStat {
    int64_t num_sequences;
    int64_t memory_for_data;
    int64_t words_per_lv2;
    int64_t aux_words_per_lv2;
  };

 public:
  BaseSequenceSortingEngine(int64_t mem, int mem_flag, int n_threads)
      : host_mem_(mem), mem_flag_(mem_flag), n_threads_(n_threads){};
  virtual ~BaseSequenceSortingEngine() { JoinPrefaultThreads(); }

  /**
   * Data that will change at every Lv1 stage
   */
 private:
  unsigned lv1_start_bucket_{}, lv1_end_bucket_{};
  int64_t lv1_num_items_{};
  std::array<uint64_t, kNumBuckets / 64> cur_lv1_buckets_{};
  struct SpecialOffset {
    int64_t item_index;
    int64_t full_offset;
  };
  std::vector<SpecialOffset> lv1_special_offsets_;
  // Lv1 locators and Lv2 sortable records have different element widths.
  // Separate aligned regions let packed-locator producers use three bytes per
  // occurrence without misaligning the uint32_t radix scratch.
  struct FreeDeleter {
    void operator()(uint32_t *p) const { std::free(p); }
  };
  struct ByteFreeDeleter {
    void operator()(uint8_t *p) const { std::free(p); }
  };
  std::unique_ptr<uint8_t[], ByteFreeDeleter> lv1_offsets_;
  int64_t lv1_offsets_capacity_items_{};
  unsigned lv1_offset_bytes_{kDefaultLv1BytesPerItem};
  std::unique_ptr<uint32_t[], FreeDeleter> lv2_items_;
  int64_t lv2_capacity_items_{};
  // Some producers already have the complete sortable record in registers
  // while scanning the input.  When the whole pass fits in memory, materialize
  // those records directly at their final bucket positions and eliminate the
  // compact-locator stream plus the random Lv2 source reread.
  std::unique_ptr<uint32_t[], FreeDeleter> lv1_direct_items_;
  int64_t lv1_direct_capacity_items_{};
  int64_t lv1_direct_words_per_item_{};
  int64_t lv1_direct_aux_words_per_item_{};
  bool lv1_direct_mode_{};
  // Speculative direct-item buffer, allocated from input metadata before the
  // input is even read so that the kernel's multi-GiB page zeroing overlaps
  // the derived class's Initialize() and the Lv0 bucket scan instead of
  // serializing into the Lv1 fill.  AdjustMemory() adopts it when the real
  // bucket totals confirm the estimate, otherwise frees it.
  std::unique_ptr<uint32_t[], FreeDeleter> prefault_direct_items_;
  int64_t prefault_capacity_items_{};
  int64_t prefault_words_per_item_{};
  std::vector<std::thread> prefault_threads_;
  std::mutex special_item_lock_;

 private:
  //  Sequence divided into n_threads_ partitions
  //  Each thread read one partition when filling Lv1 offsets, but it will fill
  //  multiple buckets
  //
  //  SeqView
  //  <--------------p0--------------><--------------p1-------------->
  //  |t0---------------------------->|t1---------------------------->
  //
  //  BucketView:
  //  <------b0------><------b1------><------b2------><------b3------>
  //  |t0---->|t1---->|t0---->|t1---->|t0---->|t1---->|t0---->|t1---->
  //  |                                       |
  //  \ bucket_begin[0] in thread_meta_[0]    \ bucket_begin[2] in
  //  thread_meta_[1]

  struct ThreadMeta {
    int64_t seq_from, seq_to;  // start and end IDs of this sequence partition
                               // (end is exclusive)
    std::array<int64_t, kNumBuckets> bucket_sizes;
    std::array<int64_t, kNumBuckets> bucket_begin;
    int64_t offset_base;  // the initial offset globals.lv1_items
  };

  std::vector<ThreadMeta> thread_meta_;
  std::array<int64_t, kNumBuckets> bucket_sizes_{};
  std::array<unsigned, kNumBuckets> bucket_real_id{};
  std::array<unsigned, kNumBuckets> bucket_rank_{};

 protected:
  using SubstrPtr = uint32_t *;
  /**
   * For derived class to fill offsets
   */
  class OffsetFiller {
   public:
    OffsetFiller(BaseSequenceSortingEngine *engine, unsigned start_bucket,
                 unsigned end_bucket, const ThreadMeta &sp)
        : engine_(engine),
          start_bucket_(start_bucket),
          thread_meta_(&sp) {
      const size_t num_buckets = end_bucket - start_bucket;
      compact_cursors_ =
          !engine_->lv1_direct_mode_ &&
          engine_->lv1_num_items_ <=
              static_cast<int64_t>(std::numeric_limits<uint32_t>::max()) &&
          engine_->Lv1CanUseCompactCursor(sp.seq_from, sp.seq_to);
      if (compact_cursors_) {
        compact_cursor_.resize(num_buckets);
        for (size_t i = 0; i < num_buckets; ++i) {
          const int64_t index = sp.bucket_begin[start_bucket + i];
          assert(index >= 0 &&
                 index <=
                     static_cast<int64_t>(std::numeric_limits<uint32_t>::max()));
          compact_cursor_[i] = static_cast<uint64_t>(index) << 32u;
        }
        write_combine_ = engine_->Lv1UseWriteCombine();
        if (write_combine_) {
          void *raw = nullptr;
          const size_t bytes_per_bucket =
              engine_->lv1_offset_bytes_ == 3
                  ? kPackedWcEntries * size_t{3}
                  : kWcLineEntries * sizeof(uint32_t);
          if (posix_memalign(&raw, 64, num_buckets * bytes_per_bucket) == 0 &&
              raw != nullptr) {
            if (engine_->lv1_offset_bytes_ == 3) {
              packed_wc_buffer_.reset(static_cast<uint8_t *>(raw));
            } else {
              wc_buffer_.reset(static_cast<uint32_t *>(raw));
            }
            wc_count_.resize(num_buckets);
            for (size_t i = 0; i < num_buckets; ++i) {
              const uint64_t index = compact_cursor_[i] >> 32u;
              wc_count_[i] =
                  engine_->lv1_offset_bytes_ == 3
                      ? uint8_t{0}
                      : ((index & (kWcLineEntries - 1u)) == 0
                             ? uint8_t{0}
                             : kWcHeadMode);
            }
          } else {
            write_combine_ = false;
          }
        }
      } else {
        if (!engine_->lv1_direct_mode_) {
          prev_full_offsets_.assign(num_buckets, sp.offset_base);
        }
        bucket_index_.assign(sp.bucket_begin.begin() + start_bucket,
                             sp.bucket_begin.begin() + end_bucket);
        if (engine_->lv1_direct_mode_) {
          const bool write_combine_enabled = engine_->Lv1UseWriteCombine();
          const int64_t item_words = engine_->lv1_direct_words_per_item_;
          const int64_t item_bytes = item_words * sizeof(uint32_t);
          int64_t gcd = 64, other = item_bytes;
          while (other != 0) {
            const int64_t r = gcd % other;
            gcd = other;
            other = r;
          }
          const int64_t block_bytes = item_bytes / gcd * 64;
          const int64_t buffer_bytes =
              static_cast<int64_t>(num_buckets) * block_bytes;
          // The staging area must stay a small, cache-friendly working set;
          // an oversized block geometry would defeat its own purpose.
          if (write_combine_enabled && block_bytes <= 1024 &&
              buffer_bytes <= (int64_t{64} << 20u)) {
            void *raw = nullptr;
            if (posix_memalign(&raw, 64, buffer_bytes) == 0 &&
                raw != nullptr) {
              wc_buffer_.reset(static_cast<uint32_t *>(raw));
              wc_block_items_ = block_bytes / item_bytes;
              wc_block_words_ = block_bytes / sizeof(uint32_t);
              wc_count_.resize(num_buckets);
              for (size_t i = 0; i < num_buckets; ++i) {
                const int64_t index = bucket_index_[i];
                wc_count_[i] = (index % wc_block_items_) == 0
                                   ? uint8_t{0}
                                   : kWcHeadMode;
              }
              write_combine_ = true;
            }
          }
        }
      }
    }

    ~OffsetFiller() {
      if (write_combine_ && compact_cursors_) {
        // Buffered tails cover destinations [index - count, index).  They are
        // written with ordinary cached stores because those lines may be
        // shared with a neighbouring thread's region.
        for (size_t i = 0; i < wc_count_.size(); ++i) {
          const uint8_t count = wc_count_[i];
          if (count == kWcHeadMode || count == 0) {
            continue;
          }
          const uint64_t next_index = compact_cursor_[i] >> 32u;
          if (engine_->lv1_offset_bytes_ == 3) {
            std::memcpy(engine_->lv1_offsets_.get() +
                            (next_index - count) * size_t{3},
                        packed_wc_buffer_.get() +
                            i * kPackedWcEntries * size_t{3},
                        count * size_t{3});
          } else {
            std::memcpy(engine_->lv1_offsets_.get() +
                            (next_index - count) * sizeof(uint32_t),
                        wc_buffer_.get() + i * kWcLineEntries,
                        count * sizeof(uint32_t));
          }
        }
      } else if (write_combine_) {
        const int64_t item_words = engine_->lv1_direct_words_per_item_;
        for (size_t i = 0; i < wc_count_.size(); ++i) {
          const uint8_t count = wc_count_[i];
          if (count == kWcHeadMode || count == 0) {
            continue;
          }
          const int64_t next_index = bucket_index_[i];
          const uint32_t *block =
              wc_buffer_.get() + i * static_cast<size_t>(wc_block_words_);
          if (count == kWcFullPending) {
            FlushBlock(block,
                       engine_->lv1_direct_items_.get() +
                           (next_index - wc_block_items_) * item_words,
                       wc_block_words_);
          } else {
            std::memcpy(engine_->lv1_direct_items_.get() +
                            (next_index - count) * item_words,
                        block, count * item_words * sizeof(uint32_t));
          }
        }
      }
      if (write_combine_) {
#if defined(__SSE2__)
        // Make the non-temporal stores globally visible before the parallel
        // region ends and other threads read the arrays.
        _mm_sfence();
#endif
      }
      if (!special_offsets_.empty()) {
        engine_->AddSpecialOffsets(std::move(special_offsets_));
      }
#ifndef NDEBUG
      for (unsigned bucket = start_bucket_;
           bucket < start_bucket_ +
                        (compact_cursors_ ? compact_cursor_.size()
                                          : bucket_index_.size());
           ++bucket) {
        const size_t local_bucket = bucket - start_bucket_;
        const int64_t final_index =
            compact_cursors_
                ? static_cast<int64_t>(compact_cursor_[local_bucket] >> 32u)
                : bucket_index_[local_bucket];
        assert(final_index == thread_meta_->bucket_begin[bucket] +
                                  thread_meta_->bucket_sizes[bucket]);
      }
#endif
    }

    void WriteNextOffset(unsigned bucket, int64_t full_offset) {
      assert(IsHandling(bucket));
      bucket = engine_->bucket_rank_[bucket] - start_bucket_;
      if (compact_cursors_) {
        assert(bucket < compact_cursor_.size());
        const uint64_t cursor = compact_cursor_[bucket];
        const uint32_t index = static_cast<uint32_t>(cursor >> 32u);
        const uint32_t previous = static_cast<uint32_t>(cursor);
        const uint64_t relative = static_cast<uint64_t>(
            full_offset - thread_meta_->offset_base);
        assert(relative <= std::numeric_limits<uint32_t>::max());
        assert(relative >= previous);
        if (write_combine_) {
          // Scattered single stores into the many-GiB offset array stall on a
          // read-for-ownership DRAM round trip almost every time: each of the
          // up-to-65536 open per-thread streams touches a different cache
          // line.  Stage entries in one 64-byte line per bucket instead and
          // flush complete destination lines with non-temporal stores, which
          // write without reading the destination.  Unaligned head/tail
          // entries take ordinary cached stores, so every byte of the offset
          // array is identical to the unbuffered path.
          const uint32_t value = EncodeOffsetValue(
              index, static_cast<int64_t>(relative - previous), full_offset);
          uint8_t count = wc_count_[bucket];
          if (engine_->lv1_offset_bytes_ == 3) {
            uint8_t *block = packed_wc_buffer_.get() +
                             bucket * kPackedWcEntries * size_t{3};
            StorePacked24(block + count * size_t{3}, value);
            if (++count == kPackedWcEntries) {
              std::memcpy(engine_->lv1_offsets_.get() +
                              (static_cast<uint64_t>(index) + 1u -
                               kPackedWcEntries) *
                                  size_t{3},
                          block, kPackedWcEntries * size_t{3});
              count = 0;
            }
            wc_count_[bucket] = count;
          } else if (count == kWcHeadMode) {
            engine_->StoreOffsetValue(index, value);
            if (((index + 1u) & (kWcLineEntries - 1u)) == 0) {
              wc_count_[bucket] = 0;
            }
          } else {
            uint32_t *line = wc_buffer_.get() + bucket * kWcLineEntries;
            line[count] = value;
            if (++count == kWcLineEntries) {
              FlushLine(line,
                        reinterpret_cast<uint32_t *>(
                            engine_->lv1_offsets_.get()) +
                            (index & ~uint32_t{kWcLineEntries - 1u}));
              count = 0;
            }
            wc_count_[bucket] = count;
          }
        } else {
          WriteOffset(index, relative - previous, full_offset);
        }
        compact_cursor_[bucket] =
            (static_cast<uint64_t>(index + 1u) << 32u) |
            static_cast<uint32_t>(relative);
      } else {
        assert(bucket < bucket_index_.size());
        const int64_t diff = full_offset - prev_full_offsets_[bucket];
        const int64_t index = bucket_index_[bucket]++;
        WriteOffset(index, diff, full_offset);
        prev_full_offsets_[bucket] = full_offset;
      }
    }

    uint32_t *ReserveNextItem(unsigned bucket) {
      assert(engine_->lv1_direct_mode_);
      const unsigned real_bucket = bucket;
      const unsigned ranked_bucket = engine_->bucket_rank_[real_bucket];
      if (UNLIKELY(ranked_bucket < start_bucket_ ||
                   ranked_bucket >= start_bucket_ + bucket_index_.size())) {
        xfatal("Direct Lv1 received bucket {} at rank {} outside pass [{}, "
               "{})\n",
               real_bucket, ranked_bucket, start_bucket_,
               start_bucket_ + bucket_index_.size());
      }
      bucket = ranked_bucket - start_bucket_;
      const int64_t index = bucket_index_[bucket]++;
      const int64_t bucket_limit =
          thread_meta_->bucket_begin[start_bucket_ + bucket] +
          thread_meta_->bucket_sizes[start_bucket_ + bucket];
      if (UNLIKELY(index < 0 ||
                   index >= engine_->lv1_direct_capacity_items_ ||
                   index >= bucket_limit)) {
        xfatal("Direct Lv1 cursor overflow: pass [{}, {}), rank {}, index "
               "{}, bucket limit {}, arena {}\n",
               start_bucket_, start_bucket_ + bucket_index_.size(),
               start_bucket_ + bucket, index, bucket_limit,
               engine_->lv1_direct_capacity_items_);
      }
      if (!write_combine_) {
        return engine_->lv1_direct_items_.get() +
               index * engine_->lv1_direct_words_per_item_;
      }
      // Software write-combining for materialized records.  Records are
      // staged in a per-bucket block whose size is the smallest common
      // multiple of the record size and the cache line, so complete blocks
      // are flushed with aligned non-temporal stores that never read the
      // multi-GiB destination.  Head/tail records use ordinary stores.
      uint8_t count = wc_count_[bucket];
      if (count == kWcHeadMode) {
        if ((static_cast<uint64_t>(index) + 1) % wc_block_items_ == 0) {
          wc_count_[bucket] = 0;
        }
        return engine_->lv1_direct_items_.get() +
               index * engine_->lv1_direct_words_per_item_;
      }
      uint32_t *block =
          wc_buffer_.get() + static_cast<size_t>(bucket) * wc_block_words_;
      if (count == kWcFullPending) {
        FlushBlock(block,
                   engine_->lv1_direct_items_.get() +
                       (index - wc_block_items_) *
                           engine_->lv1_direct_words_per_item_,
                   wc_block_words_);
        count = 0;
      }
      uint32_t *slot =
          block + static_cast<size_t>(count) *
                      engine_->lv1_direct_words_per_item_;
      wc_count_[bucket] =
          (count + 1u == wc_block_items_) ? kWcFullPending : count + 1u;
      return slot;
    }

    bool IsHandling(unsigned bucket) const {
      // Rank bounds are the authoritative description of a reusable-arena
      // pass.  The historical real-bucket bitmap is sufficient for a single
      // full pass, but it is mutable shared state and proved fragile once an
      // arena is reused over arbitrary size-ranked bucket subsets.  This
      // check is equally cheap in the partial path and is tied directly to
      // the cursor table that ReserveNextItem indexes.
      if (engine_->lv1_direct_mode_) {
        if (bucket_index_.size() == kNumBuckets) return true;
        const unsigned rank = engine_->bucket_rank_[bucket];
        return rank >= start_bucket_ &&
               rank < start_bucket_ + bucket_index_.size();
      }
      return (engine_->cur_lv1_buckets_[bucket >> 6u] >>
              (bucket & 63u)) & 1u;
    }

   private:
    static constexpr unsigned kWcLineEntries = 16;  // one 64-byte cache line
    static constexpr unsigned kPackedWcEntries = 16;
    static constexpr uint8_t kWcHeadMode = 0xFF;
    static constexpr uint8_t kWcFullPending = 0xFE;

    uint32_t EncodeOffsetValue(int64_t index, int64_t diff,
                               int64_t full_offset) {
      const uint32_t escape = engine_->OffsetEscapeValue();
      if (UNLIKELY(diff < 0)) {
        xfatal("Negative differential locator: {}\n", diff);
      }
      if (static_cast<uint64_t>(diff) < escape) {
        return static_cast<uint32_t>(diff);
      }
      special_offsets_.push_back({index, full_offset});
      return escape;
    }

    void WriteOffset(int64_t index, int64_t diff, int64_t full_offset) {
      engine_->StoreOffsetValue(
          index, EncodeOffsetValue(index, diff, full_offset));
    }

    static void StorePacked24(uint8_t *dest, uint32_t value) {
      dest[0] = static_cast<uint8_t>(value);
      dest[1] = static_cast<uint8_t>(value >> 8u);
      dest[2] = static_cast<uint8_t>(value >> 16u);
    }

    static void FlushBlock(const uint32_t *block, uint32_t *dest,
                           int64_t words) {
      assert(words % 16 == 0);
      assert(reinterpret_cast<uintptr_t>(block) % 64 == 0);
      assert(reinterpret_cast<uintptr_t>(dest) % 64 == 0);
#if defined(__SSE2__)
      const __m128i *src = reinterpret_cast<const __m128i *>(block);
      __m128i *dst = reinterpret_cast<__m128i *>(dest);
      for (int64_t i = 0; i < words / 4; i += 4) {
        _mm_stream_si128(dst + i, _mm_load_si128(src + i));
        _mm_stream_si128(dst + i + 1, _mm_load_si128(src + i + 1));
        _mm_stream_si128(dst + i + 2, _mm_load_si128(src + i + 2));
        _mm_stream_si128(dst + i + 3, _mm_load_si128(src + i + 3));
      }
#else
      std::memcpy(dest, block, words * sizeof(uint32_t));
#endif
    }

    static void FlushLine(const uint32_t *line, uint32_t *dest) {
      FlushBlock(line, dest, kWcLineEntries);
    }

    BaseSequenceSortingEngine *engine_;
    unsigned start_bucket_;
    const ThreadMeta *thread_meta_;
    // Only the current rank interval is touched in an Lv1 pass. Keeping the
    // two tables compact avoids putting 1 MiB on every worker's stack and
    // improves locality when the default multi-pass memory policy is used.
    std::vector<int64_t> prev_full_offsets_;
    std::vector<int64_t> bucket_index_;
    std::vector<uint64_t> compact_cursor_;
    // Software write-combining state: one staging block per bucket plus the
    // number of staged entries (kWcHeadMode while aligning the head,
    // kWcFullPending when a complete block awaits flushing).  The compact
    // offset path stages 16 uint32 values per 64-byte line; the direct item
    // path stages whole records in lcm(record, cache line) blocks.
    std::unique_ptr<uint32_t[], FreeDeleter> wc_buffer_;
    std::unique_ptr<uint8_t[], ByteFreeDeleter> packed_wc_buffer_;
    std::vector<uint8_t> wc_count_;
    int64_t wc_block_items_{0};
    int64_t wc_block_words_{0};
    bool compact_cursors_{false};
    bool write_combine_{false};
    std::vector<SpecialOffset> special_offsets_;
  };

  /**
   * Iterator to fetch full offset from given buckets
   */
  class OffsetFetcher {
   public:
    OffsetFetcher(BaseSequenceSortingEngine *engine, unsigned start_bucket,
                  unsigned end_bucket)
        : engine_(engine), end_bucket_(end_bucket), cur_bucket_(start_bucket) {
      cur_thread_id_ = 0;
      cur_item_index_ = 0;
      cur_full_offset_ = engine_->thread_meta_[cur_thread_id_].offset_base;
      cur_n_items_ =
          engine_->thread_meta_[cur_thread_id_].bucket_sizes[cur_bucket_];
      diff_index_ = engine_->thread_meta_[0].bucket_begin[cur_bucket_];
      special_index_ = static_cast<size_t>(
          std::lower_bound(engine_->lv1_special_offsets_.begin(),
                           engine_->lv1_special_offsets_.end(), diff_index_,
                           [](const SpecialOffset &entry, int64_t index) {
                             return entry.item_index < index;
                           }) -
          engine_->lv1_special_offsets_.begin());
      TryRefill();
    }

    bool HasNext() const { return cur_bucket_ < end_bucket_; }

    int64_t Next() {
      const uint32_t diff = engine_->LoadOffsetValue(diff_index_);
      if (diff != engine_->OffsetEscapeValue()) {
        cur_full_offset_ += diff;
      } else {
        assert(special_index_ < engine_->lv1_special_offsets_.size());
        const auto &special = engine_->lv1_special_offsets_[special_index_++];
        assert(special.item_index == diff_index_);
        cur_full_offset_ = special.full_offset;
      }

      auto ret = cur_full_offset_;
      ++diff_index_;
      ++cur_item_index_;
      TryRefill();
      return ret;
    }

   private:
    void TryRefill() {
      while (cur_item_index_ >= cur_n_items_) {
        ++cur_thread_id_;
        if (cur_thread_id_ == engine_->n_threads_) {
          cur_thread_id_ = 0;
          ++cur_bucket_;
          if (cur_bucket_ >= end_bucket_) {
            return;
          }
        }
        cur_full_offset_ = engine_->thread_meta_[cur_thread_id_].offset_base;
        cur_item_index_ = 0;
        cur_n_items_ =
            engine_->thread_meta_[cur_thread_id_].bucket_sizes[cur_bucket_];
      }
    }

   private:
    BaseSequenceSortingEngine *engine_;
    unsigned end_bucket_;
    unsigned cur_bucket_;
    unsigned cur_thread_id_;
    int64_t cur_full_offset_;
    size_t cur_item_index_;
    size_t cur_n_items_;
    int64_t diff_index_{};
    size_t special_index_{};
  };

 private:
  uint32_t OffsetEscapeValue() const {
    assert(lv1_offset_bytes_ >= 1 && lv1_offset_bytes_ <= 4);
    return lv1_offset_bytes_ == 4
               ? std::numeric_limits<uint32_t>::max()
               : (uint32_t{1} << (lv1_offset_bytes_ * 8u)) - 1u;
  }

  void StoreOffsetValue(int64_t index, uint32_t value) {
    assert(index >= 0 && index < lv1_offsets_capacity_items_);
    uint8_t *dest = lv1_offsets_.get() + index * lv1_offset_bytes_;
    switch (lv1_offset_bytes_) {
      case 4:
        dest[3] = static_cast<uint8_t>(value >> 24u);
        // fall through
      case 3:
        dest[2] = static_cast<uint8_t>(value >> 16u);
        // fall through
      case 2:
        dest[1] = static_cast<uint8_t>(value >> 8u);
        // fall through
      case 1:
        dest[0] = static_cast<uint8_t>(value);
        return;
      default:
        assert(false);
    }
  }

  uint32_t LoadOffsetValue(int64_t index) const {
    assert(index >= 0 && index < lv1_offsets_capacity_items_);
    const uint8_t *src = lv1_offsets_.get() + index * lv1_offset_bytes_;
    switch (lv1_offset_bytes_) {
      case 1:
        return src[0];
      case 2:
        return static_cast<uint32_t>(src[0]) |
               (static_cast<uint32_t>(src[1]) << 8u);
      case 3:
        return static_cast<uint32_t>(src[0]) |
               (static_cast<uint32_t>(src[1]) << 8u) |
               (static_cast<uint32_t>(src[2]) << 16u);
      case 4:
        return static_cast<uint32_t>(src[0]) |
               (static_cast<uint32_t>(src[1]) << 8u) |
               (static_cast<uint32_t>(src[2]) << 16u) |
               (static_cast<uint32_t>(src[3]) << 24u);
      default:
        assert(false);
        return 0;
    }
  }

  void AddSpecialOffsets(std::vector<SpecialOffset> offsets) {
    std::lock_guard<std::mutex> lk(special_item_lock_);
    lv1_special_offsets_.insert(lv1_special_offsets_.end(), offsets.begin(),
                                offsets.end());
  }

 private:
  /**
   * Memory related data
   */
  int64_t host_mem_{};
  int mem_flag_{};
  unsigned n_threads_{};
  MemoryStat meta_{};  // set by Initialize return
  std::function<void(uint32_t *, int64_t)>
      substr_sort_;  // set after Initialize

  /**
   * Interfaces used by `Run` and must be implemented in derived class
   */
 public:
  virtual MemoryStat Initialize() = 0;

 protected:
  virtual int64_t Lv0EncodeDiffBase(int64_t) = 0;
  virtual void Lv0CalcBucketSize(int64_t seq_from, int64_t seq_to,
                                 std::array<int64_t, kNumBuckets> *out) = 0;
  virtual void Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                              int64_t seq_to) = 0;
  virtual void Lv2ExtractSubString(OffsetFetcher &fetcher,
                                   SubstrPtr substr_ptr) = 0;
  virtual void Lv2Postprocess(int64_t start_index, int64_t end_index,
                              int thread_id, uint32_t *substr_ptr,
                              unsigned bucket_id) = 0;
  virtual void Lv0Postprocess() = 0;
  // A derived engine may replace the generic histogram -> repeated Lv1 scan
  // -> Lv2 gather loop after Initialize() has populated its immutable input
  // state.  Returning true means the complete main-loop work has been done;
  // BaseSequenceSortingEngine still runs the common Lv0Postprocess() tail and
  // timing/reporting code.  The default keeps the traditional path.
  virtual bool RunSpecializedMainLoop() { return false; }
  virtual bool Lv1SupportsDirectItems() const { return false; }
  virtual int64_t Lv1DirectWordsPerItem() const {
    return meta_.words_per_lv2;
  }
  virtual int64_t Lv1DirectAuxWordsPerItem() const {
    return meta_.aux_words_per_lv2;
  }
  virtual int64_t Lv1DirectMemoryLimit() const {
    return std::numeric_limits<int64_t>::max();
  }
  // A producer that can finish every bucket independently may reuse a
  // bounded direct-record arena across several complete input scans.  The
  // default remains all-or-nothing because some algorithms defer a global
  // operation until every sorted bucket is resident.
  virtual bool Lv1AllowsPartialDirectItems() const { return false; }
  virtual unsigned Lv1BytesPerOffset() const {
    return kDefaultLv1BytesPerItem;
  }
  // A positive value makes automatic mode use up to this much memory for the
  // compact Lv1 locator stream plus Lv2 radix scratch.  Zero retains the
  // traditional fixed-pass heuristic for producers that have not opted in.
  virtual int64_t Lv1AutoWorkspaceLimit() const { return 0; }
  virtual bool Lv1UseWriteCombine() const {
    return std::getenv("MEGAHIT_DISABLE_LV1_WRITE_COMBINE") == nullptr;
  }
  virtual int Lv2SortIgnoredLowBytes() const { return 0; }
  virtual int Lv2SortIgnoredHighBytes() const { return 2; }
  virtual bool Lv2NeedsEmptyBucketPostprocess(unsigned) const {
    return false;
  }
  // Some derived algorithms need every direct bucket to remain sorted and
  // resident before they can perform a global exact set operation.  The
  // default preserves the traditional per-bucket postprocess path.
  virtual bool Lv2DeferPostprocess() const { return false; }
  virtual void Lv2PostprocessDeferred() {}
  // Derived producers with highly variable sequence lengths can provide
  // contiguous ranges balanced by the number of records they will emit.
  // Contiguity preserves the monotone locator stream required by Lv1's
  // differential encoding.  Returning false keeps the traditional equal
  // sequence-count partition.
  virtual bool Lv0BuildBalancedRanges(
      std::vector<std::pair<int64_t, int64_t>> *ranges) const {
    (void)ranges;
    return false;
  }
  // A derived packed-locator producer can certify that all locators emitted
  // by this sequence partition are within one uint32 range.  OffsetFiller can
  // then combine its write index and previous relative locator in one cursor.
  virtual bool Lv1CanUseCompactCursor(int64_t seq_from,
                                      int64_t seq_to) const {
    (void)seq_from;
    (void)seq_to;
    return false;
  }

  bool UsingDirectLv1Items() const { return lv1_direct_mode_; }
  uint32_t *SortedDirectBucketData(unsigned real_bucket) {
    assert(real_bucket < kNumBuckets);
    if (!lv1_direct_mode_) {
      assert(bucket_sizes_[bucket_rank_[real_bucket]] == 0);
      return nullptr;
    }
    const unsigned rank = bucket_rank_[real_bucket];
    return lv1_direct_items_.get() +
           thread_meta_[0].bucket_begin[rank] *
               lv1_direct_words_per_item_;
  }
  int64_t SortedDirectBucketSize(unsigned real_bucket) const {
    assert(real_bucket < kNumBuckets);
    return bucket_sizes_[bucket_rank_[real_bucket]];
  }

  // May be called from a derived Initialize() as soon as the input metadata
  // gives an item-count estimate, before the input itself is loaded.  Applies
  // the same policy gates as AdjustMemory(); a rejected or outgrown
  // speculation is simply freed there, falling back to the normal path.
  void PreallocateDirectItems(int64_t estimated_items,
                              int64_t estimated_memory_for_data);

 private:
  void JoinPrefaultThreads();
  void AdjustMemory();
  void Lv0PrepareThreadPartition();
  void Lv0CalcBucketSizeLaunchMt();
  void Lv0ReorderBuckets();
  unsigned Lv1FindEndBuckets(unsigned start_bucket);
  void Lv1ComputeThreadBegin();
  inline void Lv1FillOffsetsLaunchMt();
  void Lv1FetchAndSortLaunchMt();

  struct Lv2ThreadStatus {
    std::vector<int64_t> thread_offset;
    std::vector<int64_t> rank;
    std::vector<double> extract_seconds;
    std::vector<double> sort_seconds;
    std::vector<double> postprocess_seconds;
    std::mutex mutex;
    int64_t acc{0};
    int seen{0};
    bool profile{false};
  };

  void Lv2Sort(Lv2ThreadStatus *thread_status, unsigned b, int tid);

 public:
  void Run();
};

#endif  // MEGAHIT_BASE_ENGINE_H
