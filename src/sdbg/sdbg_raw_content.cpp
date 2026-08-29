//
// Created by vout on 11/5/18.
//

#include "sdbg_raw_content.h"
#include "sdbg_meta.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#endif

#include "kmlib/kmbit.h"
#include "sdbg_item.h"
#include "utils/buffered_reader.h"

namespace {

struct SdbgFileBuckets {
  size_t file_id;
  std::vector<const SdbgBucketRecord *> buckets;
};

// std::vector::resize value-initializes large arrays on the calling thread.
// On a multi-socket host that silently places every hot SDBG page on one NUMA
// node even though the following file decode is parallel.  The arrays below
// are completely populated by the decode (packed zero fields deliberately do
// not need a store), so discard only their page-aligned interiors before the
// workers start.  Linux supplies fresh zero pages on first write and physical
// placement then follows the worker that owns that bucket.  Boundary pages
// are retained to avoid touching unrelated allocator objects.  On other
// platforms, or if madvise is unavailable, this is simply a no-op.
template <typename T>
void RearmForParallelFirstTouch(T *data, size_t count) {
#ifdef __linux__
  if (data == nullptr || count == 0 ||
      count > std::numeric_limits<size_t>::max() / sizeof(T)) {
    return;
  }
  const long page_size_value = ::sysconf(_SC_PAGESIZE);
  if (page_size_value <= 0) {
    return;
  }
  const uintptr_t page_size = static_cast<uintptr_t>(page_size_value);
  const uintptr_t raw_begin = reinterpret_cast<uintptr_t>(data);
  const size_t bytes = count * sizeof(T);
  if (raw_begin > std::numeric_limits<uintptr_t>::max() - bytes) {
    return;
  }
  const uintptr_t raw_end = raw_begin + bytes;
  const uintptr_t page_begin =
      raw_begin + ((page_size - raw_begin % page_size) % page_size);
  const uintptr_t page_end = raw_end - raw_end % page_size;
  if (page_begin < page_end) {
    (void)::madvise(reinterpret_cast<void *>(page_begin),
                    page_end - page_begin, MADV_DONTNEED);
  }
#else
  (void)data;
  (void)count;
#endif
}

// CompactVector::Adapter updates a packed word with a read-modify-write.  Two
// adjacent buckets may therefore race on the word containing their boundary.
// All packed vectors are zero initialized here, so boundary fields can be set
// with an atomic OR.  Words strictly inside a bucket are owned by that bucket
// and keep the cheaper non-atomic update.
template <unsigned BaseSize, typename WordType>
inline void StorePacked(kmlib::CompactVector<BaseSize, WordType> *dst,
                        size_t index, WordType value, size_t bucket_begin,
                        size_t bucket_end) {
  if (value == 0) {
    return;
  }

  const size_t word_index = index / dst->kBasesPerWord;
  const size_t first_word = bucket_begin / dst->kBasesPerWord;
  const size_t last_word = (bucket_end - 1) / dst->kBasesPerWord;
  if (word_index == first_word || word_index == last_word) {
    const WordType shifted =
        value << dst->bit_shift(index % dst->kBasesPerWord);
    __atomic_fetch_or(dst->data() + word_index, shifted, __ATOMIC_RELAXED);
  } else {
    (*dst)[index] = value;
  }
}

}  // namespace

/**
 * load SDBG raw content from disk
 * @param raw_content the pointer the raw_content to be stored to
 * @param file_prefix the prefix of the SDBG files
 */
void LoadSdbgRawContent(SdbgRawContent *raw_content,
                        const std::string &file_prefix) {
  std::ifstream is(file_prefix + ".sdbg_info");
  if (!is.is_open()) {
    throw std::runtime_error("Failed to open SDBG metadata for " +
                             file_prefix);
  }
  raw_content->meta.Deserialize(is);
  const auto &metadata = raw_content->meta;
  // clear() is needed when an SdbgRawContent object is loaded more than once:
  // StorePacked relies on initially zero packed words.
  raw_content->w.clear();
  raw_content->last.clear();
  raw_content->tip.clear();
  raw_content->w.resize(metadata.item_count());
  raw_content->last.resize(metadata.item_count());
  raw_content->tip.resize(metadata.item_count());
  raw_content->tip_lables.resize(metadata.words_per_tip_label() *
                                 metadata.tip_count());
  raw_content->small_mul.clear();
  raw_content->large_mul_bits.clear();
  raw_content->large_mul_rank.clear();
  raw_content->large_mul_values.clear();
  raw_content->full_mul.clear();
  bool use_full_mul =
      metadata.large_mul_count() >= metadata.item_count() * 0.08;
  if (use_full_mul) {
    raw_content->full_mul.resize(metadata.item_count());
  } else {
    raw_content->small_mul.resize(metadata.item_count());
    raw_content->large_mul_bits.resize(
        (metadata.item_count() + 63u) / 64u);
    raw_content->large_mul_values.resize(metadata.large_mul_count());
  }

  // Let the parallel bucket decoder, rather than this serial allocation
  // phase, establish the NUMA ownership of every long-lived hot array.
  RearmForParallelFirstTouch(raw_content->w.data(),
                             raw_content->w.word_count());
  RearmForParallelFirstTouch(raw_content->last.data(),
                             raw_content->last.word_count());
  RearmForParallelFirstTouch(raw_content->tip.data(),
                             raw_content->tip.word_count());
  RearmForParallelFirstTouch(raw_content->tip_lables.data(),
                             raw_content->tip_lables.size());
  RearmForParallelFirstTouch(raw_content->small_mul.data(),
                             raw_content->small_mul.size());
  RearmForParallelFirstTouch(raw_content->full_mul.data(),
                             raw_content->full_mul.size());
  RearmForParallelFirstTouch(raw_content->large_mul_bits.data(),
                             raw_content->large_mul_bits.size());
  RearmForParallelFirstTouch(raw_content->large_mul_values.data(),
                             raw_content->large_mul_values.size());

  std::vector<SdbgFileBuckets> files;
  for (auto bucket_it = metadata.begin_bucket();
       bucket_it != metadata.end_bucket() &&
       bucket_it->bucket_id != bucket_it->kNullID;
       ++bucket_it) {
    if (files.empty() || files.back().file_id != bucket_it->file_id) {
      files.push_back(SdbgFileBuckets{bucket_it->file_id, {}});
    }
    files.back().buckets.push_back(&*bucket_it);
  }

  std::atomic<bool> load_failed{false};

#pragma omp parallel for schedule(dynamic, 1)
  for (size_t file_index = 0; file_index < files.size(); ++file_index) {
    const auto &file = files[file_index];
    std::ifstream input((file_prefix + ".sdbg." + std::to_string(file.file_id))
                            .c_str(),
                        std::ifstream::binary | std::ifstream::in);
    if (!input.is_open()) {
      load_failed.store(true, std::memory_order_relaxed);
      continue;
    }
    input.seekg(0, std::ifstream::end);
    const std::streamoff file_size_stream = input.tellg();
    if (file_size_stream < 0) {
      load_failed.store(true, std::memory_order_relaxed);
      continue;
    }
    const size_t file_size = static_cast<size_t>(file_size_stream);
    input.seekg(0, std::ifstream::beg);

    BufferedReader in;
    in.reset(&input);
    size_t file_offset = 0;
    bool local_failed = false;

    for (const SdbgBucketRecord *bucket : file.buckets) {
      if (file_offset != bucket->starting_offset) {
        local_failed = true;
        break;
      }

      const size_t bucket_begin = bucket->accumulate_item_count;
      const size_t bucket_end = bucket_begin + bucket->num_items;
      size_t tip_label_index =
          bucket->accumulate_tip_count * metadata.words_per_tip_label();
      size_t large_mul_index = bucket->accumulate_large_mul_count;

      for (size_t i = 0; i < bucket->num_items; ++i) {
        SdbgItem item;
        size_t bytes_read = in.read(&item);
        file_offset += bytes_read;
        if (bytes_read != sizeof(item)) {
          local_failed = true;
          break;
        }

        const size_t index = bucket_begin + i;
        StorePacked(&raw_content->w, index,
                    static_cast<uint64_t>(item.w), bucket_begin, bucket_end);
        StorePacked(&raw_content->last, index,
                    static_cast<uint64_t>(item.last), bucket_begin, bucket_end);
        StorePacked(&raw_content->tip, index,
                    static_cast<uint64_t>(item.tip), bucket_begin, bucket_end);

        mul_t mul = item.mul;
        if (mul == kSmallMulSentinel) {
          bytes_read = in.read(&mul);
          file_offset += bytes_read;
          if (bytes_read != sizeof(mul)) {
            local_failed = true;
            break;
          }
        }
        if (use_full_mul) {
          raw_content->full_mul[index] = mul;
        } else if (mul <= kMaxSmallMul) {
          raw_content->small_mul[index] = mul;
        } else {
          raw_content->small_mul[index] = kSmallMulSentinel;
          __atomic_fetch_or(
              raw_content->large_mul_bits.data() + (index >> 6u),
              uint64_t{1} << (index & 63u), __ATOMIC_RELAXED);
          raw_content->large_mul_values[large_mul_index++] = mul;
        }

        if (item.tip) {
          const size_t label_bytes =
              sizeof(label_word_t) * metadata.words_per_tip_label();
          label_word_t *tip_label_ptr =
              raw_content->tip_lables.data() + tip_label_index;
          bytes_read = in.read(tip_label_ptr,
                               metadata.words_per_tip_label());
          file_offset += bytes_read;
          if (bytes_read != label_bytes) {
            local_failed = true;
            break;
          }
          for (unsigned j = 0; j < metadata.words_per_tip_label(); ++j) {
            tip_label_ptr[j] =
                kmlib::bit::Reverse<kBitsPerChar>(tip_label_ptr[j]);
          }
          tip_label_index += metadata.words_per_tip_label();
        }
      }

      if (local_failed) {
        break;
      }
      if (tip_label_index !=
          (bucket->accumulate_tip_count + bucket->num_tips) *
              metadata.words_per_tip_label()) {
        local_failed = true;
        break;
      }
      if (large_mul_index !=
          bucket->accumulate_large_mul_count + bucket->num_large_mul) {
        local_failed = true;
        break;
      }
    }

    if (file_offset != file_size) {
      local_failed = true;
    }
    if (local_failed) {
      load_failed.store(true, std::memory_order_relaxed);
    }
  }

  if (load_failed.load(std::memory_order_relaxed)) {
    throw std::runtime_error("Failed to load SDBG raw content from " +
                             file_prefix);
  }

  if (!use_full_mul) {
    const size_t block_edges = SdbgRawContent::kLargeMulRankBlockEdges;
    const size_t num_blocks =
        (metadata.item_count() + block_edges - 1u) / block_edges;
    raw_content->large_mul_rank.resize(num_blocks + 1u);
#pragma omp parallel for schedule(static)
    for (int64_t block = 0; block < static_cast<int64_t>(num_blocks);
         ++block) {
      const size_t first_word = static_cast<size_t>(block) * block_edges / 64u;
      const size_t end_edge = std::min(
          metadata.item_count(),
          (static_cast<size_t>(block) + 1u) * block_edges);
      const size_t end_word = (end_edge + 63u) / 64u;
      uint64_t count = 0;
      for (size_t word = first_word; word < end_word; ++word) {
        uint64_t bits = raw_content->large_mul_bits[word];
        if (word + 1u == end_word && (end_edge & 63u) != 0u) {
          bits &= (uint64_t{1} << (end_edge & 63u)) - 1u;
        }
        count += static_cast<uint64_t>(__builtin_popcountll(bits));
      }
      raw_content->large_mul_rank[static_cast<size_t>(block) + 1u] = count;
    }
    for (size_t block = 0; block < num_blocks; ++block) {
      raw_content->large_mul_rank[block + 1u] +=
          raw_content->large_mul_rank[block];
    }
    if (raw_content->large_mul_rank.back() !=
        raw_content->large_mul_values.size()) {
      throw std::runtime_error("Corrupt large-multiplicity rank directory");
    }
  }
}
