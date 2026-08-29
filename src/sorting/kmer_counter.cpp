#include "kmer_counter.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <limits>
#include <mutex>
#include <numeric>
#include <omp.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

#ifdef __linux__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/uio.h>
#include <unistd.h>
#endif

#include "sequence/copy_substr.h"
#include "sequence/io/edge/edge_bucket_histogram.h"
#include "sequence/io/read_chunk_index.h"
#include "sequence/io/sequence_lib.h"
#include "sequence/kmer.h"
#include "kmlib/kmsort.h"
#include "sorting/kmsort_selector.h"
#include "utils/startup_affinity.h"
#include "utils/utils.h"

// A read-only view of buildlib's native [length][packed words] stream.  Count
// only needs one producer scan when its minimizer runs carry the exact bases
// required by the consumer.  Mapping the already-packed stream therefore
// avoids creating a second 2-bit copy of every read.  The persistent chunk
// index supplies independent, validated record ranges; completed source pages
// can be dropped immediately because no later phase dereferences them.
class MappedCountReads {
 public:
  using Word = SeqPackage::TWord;
  using Vector = SeqPackage::TVector;

  class ReverseView {
   public:
    ReverseView(const Word *words, unsigned length)
        : words_(words), length_(length) {}
    unsigned length() const { return length_; }
    uint8_t base_at(unsigned position) const {
      assert(position < length_);
      return static_cast<uint8_t>(
          Vector::at(words_, static_cast<size_t>(length_ - 1u - position)));
    }

   private:
    const Word *words_;
    unsigned length_;
  };

  MappedCountReads() = default;
  ~MappedCountReads() { Close(); }
  MappedCountReads(const MappedCountReads &) = delete;
  MappedCountReads &operator=(const MappedCountReads &) = delete;

  bool Open(const std::string &path, uint64_t expected_reads,
            uint64_t expected_bases, unsigned expected_max_length) {
#ifdef __linux__
    Close();
    if (!LoadPackedReadChunkIndex(path, &index_) ||
        index_.num_reads != expected_reads ||
        index_.num_bases != expected_bases ||
        index_.max_read_len != expected_max_length) {
      return false;
    }
    fd_ = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) return false;
    struct stat status;
    if (fstat(fd_, &status) != 0 || status.st_size <= 0 ||
        status.st_size % static_cast<off_t>(sizeof(Word)) != 0) {
      Close();
      return false;
    }
    bytes_ = static_cast<size_t>(status.st_size);
    void *mapping = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (mapping == MAP_FAILED) {
      mapping_ = nullptr;
      Close();
      return false;
    }
    mapping_ = static_cast<const Word *>(mapping);
#if defined(MADV_SEQUENTIAL)
    (void)madvise(const_cast<Word *>(mapping_), bytes_, MADV_SEQUENTIAL);
#endif
    // Count's first pass already visits every read.  Record one packed-word
    // offset per small read block during that pass so the later mercy
    // candidate writer can jump over blocks with no selected reads instead
    // of re-reading the complete packed library.  The stride is a layout
    // property, not a data-dependent threshold: smaller strides trade a
    // compact O(reads / stride) exact index for less sparse source traffic.
    candidate_chunk_offsets_.resize(index_.chunks.size() + 1u, 0);
    uint64_t checkpoint_count = 0;
    for (size_t chunk_id = 0; chunk_id < index_.chunks.size(); ++chunk_id) {
      candidate_chunk_offsets_[chunk_id] = checkpoint_count;
      const PackedReadChunk &chunk = index_.chunks[chunk_id];
      checkpoint_count += DivCeiling(chunk.read_end - chunk.read_begin,
                                     kCandidateReadStride);
    }
    candidate_chunk_offsets_.back() = checkpoint_count;
    candidate_word_offsets_.reset(new uint64_t[checkpoint_count]);
    return true;
#else
    (void)path;
    (void)expected_reads;
    (void)expected_bases;
    (void)expected_max_length;
    return false;
#endif
  }

  void Close() {
#ifdef __linux__
    if (mapping_ != nullptr) {
      munmap(const_cast<Word *>(mapping_), bytes_);
      mapping_ = nullptr;
    }
    if (fd_ >= 0) {
      close(fd_);
      fd_ = -1;
    }
#endif
    bytes_ = 0;
    index_ = PackedReadChunkIndex();
    candidate_word_offsets_.reset();
    candidate_chunk_offsets_.clear();
  }

  const Word *words() const { return mapping_; }
  const PackedReadChunkIndex &index() const { return index_; }
  uint64_t num_reads() const { return index_.num_reads; }
  uint64_t num_bases() const { return index_.num_bases; }
  unsigned max_length() const { return index_.max_read_len; }
  size_t mapped_bytes() const { return bytes_; }

  void RecordCandidateCheckpoint(size_t chunk_id, uint64_t read_id,
                                 const Word *record) {
    assert(chunk_id < index_.chunks.size());
    const PackedReadChunk &chunk = index_.chunks[chunk_id];
    assert(read_id >= chunk.read_begin && read_id < chunk.read_end);
    const uint64_t local_read = read_id - chunk.read_begin;
    if (local_read % kCandidateReadStride != 0) return;
    const uint64_t checkpoint =
        candidate_chunk_offsets_[chunk_id] +
        local_read / kCandidateReadStride;
    assert(checkpoint < candidate_chunk_offsets_[chunk_id + 1u]);
    candidate_word_offsets_[checkpoint] =
        static_cast<uint64_t>(record - mapping_);
  }

  void DiscardChunk(const PackedReadChunk &chunk) const {
    DiscardMemoryPages(const_cast<Word *>(mapping_ + chunk.word_begin),
                       static_cast<size_t>(chunk.word_end - chunk.word_begin) *
                           sizeof(Word));
  }

  bool WriteSelectedReversed(std::ostream *output,
                             const std::vector<uint64_t> &selected) const {
    if (output == nullptr || mapping_ == nullptr) return false;
#if defined(__linux__) && defined(MADV_SEQUENTIAL)
    // Candidate blocks are sparse logically but are consumed in strictly
    // increasing file order.  MADV_RANDOM turned every selected 4-KiB page
    // into an independent major fault after the count spool evicted the read
    // package (millions of faults on CAMI3).  Restore sequential readahead;
    // DiscardChunk below drops each completed bounded chunk immediately, so
    // readahead no longer accumulates the whole library in RSS.
    (void)madvise(const_cast<Word *>(mapping_), bytes_, MADV_SEQUENTIAL);
    if (fd_ >= 0) {
      (void)posix_fadvise(fd_, 0, static_cast<off_t>(bytes_),
                         POSIX_FADV_SEQUENTIAL);
    }
#endif
    constexpr size_t kBitsPerWord = 64u;
    std::array<uint32_t, 16> reversed{};
    uint64_t visited_blocks = 0;
    uint64_t total_blocks = 0;
    uint64_t visited_source_bytes = 0;
    auto any_selected = [&](uint64_t begin, uint64_t end) {
      if (begin >= end || begin / kBitsPerWord >= selected.size()) {
        return false;
      }
      const size_t first_word = static_cast<size_t>(begin / kBitsPerWord);
      const size_t last_word = static_cast<size_t>((end - 1u) / kBitsPerWord);
      for (size_t word_id = first_word;
           word_id <= last_word && word_id < selected.size(); ++word_id) {
        uint64_t mask = std::numeric_limits<uint64_t>::max();
        if (word_id == first_word) {
          mask &= std::numeric_limits<uint64_t>::max()
                  << (begin % kBitsPerWord);
        }
        if (word_id == last_word && end % kBitsPerWord != 0) {
          mask &= (uint64_t{1} << (end % kBitsPerWord)) - 1u;
        }
        if ((selected[word_id] & mask) != 0) return true;
      }
      return false;
    };

    for (size_t chunk_id = 0; chunk_id < index_.chunks.size(); ++chunk_id) {
      const PackedReadChunk &chunk = index_.chunks[chunk_id];
      const uint64_t num_checkpoints =
          candidate_chunk_offsets_[chunk_id + 1u] -
          candidate_chunk_offsets_[chunk_id];
      total_blocks += num_checkpoints;
      for (uint64_t block = 0; block < num_checkpoints; ++block) {
        const uint64_t read_begin =
            chunk.read_begin + block * kCandidateReadStride;
        const uint64_t read_end = std::min<uint64_t>(
            chunk.read_end, read_begin + kCandidateReadStride);
        if (!any_selected(read_begin, read_end)) continue;
        ++visited_blocks;
        const uint64_t checkpoint =
            candidate_chunk_offsets_[chunk_id] + block;
        const Word *cursor = mapping_ + candidate_word_offsets_[checkpoint];
        if (cursor < mapping_ + chunk.word_begin ||
            cursor >= mapping_ + chunk.word_end) {
          return false;
        }
        for (uint64_t read_id = read_begin; read_id < read_end; ++read_id) {
          if (cursor >= mapping_ + chunk.word_end) return false;
          const uint32_t length = *cursor++;
          const size_t packed_words =
              DivCeiling(static_cast<size_t>(length),
                         static_cast<size_t>(SeqPackage::kBasesPerWord));
          if (packed_words >
                  static_cast<size_t>(mapping_ + chunk.word_end - cursor) ||
              packed_words > reversed.size()) {
            return false;
          }
          visited_source_bytes +=
              (packed_words + 1u) * sizeof(uint32_t);
          const bool keep =
              read_id / kBitsPerWord < selected.size() &&
              ((selected[read_id / kBitsPerWord] >>
                (read_id % kBitsPerWord)) &
               uint64_t{1}) != 0;
          if (keep) {
            reversed.fill(0);
            ReverseView view(cursor, length);
            for (unsigned position = 0; position < length; ++position) {
              const unsigned word = position / SeqPackage::kBasesPerWord;
              const unsigned in_word = position % SeqPackage::kBasesPerWord;
              reversed[word] |=
                  static_cast<uint32_t>(view.base_at(position))
                  << ((SeqPackage::kBasesPerWord - 1u - in_word) * 2u);
            }
            output->write(reinterpret_cast<const char *>(&length),
                          sizeof(length));
            output->write(reinterpret_cast<const char *>(reversed.data()),
                          packed_words * sizeof(uint32_t));
            if (!*output) return false;
          }
          cursor += packed_words;
        }
      }
      // Candidate blocks are sparse but span the complete library.  Drop the
      // completed chunk immediately; otherwise those sparse 4-KiB pages
      // accumulate until the process has effectively retained the whole
      // packed-read mapping despite visiting only a small byte subset.
      DiscardChunk(chunk);
    }
    xinfo("Sparse mercy candidate gather: {} of {} read blocks, {.3} GiB "
          "packed source visited\n",
          visited_blocks, total_blocks,
          static_cast<double>(visited_source_bytes) /
              static_cast<double>(uint64_t{1} << 30u));
    return true;
  }

 private:
  static constexpr uint64_t kCandidateReadStride = 32u;
  int fd_{-1};
  const Word *mapping_{nullptr};
  size_t bytes_{0};
  PackedReadChunkIndex index_;
  std::unique_ptr<uint64_t[]> candidate_word_offsets_;
  std::vector<uint64_t> candidate_chunk_offsets_;
};

KmerCounter::KmerCounter(const KmerCounterOption &opt)
    : BaseSequenceSortingEngine(opt.host_mem, opt.mem_flag, opt.n_threads),
      opt_(opt) {}

KmerCounter::~KmerCounter() = default;

namespace {

// A minimizer-run descriptor contains a producer-local read ID, two byte-sized
// coordinates, and a temporary shard ID.  The shard field is discarded after
// the in-place grouping step.  Runtime shard selection is capped by this
// representation rather than by a particular socket/thread count.
constexpr unsigned kCountMaxShardBits = 16;
constexpr unsigned kCountExternalMaxShardBits = 20;
constexpr unsigned kCountSegmentShardShift = 48;
constexpr unsigned kCountSegmentLengthShift = 40;
constexpr unsigned kCountSegmentStartShift = 32;

inline unsigned CountBytesForValue(uint64_t value) {
  unsigned bytes = 1;
  while (bytes < sizeof(uint32_t) && value >= (uint64_t{1} << (bytes * 8u))) {
    ++bytes;
  }
  return bytes;
}

// A vector's geometric doubling is costly when hundreds of millions of run
// descriptors are produced concurrently: the unused capacity alone can be
// several GiB.  This movable stream grows by 12.5% and uses realloc/mremap when
// the allocator can extend it in place.  Read IDs are relative to the
// producer's contiguous input range, so the common large-team case needs three
// bytes instead of four.  Correctness never depends on the statistical initial
// reserve or on a fixed number of workers.
class CountSegmentStream {
 public:
  CountSegmentStream() = default;
  ~CountSegmentStream() { std::free(data_); }
  CountSegmentStream(const CountSegmentStream &) = delete;
  CountSegmentStream &operator=(const CountSegmentStream &) = delete;
  CountSegmentStream(CountSegmentStream &&other) noexcept
      : data_(other.data_),
        size_(other.size_),
        capacity_(other.capacity_),
        read_base_(other.read_base_),
        read_id_bytes_(other.read_id_bytes_),
        shard_bytes_(other.shard_bytes_),
        record_bytes_(other.record_bytes_),
        grouped_bytes_(other.grouped_bytes_),
        pending_size_(other.pending_size_),
        grouped_compact_(other.grouped_compact_) {
    std::memcpy(pending_, other.pending_, sizeof(pending_));
    other.data_ = nullptr;
    other.size_ = 0;
    other.capacity_ = 0;
    other.record_bytes_ = 0;
    other.grouped_bytes_ = 0;
    other.pending_size_ = 0;
    other.grouped_compact_ = false;
  }
  CountSegmentStream &operator=(CountSegmentStream &&other) noexcept {
    if (this != &other) {
      std::free(data_);
      data_ = other.data_;
      size_ = other.size_;
      capacity_ = other.capacity_;
      read_base_ = other.read_base_;
      read_id_bytes_ = other.read_id_bytes_;
      shard_bytes_ = other.shard_bytes_;
      record_bytes_ = other.record_bytes_;
      grouped_bytes_ = other.grouped_bytes_;
      pending_size_ = other.pending_size_;
      std::memcpy(pending_, other.pending_, sizeof(pending_));
      grouped_compact_ = other.grouped_compact_;
      other.data_ = nullptr;
      other.size_ = 0;
      other.capacity_ = 0;
      other.record_bytes_ = 0;
      other.grouped_bytes_ = 0;
      other.pending_size_ = 0;
      other.grouped_compact_ = false;
    }
    return *this;
  }

  void Configure(uint32_t read_base, uint64_t read_count,
                 unsigned shard_bits) {
    assert(data_ == nullptr && size_ == 0 && capacity_ == 0);
    assert(shard_bits > 0 && shard_bits <= kCountMaxShardBits);
    read_base_ = read_base;
    read_id_bytes_ = static_cast<uint8_t>(
        CountBytesForValue(read_count == 0 ? 0 : read_count - 1u));
    shard_bytes_ = static_cast<uint8_t>((shard_bits + 7u) / 8u);
    grouped_bytes_ = static_cast<uint8_t>(read_id_bytes_ + 2u);
    record_bytes_ = static_cast<uint8_t>(grouped_bytes_ + shard_bytes_);
    assert(record_bytes_ <= sizeof(uint64_t));
  }

  void reserve(size_t capacity) {
    assert(!grouped_compact_ && record_bytes_ != 0);
    if (capacity > capacity_) {
      Reallocate(capacity);
    }
  }
  void push_back(uint32_t read_id, unsigned start, unsigned length,
                 unsigned shard) {
    assert(!grouped_compact_ && record_bytes_ != 0);
    assert(read_id >= read_base_);
    assert(start <= std::numeric_limits<uint8_t>::max());
    assert(length > 0 && length <= std::numeric_limits<uint8_t>::max());
    assert(shard < (1u << kCountMaxShardBits));
    if (size_ == capacity_) {
      const size_t increment =
          std::max<size_t>(capacity_ / 8u, size_t{1} << 16u);
      if (capacity_ > std::numeric_limits<size_t>::max() - increment) {
        xfatal("Count segment stream capacity overflow\n");
      }
      Reallocate(capacity_ + increment);
    }
    const uint64_t read_delta = static_cast<uint64_t>(read_id - read_base_);
    assert(read_id_bytes_ == sizeof(uint32_t) ||
           read_delta < (uint64_t{1} << (read_id_bytes_ * 8u)));
    const unsigned start_shift = read_id_bytes_ * 8u;
    const unsigned length_shift = start_shift + 8u;
    const unsigned shard_shift = length_shift + 8u;
    const uint64_t value = read_delta |
                           (static_cast<uint64_t>(start) << start_shift) |
                           (static_cast<uint64_t>(length) << length_shift) |
                           (static_cast<uint64_t>(shard) << shard_shift);
    ++size_;
    if (record_bytes_ == 7u) {
      pending_[pending_size_++] = value;
      if (pending_size_ == 8u) {
        FlushPendingSeven();
      }
    } else {
      StoreLowBytes(data_ + (size_ - 1u) * record_bytes_, value,
                    record_bytes_);
    }
  }

  void Finalize() {
    assert(!grouped_compact_);
    if (pending_size_ != 0) {
      assert(record_bytes_ == 7u);
      uint8_t *dest = data_ + (size_ - pending_size_) * record_bytes_;
      for (unsigned i = 0; i < pending_size_; ++i) {
        StoreLowBytes(dest + i * record_bytes_, pending_[i], record_bytes_);
      }
      pending_size_ = 0;
    }
  }

  unsigned ShardAt(size_t index) const {
    assert(!grouped_compact_ && pending_size_ == 0 && index < size_);
    const uint8_t *source =
        data_ + index * record_bytes_ + grouped_bytes_;
    if (shard_bytes_ == 1) {
      return source[0];
    }
    uint16_t shard;
    std::memcpy(&shard, source, sizeof(shard));
    return shard;
  }

  void Swap(size_t lhs, size_t rhs) {
    assert(!grouped_compact_ && pending_size_ == 0 && lhs < size_ &&
           rhs < size_);
    if (lhs == rhs) {
      return;
    }
    uint8_t *lhs_ptr = data_ + lhs * record_bytes_;
    uint8_t *rhs_ptr = data_ + rhs * record_bytes_;
    const uint64_t lhs_value = LoadLowBytes(lhs_ptr, record_bytes_);
    const uint64_t rhs_value = LoadLowBytes(rhs_ptr, record_bytes_);
    StoreLowBytes(lhs_ptr, rhs_value, record_bytes_);
    StoreLowBytes(rhs_ptr, lhs_value, record_bytes_);
  }

  uint8_t *data() { return data_; }
  const uint8_t *data() const { return data_; }
  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }
  size_t used_bytes() const { return size_ * record_bytes_; }
  size_t capacity_bytes() const { return capacity_ * record_bytes_; }
  size_t allocation_bytes() const {
    return capacity_ == 0 ? 0 : capacity_bytes() + sizeof(uint64_t);
  }
  size_t grouped_used_bytes() const { return size_ * grouped_bytes_; }

  // Once the shard permutation is complete, the shard ID is implicit in the
  // immutable offset table.  Pack read-delta/start/length forward in the same
  // allocation, then shrink it.  Forward packing is overlap-safe because the
  // destination record is smaller than the source record.
  void CompactGrouped() {
    assert(!grouped_compact_ && pending_size_ == 0);
    const size_t old_record_bytes = record_bytes_;
    for (size_t i = 0; i < size_; ++i) {
      const uint64_t descriptor =
          LoadLowBytes(data_ + i * old_record_bytes, old_record_bytes);
      StoreLowBytes(data_ + i * grouped_bytes_, descriptor, grouped_bytes_);
    }
    const size_t compact_bytes = size_ * grouped_bytes_;
    const size_t compact_allocation = compact_bytes + sizeof(uint64_t);
    const size_t old_allocation =
        capacity_ * old_record_bytes + sizeof(uint64_t);
    if (compact_allocation < old_allocation) {
      DiscardMemoryPages(data_ + compact_allocation,
                         old_allocation - compact_allocation);
    }
    if (compact_allocation != 0) {
      void *memory = std::realloc(data_, compact_allocation);
      if (memory != nullptr) {
        data_ = static_cast<uint8_t *>(memory);
      }
      AdviseHugePages(data_, compact_allocation);
    }
    capacity_ = size_;
    record_bytes_ = grouped_bytes_;
    grouped_compact_ = true;
  }

  uint64_t GroupedAt(size_t index) const {
    assert(grouped_compact_ && index < size_);
    const uint64_t descriptor =
        LoadLowBytes(data_ + index * record_bytes_, record_bytes_);
    const unsigned start_shift = read_id_bytes_ * 8u;
    const unsigned length_shift = start_shift + 8u;
    const uint64_t read_mask =
        read_id_bytes_ == sizeof(uint32_t)
            ? std::numeric_limits<uint32_t>::max()
            : (uint64_t{1} << start_shift) - 1u;
    const uint32_t read_id =
        read_base_ + static_cast<uint32_t>(descriptor & read_mask);
    return static_cast<uint64_t>(read_id) |
           (((descriptor >> start_shift) & 0xFFu)
            << kCountSegmentStartShift) |
           (((descriptor >> length_shift) & 0xFFu)
            << kCountSegmentLengthShift);
  }

 private:
  void FlushPendingSeven() {
    assert(record_bytes_ == 7u && pending_size_ == 8u && size_ >= 8u);
    // Eight 56-bit descriptors form exactly seven machine words.  Packing a
    // small L1-resident batch restores streaming full-width stores without the
    // overlapping writes that hurt the large working set.
    uint64_t packed[7];
    packed[0] = pending_[0] | (pending_[1] << 56u);
    packed[1] = (pending_[1] >> 8u) | (pending_[2] << 48u);
    packed[2] = (pending_[2] >> 16u) | (pending_[3] << 40u);
    packed[3] = (pending_[3] >> 24u) | (pending_[4] << 32u);
    packed[4] = (pending_[4] >> 32u) | (pending_[5] << 24u);
    packed[5] = (pending_[5] >> 40u) | (pending_[6] << 16u);
    packed[6] = (pending_[6] >> 48u) | (pending_[7] << 8u);
    std::memcpy(data_ + (size_ - 8u) * record_bytes_, packed,
                sizeof(packed));
    pending_size_ = 0;
  }

  static uint64_t LoadLowBytes(const uint8_t *source, unsigned bytes) {
    uint64_t value;
    std::memcpy(&value, source, sizeof(value));
    if (bytes != sizeof(uint64_t)) {
      value &= (uint64_t{1} << (bytes * 8u)) - 1u;
    }
    return value;
  }

  static void StoreLowBytes(uint8_t *dest, uint64_t value, unsigned bytes) {
    // These fixed stores keep the seven-byte common path branch-free inside
    // the allocator and avoid a library memcpy call in the grouping loop.
    switch (bytes) {
      case 8:
        std::memcpy(dest, &value, sizeof(value));
        break;
      case 7: {
        const uint16_t middle = static_cast<uint16_t>(value >> 32u);
        const uint32_t low = static_cast<uint32_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        std::memcpy(dest + 4u, &middle, sizeof(middle));
        dest[6] = static_cast<uint8_t>(value >> 48u);
        break;
      }
      case 6: {
        const uint16_t high = static_cast<uint16_t>(value >> 32u);
        const uint32_t low = static_cast<uint32_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        std::memcpy(dest + 4u, &high, sizeof(high));
        break;
      }
      case 5: {
        const uint32_t low = static_cast<uint32_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        dest[4] = static_cast<uint8_t>(value >> 32u);
        break;
      }
      case 4: {
        const uint32_t low = static_cast<uint32_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        break;
      }
      case 3: {
        const uint16_t low = static_cast<uint16_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        dest[2] = static_cast<uint8_t>(value >> 16u);
        break;
      }
      case 2: {
        const uint16_t low = static_cast<uint16_t>(value);
        std::memcpy(dest, &low, sizeof(low));
        break;
      }
      case 1:
        dest[0] = static_cast<uint8_t>(value);
        break;
      default:
        assert(false);
    }
  }

  void Reallocate(size_t capacity) {
    assert(record_bytes_ != 0);
    if (capacity >
        (std::numeric_limits<size_t>::max() - sizeof(uint64_t)) /
            record_bytes_) {
      xfatal("Count segment stream allocation overflow\n");
    }
    void *memory = std::realloc(
        data_, capacity * record_bytes_ + sizeof(uint64_t));
    if (memory == nullptr) {
      xfatal("Failed to allocate count segment stream\n");
    }
    data_ = static_cast<uint8_t *>(memory);
    capacity_ = capacity;
  }

  uint8_t *data_{nullptr};
  size_t size_{0};
  size_t capacity_{0};
  uint32_t read_base_{0};
  uint8_t read_id_bytes_{0};
  uint8_t shard_bytes_{0};
  uint8_t record_bytes_{0};
  uint8_t grouped_bytes_{0};
  uint8_t pending_size_{0};
  uint64_t pending_[8]{};
  bool grouped_compact_{false};
};

// Temporary form used by the 80-bit biological-key path.  One direct 16-bit
// partition makes prefix implicit; compacting away the trailing prefix then
// leaves the exact 64-bit suffix and a 48-bit locator/context payload.  The
// packed 14-byte sort record reduces every later radix swap without allocating
// a second key/index array.
struct CountPrefix64Record {
  uint64_t suffix;
  uint8_t payload[6];
  uint16_t prefix;
} __attribute__((packed));

struct CountSuffix64Record {
  uint64_t suffix;
  uint8_t payload[6];
  static const int n_bytes = sizeof(uint64_t);
  int kth_byte(int byte) const {
    return static_cast<int>((suffix >> (byte * 8u)) & 0xFFu);
  }
  bool operator<(const CountSuffix64Record &other) const {
    return suffix < other.suffix;
  }
} __attribute__((packed));

static_assert(sizeof(CountPrefix64Record) == 16,
              "prefix64 count record must stay compact");
static_assert(sizeof(CountSuffix64Record) == 14,
              "suffix64 count record must stay compact");

inline void StoreCountPayload48(uint8_t *dest, uint64_t payload) {
  assert((payload >> 48u) == 0);
  for (unsigned byte = 0; byte < 6; ++byte) {
    dest[byte] = static_cast<uint8_t>(payload >> (byte * 8u));
  }
}

inline uint64_t LoadCountPayload48(const uint8_t *source) {
  uint64_t payload = 0;
  for (unsigned byte = 0; byte < 6; ++byte) {
    payload |= static_cast<uint64_t>(source[byte]) << (byte * 8u);
  }
  return payload;
}

#ifdef __linux__
struct CountRunSpoolBlockHeader {
  uint32_t payload_bytes;
  uint32_t num_records;
  uint32_t producer_thread;
};

uint64_t CountTargetShardBytes(int num_threads) {
  const uint64_t workers =
      static_cast<uint64_t>(std::max(1, num_threads));
  const uint64_t llc_bytes =
      GetNumaTopology().total_last_level_cache_bytes();
  if (llc_bytes == 0) return uint64_t{8} << 20u;
  const uint64_t cache_share = std::max<uint64_t>(1, llc_bytes / workers);
  return std::max<uint64_t>(uint64_t{1} << 20u, cache_share * 4u);
}

unsigned SelectCountRunSpoolLowBits(
    unsigned shard_bits, int num_threads, unsigned read_id_spare_bits,
    uint64_t estimated_occurrences, unsigned edge_length,
    unsigned minimizer_window, double host_mem) {
  const unsigned max_encoded_bits = std::min(
      shard_bits, 10u + std::min(22u, read_id_spare_bits));
  const char *override_value =
      std::getenv("MEGAHIT_EXTERNAL_COUNT_LOW_BITS");
  if (override_value != nullptr) {
    char *end = nullptr;
    const unsigned long parsed =
        std::strtoul(override_value, &end, 10);
    if (end != override_value && *end == '\0' &&
        parsed <= max_encoded_bits) {
      return static_cast<unsigned>(parsed);
    }
  }

  // A coarse file is also one scheduling epoch.  Give each worker several
  // independent shards per epoch so a heavy tail cannot strand most cores,
  // while keeping the number of resident run references proportional to the
  // machine's requested parallelism rather than to a particular data set.
  const uint64_t target_tasks =
      static_cast<uint64_t>(std::max(1, num_threads)) * 4u;
  unsigned concurrency_bits = 0;
  while ((uint64_t{1} << concurrency_bits) < target_tasks &&
         concurrency_bits < max_encoded_bits) {
    ++concurrency_bits;
  }
  unsigned selected = std::min(
      max_encoded_bits, std::max(8u, concurrency_bits));

  // Parser memory is the mapped run stream plus one compact 8-byte reference
  // per run.  Increase the number of shards in an epoch only while that exact
  // dimensional estimate fits a cache-derived per-worker budget.  This keeps
  // more independent sort tasks available without keying the decision to a
  // named machine or input data set.
  const uint64_t expected_run =
      std::max<uint64_t>(1u, minimizer_window / 2u);
  const uint64_t sequence_bases = edge_length + expected_run - 1u;
  const uint64_t record_bytes =
      8u + DivCeiling(sequence_bases, SeqPackage::kBasesPerWord) *
               sizeof(uint32_t);
  const long double estimated_runs =
      static_cast<long double>(estimated_occurrences) / expected_run;
  const long double parser_bytes =
      estimated_runs * (record_bytes + sizeof(uint64_t));
  long double parser_budget =
      static_cast<long double>(CountTargetShardBytes(num_threads)) *
      std::max(1, num_threads) * 8.0L;
  parser_budget = std::max<long double>(uint64_t{512} << 20u,
                                        parser_budget);
  if (host_mem > 0) {
    parser_budget = std::min<long double>(
        parser_budget,
        std::max<long double>(uint64_t{256} << 20u, host_mem / 4.0L));
  }
  while (selected < max_encoded_bits) {
    const unsigned candidate = selected + 1u;
    const uint64_t containers =
        uint64_t{1} << (shard_bits - candidate);
    if (parser_bytes / containers > parser_budget) break;
    selected = candidate;
  }
  return selected;
}

// File-backed, coarse-sharded run stream.  Producers retain small private
// buffers, reserve disjoint file offsets atomically, and issue only contiguous
// writes.  A record is self-contained: it carries the exact read slice needed
// to regenerate every (k+1)-mer in the minimizer run.  Consequently neither
// the packed read package nor the complete descriptor stream is resident
// during sort/reduce.
class CountRunSpool {
 public:
  struct Buffer {
    std::vector<uint8_t> bytes;
    uint32_t records{0};
  };

  CountRunSpool(const std::string &prefix, unsigned shard_bits,
                unsigned low_bits, unsigned read_shard_bits,
                int num_threads, uint64_t host_mem)
      : prefix_(prefix),
        shard_bits_(shard_bits),
        low_bits_(low_bits),
        read_shard_bits_(read_shard_bits),
        num_containers_(1u << (shard_bits - low_bits_)),
        num_threads_(std::max(1, num_threads)),
        fds_(num_containers_, -1),
        offsets_(new std::atomic<uint64_t>[num_containers_]),
        buffers_(static_cast<size_t>(num_threads_) * num_containers_) {
    assert(low_bits_ <= shard_bits_ && low_bits_ <= 10u + read_shard_bits_);
    const uint64_t aggregate_buffer_budget = std::max<uint64_t>(
        uint64_t{64} << 20u,
        std::min<uint64_t>(uint64_t{1} << 30u,
                           std::max<uint64_t>(1u, host_mem) / 256u));
    const uint64_t buffer_count =
        static_cast<uint64_t>(num_threads_) * num_containers_;
    target_buffer_bytes_ = static_cast<size_t>(std::max<uint64_t>(
        uint64_t{4} << 10u,
        std::min<uint64_t>(uint64_t{256} << 10u,
                           aggregate_buffer_budget /
                               std::max<uint64_t>(1u, buffer_count))));
    target_buffer_bytes_ =
        std::max(target_buffer_bytes_, sizeof(CountRunSpoolBlockHeader) + 128u);

    for (unsigned container = 0; container < num_containers_; ++container) {
      offsets_[container].store(0, std::memory_order_relaxed);
      const std::string path = Path(container);
      fds_[container] =
          open(path.c_str(), O_CREAT | O_TRUNC | O_RDWR | O_CLOEXEC, 0664);
      if (fds_[container] < 0) {
        failed_.store(true, std::memory_order_relaxed);
        break;
      }
    }
  }

  ~CountRunSpool() {
    CloseWriters();
    for (unsigned container = 0; container < num_containers_; ++container) {
      std::remove(Path(container).c_str());
    }
  }

  CountRunSpool(const CountRunSpool &) = delete;
  CountRunSpool &operator=(const CountRunSpool &) = delete;

  bool ok() const { return !failed_.load(std::memory_order_relaxed); }
  unsigned low_bits() const { return low_bits_; }
  unsigned read_shard_bits() const { return read_shard_bits_; }
  unsigned num_containers() const { return num_containers_; }
  size_t target_buffer_bytes() const { return target_buffer_bytes_; }
  std::string Path(unsigned container) const {
    return prefix_ + ".count_runs." + std::to_string(container);
  }
  uint64_t Size(unsigned container) const {
    return offsets_[container].load(std::memory_order_relaxed);
  }

  void Append(int thread_id, unsigned shard, const uint8_t *record,
              size_t record_bytes) {
    if (!ok()) return;
    assert(thread_id >= 0 && thread_id < num_threads_);
    const unsigned container = shard >> low_bits_;
    assert(container < num_containers_);
    Buffer &buffer =
        buffers_[static_cast<size_t>(thread_id) * num_containers_ + container];
    if (buffer.bytes.empty()) {
      buffer.bytes.reserve(target_buffer_bytes_);
      buffer.bytes.resize(sizeof(CountRunSpoolBlockHeader));
    }
    if (buffer.records != 0 &&
        buffer.bytes.size() + record_bytes > target_buffer_bytes_) {
      Flush(&buffer, container, thread_id);
    }
    if (buffer.bytes.size() + record_bytes > buffer.bytes.capacity()) {
      buffer.bytes.reserve(std::max(buffer.bytes.size() + record_bytes,
                                    target_buffer_bytes_));
    }
    buffer.bytes.insert(buffer.bytes.end(), record, record + record_bytes);
    ++buffer.records;
  }

  void FlushThread(int thread_id) {
    for (unsigned container = 0; container < num_containers_; ++container) {
      Buffer &buffer = buffers_[static_cast<size_t>(thread_id) *
                                    num_containers_ +
                                container];
      if (buffer.records != 0) Flush(&buffer, container, thread_id);
      std::vector<uint8_t>().swap(buffer.bytes);
    }
  }

  bool FinalizeWriters() {
    CloseWriters();
    return ok();
  }

  void RemoveContainer(unsigned container) {
    std::remove(Path(container).c_str());
  }

 private:
  void Flush(Buffer *buffer, unsigned container, int producer_thread) {
    if (buffer->records == 0 || !ok()) return;
    CountRunSpoolBlockHeader header{
        static_cast<uint32_t>(buffer->bytes.size() -
                              sizeof(CountRunSpoolBlockHeader)),
        buffer->records,
        static_cast<uint32_t>(producer_thread)};
    std::memcpy(buffer->bytes.data(), &header, sizeof(header));
    const uint64_t offset = offsets_[container].fetch_add(
        buffer->bytes.size(), std::memory_order_relaxed);
    size_t written = 0;
    while (written < buffer->bytes.size()) {
      const ssize_t result =
          pwrite(fds_[container], buffer->bytes.data() + written,
                 buffer->bytes.size() - written,
                 static_cast<off_t>(offset + written));
      if (result <= 0) {
        failed_.store(true, std::memory_order_relaxed);
        break;
      }
      written += static_cast<size_t>(result);
    }
    buffer->bytes.resize(sizeof(CountRunSpoolBlockHeader));
    buffer->records = 0;
  }

  void CloseWriters() {
    for (unsigned container = 0; container < num_containers_; ++container) {
      if (fds_[container] >= 0) {
        const uint64_t bytes =
            offsets_[container].load(std::memory_order_relaxed);
        if (ftruncate(fds_[container], static_cast<off_t>(bytes)) != 0) {
          failed_.store(true, std::memory_order_relaxed);
        }
        close(fds_[container]);
        fds_[container] = -1;
      }
    }
  }

  std::string prefix_;
  unsigned shard_bits_;
  unsigned low_bits_;
  unsigned read_shard_bits_;
  unsigned num_containers_;
  int num_threads_;
  std::vector<int> fds_;
  std::unique_ptr<std::atomic<uint64_t>[]> offsets_;
  std::vector<Buffer> buffers_;
  size_t target_buffer_bytes_{0};
  std::atomic<bool> failed_{false};
};

// Removing a multi-GiB spool file can synchronously walk and release a large
// extent/page-cache tree.  Serialize those metadata operations on one helper
// and overlap them with the next container's radix work; this bounds IO
// pressure while taking deletion latency off the critical path.
class AsyncCountFileRemover {
 public:
  AsyncCountFileRemover() : worker_([this]() { Run(); }) {}
  ~AsyncCountFileRemover() { Finish(); }
  AsyncCountFileRemover(const AsyncCountFileRemover &) = delete;
  AsyncCountFileRemover &operator=(const AsyncCountFileRemover &) = delete;

  void Enqueue(std::string path) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      paths_.push_back(std::move(path));
    }
    ready_.notify_one();
  }

  void Finish() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (closing_) return;
      closing_ = true;
    }
    ready_.notify_one();
    if (worker_.joinable()) worker_.join();
  }

 private:
  void Run() {
    for (;;) {
      std::string path;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        ready_.wait(lock, [&]() { return closing_ || !paths_.empty(); });
        if (paths_.empty()) {
          if (closing_) return;
          continue;
        }
        path = std::move(paths_.front());
        paths_.pop_front();
      }
      std::remove(path.c_str());
    }
  }

  std::thread worker_;
  std::mutex mutex_;
  std::condition_variable ready_;
  std::deque<std::string> paths_;
  bool closing_{false};
};

struct CountExternalRunRef {
  uint64_t offset;
};

inline uint32_t LoadCountRunEncodedRead(const uint8_t *record) {
  uint32_t encoded_read;
  std::memcpy(&encoded_read, record + 1u, sizeof(encoded_read));
  return encoded_read;
}

inline unsigned DecodeCountRunLowShard(const uint8_t *record,
                                       unsigned read_shard_bits) {
  const uint32_t encoded_read = LoadCountRunEncodedRead(record);
  const unsigned read_shard_mask =
      read_shard_bits == 0 ? 0u : (1u << read_shard_bits) - 1u;
  return record[0] |
         (static_cast<unsigned>(record[7] >> 6u) << 8u) |
         ((encoded_read & read_shard_mask) << 10u);
}

inline uint32_t DecodeCountRunReadId(const uint8_t *record,
                                     unsigned read_shard_bits) {
  return LoadCountRunEncodedRead(record) >> read_shard_bits;
}

struct CountExternalWorkspace {
  ~CountExternalWorkspace() {
    std::free(record_workspace);
    std::free(suffix_workspace);
  }
  CountExternalWorkspace() = default;
  CountExternalWorkspace(const CountExternalWorkspace &) = delete;
  CountExternalWorkspace &operator=(const CountExternalWorkspace &) = delete;
  CountExternalWorkspace(CountExternalWorkspace &&other) noexcept
      : record_workspace(other.record_workspace),
        record_workspace_bytes(other.record_workspace_bytes),
        suffix_workspace(other.suffix_workspace),
        suffix_workspace_bytes(other.suffix_workspace_bytes),
        prefix_counts(std::move(other.prefix_counts)),
        prefix_cursors(std::move(other.prefix_cursors)),
        touched_prefixes(std::move(other.touched_prefixes)),
        prefix_groups(std::move(other.prefix_groups)) {
    other.record_workspace = nullptr;
    other.record_workspace_bytes = 0;
    other.suffix_workspace = nullptr;
    other.suffix_workspace_bytes = 0;
  }
  CountExternalWorkspace &operator=(CountExternalWorkspace &&other) noexcept {
    if (this != &other) {
      std::free(record_workspace);
      std::free(suffix_workspace);
      record_workspace = other.record_workspace;
      record_workspace_bytes = other.record_workspace_bytes;
      suffix_workspace = other.suffix_workspace;
      suffix_workspace_bytes = other.suffix_workspace_bytes;
      prefix_counts = std::move(other.prefix_counts);
      prefix_cursors = std::move(other.prefix_cursors);
      touched_prefixes = std::move(other.touched_prefixes);
      prefix_groups = std::move(other.prefix_groups);
      other.record_workspace = nullptr;
      other.record_workspace_bytes = 0;
      other.suffix_workspace = nullptr;
      other.suffix_workspace_bytes = 0;
    }
    return *this;
  }

  uint32_t *record_workspace{nullptr};
  size_t record_workspace_bytes{0};
  void *suffix_workspace{nullptr};
  size_t suffix_workspace_bytes{0};
  std::vector<uint64_t> prefix_counts;
  std::vector<uint64_t> prefix_cursors;
  std::vector<uint16_t> touched_prefixes;
  std::vector<std::pair<uint16_t, std::pair<uint64_t, uint64_t>>>
      prefix_groups;
};
#endif

unsigned CeilLog2Count(uint64_t value, unsigned max_bits) {
  unsigned bits = 0;
  uint64_t rounded = 1;
  while (rounded < value && bits < max_bits) {
    rounded <<= 1u;
    ++bits;
  }
  return bits;
}

unsigned SelectCountShardBits(uint64_t estimated_occurrences,
                              size_t record_bytes, int num_threads,
                              double host_mem,
                              unsigned max_shard_bits =
                                  kCountMaxShardBits) {
  const uint64_t workers = static_cast<uint64_t>(std::max(1, num_threads));

  // An MSD shard should be large enough to amortize allocation/output setup,
  // yet small enough that each worker's active partition tree is comparable to
  // its share of the last-level cache.  Four cache shares allows the top
  // partitions to stream while the increasingly small suffix partitions stay
  // resident, and leaves enough independent tasks to absorb mapped-record
  // latency variance.  If cache topology is unavailable, use a conservative
  // portable workspace target; this is a fallback, not a machine signature.
  const uint64_t target_shard_bytes = CountTargetShardBytes(num_threads);
  const long double total_record_bytes =
      static_cast<long double>(estimated_occurrences) * record_bytes;
  uint64_t cache_shards = static_cast<uint64_t>(
      std::ceil(total_record_bytes / target_shard_bytes));
  const uint64_t balance_shards = workers * 32u;
  uint64_t desired_shards =
      std::max<uint64_t>(256u, std::max(cache_shards, balance_shards));

  unsigned bits =
      std::max(1u, CeilLog2Count(desired_shards, max_shard_bits));
  // Power-of-two rounding should not double scheduling/output overhead merely
  // because the cache-derived target was missed by a few percent.  Retain the
  // lower power when its average shard is within one third of the target.
  if (bits > 1) {
    const uint64_t lower_shards = uint64_t{1} << (bits - 1u);
    const long double lower_average = total_record_bytes / lower_shards;
    if (lower_average <=
        static_cast<long double>(target_shard_bytes) * (4.0L / 3.0L) &&
        lower_shards >= balance_shards) {
      --bits;
    }
  }
  bits = std::min(bits, max_shard_bits);

  // Per-producer segment and occurrence matrices are the price of retaining a
  // single read scan.  On very large teams, reduce the shard count until this
  // metadata is a small fraction of the user-provided memory budget.  The
  // decision depends on dimensions and budget, not on any named CPU model.
  if (host_mem > 0) {
    const long double metadata_budget =
        std::max<long double>(uint64_t{64} << 20u, host_mem / 256.0);
    while (bits > 1) {
      const long double matrix_bytes =
          static_cast<long double>(workers) * (uint64_t{1} << bits) *
          sizeof(uint64_t) * 3u;
      if (matrix_bytes <= metadata_budget) {
        break;
      }
      --bits;
    }
  }

  const char *override_value = std::getenv("MEGAHIT_COUNT_SHARD_BITS");
  if (override_value != nullptr) {
    char *end = nullptr;
    const unsigned long parsed = std::strtoul(override_value, &end, 10);
    if (end != override_value && *end == '\0' && parsed >= 1 &&
        parsed <= max_shard_bits) {
      bits = static_cast<unsigned>(parsed);
    }
  }
  return bits;
}

inline uint64_t MixCountMinimizer(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30u)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27u)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31u);
}

inline uint16_t ExtractCountEdgeWindow8(const uint32_t *edge,
                                        unsigned base_offset) {
  const unsigned word = base_offset / kCharsPerEdgeWord;
  const unsigned shift =
      (base_offset % kCharsPerEdgeWord) * kBitsPerEdgeChar;
  const uint64_t window =
      (static_cast<uint64_t>(edge[word]) << kBitsPerEdgeWord) |
      edge[word + 1u];
  return static_cast<uint16_t>((window << shift) >> 48u);
}

inline uint16_t ReverseComplementCountWindow8(uint16_t bases) {
  uint16_t value = bases;
  value = static_cast<uint16_t>(((value & 0x3333u) << 2u) |
                                ((value & 0xCCCCu) >> 2u));
  value = static_cast<uint16_t>(((value & 0x0F0Fu) << 4u) |
                                ((value & 0xF0F0u) >> 4u));
  value = static_cast<uint16_t>((value << 8u) | (value >> 8u));
  return static_cast<uint16_t>(value ^ 0xFFFFu);
}

}  // namespace

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeFullBaseOffset(int64_t read_id, unsigned offset,
                                    unsigned strand, const SeqPackage &p) {
  return ((p.GetSeqView(read_id).full_offset_in_pkg() + offset) << 1) | strand;
}

struct DecodedReadLocator {
  DecodedReadLocator(SeqPackage::SeqView seq_view, unsigned offset,
                     unsigned strand)
      : seq_view(seq_view), offset(offset), strand(strand) {}

  SeqPackage::SeqView seq_view;
  unsigned offset;
  unsigned strand;
};

template <bool PackedOffsets>
struct ReadLocatorDecoder;

template <>
struct ReadLocatorDecoder<false> {
  static inline DecodedReadLocator Decode(uint64_t locator, unsigned,
                                          uint64_t,
                                          const SeqPackage &seq_pkg) {
    auto seq_view = seq_pkg.GetSeqViewByOffset(locator >> 1u);
    return DecodedReadLocator(
        seq_view,
        static_cast<unsigned>((locator >> 1u) -
                              seq_view.full_offset_in_pkg()),
        locator & 1u);
  }
};

template <>
struct ReadLocatorDecoder<true> {
  static inline DecodedReadLocator Decode(uint64_t locator,
                                          unsigned read_locator_shift,
                                          uint64_t read_offset_mask,
                                          const SeqPackage &seq_pkg) {
    return DecodedReadLocator(
        seq_pkg.GetSeqView(locator >> read_locator_shift),
        static_cast<unsigned>((locator >> 1u) & read_offset_mask),
        locator & 1u);
  }
};

template <bool MetadataInKey>
inline bool IsDifferentEdges(uint32_t *item1, uint32_t *item2, int num_words,
                             int64_t spacing, uint32_t metadata_mask) {
  for (int i = num_words - 1; i >= 0; --i) {
    uint32_t lhs = *(item1 + i * spacing);
    uint32_t rhs = *(item2 + i * spacing);
    if (MetadataInKey && i == num_words - 1) {
      lhs &= ~metadata_mask;
      rhs &= ~metadata_mask;
    }
    if (lhs != rhs) {
      return true;
    }
  }

  return false;
}

/**
 * @brief pack an edge and its multiplicity to word-aligned spaces
 */
void KmerCounter::PackEdge(uint32_t *dest, uint32_t *item, int64_t counting) {
  for (int i = 0; i < words_per_edge_ && i < words_per_substr_; ++i) {
    dest[i] = *(item + i);
  }

  int chars_in_last_word = (opt_.k + 1) % kCharsPerEdgeWord;
  int which_word = (opt_.k + 1) / kCharsPerEdgeWord;

  if (chars_in_last_word > 0) {
    dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
  } else {
    dest[which_word] = 0;
  }

  while (++which_word < words_per_edge_) {
    dest[which_word] = 0;
  }

  dest[words_per_edge_ - 1] |= std::min(int64_t(kMaxMul), counting);
}

// function pass to BaseSequenceSortingEngine

int64_t KmerCounter::Lv0EncodeDiffBase(int64_t read_id) {
  return EncodeReadOffset(read_id, 0, 0);
}

int64_t KmerCounter::EncodeReadOffset(int64_t read_id, unsigned offset,
                                      unsigned strand) const {
  if (packed_read_offsets_) {
    assert(strand <= 1);
    assert(offset <= read_offset_mask_);
    return static_cast<int64_t>(
        (static_cast<uint64_t>(read_id) << read_locator_shift_) |
        (static_cast<uint64_t>(offset) << 1u) | strand);
  }
  return EncodeFullBaseOffset(read_id, offset, strand, seq_pkg_);
}

void KmerCounter::ConfigurePackedReadOffsets(uint64_t max_read_len,
                                             uint64_t num_reads) {
  packed_read_offsets_ = false;
  read_locator_shift_ = 0;
  read_offset_mask_ = 0;

  if (std::getenv("MEGAHIT_DISABLE_PACKED_READ_OFFSETS") != nullptr) {
    xinfo("Packed read offsets disabled by environment\n");
    return;
  }

  const uint64_t max_read_offset = max_read_len == 0 ? 0 : max_read_len - 1;
  unsigned offset_bits = 0;
  for (uint64_t value = max_read_offset; value != 0; value >>= 1u) {
    ++offset_bits;
  }
  read_locator_shift_ = offset_bits + 1;  // offset plus the strand bit
  read_offset_mask_ = offset_bits == 0 ? 0 : (uint64_t{1} << offset_bits) - 1;

  // Six low bits are appended to every locator in the Lv2 item for the
  // previous/next bases.  Keeping the locator in the remaining 58 bits also
  // guarantees it is a non-negative int64_t for differential encoding.
  const uint64_t max_stored_locator =
      std::numeric_limits<uint64_t>::max() >> 6u;
  const uint64_t max_read_id = num_reads == 0 ? 0 : num_reads - 1;
  if (read_locator_shift_ >= 64 ||
      max_read_id > (max_stored_locator >> read_locator_shift_)) {
    read_locator_shift_ = 0;
    read_offset_mask_ = 0;
    xinfo("Packed read offsets unavailable: read IDs exceed locator width\n");
    return;
  }

  const uint64_t max_locator =
      (max_read_id << read_locator_shift_) | (max_read_offset << 1u) | 1u;
  if (max_locator > max_stored_locator) {
    read_locator_shift_ = 0;
    read_offset_mask_ = 0;
    xinfo("Packed read offsets unavailable: offsets exceed locator width\n");
    return;
  }

  // Since every lower field is strictly smaller than 1 << shift, locators
  // increase across read IDs.  Within one read the base offset contributes 2
  // per step and strand contributes at most +/-1, so the differential stream
  // remains strictly increasing for every bucket scan.
  assert(((max_read_offset << 1u) | 1u) <
         (uint64_t{1} << read_locator_shift_));
  assert(max_locator <= (std::numeric_limits<uint64_t>::max() >> 6u));
  packed_read_offsets_ = true;
  xinfo("Using packed read offsets: {} locator bits for reads of length <= {}\n",
        read_locator_shift_, max_read_len);
}

// Record layout (packed locators, words per record, compact single-aux-word
// eligibility) depends only on the maximum read length and the read count,
// both of which the library metadata provides before any sequence is loaded.
// Deciding the layout up front lets the direct-item buffer be allocated and
// prefaulted while the input is still being read.
void KmerCounter::ConfigureRecordLayout(uint64_t max_read_len,
                                        uint64_t num_reads) {
  ConfigurePackedReadOffsets(max_read_len, num_reads);

  words_per_substr_ =
      DivCeiling((opt_.k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
  words_per_edge_ = DivCeiling((opt_.k + 1) * kBitsPerEdgeChar + kBitsPerMul,
                               kBitsPerEdgeWord);
  xinfo("{} words per substring, {} words per edge\n", words_per_substr_,
        words_per_edge_);

  // When the final k-mer word has at least six unused low bits, keep the
  // previous/next context there and drop the second 32-bit auxiliary word.
  // Context remains the least-significant sort key, so equal biological edges
  // stay contiguous; postprocessing masks those bits when finding group ends.
  const unsigned edge_bits = (opt_.k + 1) * kBitsPerEdgeChar;
  const unsigned unused_edge_bits =
      words_per_substr_ * kBitsPerEdgeWord - edge_bits;
  uint64_t max_packed_locator = std::numeric_limits<uint64_t>::max();
  if (packed_read_offsets_) {
    const uint64_t max_read_id = num_reads == 0 ? 0 : num_reads - 1;
    const uint64_t max_read_offset = max_read_len == 0 ? 0 : max_read_len - 1;
    max_packed_locator = (max_read_id << read_locator_shift_) |
                         (max_read_offset << 1u) | 1u;
  }
  // Keep the low locator word in the single auxiliary word.  Any locator bits
  // above bit 31 fit beside the six context bits in the otherwise unused tail
  // of the packed (k+1)-mer.  This extends the compact 12-byte k=21 record
  // beyond the old ~8M-read cliff without changing the biological sort key.
  // The embedded metadata is masked while grouping equal edges and stripped
  // by PackEdge before output.
  const unsigned locator_high_capacity =
      unused_edge_bits >= 6 ? unused_edge_bits - 6 : 0;
  const unsigned compact_locator_bits = 32 + locator_high_capacity;
  const bool locator_fits_compact =
      compact_locator_bits >= 64 ||
      (compact_locator_bits > 0 &&
       (max_packed_locator >> compact_locator_bits) == 0);
  compact_items_ =
      packed_read_offsets_ && unused_edge_bits >= 6 && locator_fits_compact &&
      std::getenv("MEGAHIT_DISABLE_COMPACT_COUNT_ITEMS") == nullptr &&
      // Retain the old A/B switch name for benchmark compatibility.
      std::getenv("MEGAHIT_DISABLE_COMPACT_DIRECT_COUNT_ITEMS") == nullptr;
  compact_locator_high_mask_ = 0;
  compact_key_mask_ = 0;
  if (compact_items_) {
    unsigned locator_high_bits = 0;
    for (uint64_t high = max_packed_locator >> 32u; high != 0;
         high >>= 1u) {
      ++locator_high_bits;
    }
    assert(locator_high_bits <= locator_high_capacity);
    compact_locator_high_mask_ =
        locator_high_bits == 0
            ? 0
            : ((uint32_t{1} << locator_high_bits) - uint32_t{1}) << 6u;
    compact_key_mask_ = uint32_t{0x3F} | compact_locator_high_mask_;
    xinfo(
        "Compact count records: {} unused key bits ({} of {} used for "
        "locator-high), 1 aux word\n",
        unused_edge_bits, locator_high_bits, locator_high_capacity);
  }

  // Exact minimizer runs use a six-byte descriptor containing a 32-bit read
  // ID and byte-sized start/length fields.  This is selected from metadata,
  // not from a particular data set; long-read inputs transparently retain the
  // generic exact bucket scanner.
  segmented_count_enabled_ =
      packed_read_offsets_ && compact_items_ && opt_.n_threads > 1 &&
      num_reads <= std::numeric_limits<uint32_t>::max() &&
      max_read_len <= std::numeric_limits<uint8_t>::max() &&
      std::getenv("MEGAHIT_DISABLE_SEGMENTED_COUNT") == nullptr;
}

bool KmerCounter::Lv1SupportsDirectItems() const {
  return std::getenv("MEGAHIT_DISABLE_DIRECT_COUNT_ITEMS") == nullptr;
}

int64_t KmerCounter::Lv1DirectWordsPerItem() const {
  return words_per_substr_ + (compact_items_ ? 1 : 2);
}

int64_t KmerCounter::Lv1DirectAuxWordsPerItem() const {
  return compact_items_ ? 1 : 2;
}

int KmerCounter::Lv2SortIgnoredLowBytes() const {
  const int64_t stored_bits = words_per_substr_ * kBitsPerEdgeWord;
  const int64_t biological_bits = (opt_.k + 1) * kBitsPerEdgeChar;
  assert(stored_bits >= biological_bits);
  // Context and packed-locator metadata occupy otherwise unused low bits.
  // Postprocessing masks those bits and only requires equal biological edges
  // to be contiguous, so complete non-biological bytes need not participate
  // in the radix key at all.
  return stored_bits - biological_bits >= 16 ? 2 : 0;
}

bool KmerCounter::Lv1CanUseCompactCursor(int64_t seq_from,
                                         int64_t seq_to) const {
  if (seq_from >= seq_to) {
    return true;
  }
  const auto last_read = seq_pkg_.GetSeqView(seq_to - 1);
  const unsigned edge_length = opt_.k + 1;
  const unsigned last_offset =
      last_read.length() > edge_length ? last_read.length() - edge_length : 0;
  const uint64_t begin = static_cast<uint64_t>(
      EncodeReadOffset(seq_from, 0, 0));
  const uint64_t end = static_cast<uint64_t>(
      EncodeReadOffset(seq_to - 1, last_offset, 1));
  return end >= begin &&
         end - begin <= std::numeric_limits<uint32_t>::max();
}

int64_t KmerCounter::Lv1DirectMemoryLimit() const {
  if (std::getenv("MEGAHIT_FORCE_DIRECT_COUNT_ITEMS") != nullptr) {
    return std::numeric_limits<int64_t>::max();
  }

  // Direct materialization is a transient acceleration structure, not the
  // primary data set.  Keep it inside the same generic bounded-workspace
  // envelope as compact counting; otherwise a modest change in input size
  // causes a 60-120 GiB RSS cliff before falling back to compact locators.
  const bool write_combined =
      std::getenv("MEGAHIT_DISABLE_LV1_WRITE_COMBINE") == nullptr;
  const double budget_fraction = write_combined ? (1.0 / 16.0) : (1.0 / 32.0);
  return static_cast<int64_t>(opt_.host_mem * budget_fraction);
}

unsigned KmerCounter::Lv1BytesPerOffset() const {
  return packed_read_offsets_ &&
                 std::getenv("MEGAHIT_DISABLE_PACKED_LV1_OFFSETS") == nullptr
             ? 3u
             : kDefaultLv1BytesPerItem;
}

int64_t KmerCounter::Lv1AutoWorkspaceLimit() const {
  // Bound the transient locator+radix workspace by a fraction of the memory
  // budget, independent of the data set and thread count.  Three-byte
  // locators then spend the saved byte on fewer complete input traversals
  // while keeping total RSS near later SDBG/assembly stages.
  return static_cast<int64_t>(opt_.host_mem / 16.0);
}

KmerCounter::MemoryStat KmerCounter::Initialize() {
  bool is_reverse = true;

  SequenceLibCollection seq_collection(opt_.read_lib_file);
  const auto size_info = seq_collection.GetSizeInfo();
  int64_t num_bases = size_info.num_bases;
  int64_t num_reads = size_info.num_reads;

  // The record layout is fully determined by library metadata, so lock it in
  // now, then overlap the direct-array page faults with the sequence load and
  // the Lv0 bucket scan.  Both the item count and the data size below are
  // estimates; AdjustMemory re-validates against the real totals and simply
  // drops the speculative buffer if they disagree.
  ConfigureRecordLayout(size_info.max_read_len, num_reads);
  // Reads shorter than k contribute zero items but subtract from this
  // difference, so pad the estimate slightly; an under-sized speculation is
  // rejected at adoption time, wasting the prefault work.
  int64_t estimated_items =
      num_bases - num_reads * static_cast<int64_t>(opt_.k);
  estimated_items += estimated_items / 1024 + (int64_t{1} << 20u);
  const int64_t estimated_memory_for_data =
      num_bases / 4 + num_reads * 16 + (int64_t{256} << 20u);
  if (!segmented_count_enabled_) {
    PreallocateDirectItems(estimated_items, estimated_memory_for_data);
  }

  if (segmented_count_enabled_ &&
      std::getenv("MEGAHIT_DISABLE_EXTERNAL_SEGMENTED_COUNT") == nullptr) {
    std::unique_ptr<MappedCountReads> mapped(new MappedCountReads());
    if (mapped->Open(opt_.read_lib_file + ".bin",
                     static_cast<uint64_t>(size_info.num_reads),
                     static_cast<uint64_t>(size_info.num_bases),
                     size_info.max_read_len)) {
      mapped_count_reads_ = std::move(mapped);
      xinfo("Using indexed mapped binary reads: {} chunks, {} bytes; "
            "no resident 2-bit package\n",
            mapped_count_reads_->index().chunks.size(),
            mapped_count_reads_->mapped_bytes());
    }
  }

  if (!mapped_count_reads_) {
    seq_pkg_.ReserveSequences(num_reads);
    seq_pkg_.ReserveBases(num_bases);
    seq_collection.Read(&seq_pkg_, is_reverse);

    num_reads = seq_pkg_.seq_count();
    xinfo("{} reads, {} max read length\n", num_reads,
          seq_pkg_.max_length());
    if (static_cast<uint64_t>(num_reads) !=
            static_cast<uint64_t>(size_info.num_reads) ||
        seq_pkg_.max_length() != size_info.max_read_len) {
      // Metadata disagreed with the actual library; recompute the layout from
      // the loaded package (the speculative buffer is discarded later).
      xwarn("Library metadata mismatch; reconfiguring record layout\n");
      ConfigureRecordLayout(seq_pkg_.max_length(), num_reads);
    }
    // pos_to_id is only used to recover a read ID from a full base offset.  The
    // packed locator already contains that ID, so avoid building and retaining
    // the lookup table on the fast path.
    if (!packed_read_offsets_) {
      seq_pkg_.BuildIndex();
    } else {
      const unsigned gap_bits = seq_pkg_.CompactReadOnlyLengthIndex();
      if (gap_bits == 1) {
        xinfo("Compacted read starts to a sparse gap index ({} exceptions)\n",
              seq_pkg_.compact_length_exceptions());
      } else if (gap_bits != 0) {
        xinfo("Compacted read starts to a {}-bit cumulative-gap index\n",
              gap_bits);
      }
    }
  } else {
    num_reads = static_cast<int64_t>(mapped_count_reads_->num_reads());
    xinfo("{} reads, {} max read length\n", num_reads,
          mapped_count_reads_->max_length());
  }

  // --- malloc read first_in / last_out ---
  // Parallel 0xFF fill with NUMA-local first touch: the serial vector
  // constructor costs ~0.5 s for tens of millions of reads and would place
  // every page on the master thread's node.
  const size_t n_seq = static_cast<size_t>(num_reads);
  const unsigned input_max_length =
      mapped_count_reads_ ? mapped_count_reads_->max_length()
                          : seq_pkg_.max_length();
  compact_read_flags_ = input_max_length <= 255u;
  const size_t flag_count = std::max<size_t>(1, n_seq);
  const size_t flag_element_bytes =
      compact_read_flags_ ? sizeof(AtomicWrapper<uint8_t>)
                          : sizeof(AtomicWrapper<uint32_t>);
  const size_t flag_bytes = flag_element_bytes * flag_count;
  if (compact_read_flags_) {
    first_0_out8_.reset(
        static_cast<AtomicWrapper<uint8_t> *>(std::malloc(flag_bytes)));
    last_0_in8_.reset(
        static_cast<AtomicWrapper<uint8_t> *>(std::malloc(flag_bytes)));
    if (!first_0_out8_ || !last_0_in8_) {
      xfatal("Failed to allocate compact read flag arrays\n");
    }
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(n_seq); ++i) {
      first_0_out8_[i].v.store(uint8_t{0xFF}, std::memory_order_relaxed);
      last_0_in8_[i].v.store(uint8_t{0xFF}, std::memory_order_relaxed);
    }
    xinfo("Compact read boundary flags: {} bytes (2 x 8-bit offsets)\n",
          flag_bytes * 2u);
  } else {
    first_0_out32_.reset(
        static_cast<AtomicWrapper<uint32_t> *>(std::malloc(flag_bytes)));
    last_0_in32_.reset(
        static_cast<AtomicWrapper<uint32_t> *>(std::malloc(flag_bytes)));
    if (!first_0_out32_ || !last_0_in32_) {
      xfatal("Failed to allocate read flag arrays\n");
    }
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(n_seq); ++i) {
      first_0_out32_[i].v.store(kSentinelOffset,
                               std::memory_order_relaxed);
      last_0_in32_[i].v.store(kSentinelOffset,
                              std::memory_order_relaxed);
    }
  }

  // --- initialize stat ---
  edge_counter_.SetNumThreads(opt_.n_threads);

  // --- initialize writer ---
  edge_writer_.SetFilePrefix(opt_.output_prefix);
  edge_writer_.SetNumThreads(opt_.n_threads);
  edge_writer_.SetKmerSize(opt_.k);
  edge_writer_.SetNumBuckets(kNumBuckets);
  if (segmented_count_enabled_) {
    // Minimizer ownership is exact for equal biological keys, but its output
    // order is intentionally independent of the biological prefix buckets.
    // Advertise that fact instead of emitting misleading sorted-bucket
    // metadata; SeqToSdbg consumes these sharded files as parallel sequential
    // streams and performs its own exact radix partition.
    edge_writer_.SetUnorderedSharded();
    seq_bucket_histograms_.resize(
        static_cast<size_t>(opt_.n_threads) * kNumBuckets);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < seq_bucket_histograms_.size(); ++i) {
      seq_bucket_histograms_[i] = 0;
    }
  }
  edge_writer_.InitFiles();

  endpoint_candidate_files_.reserve(opt_.n_threads);
  for (int thread = 0; thread < opt_.n_threads; ++thread) {
    endpoint_candidate_files_.emplace_back(new std::ofstream(
        opt_.output_prefix + ".endcand." + std::to_string(thread),
        std::ios::binary | std::ios::out | std::ios::trunc));
    if (!*endpoint_candidate_files_.back()) {
      xfatal("Cannot create endpoint-candidate shard {}\n", thread);
    }
  }

  int64_t memory_for_data = seq_pkg_.size_in_byte() +
                            +n_seq * flag_element_bytes *
                                2  // first_in0 & last_out0
                            + edge_counter_.size_in_byte();  // edge_counting

  return {
      num_reads,
      memory_for_data,
      words_per_substr_ + (compact_items_ ? 1 : 2),
      compact_items_ ? 1 : 2,
  };
}

void KmerCounter::UpdateLast0In(uint32_t read_id, uint32_t offset) {
  if (compact_read_flags_) {
    assert(offset < uint32_t{0xFF});
    auto &last = last_0_in8_[read_id].v;
    uint8_t old_value = last.load(std::memory_order_relaxed);
    const uint8_t value = static_cast<uint8_t>(offset);
    while ((old_value == uint8_t{0xFF} || old_value < value) &&
           !last.compare_exchange_weak(old_value, value,
                                       std::memory_order_relaxed,
                                       std::memory_order_relaxed)) {
    }
    return;
  }
  auto &last = last_0_in32_[read_id].v;
  uint32_t old_value = last.load(std::memory_order_relaxed);
  while ((old_value == kSentinelOffset || old_value < offset) &&
         !last.compare_exchange_weak(old_value, offset,
                                     std::memory_order_relaxed,
                                     std::memory_order_relaxed)) {
  }
}

void KmerCounter::WriteEndpointCandidate(const uint32_t *packed_edge,
                                         bool missing_in,
                                         bool missing_out,
                                         uint8_t verify_in_mask,
                                         uint8_t verify_out_mask,
                                         int thread_id) {
  const unsigned edge_length = opt_.k + 1u;
  // A reverse-complement palindrome has no unique canonical strand: context
  // observations from its two physical orientations collapse into the same
  // count group and can make both side flags appear present even when one
  // oriented node is an endpoint.  Include those rare edges explicitly.  An
  // eight-base end signature rejects all but about 1/65536 ordinary edges
  // before the exact comparison, so this closes the tie case without adding
  // a full reverse-complement operation to the hot path.
  if (missing_in == false && missing_out == false &&
      (edge_length & 1u) == 0 && edge_length >= 8u) {
    const auto window8 = [&](unsigned offset) {
      const unsigned word = offset / kCharsPerEdgeWord;
      const unsigned within = offset % kCharsPerEdgeWord;
      if (within + 8u <= kCharsPerEdgeWord) {
        const unsigned shift =
            (kCharsPerEdgeWord - within - 8u) * kBitsPerEdgeChar;
        return static_cast<uint16_t>(packed_edge[word] >> shift);
      }
      return ExtractCountEdgeWindow8(packed_edge, offset);
    };
    const uint16_t first = window8(0);
    const uint16_t reverse_last = ReverseComplementCountWindow8(
        window8(edge_length - 8u));
    if (first == reverse_last) {
      bool palindrome = true;
      for (unsigned i = 8u; i < edge_length / 2u; ++i) {
        const unsigned lhs =
            (packed_edge[i / kCharsPerEdgeWord] >>
             ((kCharsPerEdgeWord - 1u - i % kCharsPerEdgeWord) *
              kBitsPerEdgeChar)) &
            kEdgeCharMask;
        const unsigned reverse_position = edge_length - 1u - i;
        const unsigned rhs =
            (packed_edge[reverse_position / kCharsPerEdgeWord] >>
             ((kCharsPerEdgeWord - 1u -
               reverse_position % kCharsPerEdgeWord) *
              kBitsPerEdgeChar)) &
            kEdgeCharMask;
        if (lhs + rhs != 3u) {
          palindrome = false;
          break;
        }
      }
      if (palindrome) missing_in = missing_out = true;
    }
  }
  if (!missing_in && !missing_out && verify_in_mask == 0u &&
      verify_out_mask == 0u) {
    return;
  }
  assert(thread_id >= 0 &&
         thread_id < static_cast<int>(endpoint_candidate_files_.size()));
  std::ofstream &output = *endpoint_candidate_files_[thread_id];
  const auto write_record = [&](uint8_t flags) {
    output.write(reinterpret_cast<const char *>(packed_edge),
                 static_cast<std::streamsize>(words_per_edge_ *
                                              sizeof(uint32_t)));
    output.put(static_cast<char>(flags));
    if (!output) {
      xfatal("Failed writing endpoint-candidate shard {}\n", thread_id);
    }
  };

  // Emit proven sides directly.  Verification is unnecessary for that same
  // side, but the opposite side can still carry sub-2T ambiguous support and
  // therefore remains in the compact verification stream.
  if (missing_in || missing_out) {
    write_record(static_cast<uint8_t>(missing_in ? 1u : 0u) |
                 static_cast<uint8_t>(missing_out ? 2u : 0u));
    if (missing_in) verify_in_mask = 0u;
    if (missing_out) verify_out_mask = 0u;
  }

  // Encode the exact context base alongside each ambiguous side.  Usually
  // each mask is a singleton and both sides fit in one 13-byte record.  If a
  // side has several sub-2T contexts, emit another record for each; an absent
  // context may add a redundant sentinel, but the exact node-group filter in
  // seq2sdbg removes it.  This preserves a strict superset while replacing up
  // to eight random edge-index probes per record with one per encoded side.
  while (verify_in_mask != 0u || verify_out_mask != 0u) {
    uint8_t flags = 0u;
    if (verify_in_mask != 0u) {
      const unsigned base =
          static_cast<unsigned>(__builtin_ctz(verify_in_mask));
      verify_in_mask &= static_cast<uint8_t>(verify_in_mask - 1u);
      flags |= static_cast<uint8_t>(4u | (base << 4u));
    }
    if (verify_out_mask != 0u) {
      const unsigned base =
          static_cast<unsigned>(__builtin_ctz(verify_out_mask));
      verify_out_mask &= static_cast<uint8_t>(verify_out_mask - 1u);
      flags |= static_cast<uint8_t>(8u | (base << 6u));
    }
    write_record(flags);
  }
}

void KmerCounter::UpdateFirst0Out(uint32_t read_id, uint32_t offset) {
  if (compact_read_flags_) {
    assert(offset < uint32_t{0xFF});
    auto &first = first_0_out8_[read_id].v;
    uint8_t old_value = first.load(std::memory_order_relaxed);
    const uint8_t value = static_cast<uint8_t>(offset);
    while (old_value > value &&
           !first.compare_exchange_weak(old_value, value,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed)) {
    }
    return;
  }
  auto &first = first_0_out32_[read_id].v;
  uint32_t old_value = first.load(std::memory_order_relaxed);
  while (old_value > offset &&
         !first.compare_exchange_weak(old_value, offset,
                                      std::memory_order_relaxed,
                                      std::memory_order_relaxed)) {
  }
}

uint32_t KmerCounter::LoadFirst0Out(uint32_t read_id) const {
  if (compact_read_flags_) {
    const uint8_t value =
        first_0_out8_[read_id].v.load(std::memory_order_relaxed);
    return value == uint8_t{0xFF} ? kSentinelOffset : value;
  }
  return first_0_out32_[read_id].v.load(std::memory_order_relaxed);
}

uint32_t KmerCounter::LoadLast0In(uint32_t read_id) const {
  if (compact_read_flags_) {
    const uint8_t value =
        last_0_in8_[read_id].v.load(std::memory_order_relaxed);
    return value == uint8_t{0xFF} ? kSentinelOffset : value;
  }
  return last_0_in32_[read_id].v.load(std::memory_order_relaxed);
}

template <unsigned NWords>
void KmerCounter::MaterializeSegment(const CompactCountSegment &segment,
                                     uint32_t *dest) const {
  const uint32_t read_id = segment.read_id;
  const unsigned segment_start = segment.start;
  const unsigned segment_length = segment.length;
  const unsigned edge_length = opt_.k + 1;
  const auto seq_view = seq_pkg_.GetSeqView(read_id);
  assert(segment_start + segment_length + opt_.k <= seq_view.length());

  Kmer<NWords, uint32_t> edge, rev_edge;
  const auto raw = seq_view.raw_address();
  edge.InitFromPtr(raw.first, raw.second + segment_start, edge_length);
  rev_edge = edge;
  rev_edge.ReverseComplement(edge_length);

  uint64_t encoded_offset =
      (static_cast<uint64_t>(read_id) << read_locator_shift_) |
      (static_cast<uint64_t>(segment_start) << 1u);
  unsigned offset = segment_start;

  for (unsigned item_index = 0; item_index < segment_length; ++item_index) {
    const bool reverse = rev_edge.cmp(edge, edge_length) < 0;
    const auto &canonical = reverse ? rev_edge : edge;
    const unsigned prev =
        offset == 0 ? kSentinelValue : seq_view.base_at(offset - 1);
    const unsigned next =
        offset + edge_length < seq_view.length()
            ? seq_view.base_at(offset + edge_length)
            : kSentinelValue;
    const uint64_t locator = encoded_offset | static_cast<unsigned>(reverse);
    uint32_t context;
    if (!reverse) {
      context = (prev << 3u) | next;
    } else {
      const unsigned canonical_prev =
          next == kSentinelValue ? kSentinelValue : 3u - next;
      const unsigned canonical_next =
          prev == kSentinelValue ? kSentinelValue : 3u - prev;
      context = (canonical_prev << 3u) | canonical_next;
    }
    const uint32_t metadata =
        (static_cast<uint32_t>(locator >> 32u) << 6u) | context;
    assert((metadata & ~compact_key_mask_) == 0);
    assert((canonical.data()[NWords - 1] & compact_key_mask_) == 0);
    std::memcpy(dest, canonical.data(),
                (NWords - 1) * sizeof(uint32_t));
    dest[NWords - 1] = canonical.data()[NWords - 1] | metadata;
    dest[NWords] = static_cast<uint32_t>(locator);
    dest += NWords + 1;

    if (item_index + 1 == segment_length) {
      break;
    }
    ++offset;
    encoded_offset += 2u;
    const int c = seq_view.base_at(offset + edge_length - 1);
    edge.ShiftAppend(c, edge_length);
    rev_edge.ShiftPreappend(3 - c, edge_length);
  }
}

template <unsigned NWords>
bool KmerCounter::RunExternalSegmentedCountFor() {
#ifndef __linux__
  return false;
#else
  assert(mapped_count_reads_ != nullptr);
  const int n_threads = opt_.n_threads;
  const uint64_t num_reads = mapped_count_reads_->num_reads();
  const unsigned edge_length = opt_.k + 1u;
  const unsigned minimizer_length =
      std::min<unsigned>(15, std::max<unsigned>(7, edge_length / 2u));
  const unsigned minimizer_window = edge_length - minimizer_length + 1u;
  const uint64_t minimizer_mask =
      minimizer_length == 32
          ? std::numeric_limits<uint64_t>::max()
          : (uint64_t{1} << (2u * minimizer_length)) - 1u;
  const unsigned rc_shift = 2u * (minimizer_length - 1u);
  const uint64_t estimated_occurrences =
      mapped_count_reads_->num_bases() >
              num_reads * static_cast<uint64_t>(opt_.k)
          ? mapped_count_reads_->num_bases() -
                num_reads * static_cast<uint64_t>(opt_.k)
          : 0;
  const unsigned shard_bits = SelectCountShardBits(
      estimated_occurrences, (NWords + 1u) * sizeof(uint32_t), n_threads,
      opt_.host_mem, kCountExternalMaxShardBits);
  const unsigned num_shards = 1u << shard_bits;

  unsigned read_id_bits = 0;
  uint64_t read_id_capacity = 1;
  while (read_id_capacity < num_reads && read_id_bits < 32u) {
    read_id_capacity <<= 1u;
    ++read_id_bits;
  }
  const unsigned read_id_spare_bits = 32u - read_id_bits;
  const unsigned low_bits = SelectCountRunSpoolLowBits(
      shard_bits, n_threads, read_id_spare_bits, estimated_occurrences,
      edge_length, minimizer_window, opt_.host_mem);
  const unsigned read_shard_bits = low_bits > 10u ? low_bits - 10u : 0u;

  xinfo("External exact segmented count: {} minimizer shards, m={}, "
        "window={}; {} shards/container, {} read-ID spare shard bits; "
        "mapped read scan + self-contained run spool\n",
        num_shards, minimizer_length, minimizer_window,
        uint64_t{1} << low_bits, read_shard_bits);

  CountRunSpool spool(opt_.output_prefix, shard_bits, low_bits,
                      read_shard_bits, n_threads,
                      opt_.host_mem > 0
                          ? static_cast<uint64_t>(opt_.host_mem)
                          : uint64_t{1});
  if (!spool.ok()) {
    xfatal("Cannot create external count run spool\n");
  }
  xinfo("External run spool: {} coarse files, {} low shard bits, {}-byte "
        "producer buffers\n",
        spool.num_containers(), spool.low_bits(), spool.target_buffer_bytes());

  const size_t matrix_items = static_cast<size_t>(n_threads) * num_shards;
  std::vector<uint64_t> thread_run_counts(matrix_items, 0);
  std::vector<uint64_t> thread_occurrence_counts(matrix_items, 0);
  DiscardMemoryPages(thread_run_counts.data(),
                     thread_run_counts.size() * sizeof(uint64_t));
  AdviseHugePages(thread_run_counts.data(),
                  thread_run_counts.size() * sizeof(uint64_t));
  DiscardMemoryPages(thread_occurrence_counts.data(),
                     thread_occurrence_counts.size() * sizeof(uint64_t));
  AdviseHugePages(thread_occurrence_counts.data(),
                  thread_occurrence_counts.size() * sizeof(uint64_t));
  std::vector<uint64_t> thread_segments(n_threads, 0);
  std::vector<uint64_t> thread_occurrences(n_threads, 0);
  std::vector<uint64_t> thread_spool_bytes(n_threads, 0);

  const auto &read_index = mapped_count_reads_->index();
  const auto *mapped_words = mapped_count_reads_->words();
  const double scan_begin = omp_get_wtime();

#pragma omp parallel num_threads(n_threads)
  {
    const int tid = omp_get_thread_num();
    uint64_t *run_counts =
        thread_run_counts.data() + static_cast<size_t>(tid) * num_shards;
    uint64_t *occurrence_counts =
        thread_occurrence_counts.data() +
        static_cast<size_t>(tid) * num_shards;
    std::fill(run_counts, run_counts + num_shards, uint64_t{0});
    std::fill(occurrence_counts, occurrence_counts + num_shards, uint64_t{0});
    std::vector<uint64_t> hashes;
    std::vector<unsigned> min_queue;
    std::vector<uint8_t> read_bases;
    hashes.reserve(mapped_count_reads_->max_length());
    min_queue.reserve(mapped_count_reads_->max_length());
    read_bases.reserve(mapped_count_reads_->max_length());
    std::array<uint8_t, 8u + 64u> record{};
    std::array<uint32_t, 16u> packed_record_bases{};
    uint64_t local_segments = 0;
    uint64_t local_occurrences = 0;
    uint64_t local_spool_bytes = 0;

#pragma omp for schedule(dynamic, 1)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(read_index.chunks.size());
         ++chunk_id) {
      const PackedReadChunk &chunk = read_index.chunks[chunk_id];
      const uint32_t *cursor = mapped_words + chunk.word_begin;
      for (uint64_t read_id64 = chunk.read_begin;
           read_id64 < chunk.read_end; ++read_id64) {
        if (cursor >= mapped_words + chunk.word_end) {
          xfatal("Mapped count chunk ended before read {}\n", read_id64);
        }
        mapped_count_reads_->RecordCandidateCheckpoint(
            static_cast<size_t>(chunk_id), read_id64, cursor);
        const uint32_t read_length = *cursor++;
        const size_t packed_words =
            DivCeiling(static_cast<size_t>(read_length),
                       static_cast<size_t>(SeqPackage::kBasesPerWord));
        if (packed_words > static_cast<size_t>(
                               mapped_words + chunk.word_end - cursor)) {
          xfatal("Invalid mapped count record for read {}\n", read_id64);
        }
        if (read_length < edge_length) {
          cursor += packed_words;
          continue;
        }
        assert(read_id64 <= std::numeric_limits<uint32_t>::max());
        const uint32_t read_id = static_cast<uint32_t>(read_id64);
        MappedCountReads::ReverseView read(cursor, read_length);
        read_bases.resize(read_length);
        for (unsigned i = 0; i < read_length; ++i) {
          read_bases[i] = read.base_at(i);
        }
        cursor += packed_words;

        const unsigned num_mmers = read_length - minimizer_length + 1u;
        hashes.resize(num_mmers);
        min_queue.resize(num_mmers);
        uint64_t fwd = 0;
        uint64_t rev = 0;
        for (unsigned i = 0; i < minimizer_length; ++i) {
          const uint64_t c = read_bases[i];
          fwd = (fwd << 2u) | c;
          rev = (rev >> 2u) | ((uint64_t{3} - c) << rc_shift);
        }
        for (unsigned pos = 0; pos < num_mmers; ++pos) {
          hashes[pos] = MixCountMinimizer(
              std::min(fwd, rev) ^
              (uint64_t{0x6a09e667f3bcc909} + minimizer_length));
          if (pos + 1u < num_mmers) {
            const uint64_t c = read_bases[pos + minimizer_length];
            fwd = ((fwd << 2u) | c) & minimizer_mask;
            rev = (rev >> 2u) | ((uint64_t{3} - c) << rc_shift);
          }
        }

        auto emit_run = [&](unsigned start, unsigned length,
                            unsigned shard) {
          assert(length > 0 && length <= 255u && start <= 255u);
          const unsigned sequence_bases = edge_length + length - 1u;
          const unsigned sequence_words =
              DivCeiling(sequence_bases, SeqPackage::kBasesPerWord);
          const unsigned sequence_bytes = sequence_words * sizeof(uint32_t);
          assert(8u + sequence_bytes <= record.size());
          const unsigned low_shard =
              shard & ((1u << spool.low_bits()) - 1u);
          record[0] = static_cast<uint8_t>(low_shard);
          assert(spool.read_shard_bits() == 0 ||
                 read_id <= (std::numeric_limits<uint32_t>::max() >>
                             spool.read_shard_bits()));
          const uint32_t encoded_read =
              (read_id << spool.read_shard_bits()) |
              (low_shard >> 10u);
          std::memcpy(record.data() + 1u, &encoded_read,
                      sizeof(encoded_read));
          record[5] = static_cast<uint8_t>(start);
          record[6] = static_cast<uint8_t>(length);
          const unsigned previous =
              start == 0 ? kSentinelValue : read_bases[start - 1u];
          const unsigned next_position = start + sequence_bases;
          const unsigned next = next_position < read_length
                                    ? read_bases[next_position]
                                    : kSentinelValue;
          // Context uses six bits; retain up to two additional local-shard
          // bits in the otherwise unused high bits without growing the run
          // record or unaligning its native packed-base words.
          record[7] = static_cast<uint8_t>(
              (previous << 3u) | next |
              (((low_shard >> 8u) & 0x3u) << 6u));
          packed_record_bases.fill(0);
          for (unsigned i = 0; i < sequence_bases; ++i) {
            packed_record_bases[i / SeqPackage::kBasesPerWord] |=
                static_cast<uint32_t>(read_bases[start + i])
                << ((SeqPackage::kBasesPerWord - 1u -
                     i % SeqPackage::kBasesPerWord) *
                    2u);
          }
          std::memcpy(record.data() + 8u, packed_record_bases.data(),
                      sequence_bytes);
          const size_t bytes = 8u + sequence_bytes;
          spool.Append(tid, shard, record.data(), bytes);
          ++run_counts[shard];
          occurrence_counts[shard] += length;
          ++local_segments;
          local_occurrences += length;
          local_spool_bytes += bytes;
        };

        unsigned q_begin = 0;
        unsigned q_end = 0;
        unsigned current_shard = 0;
        unsigned run_start = 0;
        unsigned run_length = 0;
        for (unsigned pos = 0; pos < num_mmers; ++pos) {
          while (q_end > q_begin &&
                 hashes[min_queue[q_end - 1u]] > hashes[pos]) {
            --q_end;
          }
          min_queue[q_end++] = pos;
          if (pos + 1u < minimizer_window) continue;
          const unsigned edge_start = pos + 1u - minimizer_window;
          while (q_begin < q_end && min_queue[q_begin] < edge_start) {
            ++q_begin;
          }
          assert(edge_start < num_edges);
          const unsigned shard = static_cast<unsigned>(
              hashes[min_queue[q_begin]]) & (num_shards - 1u);
          if (run_length == 0) {
            current_shard = shard;
            run_start = edge_start;
            run_length = 1;
          } else if (shard == current_shard && run_length < 255u) {
            ++run_length;
          } else {
            emit_run(run_start, run_length, current_shard);
            current_shard = shard;
            run_start = edge_start;
            run_length = 1;
          }
        }
        if (run_length != 0) {
          emit_run(run_start, run_length, current_shard);
        }
      }
      if (cursor != mapped_words + chunk.word_end) {
        xfatal("Mapped count chunk {} has trailing words\n", chunk_id);
      }
      mapped_count_reads_->DiscardChunk(chunk);
    }
    spool.FlushThread(tid);
    thread_segments[tid] = local_segments;
    thread_occurrences[tid] = local_occurrences;
    thread_spool_bytes[tid] = local_spool_bytes;
  }
  if (!spool.FinalizeWriters()) {
    xfatal("Failed while writing external count run spool\n");
  }
  const double scan_wall = omp_get_wtime() - scan_begin;

  const uint64_t total_segments =
      std::accumulate(thread_segments.begin(), thread_segments.end(),
                      uint64_t{0});
  const uint64_t total_occurrences =
      std::accumulate(thread_occurrences.begin(), thread_occurrences.end(),
                      uint64_t{0});
  const uint64_t logical_spool_bytes =
      std::accumulate(thread_spool_bytes.begin(), thread_spool_bytes.end(),
                      uint64_t{0});
  uint64_t physical_spool_bytes = 0;
  for (unsigned container = 0; container < spool.num_containers();
       ++container) {
    physical_spool_bytes += spool.Size(container);
  }
  xinfo("External count scan produced {} self-contained runs for {} exact "
        "records ({.3} records/run), {.3} GiB payload / {.3} GiB blocked "
        "spool; scan+write {.4} s\n",
        total_segments, total_occurrences,
        total_segments == 0
            ? 0.0
            : static_cast<double>(total_occurrences) / total_segments,
        static_cast<double>(logical_spool_bytes) / (uint64_t{1} << 30u),
        static_cast<double>(physical_spool_bytes) / (uint64_t{1} << 30u),
        scan_wall);

  std::vector<uint64_t> shard_runs(num_shards, 0);
  std::vector<uint64_t> shard_occurrences(num_shards, 0);
#pragma omp parallel for schedule(static) num_threads(n_threads)
  for (unsigned shard = 0; shard < num_shards; ++shard) {
    uint64_t runs = 0;
    uint64_t count = 0;
    for (int thread = 0; thread < n_threads; ++thread) {
      runs += thread_run_counts[
          static_cast<size_t>(thread) * num_shards + shard];
      count += thread_occurrence_counts[
          static_cast<size_t>(thread) * num_shards + shard];
    }
    shard_runs[shard] = runs;
    shard_occurrences[shard] = count;
  }
  const uint64_t observed_runs =
      std::accumulate(shard_runs.begin(), shard_runs.end(), uint64_t{0});
  const uint64_t observed_occurrences =
      std::accumulate(shard_occurrences.begin(), shard_occurrences.end(),
                      uint64_t{0});
  if (observed_runs != total_segments) {
    xfatal("External count run mismatch: {} != {}\n", observed_runs,
           total_segments);
  }
  if (observed_occurrences != total_occurrences) {
    xfatal("External count occurrence mismatch: {} != {}\n",
           observed_occurrences, total_occurrences);
  }

  const bool use_suffix64 =
      NWords == 3 && edge_length == 40 &&
      num_reads != 0 &&
      (((num_reads - 1u) << read_locator_shift_) |
       ((static_cast<uint64_t>(mapped_count_reads_->max_length() - 1u)
         << 1u) |
        1u)) < (uint64_t{1} << 42u) &&
      std::getenv("MEGAHIT_DISABLE_SUFFIX64_COUNT") == nullptr;
  if (use_suffix64) {
    xinfo("External count radix layout: implicit 16-bit biological prefix + "
          "14-byte suffix64/payload records\n");
  }
  const auto shard_sort = SelectSortingFunc(
      NWords, 1, Lv2SortIgnoredLowBytes(), 0);

  auto postprocess_suffix64 =
      [&](CountSuffix64Record *records, uint64_t begin, uint64_t end,
          uint16_t prefix, int thread_id) {
        uint32_t packed_edge[3];
        int64_t count_prev[5];
        int64_t count_next[5];
        for (uint64_t from = begin; from < end;) {
          uint64_t to = from + 1u;
          const uint64_t suffix = records[from].suffix;
          while (to < end && records[to].suffix == suffix) ++to;
          const int64_t count = static_cast<int64_t>(to - from);
          std::fill(count_prev, count_prev + 5, int64_t{0});
          std::fill(count_next, count_next + 5, int64_t{0});
          for (uint64_t i = from; i < to; ++i) {
            const unsigned context = static_cast<unsigned>(
                LoadCountPayload48(records[i].payload) & 0x3Fu);
            ++count_prev[context >> 3u];
            ++count_next[context & 7u];
          }
          bool has_in = false;
          bool has_out = false;
          bool guaranteed_in = false;
          bool guaranteed_out = false;
          uint8_t verify_in_mask = 0u;
          uint8_t verify_out_mask = 0u;
          const int64_t unambiguous_context =
              std::max<int64_t>(1, opt_.solid_threshold) * 2;
          for (unsigned base = 0; base < 4; ++base) {
            has_in |= count_prev[base] >= opt_.solid_threshold;
            has_out |= count_next[base] >= opt_.solid_threshold;
            // A canonical edge can collect context observations from both
            // physical orientations.  Fewer than 2*T observations may split
            // into two sub-threshold adjacent edges, so context alone cannot
            // prove that side is connected.  At 2*T, at least one of the two
            // orientations necessarily has T observations (pigeonhole), and
            // the solid-neighbour decision is exact.
            guaranteed_in |= count_prev[base] >= unambiguous_context;
            guaranteed_out |= count_next[base] >= unambiguous_context;
            if (count_prev[base] >= opt_.solid_threshold &&
                count_prev[base] < unambiguous_context) {
              verify_in_mask |= static_cast<uint8_t>(1u << base);
            }
            if (count_next[base] >= opt_.solid_threshold &&
                count_next[base] < unambiguous_context) {
              verify_out_mask |= static_cast<uint8_t>(1u << base);
            }
          }
          if ((!has_in || !has_out) && count >= opt_.solid_threshold) {
            for (uint64_t i = from; i < to; ++i) {
              const uint64_t locator =
                  LoadCountPayload48(records[i].payload) >> 6u;
              const unsigned strand = locator & 1u;
              const uint32_t offset = static_cast<uint32_t>(
                  (locator >> 1u) & read_offset_mask_);
              const uint32_t read_id =
                  static_cast<uint32_t>(locator >> read_locator_shift_);
              if ((!has_in && strand == 0) ||
                  (!has_out && strand != 0)) {
                UpdateLast0In(read_id, offset);
              }
              if ((!has_out && strand == 0) ||
                  (!has_in && strand != 0)) {
                UpdateFirst0Out(read_id, offset + 1u);
              }
            }
          }
          edge_counter_.Add(count, thread_id);
          if (count >= opt_.solid_threshold) {
            packed_edge[0] = (static_cast<uint32_t>(prefix) << 16u) |
                             static_cast<uint32_t>(suffix >> 48u);
            packed_edge[1] = static_cast<uint32_t>(suffix >> 16u);
            packed_edge[2] =
                (static_cast<uint32_t>(suffix) << 16u) |
                static_cast<uint32_t>(std::min<int64_t>(kMaxMul, count));
            edge_writer_.WriteUnordered(packed_edge, thread_id);
            WriteEndpointCandidate(packed_edge, !has_in, !has_out,
                                   guaranteed_in ? 0u : verify_in_mask,
                                   guaranteed_out ? 0u : verify_out_mask,
                                   thread_id);
            uint64_t *histogram = seq_bucket_histograms_.data() +
                                  static_cast<size_t>(thread_id) * kNumBuckets;
            const uint16_t key = ExtractCountEdgeWindow8(packed_edge, 1u);
            const unsigned reverse_window_offset =
                opt_.k - kBucketPrefixLength;
            const uint16_t reverse_key = ReverseComplementCountWindow8(
                ExtractCountEdgeWindow8(packed_edge, reverse_window_offset));
            ++histogram[key];
            ++histogram[reverse_key];
          }
          from = to;
        }
      };

  std::vector<double> parse_wall(spool.num_containers(), 0.0);
  std::vector<double> materialize_cpu(n_threads, 0.0);
  std::vector<double> prefix_cpu(n_threads, 0.0);
  std::vector<double> sort_cpu(n_threads, 0.0);
  std::vector<double> postprocess_cpu(n_threads, 0.0);
  std::atomic<uint64_t> processed{0};
  const double consume_begin = omp_get_wtime();
  const unsigned low_shards = 1u << spool.low_bits();
  const NumaTopology &topology = GetNumaTopology();
  const size_t num_domains = topology.domain_count();
  std::vector<unsigned> worker_domains(n_threads, 0);
  std::vector<unsigned> domain_workers(num_domains, 0);
  std::vector<CountExternalWorkspace> workspaces(n_threads);
#pragma omp parallel num_threads(n_threads)
  {
    const int tid = omp_get_thread_num();
    worker_domains[tid] = std::min<unsigned>(
        CurrentNumaDomain(), static_cast<unsigned>(num_domains - 1u));
    CountExternalWorkspace &workspace = workspaces[tid];
    if (use_suffix64) {
      // These tables are reused for every coarse container.  Constructing
      // them on their owning worker also gives the first-touch policy the
      // same NUMA placement as the radix workspaces used below.
      workspace.prefix_counts.assign(kNumBuckets, uint64_t{0});
      workspace.prefix_cursors.resize(kNumBuckets);
      workspace.touched_prefixes.reserve(kNumBuckets);
      workspace.prefix_groups.reserve(kNumBuckets);
    }
  }
  for (unsigned domain : worker_domains) {
    ++domain_workers[domain];
  }

  // Reuse one anonymous run buffer for all coarse containers.  Mapping and
  // unmapping a multi-GiB file per container synchronously tears down millions
  // of 4-KiB PTEs over a full count.  A persistent anonymous arena is touched
  // in parallel, can collapse to transparent huge pages, and is overwritten
  // in place by the next container, so its RSS is exactly one-container sized
  // rather than cumulative.
  uint64_t max_container_bytes = 0;
  for (unsigned container = 0; container < spool.num_containers();
       ++container) {
    max_container_bytes =
        std::max(max_container_bytes, spool.Size(container));
  }
  if (max_container_bytes > std::numeric_limits<size_t>::max()) {
    xfatal("External count container exceeds addressable memory\n");
  }
  void *run_storage = nullptr;
  if (max_container_bytes != 0) {
    run_storage = mmap(nullptr, static_cast<size_t>(max_container_bytes),
                       PROT_READ | PROT_WRITE,
                       MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (run_storage == MAP_FAILED) {
      xfatal("Cannot allocate reusable external count run buffer\n");
    }
    AdviseHugePages(run_storage, static_cast<size_t>(max_container_bytes));
  }
  AsyncCountFileRemover run_file_remover;

  for (unsigned container = 0; container < spool.num_containers();
       ++container) {
    const double parse_begin = omp_get_wtime();
    const std::string path = spool.Path(container);
    const int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) xfatal("Cannot open external count container {s}\n", path.c_str());
    const uint64_t file_bytes = spool.Size(container);
    if (file_bytes == 0) {
      close(fd);
      run_file_remover.Enqueue(spool.Path(container));
      continue;
    }
    std::atomic<bool> read_ok{true};
#pragma omp parallel num_threads(n_threads)
    {
      const int tid = omp_get_thread_num();
      const uint64_t read_begin =
          file_bytes * static_cast<uint64_t>(tid) / n_threads;
      const uint64_t read_end =
          file_bytes * static_cast<uint64_t>(tid + 1) / n_threads;
      uint64_t offset = read_begin;
      while (offset < read_end) {
        const size_t request = static_cast<size_t>(std::min<uint64_t>(
            read_end - offset, uint64_t{1} << 30u));
        const ssize_t result =
            pread(fd, static_cast<uint8_t *>(run_storage) + offset,
                  request, static_cast<off_t>(offset));
        if (result <= 0) {
          read_ok.store(false, std::memory_order_relaxed);
          break;
        }
        offset += static_cast<uint64_t>(result);
      }
    }
    close(fd);
    if (!read_ok.load(std::memory_order_relaxed)) {
      xfatal("Cannot read external count container {s}\n", path.c_str());
    }
    const uint8_t *mapped = static_cast<const uint8_t *>(run_storage);

    struct ExternalBlock {
      uint64_t payload_begin;
      uint64_t payload_end;
      uint32_t num_records;
      uint32_t producer_thread;
    };
    uint64_t total_container_runs = 0;
    std::vector<ExternalBlock> blocks;
    std::vector<std::vector<size_t>> producer_blocks(n_threads);
    for (uint64_t block = 0; block < file_bytes;) {
      if (file_bytes - block < sizeof(CountRunSpoolBlockHeader)) {
        xfatal("Truncated external count block header\n");
      }
      CountRunSpoolBlockHeader header;
      std::memcpy(&header, mapped + block, sizeof(header));
      const uint64_t block_bytes =
          sizeof(header) + static_cast<uint64_t>(header.payload_bytes);
      if (block_bytes > file_bytes - block) {
        xfatal("Truncated external count block payload\n");
      }
      if (header.producer_thread >= static_cast<uint32_t>(n_threads)) {
        xfatal("Invalid external count block producer {}\n",
               header.producer_thread);
      }
      blocks.push_back({block + sizeof(header), block + block_bytes,
                        header.num_records, header.producer_thread});
      producer_blocks[header.producer_thread].push_back(blocks.size() - 1u);
      total_container_runs += header.num_records;
      block += block_bytes;
    }

    std::vector<uint64_t> run_counts(low_shards, 0);
    std::vector<uint64_t> occurrence_counts(low_shards, 0);
    for (unsigned low = 0; low < low_shards; ++low) {
      const unsigned shard = (container << spool.low_bits()) | low;
      if (shard >= num_shards) continue;
      for (int thread = 0; thread < n_threads; ++thread) {
        run_counts[low] += thread_run_counts[
            static_cast<size_t>(thread) * num_shards + shard];
        occurrence_counts[low] += thread_occurrence_counts[
            static_cast<size_t>(thread) * num_shards + shard];
      }
    }

    std::vector<uint64_t> run_offsets(low_shards + 1u, 0);
    for (unsigned low = 0; low < low_shards; ++low) {
      run_offsets[low + 1u] = run_offsets[low] + run_counts[low];
    }
    std::vector<uint64_t> thread_cursors(
        static_cast<size_t>(n_threads) * low_shards, 0);
    for (unsigned low = 0; low < low_shards; ++low) {
      uint64_t cursor = run_offsets[low];
      const unsigned shard = (container << spool.low_bits()) | low;
      for (int thread = 0; thread < n_threads; ++thread) {
        thread_cursors[static_cast<size_t>(thread) * low_shards + low] =
            cursor;
        if (shard < num_shards) {
          cursor += thread_run_counts[
              static_cast<size_t>(thread) * num_shards + shard];
        }
      }
      assert(cursor == run_offsets[low + 1u]);
    }
    if (run_offsets.back() != total_container_runs) {
      xfatal("External count run total mismatch: {} != {}\n",
             run_offsets.back(), total_container_runs);
    }

    // Producers already retained exact per-shard run counts during the one
    // read scan.  The block header identifies its producer, so one parser pass
    // can scatter record offsets directly into disjoint producer/shard ranges;
    // no count-only pass over the variable-length 100+ GiB stream is needed.
    std::vector<CountExternalRunRef> refs(
        static_cast<size_t>(total_container_runs));
    std::atomic<bool> parse_ok{true};
#pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int producer = 0; producer < n_threads; ++producer) {
      uint64_t *local_cursors =
          thread_cursors.data() +
          static_cast<size_t>(producer) * low_shards;
      for (size_t block_id : producer_blocks[producer]) {
        const ExternalBlock &block = blocks[block_id];
        uint64_t record_offset = block.payload_begin;
        for (uint32_t record_in_block = 0;
             record_in_block < block.num_records; ++record_in_block) {
          if (block.payload_end - record_offset < 8u) {
            parse_ok.store(false, std::memory_order_relaxed);
            break;
          }
          const uint8_t *record = mapped + record_offset;
          const unsigned low = DecodeCountRunLowShard(
              record, spool.read_shard_bits());
          const unsigned run_length = record[6];
          if (low >= low_shards || run_length == 0) {
            parse_ok.store(false, std::memory_order_relaxed);
            break;
          }
          const unsigned sequence_bases = edge_length + run_length - 1u;
          const uint64_t record_bytes =
              8u + DivCeiling(sequence_bases,
                              SeqPackage::kBasesPerWord) *
                       sizeof(uint32_t);
          if (record_bytes > block.payload_end - record_offset) {
            parse_ok.store(false, std::memory_order_relaxed);
            break;
          }
          refs[local_cursors[low]++] = {record_offset};
          record_offset += record_bytes;
        }
        if (record_offset != block.payload_end) {
          parse_ok.store(false, std::memory_order_relaxed);
        }
      }
    }
    if (!parse_ok.load(std::memory_order_relaxed)) {
      xfatal("Invalid external count run container\n");
    }

    std::vector<unsigned> task_order;
    std::vector<uint64_t> task_work(low_shards, 0);
    task_order.reserve(low_shards);
    for (unsigned low = 0; low < low_shards; ++low) {
      const unsigned shard = (container << spool.low_bits()) | low;
      if (run_counts[low] != 0) {
        if (shard >= num_shards ||
            occurrence_counts[low] != shard_occurrences[shard]) {
          xfatal("External count shard {} mismatch: {} != {}\n", shard,
                 occurrence_counts[low], shard_occurrences[shard]);
        }
        // Each occurrence pays one rolling/update cost.  Each run additionally
        // initializes a k-mer pair and starts an independent mapped-record
        // access; express that fixed cost in machine-word chunks of the
        // current edge length.  The weight follows algorithmic dimensions and
        // remains independent of input coverage or a named processor.
        const uint64_t run_start_weight =
            std::max<uint64_t>(1u, edge_length / 8u);
        task_work[low] = occurrence_counts[low] +
                         run_counts[low] * run_start_weight;
        task_order.push_back(low);
      }
    }
    std::sort(task_order.begin(), task_order.end(),
              [&](unsigned lhs, unsigned rhs) {
                return task_work[lhs] > task_work[rhs];
              });

    // Each producer's blocks are parsed by the corresponding persistent
    // worker.  Keep a shard on the domain that owns most of its already-faulted
    // pages when load permits, and let workers steal only after exhausting
    // local work.  This is discovered from the active cpuset and remains valid
    // for one or many NUMA domains.
    std::vector<long double> domain_load(num_domains, 0);
    std::vector<std::vector<unsigned>> domain_tasks(num_domains);
    for (unsigned low : task_order) {
      const unsigned shard = (container << spool.low_bits()) | low;
      unsigned best_domain = 0;
      long double best_normalized =
          std::numeric_limits<long double>::max();
      uint64_t best_local = 0;
      for (unsigned domain = 0; domain < num_domains; ++domain) {
        if (domain_workers[domain] == 0) continue;
        uint64_t local = 0;
        for (int thread = 0; thread < n_threads; ++thread) {
          if (worker_domains[thread] == domain) {
            local += thread_occurrence_counts[
                static_cast<size_t>(thread) * num_shards + shard];
          }
        }
        const long double normalized =
            (domain_load[domain] + task_work[low]) /
            domain_workers[domain];
        if (normalized < best_normalized ||
            (normalized == best_normalized && local > best_local)) {
          best_domain = domain;
          best_normalized = normalized;
          best_local = local;
        }
      }
      domain_tasks[best_domain].push_back(low);
      domain_load[best_domain] += task_work[low];
    }
    std::unique_ptr<std::atomic<size_t>[]> domain_next(
        new std::atomic<size_t>[num_domains]);
    for (size_t domain = 0; domain < num_domains; ++domain) {
      domain_next[domain].store(0, std::memory_order_relaxed);
    }
    const bool use_global_task_queue =
        std::getenv("MEGAHIT_EXTERNAL_GLOBAL_TASK_QUEUE") != nullptr;
    std::atomic<size_t> global_next{0};
    parse_wall[container] = omp_get_wtime() - parse_begin;

#pragma omp parallel num_threads(n_threads)
    {
      const int tid = omp_get_thread_num();
      const unsigned local_domain = worker_domains[tid];
      auto pop_from_domain = [&](unsigned domain, unsigned *low) {
        const size_t task =
            domain_next[domain].fetch_add(1, std::memory_order_relaxed);
        if (task >= domain_tasks[domain].size()) return false;
        *low = domain_tasks[domain][task];
        return true;
      };
      CountExternalWorkspace &workspace = workspaces[tid];
      uint32_t *&record_workspace = workspace.record_workspace;
      size_t &record_workspace_bytes = workspace.record_workspace_bytes;
      void *&suffix_workspace = workspace.suffix_workspace;
      size_t &suffix_workspace_bytes = workspace.suffix_workspace_bytes;
      std::vector<uint64_t> &prefix_counts = workspace.prefix_counts;
      std::vector<uint64_t> &prefix_cursors = workspace.prefix_cursors;
      std::vector<uint16_t> &touched_prefixes =
          workspace.touched_prefixes;
      auto &prefix_groups = workspace.prefix_groups;

      for (;;) {
        unsigned low = 0;
        bool found = false;
        if (use_global_task_queue) {
          const size_t task =
              global_next.fetch_add(1, std::memory_order_relaxed);
          if (task < task_order.size()) {
            low = task_order[task];
            found = true;
          }
        } else {
          found = pop_from_domain(local_domain, &low);
          if (!found) {
            for (size_t distance = 1; distance < num_domains; ++distance) {
              const unsigned victim = static_cast<unsigned>(
                  (local_domain + distance) % num_domains);
              if (pop_from_domain(victim, &low)) {
                found = true;
                break;
              }
            }
          }
        }
        if (!found) break;
        const uint64_t n_records = occurrence_counts[low];
        const size_t workspace_bytes =
            static_cast<size_t>(n_records) * (NWords + 1u) *
            sizeof(uint32_t);
        if (workspace_bytes > record_workspace_bytes) {
          void *memory = std::realloc(record_workspace, workspace_bytes);
          if (memory == nullptr) {
            std::free(record_workspace);
            xfatal("Cannot allocate external count workspace\n");
          }
          record_workspace = static_cast<uint32_t *>(memory);
          record_workspace_bytes = workspace_bytes;
          BindMemoryPagesToNumaDomain(record_workspace,
                                      record_workspace_bytes, local_domain);
          AdviseHugePages(record_workspace, record_workspace_bytes);
        }

        const double materialize_begin = omp_get_wtime();
        uint32_t *out = record_workspace;
        CountPrefix64Record *prefix_out =
            reinterpret_cast<CountPrefix64Record *>(record_workspace);
        for (uint64_t run_index = run_offsets[low];
             run_index < run_offsets[low + 1u]; ++run_index) {
          const uint8_t *record = mapped + refs[run_index].offset;
          const uint32_t read_id = DecodeCountRunReadId(
              record, spool.read_shard_bits());
          const unsigned start = record[5];
          const unsigned run_length = record[6];
          const unsigned boundary_context = record[7] & 0x3Fu;
          const uint32_t *packed_bases =
              reinterpret_cast<const uint32_t *>(record + 8u);
          const auto base_at = [&](unsigned position) -> uint8_t {
            return static_cast<uint8_t>(
                SeqPackage::TVector::at(packed_bases, position));
          };
          Kmer<NWords, uint32_t> edge;
          edge.InitFromPtr(packed_bases, 0, edge_length);
          Kmer<NWords, uint32_t> reverse = edge;
          reverse.ReverseComplement(edge_length);
          uint64_t encoded_offset =
              (static_cast<uint64_t>(read_id) << read_locator_shift_) |
              (static_cast<uint64_t>(start) << 1u);

          for (unsigned item = 0; item < run_length; ++item) {
            const bool is_reverse =
                reverse.cmp(edge, edge_length) < 0;
            const auto &canonical = is_reverse ? reverse : edge;
            const unsigned previous =
                item == 0 ? boundary_context >> 3u : base_at(item - 1u);
            const unsigned next =
                item + edge_length < edge_length + run_length - 1u
                    ? base_at(item + edge_length)
                    : boundary_context & 7u;
            const uint64_t locator =
                encoded_offset | static_cast<unsigned>(is_reverse);
            uint64_t context;
            if (!is_reverse) {
              context = (previous << 3u) | next;
            } else {
              const unsigned canonical_previous =
                  next == kSentinelValue ? kSentinelValue : 3u - next;
              const unsigned canonical_next =
                  previous == kSentinelValue ? kSentinelValue
                                             : 3u - previous;
              context = (canonical_previous << 3u) | canonical_next;
            }

            if (use_suffix64) {
              const uint32_t *words = canonical.data();
              prefix_out->prefix =
                  static_cast<uint16_t>(words[0] >> 16u);
              prefix_out->suffix =
                  (static_cast<uint64_t>(words[0] & 0xFFFFu) << 48u) |
                  (static_cast<uint64_t>(words[1]) << 16u) |
                  (static_cast<uint64_t>(words[2]) >> 16u);
              StoreCountPayload48(prefix_out->payload,
                                  (locator << 6u) | context);
              if (prefix_counts[prefix_out->prefix]++ == 0) {
                touched_prefixes.push_back(prefix_out->prefix);
              }
              ++prefix_out;
            } else {
              const uint32_t metadata =
                  (static_cast<uint32_t>(locator >> 32u) << 6u) |
                  static_cast<uint32_t>(context);
              assert((metadata & ~compact_key_mask_) == 0);
              std::memcpy(out, canonical.data(),
                          NWords * sizeof(uint32_t));
              assert((out[NWords - 1u] & compact_key_mask_) == 0);
              out[NWords - 1u] |= metadata;
              out[NWords] = static_cast<uint32_t>(locator);
              out += NWords + 1u;
            }

            if (item + 1u != run_length) {
              encoded_offset += 2u;
              const uint8_t c = base_at(item + edge_length);
              edge.ShiftAppend(c, edge_length);
              reverse.ShiftPreappend(3u - c, edge_length);
            }
          }
        }
        materialize_cpu[tid] += omp_get_wtime() - materialize_begin;

        if (use_suffix64) {
          assert(prefix_out ==
                 reinterpret_cast<CountPrefix64Record *>(record_workspace) +
                     n_records);
          const double prefix_begin = omp_get_wtime();
          prefix_groups.clear();
          uint64_t running = 0;
          for (uint16_t prefix : touched_prefixes) {
            const uint64_t count = prefix_counts[prefix];
            prefix_cursors[prefix] = running;
            prefix_groups.push_back(
                {prefix, {running, running + count}});
            running += count;
          }
          const size_t suffix_bytes =
              static_cast<size_t>(n_records) * sizeof(CountSuffix64Record);
          if (suffix_bytes > suffix_workspace_bytes) {
            void *memory = std::realloc(suffix_workspace, suffix_bytes);
            if (memory == nullptr) {
              std::free(suffix_workspace);
              std::free(record_workspace);
              xfatal("Cannot allocate external suffix workspace\n");
            }
            suffix_workspace = memory;
            suffix_workspace_bytes = suffix_bytes;
            BindMemoryPagesToNumaDomain(suffix_workspace,
                                        suffix_workspace_bytes,
                                        local_domain);
            AdviseHugePages(suffix_workspace, suffix_workspace_bytes);
          }
          CountPrefix64Record *prefix_records =
              reinterpret_cast<CountPrefix64Record *>(record_workspace);
          CountSuffix64Record *suffix_records =
              static_cast<CountSuffix64Record *>(suffix_workspace);
          for (uint64_t i = 0; i < n_records; ++i) {
            const uint16_t prefix = prefix_records[i].prefix;
            const uint64_t destination = prefix_cursors[prefix]++;
            suffix_records[destination].suffix = prefix_records[i].suffix;
            std::memcpy(suffix_records[destination].payload,
                        prefix_records[i].payload,
                        sizeof(suffix_records[destination].payload));
          }
          for (uint16_t prefix : touched_prefixes) {
            prefix_counts[prefix] = 0;
          }
          touched_prefixes.clear();
          prefix_cpu[tid] += omp_get_wtime() - prefix_begin;

          const double sort_begin = omp_get_wtime();
          for (const auto &group : prefix_groups) {
            if (group.second.second - group.second.first > 1u) {
              kmlib::kmsort(suffix_records + group.second.first,
                            suffix_records + group.second.second);
            }
          }
          sort_cpu[tid] += omp_get_wtime() - sort_begin;

          const double post_begin = omp_get_wtime();
          for (const auto &group : prefix_groups) {
            postprocess_suffix64(suffix_records, group.second.first,
                                 group.second.second, group.first, tid);
          }
          postprocess_cpu[tid] += omp_get_wtime() - post_begin;
        } else {
          assert(out == record_workspace +
                            static_cast<size_t>(n_records) * (NWords + 1u));
          const double sort_begin = omp_get_wtime();
          shard_sort(record_workspace, static_cast<int64_t>(n_records));
          sort_cpu[tid] += omp_get_wtime() - sort_begin;
          const double post_begin = omp_get_wtime();
          if (compact_locator_high_mask_ != 0) {
            Lv2PostprocessFor<true, true, true>(
                0, static_cast<int64_t>(n_records), tid, record_workspace);
          } else {
            Lv2PostprocessFor<true, true, false>(
                0, static_cast<int64_t>(n_records), tid, record_workspace);
          }
          postprocess_cpu[tid] += omp_get_wtime() - post_begin;
        }
        processed.fetch_add(n_records, std::memory_order_relaxed);
      }
    }

    run_file_remover.Enqueue(spool.Path(container));
  }

  run_file_remover.Finish();
  if (run_storage != nullptr) {
    munmap(run_storage, static_cast<size_t>(max_container_bytes));
  }

  const double consume_wall = omp_get_wtime() - consume_begin;
  if (processed.load(std::memory_order_relaxed) != total_occurrences) {
    xfatal("External count processed {} of {} records\n",
           processed.load(std::memory_order_relaxed), total_occurrences);
  }
  double min_worker_compute = std::numeric_limits<double>::max();
  double max_worker_compute = 0;
  for (int thread = 0; thread < n_threads; ++thread) {
    const double compute = materialize_cpu[thread] + prefix_cpu[thread] +
                           sort_cpu[thread] + postprocess_cpu[thread];
    min_worker_compute = std::min(min_worker_compute, compute);
    max_worker_compute = std::max(max_worker_compute, compute);
  }
  xinfo("External segmented count profile: parse={.4} s, consume={.4} s; "
        "CPU-s materialize={.3}, prefix/compact={.3}, sort={.3}, "
        "count/output={.3}; worker compute min/max={.3}/{.3} s\n",
        std::accumulate(parse_wall.begin(), parse_wall.end(), 0.0),
        consume_wall,
        std::accumulate(materialize_cpu.begin(), materialize_cpu.end(), 0.0),
        std::accumulate(prefix_cpu.begin(), prefix_cpu.end(), 0.0),
        std::accumulate(sort_cpu.begin(), sort_cpu.end(), 0.0),
        std::accumulate(postprocess_cpu.begin(), postprocess_cpu.end(), 0.0),
        min_worker_compute, max_worker_compute);
  return true;
#endif
}

template <unsigned NWords>
bool KmerCounter::RunSegmentedCountFor() {
  const int n_threads = opt_.n_threads;
  const uint64_t num_reads = seq_pkg_.seq_count();
  const unsigned edge_length = opt_.k + 1;
  // A longer minimizer increases the number of shared windows between
  // adjacent (k+1)-mers without making the rolling code overflow uint64_t.
  const unsigned minimizer_length =
      std::min<unsigned>(15, std::max<unsigned>(7, edge_length / 2));
  const unsigned minimizer_window = edge_length - minimizer_length + 1;
  const uint64_t minimizer_mask =
      minimizer_length == 32
          ? std::numeric_limits<uint64_t>::max()
          : (uint64_t{1} << (2u * minimizer_length)) - 1u;
  const unsigned rc_shift = 2u * (minimizer_length - 1u);

  const uint64_t estimated_occurrences =
      seq_pkg_.base_count() > num_reads * static_cast<uint64_t>(opt_.k)
          ? seq_pkg_.base_count() -
                num_reads * static_cast<uint64_t>(opt_.k)
          : 0;
  const unsigned shard_bits = SelectCountShardBits(
      estimated_occurrences, (NWords + 1u) * sizeof(uint32_t), n_threads,
      opt_.host_mem);
  const unsigned num_shards = 1u << shard_bits;
  const NumaTopology &topology = GetNumaTopology();

  xinfo(
      "Exact segmented count: {} minimizer shards (runtime {} bits), m={}, "
      "window={}; one read scan; {} NUMA domains, {} MiB aggregate LLC\n",
      num_shards, shard_bits, minimizer_length, minimizer_window,
      topology.domain_count(),
      topology.total_last_level_cache_bytes() >> 20u);

  // RabbitBin-style producer ownership remains lock-free, but descriptors are
  // now one contiguous stream per producer.  This removes millions of
  // per-(thread,shard) vector growth allocations.  Rows in the two count
  // matrices are first-touched by the producer that owns them.
  std::vector<CountSegmentStream> thread_streams(n_threads);
  const size_t matrix_items = static_cast<size_t>(n_threads) * num_shards;
  std::vector<uint64_t> thread_segment_counts(matrix_items, 0);
  std::vector<uint64_t> thread_occurrence_counts(matrix_items, 0);
  DiscardMemoryPages(thread_segment_counts.data(),
                     thread_segment_counts.size() * sizeof(uint64_t));
  DiscardMemoryPages(thread_occurrence_counts.data(),
                     thread_occurrence_counts.size() * sizeof(uint64_t));
  AdviseHugePages(thread_segment_counts.data(),
                  thread_segment_counts.size() * sizeof(uint64_t));
  AdviseHugePages(thread_occurrence_counts.data(),
                  thread_occurrence_counts.size() * sizeof(uint64_t));
  std::vector<uint64_t> thread_occurrences(n_threads, 0);
  std::vector<unsigned> producer_domains(n_threads, 0);

  const double scan_wall_begin = omp_get_wtime();

#pragma omp parallel num_threads(n_threads)
  {
    const int tid = omp_get_thread_num();
    producer_domains[tid] = CurrentNumaDomain();
    uint64_t *segment_counts =
        thread_segment_counts.data() + static_cast<size_t>(tid) * num_shards;
    uint64_t *occurrence_counts = thread_occurrence_counts.data() +
                                  static_cast<size_t>(tid) * num_shards;
    std::fill(segment_counts, segment_counts + num_shards, uint64_t{0});
    std::fill(occurrence_counts, occurrence_counts + num_shards, uint64_t{0});
    const uint64_t read_from = num_reads * static_cast<uint64_t>(tid) /
                               static_cast<uint64_t>(n_threads);
    const uint64_t read_to = num_reads * static_cast<uint64_t>(tid + 1) /
                             static_cast<uint64_t>(n_threads);
    assert(read_from <= std::numeric_limits<uint32_t>::max());
    thread_streams[tid].Configure(static_cast<uint32_t>(read_from),
                                  read_to - read_from, shard_bits);
    std::vector<uint64_t> hashes;
    std::vector<unsigned> min_queue;
    hashes.reserve(seq_pkg_.max_length());
    min_queue.reserve(seq_pkg_.max_length());
    uint64_t local_occurrences = 0;

    // A random minimizer changes about twice per minimizer window.  This
    // reserve is only a capacity hint; repeats use less and adversarial input
    // grows safely in 12.5% increments.
    const uint64_t local_estimated_occurrences =
        estimated_occurrences / n_threads +
        (static_cast<uint64_t>(tid) < estimated_occurrences % n_threads);
    const long double expected_runs =
        static_cast<long double>(local_estimated_occurrences) * 2.0L /
        std::max(1u, minimizer_window);
    if (expected_runs < std::numeric_limits<size_t>::max()) {
      thread_streams[tid].reserve(
          static_cast<size_t>(expected_runs) + (size_t{1} << 16u));
    }

    auto emit_segment = [&](uint32_t read_id, unsigned start,
                            unsigned length, unsigned shard) {
      thread_streams[tid].push_back(read_id, start, length, shard);
      ++segment_counts[shard];
      occurrence_counts[shard] += length;
    };

    for (uint64_t rid = read_from; rid < read_to; ++rid) {
      const auto seq = seq_pkg_.GetSeqView(rid);
      const unsigned read_length = seq.length();
      if (read_length < edge_length) {
        continue;
      }

      const unsigned num_mmers = read_length - minimizer_length + 1;
      const unsigned num_edges = read_length - edge_length + 1;
      hashes.resize(num_mmers);
      min_queue.resize(num_mmers);

      uint64_t fwd = 0;
      uint64_t rev = 0;
      for (unsigned i = 0; i < minimizer_length; ++i) {
        const uint64_t c = seq.base_at(i);
        fwd = (fwd << 2u) | c;
        rev = (rev >> 2u) | ((uint64_t{3} - c) << rc_shift);
      }
      for (unsigned pos = 0; pos < num_mmers; ++pos) {
        const uint64_t canonical = std::min(fwd, rev);
        hashes[pos] = MixCountMinimizer(
            canonical ^ (uint64_t{0x6a09e667f3bcc909} +
                         minimizer_length));
        if (pos + 1 < num_mmers) {
          const uint64_t c = seq.base_at(pos + minimizer_length);
          fwd = ((fwd << 2u) | c) & minimizer_mask;
          rev = (rev >> 2u) | ((uint64_t{3} - c) << rc_shift);
        }
      }

      unsigned q_begin = 0;
      unsigned q_end = 0;
      unsigned current_shard = 0;
      unsigned run_start = 0;
      unsigned run_length = 0;
      for (unsigned pos = 0; pos < num_mmers; ++pos) {
        while (q_end > q_begin &&
               hashes[min_queue[q_end - 1]] > hashes[pos]) {
          --q_end;
        }
        min_queue[q_end++] = pos;
        if (pos + 1 < minimizer_window) {
          continue;
        }
        const unsigned edge_start = pos + 1 - minimizer_window;
        while (q_begin < q_end && min_queue[q_begin] < edge_start) {
          ++q_begin;
        }
        assert(edge_start < num_edges);
        const unsigned shard = static_cast<unsigned>(
            hashes[min_queue[q_begin]]) & (num_shards - 1u);
        if (run_length == 0) {
          current_shard = shard;
          run_start = edge_start;
          run_length = 1;
        } else if (shard == current_shard &&
                   run_length < std::numeric_limits<uint8_t>::max()) {
          ++run_length;
        } else {
          emit_segment(static_cast<uint32_t>(rid), run_start, run_length,
                       current_shard);
          current_shard = shard;
          run_start = edge_start;
          run_length = 1;
        }
      }
      if (run_length != 0) {
        emit_segment(static_cast<uint32_t>(rid), run_start, run_length,
                     current_shard);
      }
      local_occurrences += num_edges;
    }
    thread_occurrences[tid] = local_occurrences;
    thread_streams[tid].Finalize();
    AdviseHugePages(thread_streams[tid].data(),
                    thread_streams[tid].allocation_bytes());
  }
  const double scan_wall = omp_get_wtime() - scan_wall_begin;

  uint64_t total_segments = 0;
  uint64_t descriptor_bytes = 0;
  uint64_t descriptor_capacity_bytes = 0;
  for (const CountSegmentStream &stream : thread_streams) {
    total_segments += stream.size();
    descriptor_bytes += stream.used_bytes();
    descriptor_capacity_bytes += stream.capacity_bytes();
  }
  const uint64_t total_occurrences =
      std::accumulate(thread_occurrences.begin(), thread_occurrences.end(),
                      uint64_t{0});
  xinfo("Segmented count produced {} runs for {} exact records ({.3} "
        "records/run), {.3} GiB used / {.3} GiB capacity; scan {.4} s\n",
        total_segments, total_occurrences,
        total_segments == 0
            ? 0.0
            : static_cast<double>(total_occurrences) / total_segments,
        static_cast<double>(descriptor_bytes) / (uint64_t{1} << 30u),
        static_cast<double>(descriptor_capacity_bytes) /
            (uint64_t{1} << 30u),
        scan_wall);

  std::vector<uint64_t> shard_occurrences(num_shards, 0);
#pragma omp parallel for schedule(static) num_threads(n_threads)
  for (unsigned shard = 0; shard < num_shards; ++shard) {
    uint64_t count = 0;
    for (int producer = 0; producer < n_threads; ++producer) {
      count += thread_occurrence_counts[
          static_cast<size_t>(producer) * num_shards + shard];
    }
    shard_occurrences[shard] = count;
  }

  const uint64_t observed_occurrences =
      std::accumulate(shard_occurrences.begin(), shard_occurrences.end(),
                      uint64_t{0});
  if (observed_occurrences != total_occurrences) {
    xfatal("Segmented count record mismatch: {} != {}\n",
           observed_occurrences, total_occurrences);
  }

  // Convert each producer stream to contiguous shard ranges with an in-place
  // American-flag partition.  There is no second read scan and no descriptor
  // copy.  The offsets remain immutable while the segment-count matrix is
  // reused as the per-bucket cursor array.
  const double grouping_wall_begin = omp_get_wtime();
  std::vector<uint64_t> producer_shard_offsets(
      static_cast<size_t>(n_threads) * (num_shards + 1u), 0);
  DiscardMemoryPages(producer_shard_offsets.data(),
                     producer_shard_offsets.size() * sizeof(uint64_t));
#pragma omp parallel for schedule(static) num_threads(n_threads)
  for (int producer = 0; producer < n_threads; ++producer) {
    uint64_t *offsets = producer_shard_offsets.data() +
                        static_cast<size_t>(producer) * (num_shards + 1u);
    uint64_t *cursors = thread_segment_counts.data() +
                        static_cast<size_t>(producer) * num_shards;
    offsets[0] = 0;
    for (unsigned shard = 0; shard < num_shards; ++shard) {
      offsets[shard + 1u] = offsets[shard] + cursors[shard];
      cursors[shard] = offsets[shard];
    }
    CountSegmentStream &stream = thread_streams[producer];
    assert(offsets[num_shards] == stream.size());
    for (unsigned shard = 0; shard < num_shards; ++shard) {
      while (cursors[shard] < offsets[shard + 1u]) {
        const unsigned target =
            stream.ShardAt(cursors[shard]);
        assert(target < num_shards);
        if (target == shard) {
          ++cursors[shard];
        } else {
          assert(cursors[target] < offsets[target + 1u]);
          stream.Swap(cursors[shard], cursors[target]++);
        }
      }
    }
    stream.CompactGrouped();
  }
  std::vector<uint64_t>().swap(thread_segment_counts);
  const double grouping_wall = omp_get_wtime() - grouping_wall_begin;
  uint64_t grouped_descriptor_bytes = 0;
  for (const CountSegmentStream &stream : thread_streams) {
    grouped_descriptor_bytes += stream.grouped_used_bytes();
  }
  xinfo("Grouped run descriptors compacted from {.3} to {.3} GiB\n",
        static_cast<double>(descriptor_bytes) / (uint64_t{1} << 30u),
        static_cast<double>(grouped_descriptor_bytes) /
            (uint64_t{1} << 30u));

  // Shards are unrelated to the biological prefix buckets, so all biological
  // key bytes participate in the radix sort.  Only complete low bytes that
  // contain padding/context/locator metadata are ignored.
  const auto shard_sort = SelectSortingFunc(
      NWords, 1, Lv2SortIgnoredLowBytes(), 0);
  std::atomic<uint64_t> processed{0};

  // Assign shards to every runtime-discovered domain in proportion to the
  // workers actually placed there.  Weight-first greedy assignment balances
  // skewed minimizers; exact-load ties prefer the domain that produced more of
  // the shard's occurrences.  Workers exhaust their local queue before
  // stealing, so the policy remains useful for 2, 4, or many NUMA nodes while
  // retaining progress under asymmetric cpusets.
  const size_t num_domains = topology.domain_count();
  std::vector<unsigned> domain_workers(num_domains, 0);
  for (unsigned domain : producer_domains) {
    if (domain < num_domains) {
      ++domain_workers[domain];
    }
  }
  std::vector<unsigned> shard_order(num_shards);
  std::iota(shard_order.begin(), shard_order.end(), 0u);
  std::sort(shard_order.begin(), shard_order.end(),
            [&](unsigned lhs, unsigned rhs) {
              return shard_occurrences[lhs] > shard_occurrences[rhs];
            });
  std::vector<long double> domain_load(num_domains, 0);
  std::vector<std::vector<unsigned>> domain_shards(num_domains);
  for (unsigned shard : shard_order) {
    unsigned best_domain = 0;
    long double best_normalized = std::numeric_limits<long double>::max();
    uint64_t best_local = 0;
    for (unsigned domain = 0; domain < num_domains; ++domain) {
      if (domain_workers[domain] == 0) {
        continue;
      }
      uint64_t local = 0;
      for (int producer = 0; producer < n_threads; ++producer) {
        if (producer_domains[producer] == domain) {
          local += thread_occurrence_counts[
              static_cast<size_t>(producer) * num_shards + shard];
        }
      }
      const long double normalized =
          (domain_load[domain] + shard_occurrences[shard]) /
          domain_workers[domain];
      if (normalized < best_normalized ||
          (normalized == best_normalized && local > best_local)) {
        best_domain = domain;
        best_normalized = normalized;
        best_local = local;
      }
    }
    domain_shards[best_domain].push_back(shard);
    domain_load[best_domain] += shard_occurrences[shard];
  }
  std::vector<uint64_t>().swap(thread_occurrence_counts);
  std::unique_ptr<std::atomic<size_t>[]> domain_next(
      new std::atomic<size_t>[num_domains]);
  for (size_t domain = 0; domain < num_domains; ++domain) {
    domain_next[domain].store(0, std::memory_order_relaxed);
  }

  const uint64_t max_suffix64_locator =
      num_reads == 0
          ? 0
          : ((num_reads - 1u) << read_locator_shift_) |
                ((static_cast<uint64_t>(seq_pkg_.max_length() - 1u) << 1u) |
                 1u);
  const bool use_suffix64 =
      NWords == 3 && edge_length == 40 &&
      max_suffix64_locator < (uint64_t{1} << 42u) &&
      std::getenv("MEGAHIT_DISABLE_SUFFIX64_COUNT") == nullptr;

  if (use_suffix64) {
    xinfo("Count radix layout: implicit 16-bit biological prefix + "
          "14-byte suffix64/payload records\n");
  }

  const unsigned segment_prefetch_distance = [&]() {
    const char *value =
        std::getenv("MEGAHIT_SEGMENT_READ_PREFETCH_DISTANCE");
    if (value != nullptr) {
      char *end = nullptr;
      const unsigned long parsed = std::strtoul(value, &end, 10);
      if (end != value && *end == '\0' && parsed <= 256u) {
        return static_cast<unsigned>(parsed);
      }
    }
    // A minimizer run performs several rolling updates after one random read
    // fetch.  Keep enough independent reads in flight to cover latency, but
    // reduce the queue as the expected work per run grows.  This expression is
    // based on algorithmic run length and is re-evaluated for every k.
    const unsigned expected_run = std::max(1u, minimizer_window / 2u);
    return std::max(4u, std::min(32u, 96u / expected_run));
  }();
  if (use_suffix64 && segment_prefetch_distance != 0) {
    xinfo("Segmented read-gather pipeline: {} prepared reads per worker\n",
          segment_prefetch_distance);
  }

  auto materialize_suffix64 =
      [&](const CompactCountSegment &segment,
          const SeqPackage::SeqView &seq_view, CountPrefix64Record *dest,
          uint64_t *prefix_counts, std::vector<uint16_t> *touched_prefixes) {
        const uint32_t read_id = segment.read_id;
        const unsigned segment_start = segment.start;
        const unsigned segment_length = segment.length;
        assert(seq_view.id() == read_id);
        assert(segment_start + segment_length + opt_.k <= seq_view.length());

        Kmer<3, uint32_t> edge, rev_edge;
        const auto raw = seq_view.raw_address();
        edge.InitFromPtr(raw.first, raw.second + segment_start, edge_length);
        rev_edge = edge;
        rev_edge.ReverseComplement(edge_length);
        uint64_t encoded_offset =
            (static_cast<uint64_t>(read_id) << read_locator_shift_) |
            (static_cast<uint64_t>(segment_start) << 1u);
        unsigned offset = segment_start;

        for (unsigned item_index = 0; item_index < segment_length;
             ++item_index) {
          const bool reverse = rev_edge.cmp(edge, edge_length) < 0;
          const auto &canonical = reverse ? rev_edge : edge;
          const unsigned prev =
              offset == 0 ? kSentinelValue : seq_view.base_at(offset - 1u);
          const unsigned next =
              offset + edge_length < seq_view.length()
                  ? seq_view.base_at(offset + edge_length)
                  : kSentinelValue;
          const uint64_t locator =
              encoded_offset | static_cast<unsigned>(reverse);
          uint64_t context;
          if (!reverse) {
            context = (prev << 3u) | next;
          } else {
            const unsigned canonical_prev =
                next == kSentinelValue ? kSentinelValue : 3u - next;
            const unsigned canonical_next =
                prev == kSentinelValue ? kSentinelValue : 3u - prev;
            context = (canonical_prev << 3u) | canonical_next;
          }

          const uint32_t *words = canonical.data();
          dest->prefix = static_cast<uint16_t>(words[0] >> 16u);
          dest->suffix =
              (static_cast<uint64_t>(words[0] & 0xFFFFu) << 48u) |
              (static_cast<uint64_t>(words[1]) << 16u) |
              (static_cast<uint64_t>(words[2]) >> 16u);
          StoreCountPayload48(dest->payload, (locator << 6u) | context);
          if (prefix_counts[dest->prefix]++ == 0) {
            touched_prefixes->push_back(dest->prefix);
          }
          ++dest;

          if (item_index + 1u == segment_length) {
            break;
          }
          ++offset;
          encoded_offset += 2u;
          const int c = seq_view.base_at(offset + edge_length - 1u);
          edge.ShiftAppend(c, edge_length);
          rev_edge.ShiftPreappend(3 - c, edge_length);
        }
      };

  auto postprocess_suffix64 =
      [&](CountSuffix64Record *records, uint64_t begin, uint64_t end,
          uint16_t prefix, int thread_id) {
        uint32_t packed_edge[3];
        int64_t count_prev[5];
        int64_t count_next[5];
        const auto suffix_at = [&](uint64_t index) {
          return records[index].suffix;
        };
        const auto payload_at = [&](uint64_t index) {
          return LoadCountPayload48(records[index].payload);
        };
        for (uint64_t from = begin; from < end;) {
          uint64_t to = from + 1u;
          const uint64_t suffix = suffix_at(from);
          while (to < end && suffix_at(to) == suffix) {
            ++to;
          }
          const int64_t count = static_cast<int64_t>(to - from);
          std::fill(count_prev, count_prev + 5, int64_t{0});
          std::fill(count_next, count_next + 5, int64_t{0});
          for (uint64_t i = from; i < to; ++i) {
            const unsigned context =
                static_cast<unsigned>(payload_at(i) & 0x3Fu);
            ++count_prev[context >> 3u];
            ++count_next[context & 7u];
          }
          bool has_in = false;
          bool has_out = false;
          for (unsigned base = 0; base < 4; ++base) {
            has_in |= count_prev[base] >= opt_.solid_threshold;
            has_out |= count_next[base] >= opt_.solid_threshold;
          }

          if ((!has_in || !has_out) && count >= opt_.solid_threshold) {
            for (uint64_t i = from; i < to; ++i) {
              const uint64_t locator = payload_at(i) >> 6u;
              const unsigned strand = locator & 1u;
              const uint32_t offset = static_cast<uint32_t>(
                  (locator >> 1u) & read_offset_mask_);
              const uint32_t read_id =
                  static_cast<uint32_t>(locator >> read_locator_shift_);
              if ((!has_in && strand == 0) || (!has_out && strand != 0)) {
                UpdateLast0In(read_id, offset);
              }
              if ((!has_out && strand == 0) || (!has_in && strand != 0)) {
                const uint32_t first_offset = offset + 1u;
                UpdateFirst0Out(read_id, first_offset);
              }
            }
          }

          edge_counter_.Add(count, thread_id);
          if (count >= opt_.solid_threshold) {
            packed_edge[0] = (static_cast<uint32_t>(prefix) << 16u) |
                             static_cast<uint32_t>(suffix >> 48u);
            packed_edge[1] = static_cast<uint32_t>(suffix >> 16u);
            packed_edge[2] =
                (static_cast<uint32_t>(suffix) << 16u) |
                static_cast<uint32_t>(std::min<int64_t>(kMaxMul, count));
            edge_writer_.WriteUnordered(packed_edge, thread_id);
            uint64_t *histogram = seq_bucket_histograms_.data() +
                                  static_cast<size_t>(thread_id) * kNumBuckets;
            const uint16_t key = ExtractCountEdgeWindow8(packed_edge, 1u);
            const unsigned reverse_window_offset =
                opt_.k - kBucketPrefixLength;
            const uint16_t reverse_key = ReverseComplementCountWindow8(
                ExtractCountEdgeWindow8(packed_edge, reverse_window_offset));
            ++histogram[key];
            ++histogram[reverse_key];
          }
          from = to;
        }
      };

  std::vector<double> materialize_cpu(n_threads, 0.0);
  std::vector<double> prefix_partition_cpu(n_threads, 0.0);
  std::vector<double> sort_cpu(n_threads, 0.0);
  std::vector<double> postprocess_cpu(n_threads, 0.0);
  const double consume_wall_begin = omp_get_wtime();

#pragma omp parallel num_threads(n_threads)
  {
    const int tid = omp_get_thread_num();
    const unsigned local_domain = std::min<unsigned>(
        CurrentNumaDomain(), static_cast<unsigned>(num_domains - 1u));
    auto pop_from_domain = [&](unsigned domain, unsigned *shard) {
      const size_t index =
          domain_next[domain].fetch_add(1, std::memory_order_relaxed);
      if (index >= domain_shards[domain].size()) {
        return false;
      }
      *shard = domain_shards[domain][index];
      return true;
    };

    uint32_t *record_workspace = nullptr;
    size_t record_workspace_bytes = 0;
    void *suffix_workspace = nullptr;
    size_t suffix_workspace_bytes = 0;
    std::vector<uint64_t> prefix_ends(use_suffix64 ? kNumBuckets : 0, 0);
    std::vector<uint64_t> prefix_cursors(use_suffix64 ? kNumBuckets : 0, 0);
    std::vector<uint16_t> touched_prefixes;
    if (use_suffix64) {
      touched_prefixes.reserve(kNumBuckets);
    }
    struct PrefixGroup {
      uint16_t prefix;
      uint64_t begin;
      uint64_t end;
    };
    std::vector<PrefixGroup> prefix_groups;
    if (use_suffix64) {
      prefix_groups.reserve(kNumBuckets);
    }
    struct PreparedSegment {
      PreparedSegment(const CompactCountSegment &segment,
                      const SeqPackage::SeqView &seq_view)
          : segment(segment), seq_view(seq_view) {}
      CompactCountSegment segment;
      SeqPackage::SeqView seq_view;
    };
    std::vector<PreparedSegment> prepared_segments;
    if (use_suffix64 && segment_prefetch_distance != 0) {
      prepared_segments.reserve(segment_prefetch_distance);
    }

    for (;;) {
      unsigned shard = 0;
      bool found = pop_from_domain(local_domain, &shard);
      if (!found) {
        for (size_t distance = 1; distance < num_domains; ++distance) {
          const unsigned victim =
              static_cast<unsigned>((local_domain + distance) % num_domains);
          if (pop_from_domain(victim, &shard)) {
            found = true;
            break;
          }
        }
      }
      if (!found) {
        break;
      }

      const uint64_t n_records = shard_occurrences[shard];
      if (n_records == 0) {
        continue;
      }
      if (n_records > std::numeric_limits<size_t>::max() / (NWords + 1u)) {
        xfatal("Segmented count shard {} is too large\n", shard);
      }
      const size_t record_words =
          static_cast<size_t>(n_records) * (NWords + 1u);
      const size_t record_bytes = record_words * sizeof(uint32_t);
      if (record_bytes > record_workspace_bytes) {
        // Domain queues are weight-sorted, so each worker normally sees its
        // largest shard first.  Exact growth avoids retaining 25% slack across
        // every worker without causing a realloc ladder.
        const size_t new_capacity = record_bytes;
        void *memory = std::realloc(record_workspace, new_capacity);
        if (memory == nullptr) {
          std::free(record_workspace);
          xfatal("Failed to allocate segmented count workspace ({} bytes)\n",
                 new_capacity);
        }
        record_workspace = static_cast<uint32_t *>(memory);
        record_workspace_bytes = new_capacity;
        // The workspace lives with this bound OpenMP worker for the whole
        // consume phase.  Install the policy before its first materializing
        // write, then reuse the already-local pages across all later shards.
        BindMemoryPagesToNumaDomain(record_workspace, record_workspace_bytes,
                                    local_domain);
        AdviseHugePages(record_workspace, record_workspace_bytes);
      }
      uint32_t *records = record_workspace;

      const double materialize_begin = omp_get_wtime();
      uint32_t *out = records;
      CountPrefix64Record *prefix_out =
          reinterpret_cast<CountPrefix64Record *>(records);
      for (int producer = 0; producer < n_threads; ++producer) {
        const uint64_t *offsets = producer_shard_offsets.data() +
                                  static_cast<size_t>(producer) *
                                      (num_shards + 1u);
        const CountSegmentStream &stream = thread_streams[producer];
        const uint64_t segment_begin = offsets[shard];
        const uint64_t segment_end = offsets[shard + 1u];
        auto decode_segment = [&](uint64_t segment_index) {
          const uint64_t packed = stream.GroupedAt(segment_index);
          return CompactCountSegment{
              static_cast<uint32_t>(packed),
              static_cast<uint8_t>(packed >> kCountSegmentStartShift),
              static_cast<uint8_t>(packed >> kCountSegmentLengthShift)};
        };
        if (!use_suffix64) {
          for (uint64_t segment_index = segment_begin;
               segment_index < segment_end; ++segment_index) {
            const CompactCountSegment segment = decode_segment(segment_index);
            MaterializeSegment<NWords>(segment, out);
            out += static_cast<size_t>(segment.length) * (NWords + 1u);
          }
          continue;
        }

        auto prepare_segment = [&](uint64_t segment_index) {
          const CompactCountSegment segment = decode_segment(segment_index);
          const SeqPackage::SeqView seq_view =
              seq_pkg_.GetSeqView(segment.read_id);
          const auto raw = seq_view.raw_address(segment.start);
          __builtin_prefetch(raw.first, 0, 1);
          const unsigned touched_bases = edge_length + segment.length - 1u;
          __builtin_prefetch(
              raw.first + (raw.second + touched_bases - 1u) /
                              kCharsPerEdgeWord,
              0, 1);
          return PreparedSegment(segment, seq_view);
        };

        if (segment_prefetch_distance == 0) {
          for (uint64_t segment_index = segment_begin;
               segment_index < segment_end; ++segment_index) {
            const CompactCountSegment segment = decode_segment(segment_index);
            const SeqPackage::SeqView seq_view =
                seq_pkg_.GetSeqView(segment.read_id);
            materialize_suffix64(segment, seq_view, prefix_out,
                                 prefix_ends.data(), &touched_prefixes);
            prefix_out += segment.length;
          }
          continue;
        }

        prepared_segments.clear();
        uint64_t next_segment = segment_begin;
        while (next_segment < segment_end &&
               prepared_segments.size() < segment_prefetch_distance) {
          prepared_segments.emplace_back(prepare_segment(next_segment++));
        }
        size_t head = 0;
        size_t buffered = prepared_segments.size();
        while (buffered != 0) {
          const PreparedSegment current = prepared_segments[head];
          materialize_suffix64(current.segment, current.seq_view, prefix_out,
                               prefix_ends.data(), &touched_prefixes);
          prefix_out += current.segment.length;
          if (next_segment < segment_end) {
            prepared_segments[head] = prepare_segment(next_segment++);
          } else {
            --buffered;
          }
          if (++head == prepared_segments.size()) {
            head = 0;
          }
        }
      }
      if (use_suffix64) {
        assert(prefix_out ==
               reinterpret_cast<CountPrefix64Record *>(records) + n_records);
      } else {
        assert(out == records + record_words);
      }
      materialize_cpu[tid] += omp_get_wtime() - materialize_begin;

      if (use_suffix64) {
        const double partition_begin = omp_get_wtime();
        prefix_groups.clear();
        uint64_t running = 0;
        for (uint16_t prefix : touched_prefixes) {
          const uint64_t count = prefix_ends[prefix];
          prefix_cursors[prefix] = running;
          prefix_groups.push_back({prefix, running, running + count});
          running += count;
          prefix_ends[prefix] = running;
        }
        assert(running == n_records);

        CountPrefix64Record *prefix_records =
            reinterpret_cast<CountPrefix64Record *>(records);
        const size_t suffix_bytes =
            static_cast<size_t>(n_records) * sizeof(CountSuffix64Record);
        if (suffix_bytes > suffix_workspace_bytes) {
          const size_t new_capacity = suffix_bytes;
          void *memory = std::realloc(suffix_workspace, new_capacity);
          if (memory == nullptr) {
            std::free(suffix_workspace);
            std::free(record_workspace);
            xfatal("Failed to allocate suffix64 count workspace ({} bytes)\n",
                   new_capacity);
          }
          suffix_workspace = memory;
          suffix_workspace_bytes = new_capacity;
          BindMemoryPagesToNumaDomain(suffix_workspace,
                                      suffix_workspace_bytes, local_domain);
          AdviseHugePages(suffix_workspace, suffix_workspace_bytes);
        }

        // Histogram offsets are already known.  A single sequential read of
        // the 16-byte materialized stream and one compact scatter replaces an
        // in-place 16-byte permutation followed by a second compaction pass.
        for (uint64_t i = 0; i < n_records; ++i) {
          const uint16_t prefix = prefix_records[i].prefix;
          const uint64_t dest_index = prefix_cursors[prefix]++;
          CountSuffix64Record *compact =
              static_cast<CountSuffix64Record *>(suffix_workspace);
          compact[dest_index].suffix = prefix_records[i].suffix;
          std::memcpy(compact[dest_index].payload,
                      prefix_records[i].payload,
                      sizeof(compact[dest_index].payload));
        }
        for (uint16_t prefix : touched_prefixes) {
          prefix_ends[prefix] = 0;
        }
        touched_prefixes.clear();
        prefix_partition_cpu[tid] += omp_get_wtime() - partition_begin;

        const double sort_begin = omp_get_wtime();
        CountSuffix64Record *compact =
            static_cast<CountSuffix64Record *>(suffix_workspace);
        for (const PrefixGroup &group : prefix_groups) {
          if (group.end - group.begin > 1u) {
            kmlib::kmsort(compact + group.begin, compact + group.end);
          }
        }
        sort_cpu[tid] += omp_get_wtime() - sort_begin;

        const double postprocess_begin = omp_get_wtime();
        for (const PrefixGroup &group : prefix_groups) {
          postprocess_suffix64(
              static_cast<CountSuffix64Record *>(suffix_workspace),
              group.begin, group.end, group.prefix, tid);
        }
        postprocess_cpu[tid] += omp_get_wtime() - postprocess_begin;
      } else {
        const double sort_begin = omp_get_wtime();
        shard_sort(records, static_cast<int64_t>(n_records));
        sort_cpu[tid] += omp_get_wtime() - sort_begin;

        const double postprocess_begin = omp_get_wtime();
        Lv2Postprocess(0, static_cast<int64_t>(n_records), tid, records, shard);
        postprocess_cpu[tid] += omp_get_wtime() - postprocess_begin;
      }
      processed.fetch_add(n_records, std::memory_order_relaxed);
    }
    std::free(suffix_workspace);
    std::free(record_workspace);
  }
  const double consume_wall = omp_get_wtime() - consume_wall_begin;
  if (processed.load(std::memory_order_relaxed) != total_occurrences) {
    xfatal("Segmented count processed {} of {} records\n",
           processed.load(std::memory_order_relaxed), total_occurrences);
  }
  xinfo(
      "Segmented count profile: group={.4} s, consume={.4} s; CPU-s "
      "materialize={.3}, prefix/compact={.3}, sort={.3}, "
      "count/output={.3}\n",
      grouping_wall, consume_wall,
      std::accumulate(materialize_cpu.begin(), materialize_cpu.end(), 0.0),
      std::accumulate(prefix_partition_cpu.begin(),
                      prefix_partition_cpu.end(), 0.0),
      std::accumulate(sort_cpu.begin(), sort_cpu.end(), 0.0),
      std::accumulate(postprocess_cpu.begin(), postprocess_cpu.end(), 0.0));
  return true;
}

bool KmerCounter::RunSpecializedMainLoop() {
  if (!segmented_count_enabled_) {
    return false;
  }
#define DISPATCH_EXTERNAL_SEGMENTED_COUNT(N) \
  case N:                                    \
    return RunExternalSegmentedCountFor<N>()
  if (mapped_count_reads_) {
    switch (words_per_substr_) {
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(1);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(2);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(3);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(4);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(5);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(6);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(7);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(8);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(9);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(10);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(11);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(12);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(13);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(14);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(15);
      DISPATCH_EXTERNAL_SEGMENTED_COUNT(16);
      default:
        return false;
    }
  }
#undef DISPATCH_EXTERNAL_SEGMENTED_COUNT
#define DISPATCH_SEGMENTED_COUNT(N) \
  case N:                           \
    return RunSegmentedCountFor<N>()
  switch (words_per_substr_) {
    DISPATCH_SEGMENTED_COUNT(1);
    DISPATCH_SEGMENTED_COUNT(2);
    DISPATCH_SEGMENTED_COUNT(3);
    DISPATCH_SEGMENTED_COUNT(4);
    DISPATCH_SEGMENTED_COUNT(5);
    DISPATCH_SEGMENTED_COUNT(6);
    DISPATCH_SEGMENTED_COUNT(7);
    DISPATCH_SEGMENTED_COUNT(8);
    DISPATCH_SEGMENTED_COUNT(9);
    DISPATCH_SEGMENTED_COUNT(10);
    DISPATCH_SEGMENTED_COUNT(11);
    DISPATCH_SEGMENTED_COUNT(12);
    DISPATCH_SEGMENTED_COUNT(13);
    DISPATCH_SEGMENTED_COUNT(14);
    DISPATCH_SEGMENTED_COUNT(15);
    DISPATCH_SEGMENTED_COUNT(16);
    default:
      return false;
  }
#undef DISPATCH_SEGMENTED_COUNT
}

template <unsigned NWords>
void KmerCounter::Lv0CalcBucketSizeFor(
    int64_t seq_from, int64_t seq_to,
    std::array<int64_t, kNumBuckets> *out) {
  const unsigned edge_length = opt_.k + 1;
  const unsigned max_edge_length = Kmer<NWords, uint32_t>::max_size();
  if (UNLIKELY(edge_length == 0 || edge_length > max_edge_length)) {
    xfatal("Invalid edge length: {}\n", edge_length);
  }
  std::array<int64_t, kNumBuckets> &bucket_sizes = *out;
  std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);
  Kmer<NWords, uint32_t> edge, rev_edge;  // (k+1)-mer and its rc

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < edge_length) {
      continue;
    }

    auto start_ptr_and_offset = seq_view.raw_address();
    edge.InitFromPtr(start_ptr_and_offset.first, start_ptr_and_offset.second,
                     edge_length);
    rev_edge = edge;
    rev_edge.ReverseComplement(edge_length);

    unsigned last_char_offset = opt_.k;

    while (true) {
      if (rev_edge.cmp(edge, edge_length) < 0) {
        bucket_sizes[rev_edge.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      } else {
        bucket_sizes[edge.data()[0] >>
                     (kCharsPerEdgeWord - kBucketPrefixLength) *
                         kBitsPerEdgeChar]++;
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, edge_length);
        rev_edge.ShiftPreappend(3 - c, edge_length);
      }
    }
  }
}

void KmerCounter::Lv0CalcBucketSize(
    int64_t seq_from, int64_t seq_to,
    std::array<int64_t, kNumBuckets> *out) {
#define DISPATCH_KMER_WORDS(N)                         \
  case N:                                              \
    return Lv0CalcBucketSizeFor<N>(seq_from, seq_to, out)
  switch (words_per_substr_) {
    DISPATCH_KMER_WORDS(1);
    DISPATCH_KMER_WORDS(2);
    DISPATCH_KMER_WORDS(3);
    DISPATCH_KMER_WORDS(4);
    DISPATCH_KMER_WORDS(5);
    DISPATCH_KMER_WORDS(6);
    DISPATCH_KMER_WORDS(7);
    DISPATCH_KMER_WORDS(8);
    DISPATCH_KMER_WORDS(9);
    DISPATCH_KMER_WORDS(10);
    DISPATCH_KMER_WORDS(11);
    DISPATCH_KMER_WORDS(12);
    DISPATCH_KMER_WORDS(13);
    DISPATCH_KMER_WORDS(14);
    DISPATCH_KMER_WORDS(15);
    DISPATCH_KMER_WORDS(16);
    default:
      xfatal("Invalid number of k-mer words: {}\n", words_per_substr_);
  }
#undef DISPATCH_KMER_WORDS
}

template <unsigned NWords, bool PackedOffsets>
void KmerCounter::Lv1FillOffsetsFor(OffsetFiller &filler, int64_t seq_from,
                                    int64_t seq_to) {
  const unsigned edge_length = opt_.k + 1;
  const unsigned max_edge_length = Kmer<NWords, uint32_t>::max_size();
  if (UNLIKELY(edge_length == 0 || edge_length > max_edge_length)) {
    xfatal("Invalid edge length: {}\n", edge_length);
  }
  Kmer<NWords, uint32_t> edge, rev_edge;  // (k+1)-mer and its rc
  unsigned key;
  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    auto read_length = seq_view.length();

    if (read_length < edge_length) {
      continue;
    }

    auto ptr_and_offset = seq_view.raw_address();
    edge.InitFromPtr(ptr_and_offset.first, ptr_and_offset.second, edge_length);
    rev_edge = edge;
    rev_edge.ReverseComplement(edge_length);

    // shift the key char by char
    unsigned last_char_offset = opt_.k;
    int64_t encoded_offset = PackedOffsets
                                 ? static_cast<int64_t>(
                                       static_cast<uint64_t>(read_id)
                                       << read_locator_shift_)
                                 : static_cast<int64_t>(
                                       seq_view.full_offset_in_pkg() << 1u);

    while (true) {
      if (rev_edge.cmp(edge, edge_length) < 0) {
        key = rev_edge.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(key, encoded_offset | 1u);
        }
      } else {
        key = edge.data()[0] >>
              (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
        if (filler.IsHandling(key)) {
          filler.WriteNextOffset(key, encoded_offset);
        }
      }

      if (++last_char_offset >= read_length) {
        break;
      } else {
        encoded_offset += 2;
        int c = seq_view.base_at(last_char_offset);
        edge.ShiftAppend(c, edge_length);
        rev_edge.ShiftPreappend(3 - c, edge_length);
      }
    }
  }
}

template <unsigned NWords, bool PackedOffsets, bool CompactItems>
void KmerCounter::Lv1FillDirectItemsFor(OffsetFiller &filler,
                                        int64_t seq_from, int64_t seq_to) {
  const unsigned edge_length = opt_.k + 1;
  Kmer<NWords, uint32_t> edge, rev_edge;

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    const unsigned read_length = seq_view.length();
    if (read_length < edge_length) {
      continue;
    }

    auto raw = seq_view.raw_address();
    edge.InitFromPtr(raw.first, raw.second, edge_length);
    rev_edge = edge;
    rev_edge.ReverseComplement(edge_length);

    int64_t encoded_offset =
        PackedOffsets
            ? static_cast<int64_t>(static_cast<uint64_t>(read_id)
                                   << read_locator_shift_)
            : static_cast<int64_t>(seq_view.full_offset_in_pkg() << 1u);
    unsigned offset = 0;
    unsigned last_char_offset = opt_.k;

    while (true) {
      const bool reverse = rev_edge.cmp(edge, edge_length) < 0;
      const auto &canonical = reverse ? rev_edge : edge;
      const unsigned key =
          canonical.data()[0] >>
          ((kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar);

      if (filler.IsHandling(key)) {
        uint32_t *item = filler.ReserveNextItem(key);

        const unsigned prev =
            offset == 0 ? kSentinelValue : seq_view.base_at(offset - 1);
        const unsigned next =
            offset + edge_length < read_length
                ? seq_view.base_at(offset + edge_length)
                : kSentinelValue;
        const uint64_t locator =
            static_cast<uint64_t>(encoded_offset) | unsigned(reverse);
        uint64_t read_info;
        if (!reverse) {
          read_info = (locator << 6u) | (prev << 3u) | next;
        } else {
          const unsigned canonical_prev =
              next == kSentinelValue ? kSentinelValue : 3u - next;
          const unsigned canonical_next =
              prev == kSentinelValue ? kSentinelValue : 3u - prev;
          read_info = (locator << 6u) | (canonical_prev << 3u) |
                      canonical_next;
        }
        if (CompactItems) {
          assert(locator <= std::numeric_limits<uint32_t>::max());
          assert((canonical.data()[NWords - 1] & 0x3Fu) == 0);
          std::memcpy(item, canonical.data(),
                      (NWords - 1) * sizeof(uint32_t));
          item[NWords - 1] =
              canonical.data()[NWords - 1] |
              static_cast<uint32_t>(read_info & 0x3Fu);
          item[NWords] = static_cast<uint32_t>(locator);
        } else {
          std::memcpy(item, canonical.data(), NWords * sizeof(uint32_t));
          DecomposeUint64(item + NWords, read_info);
        }
      }

      if (++last_char_offset >= read_length) {
        break;
      }
      ++offset;
      encoded_offset += 2;
      const int c = seq_view.base_at(last_char_offset);
      edge.ShiftAppend(c, edge_length);
      rev_edge.ShiftPreappend(3 - c, edge_length);
    }
  }
}

template <unsigned NWords, bool PackedOffsets>
void KmerCounter::Lv1FillWideCompactItemsFor(OffsetFiller &filler,
                                             int64_t seq_from,
                                             int64_t seq_to) {
  const unsigned edge_length = opt_.k + 1;
  Kmer<NWords, uint32_t> edge, rev_edge;

  for (int64_t read_id = seq_from; read_id < seq_to; ++read_id) {
    auto seq_view = seq_pkg_.GetSeqView(read_id);
    const unsigned read_length = seq_view.length();
    if (read_length < edge_length) {
      continue;
    }

    auto raw = seq_view.raw_address();
    edge.InitFromPtr(raw.first, raw.second, edge_length);
    rev_edge = edge;
    rev_edge.ReverseComplement(edge_length);

    int64_t encoded_offset =
        PackedOffsets
            ? static_cast<int64_t>(static_cast<uint64_t>(read_id)
                                   << read_locator_shift_)
            : static_cast<int64_t>(seq_view.full_offset_in_pkg() << 1u);
    unsigned offset = 0;
    unsigned last_char_offset = opt_.k;

    while (true) {
      const bool reverse = rev_edge.cmp(edge, edge_length) < 0;
      const auto &canonical = reverse ? rev_edge : edge;
      const unsigned key =
          canonical.data()[0] >>
          ((kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar);

      if (filler.IsHandling(key)) {
        uint32_t *item = filler.ReserveNextItem(key);

        const unsigned prev =
            offset == 0 ? kSentinelValue : seq_view.base_at(offset - 1);
        const unsigned next =
            offset + edge_length < read_length
                ? seq_view.base_at(offset + edge_length)
                : kSentinelValue;
        const uint64_t locator =
            static_cast<uint64_t>(encoded_offset) | unsigned(reverse);
        uint32_t context;
        if (!reverse) {
          context = (prev << 3u) | next;
        } else {
          const unsigned canonical_prev =
              next == kSentinelValue ? kSentinelValue : 3u - next;
          const unsigned canonical_next =
              prev == kSentinelValue ? kSentinelValue : 3u - prev;
          context = (canonical_prev << 3u) | canonical_next;
        }
        const uint32_t metadata =
            (static_cast<uint32_t>(locator >> 32u) << 6u) | context;
        assert((metadata & ~compact_key_mask_) == 0);
        assert((canonical.data()[NWords - 1] & compact_key_mask_) == 0);
        std::memcpy(item, canonical.data(),
                    (NWords - 1) * sizeof(uint32_t));
        item[NWords - 1] = canonical.data()[NWords - 1] | metadata;
        item[NWords] = static_cast<uint32_t>(locator);
      }

      if (++last_char_offset >= read_length) {
        break;
      }
      ++offset;
      encoded_offset += 2;
      const int c = seq_view.base_at(last_char_offset);
      edge.ShiftAppend(c, edge_length);
      rev_edge.ShiftPreappend(3 - c, edge_length);
    }
  }
}

void KmerCounter::Lv1FillOffsets(OffsetFiller &filler, int64_t seq_from,
                                 int64_t seq_to) {
#define DISPATCH_KMER_WORDS(N)                                      \
  case N:                                                           \
    if (UsingDirectLv1Items()) {                                    \
      if (compact_items_) {                                         \
        if (packed_read_offsets_) {                                 \
          if (compact_locator_high_mask_ != 0) {                    \
            return Lv1FillWideCompactItemsFor<N, true>(             \
                filler, seq_from, seq_to);                          \
          }                                                         \
          return Lv1FillDirectItemsFor<N, true, true>(              \
              filler, seq_from, seq_to);                            \
        }                                                           \
        return Lv1FillDirectItemsFor<N, false, true>(               \
            filler, seq_from, seq_to);                              \
      }                                                             \
      if (packed_read_offsets_) {                                   \
        return Lv1FillDirectItemsFor<N, true, false>(               \
            filler, seq_from, seq_to);                              \
      }                                                             \
      return Lv1FillDirectItemsFor<N, false, false>(                \
          filler, seq_from, seq_to);                                \
    }                                                               \
    if (packed_read_offsets_) {                                     \
      return Lv1FillOffsetsFor<N, true>(filler, seq_from, seq_to);  \
    }                                                               \
    return Lv1FillOffsetsFor<N, false>(filler, seq_from, seq_to)
  switch (words_per_substr_) {
    DISPATCH_KMER_WORDS(1);
    DISPATCH_KMER_WORDS(2);
    DISPATCH_KMER_WORDS(3);
    DISPATCH_KMER_WORDS(4);
    DISPATCH_KMER_WORDS(5);
    DISPATCH_KMER_WORDS(6);
    DISPATCH_KMER_WORDS(7);
    DISPATCH_KMER_WORDS(8);
    DISPATCH_KMER_WORDS(9);
    DISPATCH_KMER_WORDS(10);
    DISPATCH_KMER_WORDS(11);
    DISPATCH_KMER_WORDS(12);
    DISPATCH_KMER_WORDS(13);
    DISPATCH_KMER_WORDS(14);
    DISPATCH_KMER_WORDS(15);
    DISPATCH_KMER_WORDS(16);
    default:
      xfatal("Invalid number of k-mer words: {}\n", words_per_substr_);
  }
#undef DISPATCH_KMER_WORDS
}

template <bool PackedOffsets, bool CompactItems, bool WideCompactLocator>
void KmerCounter::Lv2ExtractSubStringFor(OffsetFetcher &fetcher,
                                         SubstrPtr substr_ptr) {
  // Offsets inside one bucket are monotone, but the corresponding packed read
  // words are still far enough apart that a core otherwise serializes on the
  // first cache miss of almost every item.  Walk a private copy of the offset
  // stream ahead of the real consumer and bring the small substring footprint
  // in while the current item is being decoded/copied.  The lookahead fetcher
  // is read-only with respect to the engine, so it cannot change item order or
  // differential-offset state.  Keep the full-offset fallback untouched.
  static const unsigned kPrefetchDistance = []() {
    const char *value = std::getenv("MEGAHIT_COUNT_READ_PREFETCH_DISTANCE");
    if (value == nullptr) {
      return 64u;
    }
    char *end = nullptr;
    const unsigned long parsed = std::strtoul(value, &end, 10);
    return end != value && *end == '\0' && parsed <= 64 ?
               static_cast<unsigned>(parsed) :
               64u;
  }();
  static const bool prefetch_packed_reads =
      PackedOffsets && kPrefetchDistance != 0 &&
      std::getenv("MEGAHIT_DISABLE_COUNT_READ_PREFETCH") == nullptr;
  struct StagedRead {
    StagedRead(uint64_t value, const uint32_t *words, unsigned word_offset,
               unsigned length, unsigned read_offset, unsigned read_strand)
        : locator(value),
          read_p(words),
          start_offset(word_offset),
          read_length(length),
          offset(read_offset),
          strand(read_strand) {}
    uint64_t locator;
    const uint32_t *read_p;
    unsigned start_offset;
    unsigned read_length;
    unsigned offset;
    unsigned strand;
  };
  auto decode_read = [&](uint64_t locator) {
    const auto decoded = ReadLocatorDecoder<PackedOffsets>::Decode(
        locator, read_locator_shift_, read_offset_mask_, seq_pkg_);
    const auto raw = decoded.seq_view.raw_address();
    return StagedRead(locator, raw.first, raw.second,
                      decoded.seq_view.length(), decoded.offset,
                      decoded.strand);
  };
  auto prefetch_read_words = [](const StagedRead &read, unsigned k) {
    const unsigned begin_char = read.start_offset + read.offset;
    const unsigned first_char = begin_char - (read.offset != 0);
    const unsigned last_char = begin_char + k;
    __builtin_prefetch(read.read_p + first_char / kCharsPerEdgeWord, 0, 1);
    __builtin_prefetch(read.read_p + last_char / kCharsPerEdgeWord, 0, 1);
  };

  auto extract_one = [&](const StagedRead &read) {
    const uint64_t locator = read.locator;
    const unsigned strand = read.strand;
    const unsigned offset = read.offset;
    unsigned num_chars_to_copy = opt_.k + 1;

    const unsigned read_length = read.read_length;
    const unsigned start_offset = read.start_offset;
    const unsigned words_this_seq =
        DivCeiling(start_offset + read_length, kCharsPerEdgeWord);
    const uint32_t *read_p = read.read_p;
    const auto base_at = [read_p, start_offset](unsigned index) {
      const unsigned packed_index = start_offset + index;
      const unsigned shift =
          kBitsPerEdgeWord -
          (packed_index % kCharsPerEdgeWord + 1u) * kBitsPerEdgeChar;
      return static_cast<unsigned>(
          (read_p[packed_index / kCharsPerEdgeWord] >> shift) & 3u);
    };

    unsigned char prev, next;

    if (offset > 0) {
      prev = base_at(offset - 1);
    } else {
      prev = kSentinelValue;
    }

    if (offset + opt_.k + 1 < read_length) {
      next = base_at(offset + opt_.k + 1);
    } else {
      next = kSentinelValue;
    }
    uint64_t read_info;
    if (strand == 0) {
      CopySubstring(substr_ptr, read_p, offset + start_offset,
                    num_chars_to_copy, 1, words_this_seq, words_per_substr_);
      read_info = (locator << 6u) | (prev << 3) | next;
    } else {
      CopySubstringRC(substr_ptr, read_p, offset + start_offset,
                      num_chars_to_copy, 1, words_this_seq,
                      words_per_substr_);
      read_info =
          (locator << 6u) |
          ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3) |
          (prev == kSentinelValue ? kSentinelValue : (3 - prev));
    }
    if (CompactItems) {
      uint32_t metadata = static_cast<uint32_t>(read_info & 0x3Fu);
      if (WideCompactLocator) {
        metadata |= static_cast<uint32_t>(locator >> 32u) << 6u;
      }
      assert((metadata & ~compact_key_mask_) == 0);
      assert((substr_ptr[words_per_substr_ - 1] & compact_key_mask_) == 0);
      substr_ptr[words_per_substr_ - 1] |= metadata;
      substr_ptr[words_per_substr_] = static_cast<uint32_t>(locator);
      substr_ptr += words_per_substr_ + 1;
    } else {
      DecomposeUint64(substr_ptr + words_per_substr_, read_info);
      substr_ptr += words_per_substr_ + 2;
    }
  };

  if (prefetch_packed_reads) {
    // The old look-ahead copied OffsetFetcher and decoded every locator twice:
    // once to issue the prefetch and again when materializing the record.
    // Keep a small bounded ring of already-decoded reads instead.  This is a
    // true software pipeline: input order is unchanged, each locator and read
    // boundary is decoded exactly once, and the staged view is consumed after
    // its data has had kPrefetchDistance iterations to reach cache.
    std::vector<StagedRead> staged;
    staged.reserve(kPrefetchDistance);
    while (staged.size() < kPrefetchDistance && fetcher.HasNext()) {
      const uint64_t locator = fetcher.Next();
      const auto decoded = decode_read(locator);
      prefetch_read_words(decoded, opt_.k);
      staged.emplace_back(decoded);
    }

    size_t head = 0;
    size_t buffered = staged.size();
    while (buffered != 0) {
      // Copy before replacing this ring slot with the next look-ahead item.
      const StagedRead current = staged[head];
      extract_one(current);
      if (fetcher.HasNext()) {
        const uint64_t locator = fetcher.Next();
        const auto decoded = decode_read(locator);
        prefetch_read_words(decoded, opt_.k);
        staged[head] = decoded;
      } else {
        --buffered;
      }
      if (++head == staged.size()) {
        head = 0;
      }
    }
    return;
  }

  while (fetcher.HasNext()) {
    const uint64_t locator = fetcher.Next();
    extract_one(decode_read(locator));
  }
}

void KmerCounter::Lv2ExtractSubString(OffsetFetcher &fetcher,
                                      SubstrPtr substr_ptr) {
  if (packed_read_offsets_) {
    if (compact_items_) {
      if (compact_locator_high_mask_ != 0) {
        return Lv2ExtractSubStringFor<true, true, true>(fetcher, substr_ptr);
      }
      return Lv2ExtractSubStringFor<true, true, false>(fetcher, substr_ptr);
    }
    return Lv2ExtractSubStringFor<true, false, false>(fetcher, substr_ptr);
  }
  return Lv2ExtractSubStringFor<false, false, false>(fetcher, substr_ptr);
}

template <bool WideCompactLocator>
uint64_t KmerCounter::DecodeCompactLocator(const uint32_t *item) const {
  const uint64_t low = item[words_per_substr_];
  if (!WideCompactLocator) {
    return low;
  }
  const uint64_t high =
      (item[words_per_substr_ - 1] & compact_locator_high_mask_) >> 6u;
  return low | (high << 32u);
}

template <bool PackedOffsets, bool CompactItems, bool WideCompactLocator>
void KmerCounter::Lv2PostprocessFor(int64_t start_index, int64_t end_index,
                                    int thread_id, uint32_t *substr_ptr) {
  uint32_t packed_edge[32];
  int64_t count_prev[5], count_next[5];

  int64_t from_;
  int64_t to_;
  const int64_t item_words =
      words_per_substr_ + (CompactItems ? 1 : 2);

  EdgeWriter::Snapshot snapshot;

  for (int64_t i = start_index; i < end_index; i = to_) {
    from_ = i;
    to_ = i + 1;
    uint32_t *first_item = substr_ptr + i * item_words;

    while (to_ < end_index) {
      if (IsDifferentEdges<CompactItems>(
              first_item, substr_ptr + to_ * item_words, words_per_substr_,
              1, compact_key_mask_)) {
        break;
      }

      ++to_;
    }

    int64_t count = to_ - from_;

    // update read's first and last

    memset(count_prev, 0, sizeof(count_prev[0]) * 4);
    memset(count_next, 0, sizeof(count_next[0]) * 4);
    bool has_in = false;
    bool has_out = false;
    bool guaranteed_in = false;
    bool guaranteed_out = false;
    uint8_t verify_in_mask = 0u;
    uint8_t verify_out_mask = 0u;
    const int64_t unambiguous_context =
        std::max<int64_t>(1, opt_.solid_threshold) * 2;

    for (int64_t j = from_; j < to_; ++j) {
      auto *item = substr_ptr + j * item_words;
      const int prev_and_next =
          CompactItems
              ? static_cast<int>(item[words_per_substr_ - 1] & 0x3Fu)
              : static_cast<int>(
                    ComposeUint64(item + words_per_substr_) & 0x3Fu);
      count_prev[prev_and_next >> 3]++;
      count_next[prev_and_next & 7]++;
    }

    for (int j = 0; j < 4; ++j) {
      if (count_prev[j] >= opt_.solid_threshold) {
        has_in = true;
      }

      if (count_next[j] >= opt_.solid_threshold) {
        has_out = true;
      }
      guaranteed_in |= count_prev[j] >= unambiguous_context;
      guaranteed_out |= count_next[j] >= unambiguous_context;
      if (count_prev[j] >= opt_.solid_threshold &&
          count_prev[j] < unambiguous_context) {
        verify_in_mask |= static_cast<uint8_t>(1u << j);
      }
      if (count_next[j] >= opt_.solid_threshold &&
          count_next[j] < unambiguous_context) {
        verify_out_mask |= static_cast<uint8_t>(1u << j);
      }
    }

    // Missing incoming/outgoing support used to be handled by two separate
    // passes over the same equal-edge run.  In the common isolated-edge case
    // that decoded every read locator twice and revisited the same scattered
    // read-boundary arrays twice.  Both operations are independent numeric
    // min/max reductions, so decode once and update the requested endpoints
    // together.  The OpenMP region boundary publishes the completed arrays;
    // these atomics do not carry data-dependency ordering themselves.
    if ((!has_in || !has_out) && count >= opt_.solid_threshold) {
      for (int64_t j = from_; j < to_; ++j) {
        auto *item = substr_ptr + j * item_words;
        const uint64_t locator =
            CompactItems ? DecodeCompactLocator<WideCompactLocator>(item)
                          : ComposeUint64(item + words_per_substr_) >> 6u;
        unsigned strand;
        uint32_t offset;
        uint32_t read_id;
        if (PackedOffsets && mapped_count_reads_) {
          strand = static_cast<unsigned>(locator & 1u);
          offset = static_cast<uint32_t>((locator >> 1u) &
                                         read_offset_mask_);
          read_id = static_cast<uint32_t>(locator >> read_locator_shift_);
        } else {
          auto decoded = ReadLocatorDecoder<PackedOffsets>::Decode(
              locator, read_locator_shift_, read_offset_mask_, seq_pkg_);
          strand = decoded.strand;
          offset = decoded.offset;
          read_id = static_cast<uint32_t>(decoded.seq_view.id());
        }

        if ((!has_in && strand == 0) || (!has_out && strand != 0)) {
          UpdateLast0In(read_id, offset);
        }
        if ((!has_out && strand == 0) || (!has_in && strand != 0)) {
          const uint32_t first_offset = offset + 1;
          UpdateFirst0Out(read_id, first_offset);
        }
      }
    }
    edge_counter_.Add(count, thread_id);

    if (count >= opt_.solid_threshold) {
      PackEdge(packed_edge, first_item, count);
      WriteEndpointCandidate(packed_edge, !has_in, !has_out,
                             guaranteed_in ? 0u : verify_in_mask,
                             guaranteed_out ? 0u : verify_out_mask,
                             thread_id);
      if (segmented_count_enabled_) {
        edge_writer_.WriteUnordered(packed_edge, thread_id);
        uint64_t *histogram = seq_bucket_histograms_.data() +
                              static_cast<size_t>(thread_id) * kNumBuckets;
        const uint16_t key = ExtractCountEdgeWindow8(packed_edge, 1u);
        const unsigned reverse_window_offset =
            opt_.k - kBucketPrefixLength;
        const uint16_t reverse_key = ReverseComplementCountWindow8(
            ExtractCountEdgeWindow8(packed_edge, reverse_window_offset));
        ++histogram[key];
        ++histogram[reverse_key];
      } else {
        edge_writer_.Write(packed_edge,
                           packed_edge[0] >> (32 - 2 * kBucketPrefixLength),
                           thread_id, &snapshot);
      }
    }
  }

  edge_writer_.SaveSnapshot(snapshot);
}

void KmerCounter::Lv2Postprocess(int64_t start_index, int64_t end_index,
                                 int thread_id, uint32_t *substr_ptr,
                                 unsigned bucket_id) {
  (void)bucket_id;
  const bool compact_items = compact_items_;
  if (packed_read_offsets_) {
    if (compact_items) {
      if (compact_locator_high_mask_ != 0) {
        return Lv2PostprocessFor<true, true, true>(
            start_index, end_index, thread_id, substr_ptr);
      }
      return Lv2PostprocessFor<true, true, false>(
          start_index, end_index, thread_id, substr_ptr);
    }
    return Lv2PostprocessFor<true, false, false>(
        start_index, end_index, thread_id, substr_ptr);
  }
  return Lv2PostprocessFor<false, false, false>(
      start_index, end_index, thread_id, substr_ptr);
}

void KmerCounter::Lv0Postprocess() {
  uint64_t endpoint_candidate_bytes = 0;
  for (int thread = 0; thread < opt_.n_threads; ++thread) {
    if (thread < static_cast<int>(endpoint_candidate_files_.size()) &&
        endpoint_candidate_files_[thread]) {
      endpoint_candidate_files_[thread]->flush();
      const std::streampos end = endpoint_candidate_files_[thread]->tellp();
      if (end >= std::streampos{0}) {
        endpoint_candidate_bytes += static_cast<uint64_t>(end);
      }
      endpoint_candidate_files_[thread]->close();
    }
  }
  endpoint_candidate_files_.clear();
  xinfo("Endpoint candidate stream: {} bytes in {} shards\n",
        endpoint_candidate_bytes, opt_.n_threads);

  // --- output reads for mercy ---
  int64_t num_candidate_reads = 0;
  int64_t num_has_tips = 0;
  std::ofstream candidate_file(opt_.output_prefix + ".cand",
                               std::ofstream::binary | std::ofstream::out);

  const size_t num_reads =
      mapped_count_reads_ ? static_cast<size_t>(mapped_count_reads_->num_reads())
                          : seq_pkg_.seq_count();
  constexpr size_t kParallelCandidateThreshold = size_t{1} << 18u;
  if (mapped_count_reads_ ||
      (num_reads >= kParallelCandidateThreshold && omp_get_max_threads() > 1)) {
    // Classification is embarrassingly parallel, but the candidate stream is
    // required to remain in ascending read-ID order.  Materialize only one bit
    // per read, then visit set bits in word order.  This removes the serial
    // all-read scan without changing the byte order of the sparse output.
    constexpr size_t kBitsPerWord = 64;
    std::vector<uint64_t> candidate_words(
        DivCeiling(num_reads, kBitsPerWord), 0);
#pragma omp parallel for reduction(+ : num_candidate_reads, num_has_tips) schedule(static)
    for (size_t word_idx = 0; word_idx < candidate_words.size(); ++word_idx) {
      const size_t begin = word_idx * kBitsPerWord;
      const size_t end = std::min(begin + kBitsPerWord, num_reads);
      uint64_t word = 0;
      for (size_t read_id = begin; read_id < end; ++read_id) {
        const auto first = LoadFirst0Out(static_cast<uint32_t>(read_id));
        const auto last = LoadLast0In(static_cast<uint32_t>(read_id));
        if (first != kSentinelOffset && last != kSentinelOffset) {
          ++num_has_tips;
          if (last > first) {
            ++num_candidate_reads;
            word |= uint64_t{1} << (read_id - begin);
          }
        }
      }
      candidate_words[word_idx] = word;
    }

    if (mapped_count_reads_) {
      if (!mapped_count_reads_->WriteSelectedReversed(&candidate_file,
                                                       candidate_words)) {
        xfatal("Failed to write mapped mercy candidate reads\n");
      }
    } else {
      for (size_t word_idx = 0; word_idx < candidate_words.size();
           ++word_idx) {
        uint64_t word = candidate_words[word_idx];
        while (word != 0) {
          const unsigned bit = static_cast<unsigned>(__builtin_ctzll(word));
          const size_t read_id = word_idx * kBitsPerWord + bit;
          seq_pkg_.WriteSequences(candidate_file, read_id, read_id);
          word &= word - 1;
        }
      }
    }
  } else {
    for (size_t i = 0; i < num_reads; ++i) {
      const auto first = LoadFirst0Out(static_cast<uint32_t>(i));
      const auto last = LoadLast0In(static_cast<uint32_t>(i));

      if (first != kSentinelOffset && last != kSentinelOffset) {
        ++num_has_tips;

        if (last > first) {
          ++num_candidate_reads;
          seq_pkg_.WriteSequences(candidate_file, i, i);
        }
      }
    }
  }

  xinfo("Total number of candidate reads: {} ({})\n", num_candidate_reads,
        num_has_tips);
  xinfo("Total number of solid edges: {}\n",
        edge_counter_.GetNumSolidEdges(opt_.solid_threshold));
  std::ofstream counting_file(opt_.output_prefix + ".counting");
  edge_counter_.DumpStat(counting_file);

  // --- cleaning ---
  edge_writer_.Finalize();
  if (!seq_bucket_histograms_.empty()) {
    EdgeBucketHistogram histogram;
    histogram.kmer_size = opt_.k;
    histogram.num_files = opt_.n_threads;
    histogram.counts.swap(seq_bucket_histograms_);
    const uint64_t total_records = std::accumulate(
        histogram.counts.begin(), histogram.counts.end(), uint64_t{0});
    if ((total_records & 1u) != 0) {
      xfatal("Odd count->SDBG histogram total: {}\n", total_records);
    }
    histogram.num_edges = total_records / 2u;
    if (!histogram.Write(opt_.output_prefix)) {
      xwarn("Failed to write count->SDBG bucket histogram; the next stage "
            "will rescan edges\n");
    } else {
      xinfo("Wrote count->SDBG bucket histogram for {} edges in {} shards\n",
            histogram.num_edges, histogram.num_files);
    }
  }
}
