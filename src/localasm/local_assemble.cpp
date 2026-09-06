#include "local_assemble.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <omp.h>
#include "idba/contig_graph.h"
#include "idba/hash_graph.h"
#include "idba/sequence.h"
#include "kmlib/kmbit.h"
#include "kmlib/kmsort.h"

#include "hash_mapper.h"
#include "mapping_result_collector.h"
#include "sequence/io/contig/contig_reader.h"
#include "sequence/io/contig/contig_writer.h"
#include "sequence/io/local_candidate_index.h"
#include "sequence/io/local_seed_positions.h"
#include "sequence/io/read_chunk_index.h"
#include "sequence/io/sequence_lib.h"
#include "utils/histgram.h"
#include "utils/utils.h"

namespace {

static const int kMaxLocalRange = 650;
using TInsertSize = std::pair<double, double>;

struct PackedReadRecord {
  const uint32_t *words;
  uint32_t length;
};

/**
 * Read-only view of buildlib's indexed binary stream.  Mapping works directly
 * on word-aligned records; no all-library SeqPackage is constructed.  The
 * virtual mapping may be arbitrarily large, while completed chunk pages are
 * discarded so resident memory follows active workers rather than input size.
 */
class MappedReadFile {
 public:
  MappedReadFile() = default;
  ~MappedReadFile() { Close(); }
  MappedReadFile(const MappedReadFile &) = delete;
  MappedReadFile &operator=(const MappedReadFile &) = delete;

  bool Open(const std::string &lib_prefix,
            const SequenceLibCollection::SizeInfo &expected) {
    Close();
    path_ = lib_prefix + ".bin";
    if (!LoadPackedReadChunkIndex(path_, &index_) ||
        index_.num_reads != static_cast<uint64_t>(expected.num_reads) ||
        index_.num_bases != static_cast<uint64_t>(expected.num_bases) ||
        index_.max_read_len != expected.max_read_len) {
      return false;
    }
    fd_ = open(path_.c_str(), O_RDONLY);
    struct stat status;
    if (fd_ < 0 || fstat(fd_, &status) != 0 || status.st_size <= 0 ||
        status.st_size % static_cast<off_t>(sizeof(uint32_t)) != 0) {
      Close();
      return false;
    }
    bytes_ = static_cast<size_t>(status.st_size);
    void *address = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (address == MAP_FAILED) {
      words_ = nullptr;
      Close();
      return false;
    }
    words_ = static_cast<const uint32_t *>(address);
    return true;
  }

  void Close() {
    if (words_ != nullptr) {
      munmap(const_cast<uint32_t *>(words_), bytes_);
      words_ = nullptr;
    }
    if (fd_ >= 0) {
      close(fd_);
      fd_ = -1;
    }
    bytes_ = 0;
    index_ = PackedReadChunkIndex();
  }

  const PackedReadChunkIndex &index() const { return index_; }
  const uint32_t *words() const { return words_; }
  size_t total_words() const { return bytes_ / sizeof(uint32_t); }

  const uint32_t *LocateRead(const PackedReadChunk &chunk,
                             uint64_t read_id) const {
    assert(read_id >= chunk.read_begin && read_id < chunk.read_end);
    const uint32_t *cursor = words_ + chunk.word_begin;
    for (uint64_t id = chunk.read_begin; id < read_id; ++id) {
      const uint32_t len = *cursor++;
      cursor += DivCeiling(static_cast<size_t>(len),
                           SeqPackage::kBasesPerWord);
    }
    return cursor;
  }

  static PackedReadRecord Next(const uint32_t **cursor) {
    const uint32_t len = *(*cursor)++;
    const uint32_t *words = *cursor;
    *cursor += DivCeiling(static_cast<size_t>(len),
                          SeqPackage::kBasesPerWord);
    return {words, len};
  }

  void DropChunk(const PackedReadChunk &chunk) const {
    DiscardMemoryPages(
        const_cast<uint32_t *>(words_ + chunk.word_begin),
        (chunk.word_end - chunk.word_begin) * sizeof(uint32_t));
  }

 private:
  std::string path_;
  int fd_{-1};
  size_t bytes_{0};
  const uint32_t *words_{nullptr};
  PackedReadChunkIndex index_;
};

class MappedLocalCandidateFile {
 public:
  MappedLocalCandidateFile() = default;
  ~MappedLocalCandidateFile() { Close(); }
  MappedLocalCandidateFile(const MappedLocalCandidateFile &) = delete;
  MappedLocalCandidateFile &operator=(const MappedLocalCandidateFile &) =
      delete;

  bool Open(const std::string &path, uint64_t expected_words,
            uint64_t expected_reads) {
    Close();
    fd_ = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    struct stat status {};
    if (fd_ < 0 || fstat(fd_, &status) != 0 ||
        status.st_size <
            static_cast<off_t>(sizeof(LocalCandidateFileHeader))) {
      Close();
      return false;
    }
    bytes_ = static_cast<size_t>(status.st_size);
    void *address = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (address == MAP_FAILED) {
      mapping_ = nullptr;
      Close();
      return false;
    }
    mapping_ = static_cast<const uint8_t *>(address);
    const auto *header =
        reinterpret_cast<const LocalCandidateFileHeader *>(mapping_);
    if (!IsValidLocalCandidateFileHeader(*header) ||
        header->source_words != expected_words ||
        header->source_reads != expected_reads ||
        header->candidate_count >
            (std::numeric_limits<size_t>::max() -
             sizeof(LocalCandidateFileHeader)) /
                sizeof(uint64_t) ||
        sizeof(LocalCandidateFileHeader) +
                header->candidate_count * sizeof(uint64_t) !=
            bytes_) {
      Close();
      return false;
    }
    candidates_ = reinterpret_cast<const uint64_t *>(
        mapping_ + sizeof(LocalCandidateFileHeader));
    count_ = static_cast<size_t>(header->candidate_count);
    for (size_t i = 0; i < count_; ++i) {
      if (candidates_[i] >= expected_words ||
          (i != 0u && candidates_[i - 1u] >= candidates_[i])) {
        Close();
        return false;
      }
    }
    return true;
  }

  void Close() {
    if (mapping_ != nullptr) {
      munmap(const_cast<uint8_t *>(mapping_), bytes_);
    }
    if (fd_ >= 0) close(fd_);
    fd_ = -1;
    bytes_ = 0;
    mapping_ = nullptr;
    candidates_ = nullptr;
    count_ = 0;
  }

  std::pair<const uint64_t *, const uint64_t *> Range(
      uint64_t word_begin, uint64_t word_end) const {
    const uint64_t *begin =
        std::lower_bound(candidates_, candidates_ + count_, word_begin);
    const uint64_t *end =
        std::lower_bound(begin, candidates_ + count_, word_end);
    return std::make_pair(begin, end);
  }

  size_t size() const { return count_; }

 private:
  int fd_{-1};
  size_t bytes_{0};
  const uint8_t *mapping_{nullptr};
  const uint64_t *candidates_{nullptr};
  size_t count_{0};
};


bool PackSequence2Bit(const Sequence &sequence,
                      std::vector<uint64_t> *packed) {
  packed->assign(((sequence.size() + 31u) >> 5u) + 1u, 0u);
  uint32_t i = 0;
#ifdef USE_BMI2
  constexpr uint64_t kTwoBitsPerByte = UINT64_C(0x0303030303030303);
  constexpr uint64_t kHighBitsPerByte = UINT64_C(0xFCFCFCFCFCFCFCFC);
  const uint8_t *const bases = sequence.encoded_data();
  for (; i + 8u <= sequence.size(); i += 8u) {
    uint64_t bytes;
    std::memcpy(&bytes, bases + i, sizeof(bytes));
    if ((bytes & kHighBitsPerByte) != 0u) return false;
    const uint64_t compact = _pext_u64(bytes, kTwoBitsPerByte);
    (*packed)[i >> 5u] |= compact << ((i & 31u) << 1u);
  }
#endif
  for (; i < sequence.size(); ++i) {
    const uint8_t base = sequence[i];
    if (base >= 4u) return false;
    (*packed)[i >> 5u] |= uint64_t(base) << ((i & 31u) << 1u);
  }
  return true;
}

// Endpoint-local reads in the same two-bit form used by the global read
// package.  Every record has a zero look-ahead word, so arbitrary k-mer
// windows can be loaded without a boundary branch.  The forward stream is
// built once and shared by all inner-k rounds; reverse state is derived only
// for each round's first window and then rolled.
class LocalPackedReads {
 public:
  struct ReadView {
    const uint64_t *forward;
    const uint64_t *reverse;
    uint32_t length;
  };

  void clear() {
    forward_.clear();
    reverse_.clear();
    offsets_.clear();
    lengths_.clear();
    max_length_ = 0;
  }

  void Add(const SeqPackage::SeqView &source) {
    const auto address = source.raw_address();
    Add(address.first, address.second, source.length());
  }

  void Add(const uint32_t *source_words, uint32_t source_shift,
             uint32_t length) {
    const uint32_t num_words = (length + 31u) >> 5u;
    const uint32_t offset = static_cast<uint32_t>(forward_.size());
    offsets_.push_back(offset);
    lengths_.push_back(length);
    max_length_ = std::max(max_length_, length);
    forward_.resize(offset + num_words + 1u, 0);

    // SeqPackage stores 16 two-bit bases per uint32 in big-endian lane order;
    // the local rolling kernels consume 32 bases per uint64 in little-endian
    // lane order. Convert a complete source word at a time instead of calling
    // base_at() and issuing one read/shift/OR for every base. The source view
    // may start mid-word, so join the two adjacent words before reversing the
    // two-bit lanes. The look-ahead word allocated above remains zero.
    for (uint32_t base_begin = 0; base_begin < length; base_begin += 16u) {
      const uint32_t take = std::min<uint32_t>(16u, length - base_begin);
      const uint32_t absolute = source_shift + base_begin;
      const uint32_t word = absolute >> 4u;
      const uint32_t lane = absolute & 15u;
      uint32_t packed = source_words[word] << (lane << 1u);
      if (lane + take > 16u) {
        packed |= source_words[word + 1u] >> ((16u - lane) << 1u);
      }
      if (take != 16u) {
        packed &= UINT32_MAX << ((16u - take) << 1u);
      }
      const uint32_t little = kmlib::bit::Reverse<2>(packed);
      forward_[offset + (base_begin >> 5u)] |=
          uint64_t(little) << ((base_begin & 31u) << 1u);
    }
  }

  ReadView operator[](size_t read_id) const {
    const uint32_t offset = offsets_[read_id];
    return {forward_.data() + offset,
            reverse_.empty() ? nullptr : reverse_.data() + offset,
            lengths_[read_id]};
  }

  void BuildReverse() {
    reverse_.assign(forward_.size(), 0);
    for (size_t read_id = 0; read_id < offsets_.size(); ++read_id) {
      const uint32_t offset = offsets_[read_id];
      const uint32_t length = lengths_[read_id];
      const uint32_t words = (length + 31u) >> 5u;
      for (uint32_t i = 0; i < words; ++i) {
        uint64_t value = forward_[offset + i];
        bit_operation::ReverseComplement(value);
        reverse_[offset + words - 1u - i] = value;
      }
      const uint32_t tail = length & 31u;
      if (tail != 0u) {
        const uint32_t shift = (32u - tail) << 1u;
        for (uint32_t i = 0; i + 1u < words; ++i) {
          reverse_[offset + i] =
              (reverse_[offset + i] >> shift) |
              (reverse_[offset + i + 1u] << (64u - shift));
        }
        reverse_[offset + words - 1u] >>= shift;
      }
    }
  }

  static uint8_t PackedBase(const uint64_t *words, uint32_t index) {
    return static_cast<uint8_t>(
        (words[index >> 5u] >> ((index & 31u) << 1u)) & 3u);
  }

  size_t size() const { return offsets_.size(); }
  uint32_t max_length() const { return max_length_; }

  uint64_t EstimateKmerWork(uint32_t mink, uint32_t maxk,
                            uint32_t step) const {
    if (step == 0) return 0;
    uint64_t total = 0;
    for (uint32_t length : lengths_) {
      if (length < mink) continue;
      const uint32_t last_k = std::min(length, maxk);
      const uint64_t rounds = (last_k - mink) / step + 1u;
      total += rounds * uint64_t(length - mink + 1u) -
               uint64_t(step) * rounds * (rounds - 1u) / 2u;
    }
    return total;
  }

  // Diagnostic for deciding whether exact read multiplicity compression is
  // worth adding to the hot graph builder.  Compare the packed strings, not
  // fingerprints, so the reported opportunity is collision-free.  This is
  // deliberately opt-in: sorting tiny endpoint-local ID arrays is useful for
  // measurement but must not tax the default path when duplicates are rare.
  std::pair<uint64_t, uint64_t> ExactDuplicateWork(
      uint32_t mink, uint32_t maxk, uint32_t step) const {
    if (lengths_.size() < 2u || step == 0u) {
      return std::make_pair(uint64_t{0}, uint64_t{0});
    }
    std::vector<uint32_t> order(lengths_.size());
    std::iota(order.begin(), order.end(), uint32_t{0});
    const auto less = [this](uint32_t lhs, uint32_t rhs) {
      if (lengths_[lhs] != lengths_[rhs]) {
        return lengths_[lhs] < lengths_[rhs];
      }
      const uint32_t words = (lengths_[lhs] + 31u) >> 5u;
      const uint64_t *lhs_words = forward_.data() + offsets_[lhs];
      const uint64_t *rhs_words = forward_.data() + offsets_[rhs];
      for (uint32_t word = 0; word < words; ++word) {
        if (lhs_words[word] != rhs_words[word]) {
          return lhs_words[word] < rhs_words[word];
        }
      }
      return lhs < rhs;
    };
    std::sort(order.begin(), order.end(), less);

    uint64_t duplicates = 0;
    uint64_t duplicate_work = 0;
    for (size_t i = 1; i < order.size(); ++i) {
      const uint32_t lhs = order[i - 1u];
      const uint32_t rhs = order[i];
      const uint32_t length = lengths_[rhs];
      if (length != lengths_[lhs]) continue;
      const uint32_t words = (length + 31u) >> 5u;
      if (std::memcmp(forward_.data() + offsets_[lhs],
                      forward_.data() + offsets_[rhs],
                      static_cast<size_t>(words) * sizeof(uint64_t)) != 0) {
        continue;
      }
      ++duplicates;
      if (length >= mink) {
        const uint32_t last_k = std::min(length, maxk);
        const uint64_t rounds = (last_k - mink) / step + 1u;
        duplicate_work += rounds * uint64_t(length - mink + 1u) -
                          uint64_t(step) * rounds * (rounds - 1u) / 2u;
      }
    }
    return std::make_pair(duplicates, duplicate_work);
  }

  std::pair<uint64_t, uint64_t> InputFingerprint(
      const Sequence &endpoint) const {
    // Diagnostic only: hash the exact endpoint bytes and the complete ordered
    // local-read representation (including record boundaries).  This lets us
    // measure how much work could be resumed across outer k rounds before a
    // persistent exact cache is introduced into the default path.
    const XXH128_hash_t endpoint_hash = XXH3_128bits(
        endpoint.encoded_data(), endpoint.size() * sizeof(uint8_t));
    const XXH128_hash_t bases_hash = XXH3_128bits(
        forward_.data(), forward_.size() * sizeof(uint64_t));
    const XXH128_hash_t lengths_hash = XXH3_128bits(
        lengths_.data(), lengths_.size() * sizeof(uint32_t));
    uint64_t low = endpoint_hash.low64 ^
                   (bases_hash.low64 + UINT64_C(0x9e3779b97f4a7c15)) ^
                   (lengths_hash.high64 << 1u) ^ uint64_t(lengths_.size());
    uint64_t high = endpoint_hash.high64 ^
                    (bases_hash.high64 + UINT64_C(0xd6e8feb86659fd93)) ^
                    (lengths_hash.low64 >> 1u) ^ uint64_t(forward_.size());
    low ^= high >> 30u;
    low *= UINT64_C(0xbf58476d1ce4e5b9);
    low ^= low >> 27u;
    high ^= low >> 31u;
    high *= UINT64_C(0x94d049bb133111eb);
    high ^= high >> 29u;
    return std::make_pair(low, high);
  }

  std::pair<uint64_t, uint64_t> ReadFingerprint() const {
    const XXH128_hash_t bases_hash = XXH3_128bits(
        forward_.data(), forward_.size() * sizeof(uint64_t));
    const XXH128_hash_t lengths_hash = XXH3_128bits(
        lengths_.data(), lengths_.size() * sizeof(uint32_t));
    uint64_t low = bases_hash.low64 ^
                   (lengths_hash.high64 + UINT64_C(0x9e3779b97f4a7c15)) ^
                   uint64_t(lengths_.size());
    uint64_t high = bases_hash.high64 ^
                    (lengths_hash.low64 + UINT64_C(0xd6e8feb86659fd93)) ^
                    uint64_t(forward_.size());
    low ^= high >> 30u;
    low *= UINT64_C(0xbf58476d1ce4e5b9);
    low ^= low >> 27u;
    high ^= low >> 31u;
    high *= UINT64_C(0x94d049bb133111eb);
    high ^= high >> 29u;
    return std::make_pair(low, high);
  }

 private:
  std::vector<uint64_t> forward_;
  std::vector<uint64_t> reverse_;
  std::vector<uint32_t> offsets_;
  std::vector<uint32_t> lengths_;
  uint32_t max_length_{0};
};

// A reusable variable-order occurrence index for one endpoint.  It builds a
// generalized suffix order over each selected read and its reverse complement
// once.  Every later k obtains exact equality groups from adjacent LCP values,
// then reduces occurrences in historical read/position order.  The suffix
// order is only an equality index: canonical orientation still uses IDBA's
// exact packed-word comparison, so it cannot alter graph or tie semantics.
class MultiKReadIndex {
 public:
  struct Stats {
    double build{0};
    double finish_initial{0};
    double advance{0};
    double aggregate{0};
    double replay{0};
  };

  void Build(LocalPackedReads &reads, uint32_t min_k) {
    const double begin = omp_get_wtime();
    stats_ = Stats();
    reads.BuildReverse();
    reads_ = &reads;
    current_k_ = min_k;
    group_count_ = 0;

    oriented_offsets_.assign(reads.size() * 2u + 1u, 0);
    uint64_t total_positions = 0;
    for (size_t sequence_id = 0; sequence_id < reads.size() * 2u;
         ++sequence_id) {
      oriented_offsets_[sequence_id] =
          static_cast<uint32_t>(total_positions);
      total_positions += reads[sequence_id >> 1u].length;
    }
    if (total_positions > std::numeric_limits<uint32_t>::max()) {
      throw std::length_error("local multi-k group tape exceeds 32-bit range");
    }
    oriented_offsets_.back() = static_cast<uint32_t>(total_positions);
    groups_.assign(static_cast<size_t>(total_positions), UINT32_MAX);
    stats_.build += omp_get_wtime() - begin;
  }

  uint32_t *ForwardGroups(size_t read_id) {
    return groups_.data() + oriented_offsets_[read_id * 2u];
  }

  uint32_t *ReverseGroups(size_t read_id) {
    return groups_.data() + oriented_offsets_[read_id * 2u + 1u];
  }

  void FinishInitial(uint64_t vertex_count) {
    const double begin_time = omp_get_wtime();
    if (vertex_count > (uint64_t(UINT32_MAX) >> 1u)) {
      throw std::length_error("local multi-k group IDs exceed 32-bit range");
    }
    group_count_ = static_cast<uint32_t>(vertex_count * 2u);

    // Build the occurrence frontier once.  Later k values only refine the
    // already contiguous equivalence-class ranges; they never histogram and
    // scatter the full occurrence tape again.
    group_counts_.assign(group_count_, 0);
    group_meta_.assign(group_count_, GroupMeta());
    BuildRawOrdinalOffsets(current_k_);
    for (uint32_t sequence_id = 0; sequence_id < reads_->size() * 2u;
         ++sequence_id) {
      const LocalPackedReads::ReadView read = (*reads_)[sequence_id >> 1u];
      if (read.length < current_k_) continue;
      const uint32_t offset = oriented_offsets_[sequence_id];
      const uint32_t starts = read.length - current_k_ + 1u;
      for (uint32_t start = 0; start < starts; ++start) {
        const uint32_t group = groups_[offset + start];
        if (group >= group_count_) {
          throw std::logic_error("invalid initial local multi-k group");
        }
        UpdateGroupMeta(group, sequence_id, offset + start, current_k_);
      }
    }
    live_groups_.clear();
    live_groups_.reserve(group_count_);
    for (uint32_t group = 0; group < group_count_; ++group) {
      group_counts_[group] = group_meta_[group].count;
      if (group_counts_[group] != 0u) live_groups_.push_back(group);
    }

    // Singleton classes can never split at a larger k.  Exclude them from the
    // persistent frontier permanently; their group ID remains valid in the
    // occurrence tape without any more work.
    group_offsets_.resize(group_count_ + 1u);
    group_offsets_[0] = 0;
    for (uint32_t group = 0; group < group_count_; ++group) {
      group_offsets_[group + 1u] =
          group_offsets_[group] +
          (group_counts_[group] > 1u ? group_counts_[group] : 0u);
    }
    group_cursors_ = group_offsets_;
    refine_records_.resize(group_offsets_.back());
    active_groups_.clear();
    active_groups_.reserve(group_count_);
    for (uint32_t group = 0; group < group_count_; ++group) {
      if (group_counts_[group] > 1u) {
        active_groups_.push_back(
            {group, group_offsets_[group], group_offsets_[group + 1u]});
      }
    }
    for (uint32_t sequence_id = 0; sequence_id < reads_->size() * 2u;
         ++sequence_id) {
      const LocalPackedReads::ReadView read = (*reads_)[sequence_id >> 1u];
      if (read.length < current_k_) continue;
      const uint32_t offset = oriented_offsets_[sequence_id];
      const uint32_t starts = read.length - current_k_ + 1u;
      for (uint32_t start = 0; start < starts; ++start) {
        const uint32_t flat = offset + start;
        const uint32_t group = groups_[flat];
        if (group_counts_[group] > 1u) {
          refine_records_[group_cursors_[group]++] = {0, flat, sequence_id};
        }
      }
    }
    UpdatePartners(current_k_);
    stats_.finish_initial += omp_get_wtime() - begin_time;
  }

  uint64_t BuildGraphAggregated(HashGraph *graph, uint32_t kmer_size) {
    double phase_begin = omp_get_wtime();
    AdvanceTo(kmer_size);
    stats_.advance += omp_get_wtime() - phase_begin;
    phase_begin = omp_get_wtime();
    graph->clear();
    graph->set_kmer_size(kmer_size);

    if (owner_group_.size() < group_count_) {
      owner_group_.resize(group_count_, UINT32_MAX);
      vertex_by_group_.resize(group_count_, UINT32_MAX);
    }
    canonical_order_.clear();
    canonical_order_.reserve(live_groups_.size() / 2u + 1u);

    for (uint32_t group : live_groups_) {
      const uint32_t partner = group_meta_[group].partner;
      if (partner >= group_count_ || group >= partner) continue;
      const uint32_t canonical = GroupLess(group, partner, kmer_size)
                                     ? group
                                     : partner;
      owner_group_[group] = canonical;
      owner_group_[partner] = canonical;
      canonical_order_.push_back(
          (uint64_t(group_meta_[canonical].first_ordinal) << 32u) |
          canonical);
    }

    // Historical graph insertion is ordered by first raw occurrence.  Pack
    // that semantic rank with the group ID and radix-sort a contiguous scalar
    // stream instead of repeatedly gathering GroupMeta through comparison
    // sort.  The rank is unique for distinct canonical k-mers.
    kmlib::kmsort(canonical_order_.begin(), canonical_order_.end());

    uint32_t last_first_ordinal = 0u;
    for (uint64_t ordered_group : canonical_order_) {
      const uint32_t group = static_cast<uint32_t>(ordered_group);
      const GroupMeta &meta = group_meta_[group];
      const LocalPackedReads::ReadView read =
          (*reads_)[meta.sequence_id >> 1u];
      const uint64_t *words = (meta.sequence_id & 1u) == 0u
                                  ? read.forward
                                  : read.reverse;
      const uint32_t start =
          meta.flat - oriented_offsets_[meta.sequence_id];
      IdbaKmer key;
      key.AssignPackedBases(words, start, kmer_size);
      uint32_t vertex_index = UINT32_MAX;
      graph->InsertAggregate(key, meta.count, meta.in_edges, meta.out_edges,
                             &vertex_index);
      vertex_by_group_[group] = vertex_index;
      vertex_by_group_[meta.partner] = vertex_index;
      last_first_ordinal = static_cast<uint32_t>(ordered_group >> 32u);
    }
    graph->FinishBulkOccurrences(
        !canonical_order_.empty() &&
        uint64_t(last_first_ordinal) + 1u < total_raw_occurrences_);

    for (uint64_t ordered_group : canonical_order_) {
      const uint32_t canonical = static_cast<uint32_t>(ordered_group);
      const uint32_t vertex_index = vertex_by_group_[canonical];
      const uint32_t oriented_groups[2] = {
          canonical, group_meta_[canonical].partner};
      for (uint32_t strand = 0; strand < 2u; ++strand) {
        const uint32_t oriented = oriented_groups[strand];
        const GroupMeta &meta = group_meta_[oriented];
        const uint8_t edges = meta.out_edges;
        if (edges == 0u) {
          continue;
        }
        if ((edges & (edges - 1u)) != 0u) {
          // Preserve the historical three-state cache machine.  A branching
          // vertex must be explicitly uncacheable; leaving it in the initial
          // "not observed" state would let a later anchor insertion install a
          // spurious single successor.
          graph->SetAggregateNeighbor(vertex_index, strand != 0u, vertex_index,
                                      false);
          continue;
        }
        if (meta.sole_next_flat == UINT32_MAX) continue;
        const uint32_t next_group = groups_[meta.sole_next_flat + 1u];
        if (next_group >= group_count_) continue;
        const uint32_t next_owner = owner_group_[next_group];
        if (next_owner >= group_count_) continue;
        graph->SetAggregateNeighbor(
            vertex_index, strand != 0u, vertex_by_group_[next_owner],
            next_group != next_owner);
      }
    }

    for (uint32_t group : live_groups_) {
      owner_group_[group] = UINT32_MAX;
      vertex_by_group_[group] = UINT32_MAX;
    }
    stats_.aggregate += omp_get_wtime() - phase_begin;
    return total_raw_occurrences_;
  }

  uint64_t BuildGraph(HashGraph *graph, uint32_t kmer_size) {
    double phase_begin = omp_get_wtime();
    AdvanceTo(kmer_size);
    stats_.advance += omp_get_wtime() - phase_begin;
    phase_begin = omp_get_wtime();
    graph->clear();
    graph->set_kmer_size(kmer_size);

    // Only groups touched by this k acquire an aggregate.  Keep the lookup
    // workspace initialized and clear those touched slots while replaying
    // below; memset(group_count_) on every round increasingly scans dead
    // historical groups created by earlier refinements.
    if (aggregate_of_group_.size() < group_count_) {
      aggregate_of_group_.resize(group_count_, UINT32_MAX);
    }
    aggregates_.clear();
    uint64_t ordinal = 0;
    for (uint32_t read_id = 0; read_id < reads_->size(); ++read_id) {
      const LocalPackedReads::ReadView read = (*reads_)[read_id];
      if (read.length < kmer_size) continue;
      uint32_t previous_aggregate = UINT32_MAX;
      bool previous_reverse = false;
      const uint32_t num_kmers = read.length - kmer_size + 1u;
      for (uint32_t start = 0; start < num_kmers; ++start, ++ordinal) {
        const uint32_t reverse_start = read.length - start - kmer_size;
        const uint32_t forward_group =
            groups_[oriented_offsets_[read_id * 2u] + start];
        const uint32_t reverse_group =
            groups_[oriented_offsets_[read_id * 2u + 1u] + reverse_start];
        if (forward_group >= group_count_ || reverse_group >= group_count_) {
          throw std::logic_error("invalid local multi-k equivalence class");
        }

        // An exact oriented string belongs to exactly one reverse-complement
        // pair, so the smaller (arbitrary) group ID is a collision-free pair
        // representative.  Determine lexical canonical orientation once per
        // unique k-mer instead of doing a multiword comparison per occurrence.
        const bool representative_reverse = reverse_group < forward_group;
        const uint32_t representative_group =
            representative_reverse ? reverse_group : forward_group;
        uint32_t aggregate_id = aggregate_of_group_[representative_group];
        bool is_reverse;
        if (aggregate_id == UINT32_MAX) {
          is_reverse =
              PackedLess(read.reverse, reverse_start, read.forward, start,
                         kmer_size);
          aggregate_id = static_cast<uint32_t>(aggregates_.size());
          aggregate_of_group_[representative_group] = aggregate_id;
          aggregates_.emplace_back(representative_group, ordinal,
                                   read_id * 2u + uint32_t(is_reverse),
                                   is_reverse ? reverse_start : start,
                                   is_reverse != representative_reverse);
        } else {
          is_reverse = representative_reverse !=
                       aggregates_[aggregate_id].representative_flip;
        }
        Aggregate &aggregate = aggregates_[aggregate_id];
        ++aggregate.count;

        if (start != 0) {
          const uint8_t edge = static_cast<uint8_t>(
              1u << (3u - LocalPackedReads::PackedBase(read.forward,
                                                       start - 1u)));
          if (is_reverse)
            aggregate.out_edges |= edge;
          else
            aggregate.in_edges |= edge;
        }
        if (start + kmer_size < read.length) {
          const uint8_t edge = static_cast<uint8_t>(
              1u << LocalPackedReads::PackedBase(read.forward,
                                                 start + kmer_size));
          if (is_reverse)
            aggregate.in_edges |= edge;
          else
            aggregate.out_edges |= edge;
        }

        if (previous_aggregate != UINT32_MAX) {
          aggregates_[previous_aggregate].AddNeighbor(
              previous_reverse,
              (uint64_t(aggregate_id) << 1u) | uint64_t(is_reverse));
          aggregate.AddNeighbor(
              !is_reverse,
              (uint64_t(previous_aggregate) << 1u) |
                  uint64_t(!previous_reverse));
        }
        previous_aggregate = aggregate_id;
        previous_reverse = is_reverse;
      }
    }

    stats_.aggregate += omp_get_wtime() - phase_begin;
    phase_begin = omp_get_wtime();

    vertex_indices_.resize(aggregates_.size());
    for (uint32_t aggregate_id = 0; aggregate_id < aggregates_.size();
         ++aggregate_id) {
      Aggregate &aggregate = aggregates_[aggregate_id];
      const LocalPackedReads::ReadView read =
          (*reads_)[aggregate.sequence_id >> 1u];
      const uint64_t *words =
          (aggregate.sequence_id & 1u) == 0 ? read.forward : read.reverse;
      IdbaKmer key;
      key.AssignPackedBases(words, aggregate.start, kmer_size);
      graph->InsertAggregate(key, aggregate.count, aggregate.in_edges,
                             aggregate.out_edges,
                             &vertex_indices_[aggregate_id]);
    }
    graph->FinishBulkOccurrences(
        !aggregates_.empty() && aggregates_.back().first_ordinal + 1u < ordinal);

    for (uint32_t aggregate_id = 0; aggregate_id < aggregates_.size();
         ++aggregate_id) {
      for (uint32_t strand = 0; strand < 2u; ++strand) {
        const uint64_t candidate = aggregates_[aggregate_id].next[strand];
        if (candidate >= Aggregate::kAmbiguous) continue;
        const uint32_t neighbor = static_cast<uint32_t>(candidate >> 1u);
        graph->SetAggregateNeighbor(vertex_indices_[aggregate_id], strand != 0,
                                    vertex_indices_[neighbor],
                                    (candidate & 1u) != 0);
      }
      aggregate_of_group_[aggregates_[aggregate_id].group] = UINT32_MAX;
    }
    stats_.replay += omp_get_wtime() - phase_begin;
    return ordinal;
  }

  const Stats &stats() const { return stats_; }

 private:
  struct GroupMeta {
    uint32_t count{0};
    uint32_t first_ordinal{UINT32_MAX};
    uint32_t sequence_id{UINT32_MAX};
    uint32_t flat{UINT32_MAX};
    uint32_t partner{UINT32_MAX};
    uint32_t sole_next_flat{UINT32_MAX};
    uint8_t in_edges{0};
    uint8_t out_edges{0};

    void Reset() {
      count = 0;
      first_ordinal = UINT32_MAX;
      sequence_id = UINT32_MAX;
      flat = UINT32_MAX;
      partner = UINT32_MAX;
      sole_next_flat = UINT32_MAX;
      in_edges = 0;
      out_edges = 0;
    }
  };

  struct RefineRecord {
    uint64_t extension;
    uint32_t flat;
    uint32_t sequence_id;
  };

  struct ActiveGroup {
    uint32_t id;
    uint32_t begin;
    uint32_t end;
  };

  struct Aggregate {
    static constexpr uint64_t kNone = UINT64_MAX;
    static constexpr uint64_t kAmbiguous = UINT64_MAX - 1u;
    Aggregate(uint32_t group_id, uint64_t ordinal, uint32_t sequence,
              uint32_t sequence_start, bool flip)
        : group(group_id),
          first_ordinal(ordinal),
          sequence_id(sequence),
          start(sequence_start),
          representative_flip(flip) {
      next[0] = kNone;
      next[1] = kNone;
    }

    void AddNeighbor(bool strand, uint64_t candidate) {
      uint64_t &slot = next[strand ? 1u : 0u];
      if (slot == kNone)
        slot = candidate;
      else if (slot != candidate)
        slot = kAmbiguous;
    }

    uint32_t group;
    uint32_t count{0};
    uint8_t in_edges{0};
    uint8_t out_edges{0};
    uint64_t first_ordinal;
    uint32_t sequence_id;
    uint32_t start;
    bool representative_flip;
    uint64_t next[2];
  };

  void AdvanceTo(uint32_t target_k) {
    while (current_k_ < target_k) {
      const uint32_t chunk =
          std::min<uint32_t>(32u, target_k - current_k_);
      const uint32_t new_k = current_k_ + chunk;
      BuildRawOrdinalOffsets(new_k);
      next_live_groups_.clear();
      next_live_groups_.reserve(live_groups_.size());
      next_active_groups_.clear();
      next_active_groups_.reserve(active_groups_.size());
      uint32_t compact_write = 0;
      size_t active_cursor = 0;

      for (uint32_t old_group : live_groups_) {
        const GroupMeta old_meta = group_meta_[old_group];
        if (old_meta.count <= 1u) {
          const uint32_t start =
              old_meta.flat - oriented_offsets_[old_meta.sequence_id];
          const uint32_t length =
              (*reads_)[old_meta.sequence_id >> 1u].length;
          if (start + new_k > length) continue;
          group_meta_[old_group].Reset();
          UpdateGroupMeta(old_group, old_meta.sequence_id, old_meta.flat,
                          new_k);
          next_live_groups_.push_back(old_group);
          continue;
        }

        if (active_cursor >= active_groups_.size() ||
            active_groups_[active_cursor].id != old_group) {
          throw std::logic_error(
              "invalid variable-order duplicate-group frontier");
        }
        const ActiveGroup active = active_groups_[active_cursor++];
        uint32_t valid_end = active.begin;
        uint64_t first_extension = UINT64_MAX;
        bool extensions_differ = false;
        for (uint32_t i = active.begin; i < active.end; ++i) {
          RefineRecord record = refine_records_[i];
          const LocalPackedReads::ReadView read =
              (*reads_)[record.sequence_id >> 1u];
          const uint32_t start =
              record.flat - oriented_offsets_[record.sequence_id];
          if (start + new_k > read.length) continue;
          const uint64_t *words = (record.sequence_id & 1u) == 0u
                                      ? read.forward
                                      : read.reverse;
          record.extension = PackedWord(words, start + current_k_);
          if (chunk != 32u) {
            record.extension &=
                (uint64_t{1} << (chunk << 1u)) - 1u;
          }
          if (valid_end == active.begin)
            first_extension = record.extension;
          else
            extensions_differ |= record.extension != first_extension;
          refine_records_[valid_end++] = record;
        }
        if (valid_end == active.begin) continue;
        if (extensions_differ) {
          std::sort(refine_records_.begin() + active.begin,
                    refine_records_.begin() + valid_end,
                    [](const RefineRecord &lhs, const RefineRecord &rhs) {
                      return lhs.extension < rhs.extension;
                    });
        }

        uint32_t run_begin = active.begin;
        bool first_run = true;
        while (run_begin < valid_end) {
          uint32_t run_end = run_begin + 1u;
          while (run_end < valid_end &&
                 refine_records_[run_end].extension ==
                     refine_records_[run_begin].extension) {
            ++run_end;
          }
          uint32_t new_group = old_group;
          if (!first_run) {
            if (group_count_ == UINT32_MAX) {
              throw std::length_error(
                  "local variable-order group IDs exceed 32-bit range");
            }
            new_group = group_count_++;
            if (group_meta_.size() < group_count_) {
              group_meta_.resize(group_count_);
            }
          }
          first_run = false;
          group_meta_[new_group].Reset();
          const uint32_t compact_begin = compact_write;
          for (uint32_t i = run_begin; i < run_end; ++i) {
            const RefineRecord record = refine_records_[i];
            groups_[record.flat] = new_group;
            UpdateGroupMeta(new_group, record.sequence_id, record.flat,
                            new_k);
            if (run_end - run_begin > 1u) {
              refine_records_[compact_write++] = record;
            }
          }
          next_live_groups_.push_back(new_group);
          if (run_end - run_begin > 1u) {
            next_active_groups_.push_back(
                {new_group, compact_begin, compact_write});
          }
          run_begin = run_end;
        }
      }
      if (active_cursor != active_groups_.size()) {
        throw std::logic_error(
            "incomplete variable-order duplicate-group frontier");
      }
      refine_records_.resize(compact_write);
      active_groups_.swap(next_active_groups_);
      live_groups_.swap(next_live_groups_);
      current_k_ = new_k;
      UpdatePartners(current_k_);
    }
  }

  void AdvanceToLegacy(uint32_t target_k) {
    while (current_k_ < target_k) {
      const uint32_t chunk = std::min<uint32_t>(32u, target_k - current_k_);
      const uint32_t new_k = current_k_ + chunk;
      uint32_t compact_write = 0;
      next_active_groups_.clear();

      // For a small extension alphabet, grouping through a dense directory
      // moves far fewer bytes than sorting every 16-byte occurrence record.
      // Enable it only when its three compact uint32 directories are no
      // larger than the active tape they replace; this is a work/space cost
      // decision, not a k- or machine-specific threshold.
      uint64_t extension_space = 0;
      if (chunk < 32u) extension_space = uint64_t{1} << (chunk << 1u);
      const bool use_dense_extension =
          extension_space != 0u &&
          extension_space <= std::numeric_limits<size_t>::max() /
                                 (3u * sizeof(uint32_t)) &&
          extension_space * (3u * sizeof(uint32_t)) <=
              refine_records_.size() * sizeof(RefineRecord);
      if (use_dense_extension) {
        const size_t directory_size = static_cast<size_t>(extension_space);
        if (extension_counts_.size() < directory_size) {
          extension_counts_.resize(directory_size, 0u);
          extension_groups_.resize(directory_size);
          extension_cursors_.resize(directory_size);
        }
      }
      for (const ActiveGroup &old_group : active_groups_) {
        const uint32_t begin = old_group.begin;
        const uint32_t end = old_group.end;

        if (use_dense_extension) {
          refine_group_workspace_.clear();
          touched_extensions_.clear();
          refine_group_workspace_.reserve(end - begin);

          for (uint32_t i = begin; i < end; ++i) {
            RefineRecord record = refine_records_[i];
            const LocalPackedReads::ReadView read =
                (*reads_)[record.sequence_id >> 1u];
            const uint32_t start =
                record.flat - oriented_offsets_[record.sequence_id];
            if (start + new_k > read.length) continue;

            const uint64_t *words = (record.sequence_id & 1u) == 0
                                        ? read.forward
                                        : read.reverse;
            record.extension = PackedWord(words, start + current_k_);
            if (chunk != 32u) {
              record.extension &=
                  (uint64_t{1} << (chunk << 1u)) - 1u;
            }
            const uint32_t extension =
                static_cast<uint32_t>(record.extension);
            if (extension_counts_[extension]++ == 0u) {
              touched_extensions_.push_back(extension);
            }
            refine_group_workspace_.push_back(record);
          }
          if (refine_group_workspace_.empty()) continue;

          // Assign groups in the same ascending-extension order produced by
          // the historical record sort. This keeps even the internal group
          // IDs stable while sorting only distinct 32-bit extension values.
          std::sort(touched_extensions_.begin(), touched_extensions_.end());
          bool first_extension = true;
          for (uint32_t extension : touched_extensions_) {
            uint32_t new_group = old_group.id;
            if (!first_extension) {
              if (group_count_ == UINT32_MAX) {
                throw std::length_error(
                    "local multi-k refined group IDs exceed 32-bit range");
              }
              new_group = group_count_++;
            }
            first_extension = false;
            extension_groups_[extension] = new_group;

            const uint32_t run_size = extension_counts_[extension];
            if (run_size > 1u) {
              const uint32_t compact_begin = compact_write;
              extension_cursors_[extension] = compact_begin;
              compact_write += run_size;
              next_active_groups_.push_back(
                  {new_group, compact_begin, compact_write});
            }
          }

          for (const RefineRecord &record : refine_group_workspace_) {
            const uint32_t extension =
                static_cast<uint32_t>(record.extension);
            groups_[record.flat] = extension_groups_[extension];
            if (extension_counts_[extension] > 1u) {
              refine_records_[extension_cursors_[extension]++] = record;
            }
          }
          for (uint32_t extension : touched_extensions_) {
            extension_counts_[extension] = 0u;
          }
          continue;
        }

        // Drop suffix positions invalidated by increasing k and extract the
        // newly exposed bases in the same pass.  The previous implementation
        // walked every live occurrence twice and repeated its read/offset
        // gather; refinement is bandwidth-bound, so keep that state hot.
        uint32_t valid_end = begin;
        uint64_t first_extension = UINT64_MAX;
        bool extensions_differ = false;
        for (uint32_t i = begin; i < end; ++i) {
          RefineRecord record = refine_records_[i];
          const LocalPackedReads::ReadView read =
              (*reads_)[record.sequence_id >> 1u];
          const uint32_t start =
              record.flat - oriented_offsets_[record.sequence_id];
          if (start + new_k <= read.length) {
            const uint64_t *words = (record.sequence_id & 1u) == 0
                                        ? read.forward
                                        : read.reverse;
            record.extension = PackedWord(words, start + current_k_);
            if (chunk != 32u) {
              record.extension &=
                  (uint64_t{1} << (chunk << 1u)) - 1u;
            }
            if (valid_end == begin)
              first_extension = record.extension;
            else
              extensions_differ |= record.extension != first_extension;
            refine_records_[valid_end++] = record;
          }
        }
        if (valid_end == begin) continue;

        if (valid_end - begin > 1u) {
          if (extensions_differ) {
            std::sort(refine_records_.begin() + begin,
                      refine_records_.begin() + valid_end,
                      [](const RefineRecord &lhs, const RefineRecord &rhs) {
                        return lhs.extension < rhs.extension;
                      });
          }
        } else {
          refine_records_[begin].extension = 0;
        }

        bool first_run = true;
        uint32_t run_begin = begin;
        while (run_begin < valid_end) {
          uint32_t run_end = run_begin + 1u;
          while (run_end < valid_end &&
                 refine_records_[run_end].extension ==
                     refine_records_[run_begin].extension) {
            ++run_end;
          }
          uint32_t new_group = old_group.id;
          if (!first_run) {
            if (group_count_ == UINT32_MAX) {
              throw std::length_error(
                  "local multi-k refined group IDs exceed 32-bit range");
            }
            new_group = group_count_++;
          }
          for (uint32_t i = run_begin; i < run_end; ++i) {
            groups_[refine_records_[i].flat] = new_group;
          }
          const uint32_t run_size = run_end - run_begin;
          if (run_size > 1u) {
            const uint32_t compact_begin = compact_write;
            for (uint32_t i = run_begin; i < run_end; ++i) {
              refine_records_[compact_write++] = refine_records_[i];
            }
            next_active_groups_.push_back(
                {new_group, compact_begin, compact_write});
          }
          first_run = false;
          run_begin = run_end;
        }
      }
      refine_records_.resize(compact_write);
      active_groups_.swap(next_active_groups_);
      current_k_ = new_k;
    }
  }

  static uint64_t PackedWord(const uint64_t *words, uint32_t base_offset) {
    const uint32_t source_word = base_offset >> 5u;
    const uint32_t shift = (base_offset & 31u) << 1u;
    uint64_t value = words[source_word] >> shift;
    if (shift != 0) value |= words[source_word + 1u] << (64u - shift);
    return value;
  }

  static bool PackedLess(const uint64_t *lhs, uint32_t lhs_start,
                         const uint64_t *rhs, uint32_t rhs_start,
                         uint32_t kmer_size) {
    const uint32_t words = (kmer_size + 31u) >> 5u;
    const uint32_t tail = kmer_size & 31u;
    for (int32_t word = static_cast<int32_t>(words) - 1; word >= 0; --word) {
      uint64_t lhs_word = PackedWord(lhs, lhs_start + uint32_t(word) * 32u);
      uint64_t rhs_word = PackedWord(rhs, rhs_start + uint32_t(word) * 32u);
      if (word == static_cast<int32_t>(words) - 1 && tail != 0) {
        const uint64_t mask = (uint64_t{1} << (tail << 1u)) - 1u;
        lhs_word &= mask;
        rhs_word &= mask;
      }
      if (lhs_word != rhs_word) return lhs_word < rhs_word;
    }
    return false;
  }

  void BuildRawOrdinalOffsets(uint32_t kmer_size) {
    raw_ordinal_offsets_.resize(reads_->size() + 1u);
    uint64_t total = 0;
    for (uint32_t read_id = 0; read_id < reads_->size(); ++read_id) {
      if (total > UINT32_MAX) {
        throw std::length_error(
            "local variable-order ordinal tape exceeds 32-bit range");
      }
      raw_ordinal_offsets_[read_id] = static_cast<uint32_t>(total);
      const uint32_t length = (*reads_)[read_id].length;
      if (length >= kmer_size) total += length - kmer_size + 1u;
    }
    if (total > UINT32_MAX) {
      throw std::length_error(
          "local variable-order ordinal tape exceeds 32-bit range");
    }
    raw_ordinal_offsets_.back() = static_cast<uint32_t>(total);
    total_raw_occurrences_ = total;
  }

  uint32_t RawOrdinal(uint32_t sequence_id, uint32_t start,
                      uint32_t kmer_size) const {
    const uint32_t read_id = sequence_id >> 1u;
    const uint32_t length = (*reads_)[read_id].length;
    const uint32_t raw_start = (sequence_id & 1u) == 0u
                                   ? start
                                   : length - start - kmer_size;
    return raw_ordinal_offsets_[read_id] + raw_start;
  }

  void UpdateGroupMeta(uint32_t group, uint32_t sequence_id, uint32_t flat,
                       uint32_t kmer_size) {
    if (group >= group_meta_.size()) group_meta_.resize(group + 1u);
    GroupMeta &meta = group_meta_[group];
    const LocalPackedReads::ReadView read = (*reads_)[sequence_id >> 1u];
    const uint64_t *words =
        (sequence_id & 1u) == 0u ? read.forward : read.reverse;
    const uint32_t start = flat - oriented_offsets_[sequence_id];
    if (meta.count == 0u) {
      meta.sequence_id = sequence_id;
      meta.flat = flat;
    }
    ++meta.count;
    meta.first_ordinal =
        std::min(meta.first_ordinal,
                 RawOrdinal(sequence_id, start, kmer_size));
    if (start != 0u) {
      meta.in_edges |= static_cast<uint8_t>(
          1u << (3u - LocalPackedReads::PackedBase(words, start - 1u)));
    }
    if (start + kmer_size < read.length) {
      const uint8_t base =
          LocalPackedReads::PackedBase(words, start + kmer_size);
      const uint8_t edge = static_cast<uint8_t>(1u << base);
      if (meta.out_edges == 0u) {
        meta.sole_next_flat = flat;
      } else if ((meta.out_edges & edge) == 0u) {
        meta.sole_next_flat = UINT32_MAX;
      }
      meta.out_edges |= edge;
    }
  }

  void UpdatePartners(uint32_t kmer_size) {
    for (uint32_t group : live_groups_) {
      GroupMeta &meta = group_meta_[group];
      const uint32_t reverse_sequence = meta.sequence_id ^ 1u;
      const uint32_t start =
          meta.flat - oriented_offsets_[meta.sequence_id];
      const uint32_t length = (*reads_)[meta.sequence_id >> 1u].length;
      const uint32_t reverse_start = length - start - kmer_size;
      meta.partner =
          groups_[oriented_offsets_[reverse_sequence] + reverse_start];
    }
  }

  bool GroupLess(uint32_t lhs, uint32_t rhs, uint32_t kmer_size) const {
    const GroupMeta &left = group_meta_[lhs];
    const GroupMeta &right = group_meta_[rhs];
    const LocalPackedReads::ReadView left_read =
        (*reads_)[left.sequence_id >> 1u];
    const LocalPackedReads::ReadView right_read =
        (*reads_)[right.sequence_id >> 1u];
    const uint64_t *left_words = (left.sequence_id & 1u) == 0u
                                     ? left_read.forward
                                     : left_read.reverse;
    const uint64_t *right_words = (right.sequence_id & 1u) == 0u
                                      ? right_read.forward
                                      : right_read.reverse;
    return PackedLess(left_words,
                      left.flat - oriented_offsets_[left.sequence_id],
                      right_words,
                      right.flat - oriented_offsets_[right.sequence_id],
                      kmer_size);
  }

  LocalPackedReads *reads_{nullptr};
  uint32_t current_k_{0};
  uint32_t group_count_{0};
  std::vector<uint32_t> oriented_offsets_;
  std::vector<uint32_t> groups_;
  std::vector<uint32_t> group_counts_;
  std::vector<uint32_t> group_offsets_;
  std::vector<uint32_t> group_cursors_;
  std::vector<RefineRecord> refine_records_;
  std::vector<ActiveGroup> active_groups_;
  std::vector<ActiveGroup> next_active_groups_;
  std::vector<RefineRecord> refine_group_workspace_;
  std::vector<uint32_t> touched_extensions_;
  std::vector<uint32_t> extension_counts_;
  std::vector<uint32_t> extension_groups_;
  std::vector<uint32_t> extension_cursors_;
  std::vector<uint32_t> aggregate_of_group_;
  std::vector<Aggregate> aggregates_;
  std::vector<uint32_t> vertex_indices_;
  std::vector<GroupMeta> group_meta_;
  std::vector<uint32_t> live_groups_;
  std::vector<uint32_t> next_live_groups_;
  std::vector<uint32_t> raw_ordinal_offsets_;
  std::vector<uint64_t> canonical_order_;
  std::vector<uint32_t> owner_group_;
  std::vector<uint32_t> vertex_by_group_;
  uint64_t total_raw_occurrences_{0};
  Stats stats_;
};

// Exact cross-k occurrence refinement for the read-derived graph.  When the
// increment s does not exceed the current k, a (k+s)-mer is uniquely
// identified by its length-k prefix and suffix.  The previous round already
// assigned an exact oriented vertex code to both, so later rounds operate on
// one compact 64-bit signature instead of rescanning bases, reverse-
// complementing a full key, and probing the k-mer table for every occurrence.
// Full DNA keys are materialized only once per distinct signature, and all
// counts/edges/transitions are replayed in historical read order.
class IncrementalKReadGraph {
 public:
  struct Stats {
    uint64_t occurrences{0};
    uint64_t path_replays{0};
    uint64_t primary_inserts{0};
    uint64_t primary_hits{0};
    uint64_t fork_lookups{0};
  };

  void Initialize(LocalPackedReads *reads, uint32_t initial_k) {
    stats_ = Stats();
    reads_ = reads;
    current_k_ = initial_k;
    offsets_.resize(reads_->size() + 1u);
    uint64_t total = 0;
    for (size_t read_id = 0; read_id < reads_->size(); ++read_id) {
      offsets_[read_id] = static_cast<uint32_t>(total);
      const uint32_t length = (*reads_)[read_id].length;
      if (length >= initial_k) total += length - initial_k + 1u;
    }
    if (total > UINT32_MAX) {
      throw std::length_error(
          "local incremental-k occurrence tape exceeds 32-bit range");
    }
    offsets_.back() = static_cast<uint32_t>(total);
    codes_.resize(static_cast<size_t>(total));
  }

  uint32_t *InitialCodes(size_t read_id) {
    return codes_.data() + offsets_[read_id];
  }

  void FinishInitial(uint64_t vertex_count) {
    if (vertex_count >= (uint64_t{1} << 31u)) {
      throw std::length_error(
          "local incremental-k vertex code exceeds 31 bits");
    }
    previous_code_space_ = static_cast<uint32_t>(vertex_count * 2u);
  }

  bool CanAdvance(uint32_t target_k) const {
    return reads_ != nullptr && target_k > current_k_ &&
           target_k - current_k_ <= current_k_ &&
           (current_k_ & 1u) != 0u && (target_k & 1u) != 0u;
  }

  uint64_t BuildGraph(HashGraph *graph, uint32_t target_k) {
    if (!CanAdvance(target_k)) {
      throw std::logic_error("unsupported exact incremental-k transition");
    }
    const uint32_t delta = target_k - current_k_;
    uint64_t total64 = 0;
    for (size_t read_id = 0; read_id < reads_->size(); ++read_id) {
      const uint32_t length = (*reads_)[read_id].length;
      if (length >= target_k) total64 += length - target_k + 1u;
    }
    if (total64 > UINT32_MAX) {
      throw std::length_error(
          "local incremental-k round exceeds 32-bit occurrence range");
    }
    const uint32_t total = static_cast<uint32_t>(total64);

    // A target string is uniquely identified by its overlapping oriented
    // prefix/suffix vertex codes from the preceding k.  Most oriented k-mers
    // have one exact delta-base continuation, so keep that primary signature
    // in a direct array and hash only the exceptional forks.
    //
    // Previous-k prefix/suffix codes resolve a cold target signature without
    // hashing its full DNA. Once a current-k transition has been observed,
    // replay it directly in HashGraph. Its oriented vertex code is both the
    // exact target identity and the code required by the next k round, so the
    // hot path needs no second randomly accessed aggregate graph.
    for (uint32_t prefix : touched_prefixes_) {
      primary_suffix_[prefix] = UINT32_MAX;
    }
    touched_prefixes_.clear();
    const size_t prefix_slots = size_t(previous_code_space_);
    if (primary_suffix_.size() < prefix_slots) {
      primary_suffix_.resize(prefix_slots, UINT32_MAX);
      primary_state_.resize(prefix_slots);
    }
    signature_table_.Reset(16u);
    signature_states_.clear();
    signature_states_.reserve(
        std::min<size_t>(total, std::max<size_t>(64u, prefix_slots)));

    uint32_t ordinal = 0;
    for (size_t read_id = 0; read_id < reads_->size(); ++read_id) {
      const LocalPackedReads::ReadView read = (*reads_)[read_id];
      if (read.length < target_k) continue;
      const uint32_t offset = offsets_[read_id];
      const uint32_t starts = read.length - target_k + 1u;
      uint32_t previous_code = UINT32_MAX;
      for (uint32_t start = 0; start < starts; ++start, ++ordinal) {
        uint32_t code = UINT32_MAX;
        if (previous_code != UINT32_MAX) {
          const uint8_t transition_base = LocalPackedReads::PackedBase(
              read.forward, start + target_k - 1u);
          if (graph->ReplayResolvedTransition(previous_code,
                                              transition_base, &code)) {
            ++stats_.path_replays;
          }
        }

        if (code == UINT32_MAX) {
          // start+delta is ahead of the slot overwritten below, so both codes
          // still belong to the preceding generation.
          const uint32_t prefix_code = codes_[offset + start];
          const uint32_t suffix_code = codes_[offset + start + delta];
          if (prefix_code >= previous_code_space_ ||
              suffix_code >= previous_code_space_) {
            throw std::logic_error(
                "invalid local incremental-k oriented vertex code");
          }
          const uint64_t signature =
              (uint64_t(prefix_code) << 32u) | suffix_code;
          bool inserted = false;
          uint32_t state_id;
          if (primary_suffix_[prefix_code] == UINT32_MAX) {
            primary_suffix_[prefix_code] = suffix_code;
            state_id = static_cast<uint32_t>(signature_states_.size());
            primary_state_[prefix_code] = state_id;
            touched_prefixes_.push_back(prefix_code);
            inserted = true;
            ++stats_.primary_inserts;
          } else if (primary_suffix_[prefix_code] == suffix_code) {
            state_id = primary_state_[prefix_code];
            ++stats_.primary_hits;
          } else {
            const uint32_t proposed =
                static_cast<uint32_t>(signature_states_.size());
            state_id = signature_table_.LookupOrInsert(
                signature, proposed, &inserted);
            ++stats_.fork_lookups;
          }

          uint8_t read_in_edges = 0u;
          uint8_t read_out_edges = 0u;
          if (start != 0u) {
            read_in_edges = static_cast<uint8_t>(
                1u << (3u - LocalPackedReads::PackedBase(read.forward,
                                                         start - 1u)));
          }
          if (start + target_k < read.length) {
            read_out_edges = static_cast<uint8_t>(
                1u << LocalPackedReads::PackedBase(read.forward,
                                                   start + target_k));
          }

          if (inserted) {
            IdbaKmer forward;
            forward.AssignPackedBases(read.forward, start, target_k);
            IdbaKmer reverse = forward;
            reverse.ReverseComplement();
            const bool is_reverse = reverse < forward;
            uint32_t vertex_index = UINT32_MAX;
            graph->InsertAggregate(
                is_reverse ? reverse : forward, 1u,
                is_reverse ? read_out_edges : read_in_edges,
                is_reverse ? read_in_edges : read_out_edges, &vertex_index);
            if (vertex_index >= (UINT32_MAX >> 1u)) {
              throw std::length_error(
                  "local incremental-k vertex code exceeds 31 bits");
            }
            code = (vertex_index << 1u) | uint32_t(is_reverse);
            SignatureState state;
            state.code = code;
            signature_states_.push_back(state);
          } else {
            code = signature_states_[state_id].code;
            const bool is_reverse = (code & 1u) != 0u;
            graph->AccumulateResolvedOccurrence(
                code >> 1u,
                is_reverse ? read_out_edges : read_in_edges,
                is_reverse ? read_in_edges : read_out_edges);
          }
          if (previous_code != UINT32_MAX) {
            graph->ObserveTransition(previous_code, code);
          }
        }

        ++stats_.occurrences;
        previous_code = code;
        // All still-unread prefix/suffix codes are at larger positions. Store
        // the current code in place for the following k round.
        codes_[offset + start] = code;
      }
    }
    assert(ordinal == total);
    if (graph->num_vertices() >= (uint64_t{1} << 31u)) {
      throw std::length_error(
          "local incremental-k oriented code space exceeds 32 bits");
    }
    previous_code_space_ = static_cast<uint32_t>(graph->num_vertices() * 2u);
    current_k_ = target_k;
    return total;
  }

  const Stats &stats() const { return stats_; }

 private:
  struct SignatureState {
    uint32_t code{UINT32_MAX};
  };

  static uint64_t SignatureHash(uint64_t value) {
    value ^= value >> 30u;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27u;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31u;
    return value;
  }

  class SignatureTable {
   public:
    void Reset(size_t expected_entries) {
      size_t desired = 16u;
      while ((desired * 7u) / 10u < expected_entries) desired <<= 1u;
      if (desired > values_.size()) {
        keys_.resize(desired);
        values_.assign(desired, UINT32_MAX);
      } else {
        for (uint32_t slot : touched_) values_[slot] = UINT32_MAX;
      }
      touched_.clear();
      size_ = 0u;
      mask_ = values_.size() - 1u;
    }

    uint32_t LookupOrInsert(uint64_t key, uint32_t value, bool *inserted) {
      if ((size_ + 1u) * 10u > values_.size() * 7u) Grow();
      size_t slot = static_cast<size_t>(SignatureHash(key)) & mask_;
      while (values_[slot] != UINT32_MAX) {
        if (keys_[slot] == key) {
          *inserted = false;
          return values_[slot];
        }
        slot = (slot + 1u) & mask_;
      }
      keys_[slot] = key;
      values_[slot] = value;
      touched_.push_back(static_cast<uint32_t>(slot));
      ++size_;
      *inserted = true;
      return value;
    }

   private:
    void Grow() {
      const size_t new_capacity = values_.empty() ? 16u : values_.size() * 2u;
      std::vector<uint64_t> old_keys;
      std::vector<uint32_t> old_values;
      std::vector<uint32_t> old_touched;
      old_keys.swap(keys_);
      old_values.swap(values_);
      old_touched.swap(touched_);
      keys_.resize(new_capacity);
      values_.assign(new_capacity, UINT32_MAX);
      mask_ = new_capacity - 1u;
      touched_.reserve(old_touched.size());
      for (uint32_t old_slot : old_touched) {
        size_t slot = static_cast<size_t>(SignatureHash(old_keys[old_slot])) &
                      mask_;
        while (values_[slot] != UINT32_MAX) slot = (slot + 1u) & mask_;
        keys_[slot] = old_keys[old_slot];
        values_[slot] = old_values[old_slot];
        touched_.push_back(static_cast<uint32_t>(slot));
      }
    }

    std::vector<uint64_t> keys_;
    std::vector<uint32_t> values_;
    std::vector<uint32_t> touched_;
    size_t size_{0};
    size_t mask_{0};
  };

  LocalPackedReads *reads_{nullptr};
  uint32_t current_k_{0};
  std::vector<uint32_t> offsets_;
  std::vector<uint32_t> codes_;
  std::vector<SignatureState> signature_states_;
  std::vector<uint32_t> primary_suffix_;
  std::vector<uint32_t> primary_state_;
  std::vector<uint32_t> touched_prefixes_;
  SignatureTable signature_table_;
  uint32_t previous_code_space_{0};
  Stats stats_;
};

struct alignas(64) LocalAssemblyProfile {
  double read_graph_build{0};
  double coverage{0};
  double anchor_insert{0};
  double hash_assemble{0};
  double hash_walk{0};
  double hash_adjacency{0};
  double contig_initialize{0};
  double dead_end{0};
  double bubble{0};
  double coverage_clean{0};
  double contig_assemble{0};
  double multik_build{0};
  double multik_finish{0};
  double multik_advance{0};
  double multik_aggregate{0};
  double multik_replay{0};
  uint64_t k_rounds{0};
  uint64_t raw_kmers{0};
  uint64_t raw_vertices{0};
  uint64_t single_hash_unitigs{0};
  uint64_t long_single_hash_unitigs{0};
  uint64_t isolated_long_single_hash_unitigs{0};
  uint64_t cleaning_rounds{0};
  uint64_t deadend_changed_rounds{0};
  uint64_t bubble_changed_rounds{0};
  uint64_t coverage_changed_rounds{0};
  uint64_t branch_transition_hits{0};
  uint64_t branch_transition_misses{0};
  uint64_t incremental_occurrences{0};
  uint64_t incremental_path_replays{0};
  uint64_t incremental_primary_inserts{0};
  uint64_t incremental_primary_hits{0};
  uint64_t incremental_fork_lookups{0};
};

void LaunchIDBA(LocalPackedReads &reads,
                const Sequence &contig_end,
                std::vector<Sequence> &out_contigs,
                std::vector<ContigInfo> &out_contig_infos, uint32_t mink,
                uint32_t maxk, uint32_t step, HashGraph &hash_graph,
                ContigGraph &contig_graph,
                std::vector<ContigGraphVertex> &hash_unitigs,
                std::vector<uint32_t> &hash_unitig_neighbors,
                MultiKReadIndex &multi_k_index,
                IncrementalKReadGraph &incremental_read_graph,
                LocalAssemblyProfile *profile,
                std::vector<uint64_t> *k_round_survival,
                std::vector<double> *k_read_graph_seconds,
                std::vector<uint64_t> *k_branch_hits,
                std::vector<uint64_t> *k_branch_misses,
                uint64_t estimated_read_work) {
  int local_range = contig_end.size();
  hash_graph.reset_for_endpoint();
  out_contigs.clear();
  out_contig_infos.clear();

  // The endpoint anchor is identical in every inner-k round.  Pack it once
  // into the same two-bit layout as local reads so every round can use the
  // active-word rolling kernel instead of rebuilding fixed-capacity
  // IdbaKmers from byte bases.  Ambiguous bases retain the historical path.
  std::vector<uint64_t> packed_contig_end;
  const bool packed_contig_end_valid =
      PackSequence2Bit(contig_end, &packed_contig_end);
  std::vector<uint64_t> packed_out_contig;

  // Keep the exact incremental builder opt-in until its graph, traversal and
  // output order have passed strict A/B validation. Its choice is an
  // algorithmic reuse decision and does not depend on a particular k, input
  // size, thread count or machine topology.
  static const bool use_incremental_read_graph =
      std::getenv("MEGAHIT_EXPERIMENTAL_LOCAL_INCREMENTAL") != nullptr;
  static const bool use_multi_k_index =
      std::getenv("MEGAHIT_EXPERIMENTAL_LOCAL_MULTIK") != nullptr;
  static const bool validate_incremental_read_graph =
      std::getenv("MEGAHIT_VALIDATE_LOCAL_INCREMENTAL") != nullptr;
  static const bool validate_multi_k_index =
      std::getenv("MEGAHIT_VALIDATE_LOCAL_MULTIK") != nullptr;
  if (use_multi_k_index) {
    multi_k_index.Build(reads, mink);
  } else if (use_incremental_read_graph) {
    incremental_read_graph.Initialize(&reads, mink);
  }
  (void)estimated_read_work;

  const uint32_t max_read_len = reads.max_length();
  for (uint32_t kmer_size = mink; kmer_size <= std::min(maxk, max_read_len);
       kmer_size += step) {
    ++profile->k_rounds;
    const size_t k_index = (kmer_size - mink) / step;
    ++(*k_round_survival)[k_index];
    const uint64_t branch_hits_before =
        hash_graph.DebugBranchTransitionHits();
    const uint64_t branch_misses_before =
        hash_graph.DebugBranchTransitionMisses();
    double phase_begin = omp_get_wtime();
    hash_graph.clear();
    hash_graph.set_kmer_size(kmer_size);
    // The worker-local table retains its bucket capacity between inner-k
    // rounds.  A validation replay must start from that same capacity: using
    // the optimized graph's *final* capacity skips its real rehash sequence
    // and can manufacture a different bucket-chain traversal order even
    // when first-occurrence insertion order is identical.
    const size_t validation_initial_bucket_count =
        hash_graph.DebugBucketCount();
    if (use_multi_k_index && kmer_size != mink) {
      profile->raw_kmers +=
          multi_k_index.BuildGraphAggregated(&hash_graph, kmer_size);
      if (validate_multi_k_index) {
        HashGraph reference_graph(kmer_size);
        reference_graph.DebugResetEmptyBucketCount(
            validation_initial_bucket_count);
        for (size_t read_id = 0; read_id < reads.size(); ++read_id) {
          const LocalPackedReads::ReadView read = reads[read_id];
          if (read.length >= kmer_size) {
            reference_graph.InsertPackedKmers(read.forward, read.length);
          }
        }
        std::string difference;
        if (!hash_graph.SameGraphState(reference_graph, &difference) ||
            !hash_graph.SameTraversalState(reference_graph, &difference) ||
            !hash_graph.SameTransitionCache(reference_graph, &difference)) {
          std::fprintf(stderr,
                       "multi-k local graph differs at k=%u: %s\n",
                       kmer_size, difference.c_str());
          std::abort();
        }
      }
    } else if (use_incremental_read_graph && kmer_size != mink &&
        incremental_read_graph.CanAdvance(kmer_size)) {
      profile->raw_kmers +=
          incremental_read_graph.BuildGraph(&hash_graph, kmer_size);
      if (validate_incremental_read_graph) {
        HashGraph reference_graph(kmer_size);
        reference_graph.DebugResetEmptyBucketCount(
            validation_initial_bucket_count);
        for (size_t read_id = 0; read_id < reads.size(); ++read_id) {
          const LocalPackedReads::ReadView read = reads[read_id];
          if (read.length >= kmer_size) {
            reference_graph.InsertPackedKmers(read.forward, read.length);
          }
        }
        std::string difference;
        if (!hash_graph.SameGraphState(reference_graph, &difference)) {
          std::fprintf(stderr,
                       "incremental local graph state differs at k=%u: %s\n",
                       kmer_size, difference.c_str());
          std::abort();
        }
        if (!hash_graph.SameTraversalState(reference_graph, &difference)) {
          std::fprintf(stderr,
                       "incremental local traversal differs at k=%u: %s\n",
                       kmer_size, difference.c_str());
          std::abort();
        }
        if (!hash_graph.SameTransitionCache(reference_graph, &difference)) {
          std::fprintf(stderr,
                       "incremental local transition cache differs at k=%u: %s\n",
                       kmer_size, difference.c_str());
          std::abort();
        }
      }
    } else {
      for (size_t read_id = 0; read_id < reads.size(); ++read_id) {
        const LocalPackedReads::ReadView read = reads[read_id];
        if (read.length < kmer_size) continue;
        if (use_multi_k_index && kmer_size == mink) {
          profile->raw_kmers += hash_graph.InsertPackedKmersIndexed(
              read.forward, read.length,
              multi_k_index.ForwardGroups(read_id),
              multi_k_index.ReverseGroups(read_id));
        } else if (use_incremental_read_graph && kmer_size == mink) {
          profile->raw_kmers += hash_graph.InsertPackedKmersIndexed(
              read.forward, read.length,
              incremental_read_graph.InitialCodes(read_id), nullptr);
        } else {
          profile->raw_kmers +=
              hash_graph.InsertPackedKmers(read.forward, read.length);
        }
      }
      if (use_multi_k_index && kmer_size == mink) {
        multi_k_index.FinishInitial(hash_graph.num_vertices());
      } else if (use_incremental_read_graph && kmer_size == mink) {
        incremental_read_graph.FinishInitial(hash_graph.num_vertices());
      }
    }
    profile->raw_vertices += hash_graph.num_vertices();
    const double read_graph_seconds = omp_get_wtime() - phase_begin;
    profile->read_graph_build += read_graph_seconds;
    (*k_read_graph_seconds)[k_index] += read_graph_seconds;
    (*k_branch_hits)[k_index] +=
        hash_graph.DebugBranchTransitionHits() - branch_hits_before;
    (*k_branch_misses)[k_index] +=
        hash_graph.DebugBranchTransitionMisses() - branch_misses_before;

    phase_begin = omp_get_wtime();
    double mean = hash_graph.coverage_percentile(
        1 - 1.0 * local_range / hash_graph.num_vertices());
    double threshold = mean;
    profile->coverage += omp_get_wtime() - phase_begin;

    phase_begin = omp_get_wtime();
    if (packed_contig_end_valid) {
      hash_graph.InsertPackedKmers(packed_contig_end.data(),
                                   contig_end.size());
    } else {
      hash_graph.InsertKmers(contig_end);
    }

    for (const auto &out_contig : out_contigs) {
      if (PackSequence2Bit(out_contig, &packed_out_contig)) {
        hash_graph.InsertPackedUncountKmers(packed_out_contig.data(),
                                            out_contig.size());
      } else {
        hash_graph.InsertUncountKmers(out_contig);
      }
    }
    profile->anchor_insert += omp_get_wtime() - phase_begin;

    phase_begin = omp_get_wtime();
    hash_graph.Assemble(hash_unitigs);
    const double hash_walk_seconds = omp_get_wtime() - phase_begin;
    profile->hash_walk += hash_walk_seconds;

    // A single isolated unitig with at least 2*k k-mers is outside the
    // strict length predicate used by both tip and low-coverage cleaning;
    // with no edges it also cannot contain a bubble.  The three cleaning
    // stages and ContigGraph reconstruction are therefore exact no-ops, and
    // the historical pipeline would immediately break after emitting this
    // same sequence.
    if (hash_unitigs.size() == 1u) {
      ++profile->single_hash_unitigs;
      if (hash_unitigs.front().contig_size() >= 3u * kmer_size - 1u) {
        ++profile->long_single_hash_unitigs;
        if (hash_unitigs.front().in_edges().empty() &&
            hash_unitigs.front().out_edges().empty()) {
          ++profile->isolated_long_single_hash_unitigs;
          out_contigs.clear();
          out_contig_infos.clear();
          out_contigs.push_back(hash_unitigs.front().contig());
          out_contig_infos.push_back(hash_unitigs.front().contig_info());
          break;
        }
      }
    }

    phase_begin = omp_get_wtime();
    uint64_t hash_unitig_edges = 0;
    bool has_explicit_hash_adjacency =
        hash_graph.BuildUnitigAdjacency(hash_unitigs, hash_unitig_neighbors,
                                        &hash_unitig_edges);
    const double hash_adjacency_seconds = omp_get_wtime() - phase_begin;
    profile->hash_adjacency += hash_adjacency_seconds;
    profile->hash_assemble += hash_walk_seconds + hash_adjacency_seconds;

    phase_begin = omp_get_wtime();
    contig_graph.clear();
    contig_graph.set_kmer_size(kmer_size);
    if (has_explicit_hash_adjacency) {
      contig_graph.InitializeWithAdjacency(hash_unitigs,
                                           hash_unitig_neighbors,
                                           hash_unitig_edges);
    } else {
      contig_graph.Initialize(hash_unitigs);
    }
    profile->contig_initialize += omp_get_wtime() - phase_begin;

    phase_begin = omp_get_wtime();
    const int64_t removed_deadends =
        contig_graph.RemoveDeadEnd(kmer_size * 2);
    profile->dead_end += omp_get_wtime() - phase_begin;

    phase_begin = omp_get_wtime();
    const int64_t removed_bubbles = contig_graph.RemoveBubble();
    profile->bubble += omp_get_wtime() - phase_begin;

    phase_begin = omp_get_wtime();
    bool coverage_changed = false;
    contig_graph.IterateCoverage(kmer_size * 2, 1, threshold, 1.1,
                                 &coverage_changed);
    profile->coverage_clean += omp_get_wtime() - phase_begin;
    ++profile->cleaning_rounds;
    profile->deadend_changed_rounds += removed_deadends != 0;
    profile->bubble_changed_rounds += removed_bubbles != 0;
    profile->coverage_changed_rounds += coverage_changed;

    phase_begin = omp_get_wtime();
    contig_graph.Assemble(out_contigs, out_contig_infos);
    profile->contig_assemble += omp_get_wtime() - phase_begin;

    if (out_contigs.size() == 1) {
      break;
    }
  }

  if (use_multi_k_index) {
    const MultiKReadIndex::Stats &stats = multi_k_index.stats();
    profile->multik_build += stats.build;
    profile->multik_finish += stats.finish_initial;
    profile->multik_advance += stats.advance;
    profile->multik_aggregate += stats.aggregate;
    profile->multik_replay += stats.replay;
  }
  if (use_incremental_read_graph) {
    const IncrementalKReadGraph::Stats &stats =
        incremental_read_graph.stats();
    profile->incremental_occurrences += stats.occurrences;
    profile->incremental_path_replays += stats.path_replays;
    profile->incremental_primary_inserts += stats.primary_inserts;
    profile->incremental_primary_hits += stats.primary_hits;
    profile->incremental_fork_lookups += stats.fork_lookups;
  }

  profile->branch_transition_hits += hash_graph.DebugBranchTransitionHits();
  profile->branch_transition_misses +=
      hash_graph.DebugBranchTransitionMisses();

}

std::vector<TInsertSize> EstimateInsertSize(
    const HashMapper &mapper, const SequenceLibCollection &lib_collection) {
  std::vector<TInsertSize> insert_sizes(lib_collection.size());
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto lib = lib_collection.GetLib(lib_id);

    if (!lib.IsPaired()) {
      continue;
    }

    Histgram<int> insert_hist;
    std::vector<Histgram<int, NullMutex>> thread_insert_hist(
        omp_get_max_threads());
    const size_t min_hist_size_for_estimation = 1u << 18;
    size_t processed_reads = 0;

    while (insert_hist.size() < min_hist_size_for_estimation &&
           processed_reads < lib.seq_count()) {
      size_t start_read_id = processed_reads;
      processed_reads = std::min(min_hist_size_for_estimation + start_read_id,
                                 lib.seq_count());

      for (auto &local_hist : thread_insert_hist) {
        local_hist.clear();
      }

#pragma omp parallel for
      for (size_t i = start_read_id; i < processed_reads; i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        HashMapper::EndpointSeedWitness witness1;
        HashMapper::EndpointSeedWitness witness2;
        const bool endpoint1 = mapper.MayMapToEndpoint(seq1, &witness1);
        const bool endpoint2 =
            endpoint1 ? false
                      : mapper.MayMapToEndpoint(seq2, &witness2);
        if (!endpoint1 && !endpoint2) {
          continue;
        }
        auto rec1 = mapper.TryMap(seq1, witness1.valid ? &witness1 : nullptr);
        auto rec2 = mapper.TryMap(seq2, witness2.valid ? &witness2 : nullptr);
        if (rec1.valid && rec2.valid) {
          if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
            int insert_size;

            if (rec1.strand == 0) {
              insert_size = rec2.contig_to + seq2.length() - rec2.query_to -
                            (rec1.contig_from - rec1.query_from);
            } else {
              insert_size = rec1.contig_to + seq1.length() - rec1.query_to -
                            (rec2.contig_from - rec2.query_from);
            }

            if (insert_size >= (int)seq1.length() &&
                insert_size >= (int)seq2.length()) {
              thread_insert_hist[omp_get_thread_num()].insert(insert_size);
            }
          }
        }
      }

      for (const auto &local_hist : thread_insert_hist) {
        insert_hist.MergeFrom(local_hist);
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes[lib_id] = TInsertSize(insert_hist.mean(), insert_hist.sd());

    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id,
          insert_sizes[lib_id].first, insert_sizes[lib_id].second);
  }

  return insert_sizes;
}

std::vector<TInsertSize> EstimateInsertSizeMapped(
    const HashMapper &mapper, const SequenceLibCollection &lib_collection,
    const MappedReadFile &reads) {
  const size_t num_libs = lib_collection.size();
  const int num_threads = std::max(1, omp_get_max_threads());
  std::vector<TInsertSize> insert_sizes(num_libs);
  const auto &chunks = reads.index().chunks;
  const uint64_t sample_reads_per_wave = uint64_t{1} << 18u;
  std::vector<Histgram<int>> insert_hists(num_libs);
  std::vector<Histgram<int, NullMutex>> thread_insert_hists(
      static_cast<size_t>(num_threads) * num_libs);
  std::vector<uint64_t> processed_reads(num_libs, 0u);
  std::vector<uint64_t> range_begins(num_libs, 0u);
  std::vector<uint64_t> range_ends(num_libs, 0u);
  std::vector<std::vector<uint64_t>> pair_word_offsets(num_libs);
  std::vector<uint64_t> pair_prefix(num_libs + 1u, 0u);
  std::vector<uint8_t> touched_chunks(chunks.size(), 0u);

  while (true) {
    bool any_active = false;
    pair_prefix[0] = 0u;
    for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
      const auto &lib = lib_collection.GetLib(lib_id);
      const bool active =
          lib.IsPaired() &&
          insert_hists[lib_id].size() < sample_reads_per_wave &&
          processed_reads[lib_id] < lib.seq_count();
      if (active) {
        any_active = true;
        const uint64_t local_begin = processed_reads[lib_id];
        processed_reads[lib_id] = std::min<uint64_t>(
            local_begin + sample_reads_per_wave, lib.seq_count());
        range_begins[lib_id] = lib.global_begin() + local_begin;
        range_ends[lib_id] = lib.global_begin() + processed_reads[lib_id];
        const size_t num_pairs = static_cast<size_t>(
            (range_ends[lib_id] - range_begins[lib_id]) / 2u);
        pair_word_offsets[lib_id].resize(num_pairs);
      } else {
        range_begins[lib_id] = range_ends[lib_id] = 0u;
        pair_word_offsets[lib_id].clear();
      }
      pair_prefix[lib_id + 1u] =
          pair_prefix[lib_id] + pair_word_offsets[lib_id].size();
    }
    if (!any_active) break;

    std::fill(touched_chunks.begin(), touched_chunks.end(), uint8_t{0});
    for (auto &hist : thread_insert_hists) hist.clear();

    // Locate only the first word of each pair. The mate is the immediately
    // following packed record, halving the offset tape versus the previous
    // per-read representation.
#pragma omp parallel for schedule(static)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
      const PackedReadChunk &chunk = chunks[chunk_id];
      bool touched = false;
      for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
        const uint64_t range_begin = range_begins[lib_id];
        const uint64_t range_end = range_ends[lib_id];
        if (range_begin == range_end) continue;
        uint64_t first = std::max(chunk.read_begin, range_begin);
        const uint64_t end = std::min(chunk.read_end, range_end);
        const uint64_t lib_begin =
            lib_collection.GetLib(lib_id).global_begin();
        if (((first - lib_begin) & 1u) != 0u) ++first;
        if (first >= end) continue;
        touched = true;
        const uint32_t *cursor = reads.LocateRead(chunk, first);
        for (uint64_t read_id = first; read_id < end; read_id += 2u) {
          pair_word_offsets[lib_id][(read_id - range_begin) / 2u] =
              static_cast<uint64_t>(cursor - reads.words());
          MappedReadFile::Next(&cursor);
          MappedReadFile::Next(&cursor);
        }
      }
      touched_chunks[chunk_id] = static_cast<uint8_t>(touched);
    }

    const uint64_t total_pairs = pair_prefix.back();
#pragma omp parallel for schedule(static)
    for (int64_t flat_pair = 0;
         flat_pair < static_cast<int64_t>(total_pairs); ++flat_pair) {
      const size_t lib_id = static_cast<size_t>(
          std::upper_bound(pair_prefix.begin(), pair_prefix.end(),
                           static_cast<uint64_t>(flat_pair)) -
          pair_prefix.begin() - 1u);
      const uint64_t local_pair =
          static_cast<uint64_t>(flat_pair) - pair_prefix[lib_id];
      const uint64_t read_id = range_begins[lib_id] + local_pair * 2u;
      const uint32_t *cursor =
          reads.words() + pair_word_offsets[lib_id][local_pair];
      const PackedReadRecord seq1 = MappedReadFile::Next(&cursor);
      const PackedReadRecord seq2 = MappedReadFile::Next(&cursor);
      const MappingRecord rec1 =
          mapper.TryMap(seq1.words, seq1.length, read_id);
      const MappingRecord rec2 =
          mapper.TryMap(seq2.words, seq2.length, read_id + 1u);
      if (rec1.valid && rec2.valid && rec1.contig_id == rec2.contig_id &&
          rec1.strand != rec2.strand) {
        int insert_size;
        if (rec1.strand == 0) {
          insert_size = rec2.contig_to + seq2.length - rec2.query_to -
                        (rec1.contig_from - rec1.query_from);
        } else {
          insert_size = rec1.contig_to + seq1.length - rec1.query_to -
                        (rec2.contig_from - rec2.query_from);
        }
        if (insert_size >= static_cast<int>(seq1.length) &&
            insert_size >= static_cast<int>(seq2.length)) {
          thread_insert_hists[static_cast<size_t>(omp_get_thread_num()) *
                                  num_libs +
                              lib_id]
              .insert(insert_size);
        }
      }
    }

    for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
      for (int thread = 0; thread < num_threads; ++thread) {
        insert_hists[lib_id].MergeFrom(
            thread_insert_hists[static_cast<size_t>(thread) * num_libs +
                                lib_id]);
      }
    }

#pragma omp parallel for schedule(static)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
      if (touched_chunks[chunk_id] != 0u) reads.DropChunk(chunks[chunk_id]);
    }
  }

  for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
    const auto &lib = lib_collection.GetLib(lib_id);
    if (!lib.IsPaired()) continue;
    insert_hists[lib_id].Trim(0.01);
    insert_sizes[lib_id] =
        TInsertSize(insert_hists[lib_id].mean(), insert_hists[lib_id].sd());
    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id,
          insert_sizes[lib_id].first, insert_sizes[lib_id].second);
  }
  return insert_sizes;
}

int32_t LocalRange(const SequenceLib &lib, const TInsertSize &insert_size) {
  int32_t local_range = lib.GetMaxLength() - 1;

  if (lib.IsPaired() && insert_size.first >= lib.GetMaxLength()) {
    local_range = std::min(2 * insert_size.first,
                           insert_size.first + 3 * insert_size.second);
  }

  if (local_range > kMaxLocalRange) {
    local_range = kMaxLocalRange;
  }

  return local_range;
}

int32_t GetMaxLocalRange(const SequenceLibCollection &lib_collection,
                         const std::vector<TInsertSize> &insert_sizes) {
  int32_t max_local_range = 0;
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    max_local_range =
        std::max(max_local_range, LocalRange(lib, insert_sizes[lib_id]));
  }
  return max_local_range;
}

unsigned LegacySingleEndMappingCopies() {
  // Upstream MEGAHIT's single-end loop used `omp parallel` without `for`.
  // Every worker therefore inserted the same mapping, after which local
  // assembly retained at most the first three entries at a contig position.
  // Replay only that observable multiplicity while continuing to divide the
  // read scan among the actual OpenMP team.
  return static_cast<unsigned>(std::min(omp_get_num_threads(), 3));
}

void MapToContigs(const HashMapper &mapper,
                  const SequenceLibCollection &lib_collection,
                  const std::vector<TInsertSize> &insert_sizes,
                  MappingResultCollector *collector) {
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    int32_t local_range = LocalRange(lib, insert_sizes[lib_id]);
    bool is_paired = lib.IsPaired();

    size_t num_added = 0, num_mapped = 0;

    if (is_paired) {
#pragma omp parallel for reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        HashMapper::EndpointSeedWitness witness1;
        HashMapper::EndpointSeedWitness witness2;
        const bool endpoint1 = mapper.MayMapToEndpoint(seq1, &witness1);
        const bool endpoint2 =
            endpoint1 ? false
                      : mapper.MayMapToEndpoint(seq2, &witness2);
        if (!endpoint1 && !endpoint2) {
          continue;
        }
        auto rec1 = mapper.TryMap(seq1, witness1.valid ? &witness1 : nullptr);
        auto rec2 = mapper.TryMap(seq2, witness2.valid ? &witness2 : nullptr);

        if (rec1.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec1.contig_id).length();
          num_added += collector->AddSingle(rec1, contig_len, seq1.length(),
                                            local_range);
          num_added += collector->AddMate(rec1, rec2, contig_len, seq2.id(),
                                          local_range);
        }

        if (rec2.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec2.contig_id).length();
          num_added += collector->AddSingle(rec2, contig_len, seq2.length(),
                                            local_range);
          num_added += collector->AddMate(rec2, rec1, contig_len, seq1.id(),
                                          local_range);
        }
      }
    } else {
#pragma omp parallel reduction(+ : num_added, num_mapped)
      {
        const unsigned mapping_copies = LegacySingleEndMappingCopies();
#pragma omp for
        for (size_t i = 0; i < lib.seq_count(); ++i) {
          auto seq = lib.GetSequenceView(i);
          HashMapper::EndpointSeedWitness witness;
          if (!mapper.MayMapToEndpoint(seq, &witness)) continue;
          auto rec = mapper.TryMap(seq, &witness);

          if (rec.valid) {
            ++num_mapped;
            const int32_t contig_len =
                mapper.refseq().GetSeqView(rec.contig_id).length();
            for (unsigned copy = 0; copy < mapping_copies; ++copy) {
              num_added += collector->AddSingle(rec, contig_len, seq.length(),
                                                local_range);
            }
          }
        }
      }
    }

    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads to local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
  }
}

void MapToContigsMapped(const HashMapper &mapper,
                        const SequenceLibCollection &lib_collection,
                        const std::vector<TInsertSize> &insert_sizes,
                        const MappedReadFile &reads,
                        const MappedLocalCandidateFile *candidate_file,
                        const LocalSeedPositions *seed_positions,
                        bool mapped_assembly,
                        MappingResultCollector *collector) {
  const auto &chunks = reads.index().chunks;
  const size_t num_libs = lib_collection.size();
  std::vector<uint64_t> lib_begins(num_libs);
  std::vector<uint64_t> lib_ends(num_libs);
  std::vector<int32_t> local_ranges(num_libs);
  std::vector<uint8_t> paired(num_libs);
  for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
    const auto &lib = lib_collection.GetLib(lib_id);
    lib_begins[lib_id] = lib.global_begin();
    lib_ends[lib_id] = lib.global_end();
    local_ranges[lib_id] = LocalRange(lib, insert_sizes[lib_id]);
    paired[lib_id] = static_cast<uint8_t>(lib.IsPaired());
  }

  const size_t thread_stride = std::max<size_t>(1u, num_libs);
  std::vector<uint64_t> thread_added(
      static_cast<size_t>(omp_get_max_threads()) * thread_stride, 0u);
  std::vector<uint64_t> thread_mapped(
      static_cast<size_t>(omp_get_max_threads()) * thread_stride, 0u);
  std::vector<uint64_t> thread_mapping_attempts(
      static_cast<size_t>(omp_get_max_threads()) * thread_stride, 0u);
  std::vector<uint64_t> thread_selected_pairs(
      static_cast<size_t>(omp_get_max_threads()) * thread_stride, 0u);
  const bool validate_gate =
      std::getenv("MEGAHIT_VALIDATE_LOCAL_MINIMIZER_GATE") != nullptr;
  const auto endpoint_gate = [&](const PackedReadRecord &read, uint64_t read_id,
                                  HashMapper::EndpointSeedWitness *witness) {
    if (seed_positions == nullptr)
      return mapper.MayMapToEndpoint(read.words, read.length, witness);
    const bool result = mapper.MayMapToEndpointAnchors(
        read.words, read.length, seed_positions->Get(read_id), witness);
    if (validate_gate) {
      HashMapper::EndpointSeedWitness reference;
      const bool expected =
          mapper.MayMapToEndpoint(read.words, read.length, &reference);
      if (expected != result ||
          (result && (witness->index_value != reference.index_value ||
                      witness->end_position != reference.end_position ||
                      witness->query_strand != reference.query_strand))) {
        xfatal("Local minimizer gate differs at read {}\n", read_id);
      }
    }
    return result;
  };

#pragma omp parallel
  {
    const size_t thread_offset =
        static_cast<size_t>(omp_get_thread_num()) * thread_stride;
    const unsigned single_end_mapping_copies =
        LegacySingleEndMappingCopies();
#pragma omp for schedule(dynamic, 1)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
      const PackedReadChunk &chunk = chunks[chunk_id];
      for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
        const uint64_t lib_begin = lib_begins[lib_id];
        const uint64_t lib_end = lib_ends[lib_id];
        uint64_t first = std::max(chunk.read_begin, lib_begin);
        const uint64_t end = std::min(chunk.read_end, lib_end);
        if (first >= end) continue;
        uint64_t &num_added = thread_added[thread_offset + lib_id];
        uint64_t &num_mapped = thread_mapped[thread_offset + lib_id];
        uint64_t &mapping_attempts =
            thread_mapping_attempts[thread_offset + lib_id];
        uint64_t &selected_pairs =
            thread_selected_pairs[thread_offset + lib_id];
        const int32_t local_range = local_ranges[lib_id];

        if (paired[lib_id] != 0u) {
          if (((first - lib_begin) & 1u) != 0u) ++first;
          if (first >= end) continue;
          const uint32_t *cursor = reads.LocateRead(chunk, first);
          auto candidate_range =
              candidate_file == nullptr
                  ? std::make_pair(static_cast<const uint64_t *>(nullptr),
                                   static_cast<const uint64_t *>(nullptr))
                  : candidate_file->Range(
                        static_cast<uint64_t>(cursor - reads.words()),
                        chunk.word_end);
          auto candidate = candidate_range.first;
          const auto candidate_end = candidate_range.second;
          auto is_candidate = [&](uint64_t word_offset) mutable {
            if (candidate_file == nullptr) return true;
            while (candidate != candidate_end && *candidate < word_offset) {
              ++candidate;
            }
            if (candidate == candidate_end || *candidate != word_offset) {
              return false;
            }
            ++candidate;
            return true;
          };
          auto map_pair = [&](const PackedReadRecord &seq1,
                              const PackedReadRecord &seq2,
                              uint64_t read_id,
                              const HashMapper::EndpointSeedWitness &witness1,
                              const HashMapper::EndpointSeedWitness &witness2) {
            ++selected_pairs;
            mapping_attempts += 2u;
            const uint64_t query_id1 = mapped_assembly
                ? static_cast<uint64_t>(seq1.words - reads.words() - 1u)
                : read_id;
            const uint64_t query_id2 = mapped_assembly
                ? static_cast<uint64_t>(seq2.words - reads.words() - 1u)
                : read_id + 1u;
            const MappingRecord rec1 =
                mapper.TryMap(seq1.words, seq1.length, query_id1,
                              witness1.valid ? &witness1 : nullptr);
            const MappingRecord rec2 =
                mapper.TryMap(seq2.words, seq2.length, query_id2,
                              witness2.valid ? &witness2 : nullptr);

            if (rec1.valid) {
              ++num_mapped;
              const int32_t contig_len =
                  mapper.refseq().GetSeqView(rec1.contig_id).length();
              num_added += collector->AddSingle(
                  rec1, contig_len, seq1.length, local_range);
              num_added += collector->AddMate(
                  rec1, rec2, contig_len, query_id2, local_range);
            }
            if (rec2.valid) {
              ++num_mapped;
              const int32_t contig_len =
                  mapper.refseq().GetSeqView(rec2.contig_id).length();
              num_added += collector->AddSingle(
                  rec2, contig_len, seq2.length, local_range);
              num_added += collector->AddMate(
                  rec2, rec1, contig_len, query_id1, local_range);
            }
          };
          for (uint64_t read_id = first; read_id < end; read_id += 2u) {
            const uint64_t seq1_word_offset =
                static_cast<uint64_t>(cursor - reads.words());
            const PackedReadRecord seq1 = MappedReadFile::Next(&cursor);
            const uint64_t seq2_word_offset =
                static_cast<uint64_t>(cursor - reads.words());
            const PackedReadRecord seq2 = MappedReadFile::Next(&cursor);
            const bool select1 = is_candidate(seq1_word_offset);
            const bool select2 = is_candidate(seq2_word_offset);
            if (!select1 && !select2) continue;
            HashMapper::EndpointSeedWitness witness1;
            HashMapper::EndpointSeedWitness witness2;
            const bool endpoint1 = endpoint_gate(seq1, read_id, &witness1);
            const bool endpoint2 =
                endpoint1
                    ? false
                    : endpoint_gate(seq2, read_id + 1u, &witness2);
            if (!endpoint1 && !endpoint2) continue;
            map_pair(seq1, seq2, read_id, witness1, witness2);
          }
        } else {
          const uint32_t *cursor = reads.LocateRead(chunk, first);
          auto candidate_range =
              candidate_file == nullptr
                  ? std::make_pair(static_cast<const uint64_t *>(nullptr),
                                   static_cast<const uint64_t *>(nullptr))
                  : candidate_file->Range(
                        static_cast<uint64_t>(cursor - reads.words()),
                        chunk.word_end);
          auto candidate = candidate_range.first;
          const auto candidate_end = candidate_range.second;
          auto is_candidate = [&](uint64_t word_offset) mutable {
            if (candidate_file == nullptr) return true;
            while (candidate != candidate_end && *candidate < word_offset) {
              ++candidate;
            }
            if (candidate == candidate_end || *candidate != word_offset) {
              return false;
            }
            ++candidate;
            return true;
          };
          for (uint64_t read_id = first; read_id < end; ++read_id) {
            const uint64_t word_offset =
                static_cast<uint64_t>(cursor - reads.words());
            const PackedReadRecord seq = MappedReadFile::Next(&cursor);
            if (!is_candidate(word_offset)) continue;
            HashMapper::EndpointSeedWitness witness;
            if (!endpoint_gate(seq, read_id, &witness)) {
              continue;
            }
            ++mapping_attempts;
            const MappingRecord rec =
                mapper.TryMap(seq.words, seq.length,
                                mapped_assembly ? word_offset : read_id,
                                &witness);
            if (rec.valid) {
              ++num_mapped;
              const int32_t contig_len =
                  mapper.refseq().GetSeqView(rec.contig_id).length();
              for (unsigned copy = 0; copy < single_end_mapping_copies;
                   ++copy) {
                num_added += collector->AddSingle(
                    rec, contig_len, seq.length, local_range);
              }
            }
          }
        }
      }
      if (!mapped_assembly) reads.DropChunk(chunk);
      if (seed_positions != nullptr)
        seed_positions->Drop(chunk.read_begin, chunk.read_end - chunk.read_begin);
    }
  }

  for (size_t lib_id = 0; lib_id < num_libs; ++lib_id) {
    uint64_t num_added = 0;
    uint64_t num_mapped = 0;
    uint64_t mapping_attempts = 0;
    uint64_t selected_pairs = 0;
    for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
      const size_t index = static_cast<size_t>(thread) * thread_stride + lib_id;
      num_added += thread_added[index];
      num_mapped += thread_mapped[index];
      mapping_attempts += thread_mapping_attempts[index];
      selected_pairs += thread_selected_pairs[index];
    }
    const auto &lib = lib_collection.GetLib(lib_id);
    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads to local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
    xinfo("Lib {} endpoint gate selected {} pairs; exact mapper attempted {} "
          "reads ({.2} pct of the library)\n",
          lib_id, selected_pairs, mapping_attempts,
          lib.seq_count() == 0
              ? 0.0
              : 100.0 * mapping_attempts / lib.seq_count());
  }
}
void AssembleAndOutput(const HashMapper &mapper, const SeqPackage &read_pkg,
                       MappingResultCollector &result_collector,
                       const std::string &output_file,
                       const int32_t local_range,
                       const LocalAsmOption &opt,
                       const uint32_t *mapped_words = nullptr,
                       unsigned mapped_max_read_length = 0u) {
  const unsigned max_read_length = mapped_words != nullptr
      ? mapped_max_read_length : read_pkg.max_length();
  const size_t min_num_reads = max_read_length > 0 ?
      local_range / max_read_length : 1;
  xinfo("Minimum number of reads to do local assembly: {}\n", min_num_reads);

  Sequence contig_end;
  // Both graph workspaces are private to an OpenMP worker and persist across
  // dynamically scheduled endpoints.  clear() resets logical contents while
  // retaining flat-array capacity, eliminating hundreds of thousands of
  // allocator/zero-initialization cycles without sharing graph state.
  HashGraph hash_graph;
  MultiKReadIndex multi_k_index;
  IncrementalKReadGraph incremental_read_graph;
  ContigGraph contig_graph;
  std::vector<ContigGraphVertex> hash_unitigs;
  std::vector<uint32_t> hash_unitig_neighbors;
  // Private to an OpenMP worker and retained across dynamic endpoint tasks;
  // clear() keeps the backing capacity while replacing millions of individual
  // small string allocations with four contiguous arrays.
  LocalPackedReads reads;
  std::vector<Sequence> out_contigs;
  std::vector<ContigInfo> out_contig_infos;

  double endpoint_gather_cpu_seconds = 0.0;
  double endpoint_assembly_cpu_seconds = 0.0;
  double endpoint_gather_max_seconds = 0.0;
  uint64_t endpoint_selected_reads = 0;
  uint64_t endpoint_selected_bases = 0;
  uint64_t endpoint_duplicate_reads = 0;
  uint64_t endpoint_duplicate_kmer_work = 0;
  const bool diagnose_read_duplicates =
      std::getenv("MEGAHIT_DIAGNOSE_LOCAL_READ_DUPLICATES") != nullptr;
  std::vector<LocalAssemblyProfile> thread_profiles(omp_get_max_threads());
  const size_t num_k_values =
      opt.step == 0 ? 0 : (opt.kmax - opt.kmin) / opt.step + 1u;
  std::vector<std::vector<uint64_t>> thread_k_round_survival(
      omp_get_max_threads(), std::vector<uint64_t>(num_k_values, 0));
  std::vector<std::vector<double>> thread_k_read_graph_seconds(
      omp_get_max_threads(), std::vector<double>(num_k_values, 0.0));
  std::vector<std::vector<uint64_t>> thread_k_branch_hits(
      omp_get_max_threads(), std::vector<uint64_t>(num_k_values, 0));
  std::vector<std::vector<uint64_t>> thread_k_branch_misses(
      omp_get_max_threads(), std::vector<uint64_t>(num_k_values, 0));

  ContigWriter local_contig_writer(output_file);

  struct LocalAssemblyTask {
    uint64_t contig_id;
    uint64_t mapping_count;
    uint8_t strand;
  };
  struct LocalInputFingerprint {
    uint64_t low{0};
    uint64_t high{0};
    uint64_t reads_low{0};
    uint64_t reads_high{0};
    uint64_t potential_kmers{0};
    uint64_t actual_kmers{0};
    uint32_t read_count{0};
    uint32_t endpoint_length{0};
    uint16_t completed_rounds{0};
  };
  std::vector<LocalAssemblyTask> tasks;
  tasks.reserve(mapper.refseq().seq_count());
  for (uint64_t cid = 0; cid < mapper.refseq().seq_count(); ++cid) {
    for (uint8_t strand = 0; strand < 2; ++strand) {
      const size_t mapping_count =
          result_collector.GetMappingResults(cid, strand).size();
      if (mapping_count > min_num_reads) {
        tasks.push_back({cid, mapping_count, strand});
      }
    }
  }
  // Each contig end is independent.  Scheduling the largest read sets first
  // prevents a handful of expensive endpoints from forming a long serial
  // tail after the graph has otherwise drained, and splitting the two ends
  // doubles useful task granularity without altering an endpoint's inputs.
  std::sort(tasks.begin(), tasks.end(),
            [](const LocalAssemblyTask &lhs, const LocalAssemblyTask &rhs) {
              if (lhs.mapping_count != rhs.mapping_count) {
                return lhs.mapping_count > rhs.mapping_count;
              }
              if (lhs.contig_id != rhs.contig_id) {
                return lhs.contig_id < rhs.contig_id;
              }
              return lhs.strand < rhs.strand;
            });
  xinfo("Scheduling {} independent local-assembly endpoint tasks\n",
        tasks.size());

  // These arrays are diagnostic and scheduler inputs for the next stage.
  // Keeping one compact scalar per endpoint lets us distinguish a genuinely
  // serial heavy tail from low utilization inside otherwise plentiful tasks.
  std::vector<double> task_begin_seconds(tasks.size());
  std::vector<double> task_end_seconds(tasks.size());
  std::vector<double> task_elapsed_seconds(tasks.size());
  std::vector<uint64_t> task_potential_kmers(tasks.size());
  std::vector<uint64_t> task_actual_kmers(tasks.size());
  std::vector<uint64_t> task_unique_kmers(tasks.size());
  std::vector<uint16_t> task_k_rounds(tasks.size());
  const char *const fingerprint_path =
      std::getenv("MEGAHIT_LOCAL_INPUT_FINGERPRINT");
  std::vector<LocalInputFingerprint> task_fingerprints;
  if (fingerprint_path != nullptr && *fingerprint_path != '\0') {
    task_fingerprints.resize(tasks.size());
  }
  const double task_phase_begin = omp_get_wtime();

#pragma omp parallel for private(hash_graph, multi_k_index,                    \
                                 incremental_read_graph,                       \
                                 contig_graph, hash_unitigs,                    \
                                 hash_unitig_neighbors,                        \
                                 contig_end, reads, out_contigs,               \
                                 out_contig_infos)                              \
    reduction(+ : endpoint_gather_cpu_seconds,                             \
              endpoint_assembly_cpu_seconds, endpoint_selected_reads,      \
              endpoint_selected_bases, endpoint_duplicate_reads,           \
              endpoint_duplicate_kmer_work)                                \
    reduction(max : endpoint_gather_max_seconds) schedule(dynamic, 1)
  for (int64_t task_id = 0; task_id < static_cast<int64_t>(tasks.size());
       ++task_id) {
    task_begin_seconds[task_id] = omp_get_wtime() - task_phase_begin;
    const uint64_t cid = tasks[task_id].contig_id;
    const uint8_t strand = tasks[task_id].strand;
    auto contig_view = mapper.refseq().GetSeqView(cid);
    int cl = contig_view.length();

    auto mapping_rslts = result_collector.GetMappingResults(cid, strand);

    const double gather_begin = omp_get_wtime();

    // Collect local reads in historical mapping order.  Packing changes only
    // representation; the subsequent graph sees the same read/position order.
    reads.clear();
    uint64_t last_mapping_pos = -1;
    int pos_count = 0;

    for (const auto &encoded_rslt : mapping_rslts) {
      uint64_t pos = MappingResultCollector::GetContigAbsPos(encoded_rslt);
      pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
      last_mapping_pos = pos;

      if (pos_count <= 3) {
        const uint64_t read_id = MappingResultCollector::GetReadId(encoded_rslt);
        if (mapped_words != nullptr) {
          const uint32_t *record = mapped_words + read_id;
          endpoint_selected_bases += *record;
          reads.Add(record + 1u, 0u, *record);
        } else {
          auto read_view = read_pkg.GetSeqView(read_id);
          endpoint_selected_bases += read_view.length();
          reads.Add(read_view);
        }
      }
    }
    endpoint_selected_reads += reads.size();
    task_potential_kmers[task_id] =
        reads.EstimateKmerWork(opt.kmin, opt.kmax, opt.step);
    if (diagnose_read_duplicates) {
      const auto duplicates =
          reads.ExactDuplicateWork(opt.kmin, opt.kmax, opt.step);
      endpoint_duplicate_reads += duplicates.first;
      endpoint_duplicate_kmer_work += duplicates.second;
    }
    const double gather_seconds = omp_get_wtime() - gather_begin;
    endpoint_gather_cpu_seconds += gather_seconds;
    endpoint_gather_max_seconds =
        std::max(endpoint_gather_max_seconds, gather_seconds);

    if (strand == 0) {
      const int end = std::min(local_range, cl);
      contig_end.resize(end);
      for (int j = 0; j < end; ++j) {
        contig_end[j] = contig_view.base_at(j);
      }
    } else {
      const int begin = std::max(0, cl - local_range);
      contig_end.resize(cl - begin);
      for (int j = begin; j < cl; ++j) {
        contig_end[j - begin] = contig_view.base_at(j);
      }
    }

    if (!task_fingerprints.empty()) {
      const auto fingerprint = reads.InputFingerprint(contig_end);
      const auto reads_fingerprint = reads.ReadFingerprint();
      LocalInputFingerprint &record = task_fingerprints[task_id];
      record.low = fingerprint.first;
      record.high = fingerprint.second;
      record.reads_low = reads_fingerprint.first;
      record.reads_high = reads_fingerprint.second;
      record.potential_kmers = task_potential_kmers[task_id];
      record.read_count = static_cast<uint32_t>(reads.size());
      record.endpoint_length = static_cast<uint32_t>(contig_end.size());
    }

    out_contigs.clear();
    const double assembly_begin = omp_get_wtime();
    LocalAssemblyProfile &local_profile =
        thread_profiles[omp_get_thread_num()];
    const uint64_t rounds_before = local_profile.k_rounds;
    const uint64_t raw_kmers_before = local_profile.raw_kmers;
    const uint64_t raw_vertices_before = local_profile.raw_vertices;
    LaunchIDBA(reads, contig_end, out_contigs,
               out_contig_infos, opt.kmin, opt.kmax, opt.step, hash_graph,
               contig_graph, hash_unitigs, hash_unitig_neighbors,
               multi_k_index, incremental_read_graph,
               &local_profile, &thread_k_round_survival[omp_get_thread_num()],
               &thread_k_read_graph_seconds[omp_get_thread_num()],
               &thread_k_branch_hits[omp_get_thread_num()],
               &thread_k_branch_misses[omp_get_thread_num()],
               task_potential_kmers[task_id]);
    const double assembly_seconds = omp_get_wtime() - assembly_begin;
    endpoint_assembly_cpu_seconds += assembly_seconds;
    task_elapsed_seconds[task_id] = assembly_seconds;
    task_k_rounds[task_id] = static_cast<uint16_t>(
        local_profile.k_rounds - rounds_before);
    task_actual_kmers[task_id] = local_profile.raw_kmers - raw_kmers_before;
    task_unique_kmers[task_id] =
        local_profile.raw_vertices - raw_vertices_before;
    if (!task_fingerprints.empty()) {
      task_fingerprints[task_id].actual_kmers = task_actual_kmers[task_id];
      task_fingerprints[task_id].completed_rounds = task_k_rounds[task_id];
    }

    for (uint64_t j = 0; j < out_contigs.size(); ++j) {
      if (out_contigs[j].size() > opt.min_contig_len &&
          out_contigs[j].size() > opt.kmax) {
        auto str = out_contigs[j].str();
        local_contig_writer.WriteLocalContig(str, cid, strand, j);
      }
    }
    task_end_seconds[task_id] = omp_get_wtime() - task_phase_begin;
  }

  const double task_phase_elapsed = omp_get_wtime() - task_phase_begin;

  if (!task_fingerprints.empty()) {
    std::ofstream output(fingerprint_path, std::ios::out | std::ios::trunc);
    if (!output) {
      throw std::runtime_error("cannot write local input fingerprint file");
    }
    output << "low\thigh\treads_low\treads_high\tpotential\tactual\treads"
              "\tendpoint\trounds\n";
    output << std::hex;
    for (const LocalInputFingerprint &record : task_fingerprints) {
      output << record.low << '\t' << record.high << '\t'
             << record.reads_low << '\t' << record.reads_high << std::dec
             << '\t' << record.potential_kmers << '\t'
             << record.actual_kmers << '\t'
             << record.read_count << '\t' << record.endpoint_length << '\t'
             << record.completed_rounds << '\n' << std::hex;
    }
    output.close();
    xinfo("Wrote {} exact-input diagnostic fingerprints to {s}\n",
          task_fingerprints.size(), fingerprint_path);
  }

  xinfo("Endpoint read gather: {.6}s CPU, max {.6}s, {} reads, {} bases; "
        "endpoint graph work: {.6}s CPU\n",
        endpoint_gather_cpu_seconds, endpoint_gather_max_seconds,
        endpoint_selected_reads, endpoint_selected_bases,
        endpoint_assembly_cpu_seconds);
  if (diagnose_read_duplicates) {
    const uint64_t total_potential_work =
        std::accumulate(task_potential_kmers.begin(),
                        task_potential_kmers.end(), uint64_t{0});
    xinfo("Exact duplicate local reads: {} / {} (percent {.4}); removable "
          "inner-k occurrences {} / {} (percent {.4})\n",
          endpoint_duplicate_reads, endpoint_selected_reads,
          endpoint_selected_reads == 0u
              ? 0.0
              : 100.0 * endpoint_duplicate_reads / endpoint_selected_reads,
          endpoint_duplicate_kmer_work, total_potential_work,
          total_potential_work == 0u
              ? 0.0
              : 100.0 * endpoint_duplicate_kmer_work /
                    total_potential_work);
  }
  LocalAssemblyProfile profile;
  for (const auto &thread_profile : thread_profiles) {
    profile.read_graph_build += thread_profile.read_graph_build;
    profile.coverage += thread_profile.coverage;
    profile.anchor_insert += thread_profile.anchor_insert;
    profile.hash_assemble += thread_profile.hash_assemble;
    profile.hash_walk += thread_profile.hash_walk;
    profile.hash_adjacency += thread_profile.hash_adjacency;
    profile.contig_initialize += thread_profile.contig_initialize;
    profile.dead_end += thread_profile.dead_end;
    profile.bubble += thread_profile.bubble;
    profile.coverage_clean += thread_profile.coverage_clean;
    profile.contig_assemble += thread_profile.contig_assemble;
    profile.multik_build += thread_profile.multik_build;
    profile.multik_finish += thread_profile.multik_finish;
    profile.multik_advance += thread_profile.multik_advance;
    profile.multik_aggregate += thread_profile.multik_aggregate;
    profile.multik_replay += thread_profile.multik_replay;
    profile.k_rounds += thread_profile.k_rounds;
    profile.raw_kmers += thread_profile.raw_kmers;
    profile.raw_vertices += thread_profile.raw_vertices;
    profile.single_hash_unitigs += thread_profile.single_hash_unitigs;
    profile.long_single_hash_unitigs +=
        thread_profile.long_single_hash_unitigs;
    profile.isolated_long_single_hash_unitigs +=
        thread_profile.isolated_long_single_hash_unitigs;
    profile.cleaning_rounds += thread_profile.cleaning_rounds;
    profile.deadend_changed_rounds +=
        thread_profile.deadend_changed_rounds;
    profile.bubble_changed_rounds +=
        thread_profile.bubble_changed_rounds;
    profile.coverage_changed_rounds +=
        thread_profile.coverage_changed_rounds;
    profile.branch_transition_hits +=
        thread_profile.branch_transition_hits;
    profile.branch_transition_misses +=
        thread_profile.branch_transition_misses;
    profile.incremental_occurrences +=
        thread_profile.incremental_occurrences;
    profile.incremental_path_replays +=
        thread_profile.incremental_path_replays;
    profile.incremental_primary_inserts +=
        thread_profile.incremental_primary_inserts;
    profile.incremental_primary_hits +=
        thread_profile.incremental_primary_hits;
    profile.incremental_fork_lookups +=
        thread_profile.incremental_fork_lookups;
  }
  xinfo("Inner-k CPU seconds ({} rounds, {} raw k-mers, {} raw vertices): "
        "reads {.6}, coverage {.6}, "
        "anchors {.6}, hash-assemble {.6} (walk {.6}, adjacency {.6}), "
        "contig-init {.6}, dead-end {.6}, "
        "bubble {.6}, coverage-clean {.6}, contig-assemble {.6}\n",
        profile.k_rounds, profile.raw_kmers, profile.raw_vertices,
        profile.read_graph_build, profile.coverage,
        profile.anchor_insert, profile.hash_assemble, profile.hash_walk,
        profile.hash_adjacency,
        profile.contig_initialize, profile.dead_end, profile.bubble,
        profile.coverage_clean, profile.contig_assemble);
  if (profile.multik_build != 0.0 || profile.multik_finish != 0.0 ||
      profile.multik_advance != 0.0 || profile.multik_aggregate != 0.0 ||
      profile.multik_replay != 0.0) {
    xinfo("Cross-k CPU seconds: build {.6}, initial frontier {.6}, refine "
          "{.6}, aggregate {.6}, replay {.6}\n",
          profile.multik_build, profile.multik_finish,
          profile.multik_advance, profile.multik_aggregate,
          profile.multik_replay);
  }
  if (profile.incremental_occurrences != 0u) {
    xinfo("Incremental-k transitions: {} occurrences, {} path replays, "
          "{} primary inserts, {} primary hits, {} fork lookups\n",
          profile.incremental_occurrences,
          profile.incremental_path_replays,
          profile.incremental_primary_inserts,
          profile.incremental_primary_hits,
          profile.incremental_fork_lookups);
  }

  std::vector<uint64_t> k_round_survival(num_k_values, 0);
  std::vector<double> k_read_graph_seconds(num_k_values, 0.0);
  std::vector<uint64_t> k_branch_hits(num_k_values, 0);
  std::vector<uint64_t> k_branch_misses(num_k_values, 0);
  for (const auto &thread_counts : thread_k_round_survival) {
    for (size_t i = 0; i < num_k_values; ++i) {
      k_round_survival[i] += thread_counts[i];
    }
  }
  for (size_t thread = 0; thread < thread_k_read_graph_seconds.size();
       ++thread) {
    for (size_t i = 0; i < num_k_values; ++i) {
      k_read_graph_seconds[i] += thread_k_read_graph_seconds[thread][i];
      k_branch_hits[i] += thread_k_branch_hits[thread][i];
      k_branch_misses[i] += thread_k_branch_misses[thread][i];
    }
  }
  std::ostringstream survival_log;
  survival_log << "Inner-k survivors:";
  for (size_t i = 0; i < num_k_values; ++i) {
    survival_log << " k" << (opt.kmin + i * opt.step) << '='
                 << k_round_survival[i];
  }
  const std::string survival_text = survival_log.str();
  xinfo("{}\n", survival_text.c_str());
  std::ostringstream branch_by_k_log;
  branch_by_k_log << "Branch cache by inner k:";
  for (size_t i = 0; i < num_k_values; ++i) {
    branch_by_k_log << " k" << (opt.kmin + i * opt.step) << '='
                    << k_read_graph_seconds[i] << 's' << '/'
                    << k_branch_hits[i] << '/' << k_branch_misses[i];
  }
  const std::string branch_by_k_text = branch_by_k_log.str();
  xinfo("{}\n", branch_by_k_text.c_str());
  xinfo("Single HashGraph unitig rounds: {}, long: {}, isolated-long: {}\n",
        profile.single_hash_unitigs, profile.long_single_hash_unitigs,
        profile.isolated_long_single_hash_unitigs);
  xinfo("Cleaning mutations over {} rounds: dead-end {}, bubble {}, "
        "coverage {}\n",
        profile.cleaning_rounds, profile.deadend_changed_rounds,
        profile.bubble_changed_rounds, profile.coverage_changed_rounds);
  xinfo("Branch-transition cache: {} hits, {} cold misses\n",
        profile.branch_transition_hits, profile.branch_transition_misses);
  auto percentile_double = [](std::vector<double> values, double p) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t index = static_cast<size_t>(
        std::floor(p * static_cast<double>(values.size() - 1u)));
    return values[index];
  };
  auto percentile_u64 = [](std::vector<uint64_t> values, double p) {
    if (values.empty()) return uint64_t{0};
    std::sort(values.begin(), values.end());
    const size_t index = static_cast<size_t>(
        std::floor(p * static_cast<double>(values.size() - 1u)));
    return values[index];
  };
  xinfo("Endpoint graph seconds P50 {.6}, P95 {.6}, P99 {.6}, P99.9 {.6}, "
        "max {.6}\n",
        percentile_double(task_elapsed_seconds, 0.50),
        percentile_double(task_elapsed_seconds, 0.95),
        percentile_double(task_elapsed_seconds, 0.99),
        percentile_double(task_elapsed_seconds, 0.999),
        percentile_double(task_elapsed_seconds, 1.0));
  xinfo("Endpoint potential k-mers P50 {}, P95 {}, P99 {}, P99.9 {}, max {}; "
        "actual max {}, unique max {}\n",
        percentile_u64(task_potential_kmers, 0.50),
        percentile_u64(task_potential_kmers, 0.95),
        percentile_u64(task_potential_kmers, 0.99),
        percentile_u64(task_potential_kmers, 0.999),
        percentile_u64(task_potential_kmers, 1.0),
        percentile_u64(task_actual_kmers, 1.0),
        percentile_u64(task_unique_kmers, 1.0));

  if (task_phase_elapsed > 0) {
    const double tail_begin = task_phase_elapsed * 0.9;
    double occupied_tail_seconds = 0;
    for (size_t i = 0; i < tasks.size(); ++i) {
      const double overlap_begin = std::max(task_begin_seconds[i], tail_begin);
      const double overlap_end =
          std::min(task_end_seconds[i], task_phase_elapsed);
      if (overlap_end > overlap_begin) {
        occupied_tail_seconds += overlap_end - overlap_begin;
      }
    }
    const double average_active_endpoints =
        occupied_tail_seconds / (task_phase_elapsed - tail_begin);
    xinfo("Endpoint phase {.6}s; final-tenth average active endpoint tasks "
          "{.2}\n",
          task_phase_elapsed, average_active_endpoints);
  }
}

}  // namespace

void RunLocalAssembly(const LocalAsmOption &opt) {
  SimpleTimer timer;
  timer.reset();
  timer.start();
  HashMapper mapper;
  mapper.LoadAndBuild(opt.contig_file, opt.min_contig_len, opt.seed_kmer,
                      opt.sparsity);
  mapper.SetMappingThreshold(opt.min_mapping_len, opt.similarity);
  timer.stop();
  xinfo("Hash mapper construction time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  SequenceLibCollection lib_collection;
  SeqPackage read_pkg;
  lib_collection.SetPath(opt.lib_file_prefix);
  const SequenceLibCollection::SizeInfo read_size =
      lib_collection.GetSizeInfo();
  MappedReadFile mapped_reads;
  MappedLocalCandidateFile mapped_candidates;
  const bool use_mapped_stream =
      mapped_reads.Open(opt.lib_file_prefix, read_size);
  if (use_mapped_stream) {
    lib_collection.ReadMetadata();
    xinfo("Using {} indexed binary-read chunks; full read package is not "
          "materialized\n",
          mapped_reads.index().chunks.size());
  } else {
    lib_collection.Read(&read_pkg);
  }
  if (!opt.candidate_file.empty()) {
    if (!use_mapped_stream ||
        !mapped_candidates.Open(
            opt.candidate_file, mapped_reads.total_words(),
            static_cast<uint64_t>(read_size.num_reads))) {
      xfatal("Invalid local candidate file or packed-read source mismatch: "
             "{s}",
             opt.candidate_file.c_str());
    }
    xinfo("Using {} exact persistent-index candidate reads\n",
          mapped_candidates.size());
  }
  timer.stop();
  xinfo("Read lib time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  const auto insert_sizes =
      use_mapped_stream
          ? EstimateInsertSizeMapped(mapper, lib_collection, mapped_reads)
          : EstimateInsertSize(mapper, lib_collection);
  timer.stop();
  xinfo("Insert size estimation time elapsed: {}\n", timer.elapsed());

  const int32_t max_local_range =
      GetMaxLocalRange(lib_collection, insert_sizes);
  timer.reset();
  timer.start();
  mapper.BuildEndpointSeedFilter(max_local_range);
  LocalSeedPositions seed_positions;
  const LocalSeedPositions *seed_position_pointer = nullptr;
  if (use_mapped_stream &&
      std::getenv("MEGAHIT_EXPERIMENTAL_LOCAL_MINIMIZER_GATE") != nullptr &&
      mapper.BuildEndpointMinimizerGate()) {
    const char *configured = std::getenv("MEGAHIT_LOCAL_SEED_POSITIONS");
    const std::string position_path = configured != nullptr
        ? configured : opt.lib_file_prefix + ".bin.local_seed_positions";
    const double position_begin = omp_get_wtime();
    try {
      if (seed_positions.OpenOrBuild(position_path, opt.lib_file_prefix + ".bin",
                                    mapped_reads.words(), mapped_reads.index())) {
        seed_position_pointer = &seed_positions;
        xinfo("Local seed positions: {} bytes, {}, {.4} seconds\n",
              seed_positions.bytes(), seed_positions.built() ? "built" : "reused",
              omp_get_wtime() - position_begin);
      }
    } catch (const std::runtime_error &error) {
      xwarn("Local seed position cache unavailable: {s}; using the exact "
            "streaming endpoint gate\n", error.what());
    }
    if (seed_position_pointer == nullptr) mapper.ReleaseEndpointMinimizerGate();
  }
  timer.stop();
  xinfo("Endpoint seed filter construction time elapsed: {}\n",
        timer.elapsed());

  timer.reset();
  timer.start();

  MappingResultCollector collector(mapper.refseq().seq_count());
  // File word offsets increase strictly with read IDs. Substituting them in
  // the least-significant tie field preserves the complete historical order,
  // including the first-three-reads-at-a-position selection. Endpoint assembly
  // can then read the original record directly without ranking, gathering and
  // copying the selected library in every outer-k round.
  const uint64_t local_memory_budget =
      GetRuntimeResourcePolicy().memory_budget_per_job;
  const bool mapped_assembly = use_mapped_stream &&
      mapped_reads.total_words() < (uint64_t{1} << 44u) &&
      (local_memory_budget == 0u ||
       mapped_reads.total_words() <= local_memory_budget / 8u) &&
      std::getenv("MEGAHIT_EXPERIMENTAL_LOCAL_MAPPED_ASSEMBLY") != nullptr;
  if (use_mapped_stream) {
    MapToContigsMapped(mapper, lib_collection, insert_sizes, mapped_reads,
                       opt.candidate_file.empty() ? nullptr
                                                  : &mapped_candidates,
                       seed_position_pointer,
                       mapped_assembly,
                       &collector);
  } else {
    MapToContigs(mapper, lib_collection, insert_sizes, &collector);
  }
  timer.stop();
  xinfo("Mapping time elapsed: {}\n", timer.elapsed());
  mapper.ReportPhaseCertificateStats();

  timer.reset();
  timer.start();
  collector.Finalize();
  timer.stop();
  xinfo("Mapping result collation time elapsed: {}, retained {} mappings\n",
        timer.elapsed(), collector.size());

  mapper.ReleaseIndex();
  seed_positions.Release();
  mapped_candidates.Close();

  if (mapped_assembly) {
    const double prepare_begin = omp_get_wtime();
    const unsigned maximum = collector.MaxReferencedMappedReadLength(
        mapped_reads.words(), mapped_reads.total_words());
    xinfo("Mapped endpoint reads: {} exact records, maximum length {}, "
          "preparation {.4} seconds; selected reads are not repacked\n",
          collector.size(), maximum, omp_get_wtime() - prepare_begin);
    timer.reset();
    timer.start();
    AssembleAndOutput(mapper, read_pkg, collector, opt.output_file,
                        max_local_range, opt, mapped_reads.words(), maximum);
    timer.stop();
    xinfo("Local assembly time elapsed: {}\n", timer.elapsed());
    return;
  }

  timer.reset();
  timer.start();
  const size_t original_read_count = static_cast<size_t>(read_size.num_reads);
  const size_t original_read_bytes =
      use_mapped_stream
          ? DivCeiling(static_cast<size_t>(read_size.num_bases),
                       SeqPackage::kBasesPerWord) *
                sizeof(SeqPackage::TWord)
          : read_pkg.size_in_byte();
  std::vector<uint64_t> selected_read_ids =
      collector.CompactReadIds(original_read_count);
  SeqPackage selected_read_pkg;
  if (use_mapped_stream) {
    if (!selected_read_pkg.AssignSelectedMappedBinaryRecords(
            mapped_reads.words(), mapped_reads.total_words(),
            original_read_count, selected_read_ids,
            mapped_reads.index().chunks, omp_get_max_threads())) {
      xfatal("Indexed binary-read stream changed during local mapping\n");
    }
    mapped_reads.Close();
  } else {
    selected_read_pkg.ReserveSequences(selected_read_ids.size());
    if (read_pkg.max_length() != 0 &&
        selected_read_ids.size() <=
            std::numeric_limits<size_t>::max() / read_pkg.max_length()) {
      selected_read_pkg.ReserveBases(selected_read_ids.size() *
                                     read_pkg.max_length());
    }
    for (const uint64_t read_id : selected_read_ids) {
      selected_read_pkg.AppendSequenceView(read_pkg.GetSeqView(read_id));
    }
  }
  const unsigned selected_gap_bits =
      selected_read_pkg.CompactReadOnlyLengthIndex();
  const size_t selected_read_bytes = selected_read_pkg.size_in_byte();
  read_pkg.ReleaseStorage();
  read_pkg = std::move(selected_read_pkg);
  std::vector<uint64_t>().swap(selected_read_ids);
  timer.stop();
  xinfo("Compacted referenced reads: {} / {} reads, {} -> {} bytes "
        "(length gap index {} bits), elapsed {}\n",
        read_pkg.seq_count(), original_read_count, original_read_bytes,
        selected_read_bytes, selected_gap_bits, timer.elapsed());

  timer.reset();
  timer.start();
  AssembleAndOutput(mapper, read_pkg, collector, opt.output_file,
                    max_local_range, opt);
  timer.stop();
  xinfo("Local assembly time elapsed: {}\n", timer.elapsed());
}
