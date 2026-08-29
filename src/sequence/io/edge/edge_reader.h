//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_EDGE_READER_H
#define MEGAHIT_EDGE_READER_H

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <fcntl.h>
#include <limits>
#include <new>
#include <sys/uio.h>
#include <unistd.h>

#include "definitions.h"
#include "edge_io_meta.h"
#include "sequence/io/base_reader.h"

class EdgeReader : public BaseSequenceReader {
 public:
  EdgeReader(const std::string &file_prefix) {
    file_prefix_ = file_prefix;
    std::ifstream is(file_prefix + ".edges.info");
    metadata_.Deserialize(is);
    InitFiles();
  }

  EdgeReader *SetMultiplicityVec(std::vector<mul_t> *mul) {
    mul_ = mul;
    return this;
  }

  int64_t ReadUnsorted(SeqPackage *pkg, std::vector<mul_t> *mul,
                       int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = NextUnsortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, metadata_.kmer_size + 1);
      if (mul) {
        mul->push_back(next_edge[metadata_.words_per_edge - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  int64_t ReadSorted(SeqPackage *pkg, std::vector<mul_t> *mul,
                     int64_t max_num) {
    for (int64_t i = 0; i < max_num; ++i) {
      uint32_t *next_edge = NextSortedEdge();
      if (next_edge == nullptr) {
        return i;
      }
      pkg->AppendCompactSequence(next_edge, metadata_.kmer_size + 1);
      if (mul) {
        mul->push_back(next_edge[metadata_.words_per_edge - 1] & kMaxMul);
      }
    }
    return max_num;
  }

  /**
   * Load a complete sorted fixed-length edge set with large sequential
   * scatter reads, then pack it into SequencePackage in parallel.  Returns
   * -1 when the fast path is inapplicable or its temporary buffer does not fit
   * the caller's memory budget; in that case reader state and outputs are left
   * untouched and ReadSorted remains a safe fallback.
   */
  int64_t ReadSortedBulk(SeqPackage *pkg, std::vector<mul_t> *mul,
                         int num_threads, int64_t host_mem) {
    if (!metadata_.is_sorted || std::getenv("MEGAHIT_DISABLE_BULK_EDGE_LOAD") ||
        pkg == nullptr || pkg->seq_count() != 0 || pkg->base_count() != 0 ||
        metadata_.num_edges < 0 || host_mem <= 0) {
      return -1;
    }

    const uint64_t num_edges = static_cast<uint64_t>(metadata_.num_edges);
    const uint64_t words_per_edge = metadata_.words_per_edge;
    const uint64_t edge_length =
        static_cast<uint64_t>(metadata_.kmer_size) + 1u;
    const uint64_t expected_words_per_edge = DivCeiling(
        edge_length * kBitsPerEdgeChar + kBitsPerMul, kBitsPerEdgeWord);
    if (words_per_edge == 0 || words_per_edge != expected_words_per_edge ||
        num_edges > std::numeric_limits<size_t>::max() / words_per_edge ||
        num_edges > std::numeric_limits<size_t>::max() / edge_length) {
      return -1;
    }
    const uint64_t raw_words = num_edges * words_per_edge;
    if (raw_words > std::numeric_limits<uint64_t>::max() / sizeof(uint32_t)) {
      return -1;
    }
    const uint64_t raw_bytes = raw_words * sizeof(uint32_t);
    if (num_edges == 0) {
      return 0;
    }
    const uint64_t mul_capacity_bytes =
        mul == nullptr ? 0 : mul->capacity() * sizeof(mul_t);
    if (pkg->size_in_byte() >
        std::numeric_limits<uint64_t>::max() - mul_capacity_bytes) {
      return -1;
    }
    const uint64_t resident_capacity =
        pkg->size_in_byte() + mul_capacity_bytes;
    const uint64_t memory_budget = static_cast<uint64_t>(host_mem);
    const uint64_t headroom =
        std::max<uint64_t>(uint64_t{256} << 20u, memory_budget / 100u);
    if (resident_capacity > memory_budget || headroom > memory_budget - resident_capacity ||
        raw_bytes > memory_budget - resident_capacity - headroom) {
      xinfo("Bulk edge load skipped: {} temporary bytes exceed budget\n",
            raw_bytes);
      return -1;
    }

    std::unique_ptr<uint32_t[]> raw(
        new (std::nothrow) uint32_t[static_cast<size_t>(raw_words)]);
    if (!raw) {
      xwarn("Bulk edge allocation failed; using buffered edge reader\n");
      return -1;
    }

    struct Segment {
      uint64_t source_edge;
      uint64_t destination_edge;
      uint64_t num_edges;
    };
    const uint64_t record_bytes = words_per_edge * sizeof(uint32_t);
    const uint64_t max_batch_bytes = uint64_t{64} << 20u;
    const uint64_t max_piece_edges =
        std::max<uint64_t>(1, max_batch_bytes / record_bytes);
    const uint64_t max_file_offset =
        static_cast<uint64_t>(std::numeric_limits<off_t>::max());
    std::vector<std::vector<Segment>> file_segments(metadata_.num_files);
    uint64_t destination_edge = 0;
    for (const auto &bucket : metadata_.buckets) {
      if (bucket.total_number < 0 || bucket.file_offset < 0 ||
          bucket.file_id >= static_cast<int>(metadata_.num_files)) {
        return -1;
      }
      if (bucket.total_number == 0) {
        continue;
      }
      const uint64_t bucket_edges =
          static_cast<uint64_t>(bucket.total_number);
      if (bucket.file_id < 0 || destination_edge > num_edges ||
          bucket_edges > num_edges - destination_edge) {
        return -1;
      }
      const uint64_t source_edge = bucket.file_offset;
      if (source_edge > std::numeric_limits<uint64_t>::max() - bucket_edges ||
          (source_edge + bucket_edges) > max_file_offset / record_bytes) {
        return -1;
      }
      for (uint64_t done = 0; done < bucket_edges;) {
        const uint64_t piece_edges =
            std::min(max_piece_edges, bucket_edges - done);
        file_segments[bucket.file_id].push_back(
            {source_edge + done, destination_edge + done, piece_edges});
        done += piece_edges;
      }
      destination_edge += bucket_edges;
    }
    if (destination_edge != num_edges) {
      return -1;
    }
    for (auto &segments : file_segments) {
      std::sort(segments.begin(), segments.end(),
                [](const Segment &lhs, const Segment &rhs) {
                  return lhs.source_edge < rhs.source_edge;
                });
    }

    std::vector<int> file_descriptors(metadata_.num_files, -1);
    bool files_opened = true;
    for (unsigned file_id = 0; file_id < metadata_.num_files; ++file_id) {
      const std::string path =
          file_prefix_ + ".edges." + std::to_string(file_id);
      file_descriptors[file_id] = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
      if (file_descriptors[file_id] < 0) {
        files_opened = false;
        break;
      }
    }
    if (!files_opened) {
      for (int fd : file_descriptors) {
        if (fd >= 0) {
          ::close(fd);
        }
      }
      xwarn("Bulk edge file open failed; using buffered edge reader\n");
      return -1;
    }

    const auto io_start = std::chrono::steady_clock::now();
    std::atomic<bool> io_ok(true);
    num_threads = std::max(1, num_threads);
    const size_t max_iov = static_cast<size_t>(
        std::max<long>(1, std::min<long>(1024, ::sysconf(_SC_IOV_MAX))));
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (int file_id = 0; file_id < static_cast<int>(metadata_.num_files);
         ++file_id) {
      const int fd = file_descriptors[file_id];
      (void)::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
      const auto &segments = file_segments[file_id];
      std::vector<struct iovec> batch;
      batch.reserve(max_iov);

      for (size_t first = 0; first < segments.size() && io_ok.load(); ) {
        batch.clear();
        const uint64_t first_source_byte =
            segments[first].source_edge * record_bytes;
        uint64_t expected_source_byte = first_source_byte;
        uint64_t batch_bytes = 0;
        size_t next = first;
        while (next < segments.size() && batch.size() < max_iov) {
          const auto &segment = segments[next];
          const uint64_t source_byte = segment.source_edge * record_bytes;
          const uint64_t segment_bytes = segment.num_edges * record_bytes;
          if (source_byte != expected_source_byte ||
              (!batch.empty() && batch_bytes + segment_bytes > max_batch_bytes)) {
            break;
          }
          struct iovec part;
          part.iov_base = raw.get() + segment.destination_edge * words_per_edge;
          part.iov_len = static_cast<size_t>(segment_bytes);
          batch.push_back(part);
          batch_bytes += segment_bytes;
          expected_source_byte += segment_bytes;
          ++next;
        }

        if (!ReadvFully(fd, batch.data(), batch.size(), first_source_byte)) {
          io_ok.store(false);
          break;
        }
        first = next;
      }
    }
    for (int fd : file_descriptors) {
      ::close(fd);
    }
    if (!io_ok.load()) {
      xwarn("Bulk edge read failed; using buffered edge reader\n");
      return -1;
    }
    const auto io_end = std::chrono::steady_clock::now();

    const auto pack_start = io_end;
    pkg->AssignFixedLengthCompactSequences(
        raw.get(), static_cast<size_t>(num_edges),
        static_cast<unsigned>(edge_length),
        static_cast<unsigned>(words_per_edge), num_threads);
    if (mul != nullptr) {
      const size_t old_size = mul->size();
      mul->resize(old_size + static_cast<size_t>(num_edges));
#pragma omp parallel for schedule(static) num_threads(num_threads)
      for (int64_t edge_id = 0; edge_id < metadata_.num_edges; ++edge_id) {
        (*mul)[old_size + edge_id] =
            raw[static_cast<uint64_t>(edge_id) * words_per_edge +
                words_per_edge - 1] &
            kMaxMul;
      }
    }
    const auto pack_end = std::chrono::steady_clock::now();
    const double io_seconds =
        std::chrono::duration<double>(io_end - io_start).count();
    const double pack_seconds =
        std::chrono::duration<double>(pack_end - pack_start).count();
    xinfo("Bulk edge load: {} records, {} raw bytes, IO {.4}, pack {.4}\n",
          num_edges, raw_bytes, io_seconds, pack_seconds);
    return metadata_.num_edges;
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse = false) override {
    if (metadata_.is_sorted) {
      return ReadSorted(pkg, mul_, max_num);
    } else {
      return ReadUnsorted(pkg, mul_, max_num);
    }
  }

 private:
  static bool ReadvFully(int fd, struct iovec *iov, size_t iov_count,
                         uint64_t file_offset) {
    size_t first = 0;
    while (first < iov_count) {
      const ssize_t bytes =
          ::preadv(fd, iov + first, static_cast<int>(iov_count - first),
                   static_cast<off_t>(file_offset));
      if (bytes < 0) {
        if (errno == EINTR) {
          continue;
        }
        return false;
      }
      if (bytes == 0) {
        return false;
      }

      size_t consumed = static_cast<size_t>(bytes);
      file_offset += consumed;
      while (first < iov_count && consumed >= iov[first].iov_len) {
        consumed -= iov[first].iov_len;
        ++first;
      }
      if (consumed > 0) {
        iov[first].iov_base =
            static_cast<char *>(iov[first].iov_base) + consumed;
        iov[first].iov_len -= consumed;
      }
    }
    return true;
  }

  std::vector<mul_t> *mul_{nullptr};

  std::string file_prefix_;
  std::vector<std::unique_ptr<std::ifstream>> in_streams_;
  BufferedReader cur_reader_;
  std::vector<uint32_t> buffer_;

  int cur_bucket_{};
  int64_t cur_cnt_{};
  int64_t cur_vol_{};
  bool is_opened_{false};

  EdgeIoMetadata metadata_;

 private:
  void InitFiles() {
    assert(!is_opened_);
    buffer_.resize(metadata_.words_per_edge);

    for (unsigned i = 0; i < metadata_.num_files; ++i) {
      in_streams_.emplace_back(
          new std::ifstream(file_prefix_ + ".edges." + std::to_string(i),
                            std::ifstream::binary | std::ifstream::in));
    }

    cur_cnt_ = 0;
    cur_vol_ = 0;
    cur_bucket_ = -1;

    if (!metadata_.is_sorted) {
      cur_reader_.reset(in_streams_[0].get());
    }

    is_opened_ = true;
  }

 public:
  const EdgeIoMetadata &GetMetadata() const { return metadata_; }

 private:
  uint32_t *NextSortedEdge() {
    if (cur_bucket_ >= static_cast<int>(metadata_.buckets.size())) {
      return nullptr;
    }

    while (cur_cnt_ >= cur_vol_) {
      ++cur_bucket_;

      while (cur_bucket_ < static_cast<int>(metadata_.buckets.size()) &&
             metadata_.buckets[cur_bucket_].file_id < 0) {
        ++cur_bucket_;
      }

      if (cur_bucket_ >= static_cast<int>(metadata_.buckets.size())) {
        return nullptr;
      }

      const auto &bucket = metadata_.buckets[cur_bucket_];
      cur_cnt_ = 0;
      cur_vol_ = bucket.total_number;
      auto is = in_streams_[bucket.file_id].get();
      is->clear();
      is->seekg(bucket.file_offset * sizeof(uint32_t) *
                metadata_.words_per_edge);
      cur_reader_.reset(is, bucket.total_number * sizeof(uint32_t) *
                                metadata_.words_per_edge);
    }

    ++cur_cnt_;
    auto n_read = cur_reader_.read(buffer_.data(), metadata_.words_per_edge);
    assert(n_read == metadata_.words_per_edge * sizeof(uint32_t));
    (void)n_read;
    return buffer_.data();
  }

  uint32_t *NextUnsortedEdge() {
    if (cur_cnt_ >= metadata_.num_edges) {
      return nullptr;
    }

    ++cur_cnt_;
    auto n_read = cur_reader_.read(buffer_.data(), metadata_.words_per_edge);
    assert(n_read == metadata_.words_per_edge * sizeof(uint32_t));
    (void)n_read;
    return buffer_.data();
  }
};

#endif  // MEGAHIT_EDGE_READER_H
