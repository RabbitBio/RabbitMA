#ifndef MEGAHIT_SEQUENCE_IO_LOCAL_SEED_POSITIONS_H_
#define MEGAHIT_SEQUENCE_IO_LOCAL_SEED_POSITIONS_H_

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>
#include "sequence/io/read_anchor_positions.h"
#include "sequence/io/read_chunk_index.h"
#include "utils/utils.h"

// Immutable positions of the rightmost minimum forward 19-mer in every
// 31-base read window. The bitmap is independent of contigs and outer k.
// Endpoint lookup indexes BOTH orientations of each original 31-mer, so
// matching windows choose the same anchor and offset, including tied minima.
class LocalSeedPositions {
 public:
  LocalSeedPositions() = default;
  ~LocalSeedPositions() { Close(); }
  LocalSeedPositions(const LocalSeedPositions &) = delete;
  LocalSeedPositions &operator=(const LocalSeedPositions &) = delete;

  bool OpenOrBuild(const std::string &path, const std::string &source_path,
                   const uint32_t *source, const PackedReadChunkIndex &index) {
    Close();
    if (index.max_read_len > 512u || index.max_read_len < 31u) return false;
    struct stat source_status {};
    if (stat(source_path.c_str(), &source_status) != 0) return false;
    Header expected;
    expected.source_bytes = source_status.st_size;
    expected.source_inode = source_status.st_ino;
    expected.source_device = source_status.st_dev;
    expected.mtime_seconds = PackedReadMtimeSeconds(source_status);
    expected.mtime_nanoseconds = PackedReadMtimeNanoseconds(source_status);
    expected.read_count = index.num_reads;
    expected.max_read_length = index.max_read_len;
    expected.stride = (index.max_read_len - 19u + 1u + 7u) / 8u;
    if (index.num_reads >
        (std::numeric_limits<size_t>::max() - sizeof(Header)) /
            expected.stride) return false;
    const size_t bytes = sizeof(Header) + index.num_reads * expected.stride;
    if (OpenExisting(path, expected, bytes)) return true;

    const std::string temporary = path + ".tmp." + std::to_string(getpid());
    const int fd = open(temporary.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    if (fd < 0) throw std::runtime_error("cannot create local seed positions");
    const int allocated = posix_fallocate(fd, 0, bytes);
    if (allocated != 0) {
      close(fd);
      unlink(temporary.c_str());
      throw std::runtime_error("cannot allocate local seed positions: " +
                               std::string(std::strerror(allocated)));
    }
    void *mapping = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED,
                          fd, 0);
    if (mapping == MAP_FAILED) {
      close(fd);
      unlink(temporary.c_str());
      throw std::runtime_error("cannot map local seed positions");
    }
    auto *payload = static_cast<uint8_t *>(mapping) + sizeof(Header);
    int malformed = 0;
#pragma omp parallel reduction(| : malformed)
    {
      std::vector<uint32_t> queue_pos;
      std::vector<uint64_t> queue_hash, queue_key;
#pragma omp for schedule(dynamic, 1)
      for (int64_t chunk_id = 0;
           chunk_id < static_cast<int64_t>(index.chunks.size()); ++chunk_id) {
        const auto &chunk = index.chunks[chunk_id];
        const uint32_t *cursor = source + chunk.word_begin;
        for (uint64_t read = chunk.read_begin; read < chunk.read_end; ++read) {
          if (cursor >= source + chunk.word_end) { malformed = 1; break; }
          const unsigned length = *cursor++;
          const size_t words = (static_cast<size_t>(length) + 15u) / 16u;
          if (length > index.max_read_len ||
              words > static_cast<size_t>(source + chunk.word_end - cursor)) {
            malformed = 1;
            break;
          }
          uint8_t *bits = payload + read * expected.stride;
          std::memset(bits, 0, expected.stride);
          if (length >= 50u) {
            ForEachPackedReadAnchor(cursor, length, 19u, 31u, &queue_pos,
                                     &queue_hash, &queue_key,
                [&](uint64_t, uint32_t position) {
                  bits[position >> 3u] |= uint8_t{1} << (position & 7u);
                });
          }
          cursor += words;
        }
        if (cursor != source + chunk.word_end) malformed = 1;
        DiscardMemoryPages(const_cast<uint32_t *>(source + chunk.word_begin),
                            (chunk.word_end - chunk.word_begin) * 4u);
      }
    }
    // Publish a valid header only after every disjoint read range is complete.
    if (!malformed) std::memcpy(mapping, &expected, sizeof(expected));
    munmap(mapping, bytes);
    close(fd);
    if (malformed || rename(temporary.c_str(), path.c_str()) != 0) {
      unlink(temporary.c_str());
      throw std::runtime_error("cannot publish exact local seed positions");
    }
    if (!OpenExisting(path, expected, bytes))
      throw std::runtime_error("cannot reopen local seed positions");
    built_ = true;
    return true;
  }

  const uint8_t *Get(uint64_t read_id) const {
    return payload_ + read_id * stride_;
  }
  size_t bytes() const { return bytes_; }
  bool built() const { return built_; }
  void Release() { Close(); }
  void Drop(uint64_t first, uint64_t count) const {
    DiscardMemoryPages(const_cast<uint8_t *>(Get(first)), count * stride_);
  }

 private:
  struct Header {
    uint64_t magic{UINT64_C(0x4c5345504f533031)};
    uint64_t source_bytes{0}, source_inode{0}, source_device{0};
    int64_t mtime_seconds{0}, mtime_nanoseconds{0};
    uint64_t read_count{0};
    uint32_t max_read_length{0}, stride{0};
  };
  static_assert(sizeof(Header) == 64u, "local position header has no padding");

  bool OpenExisting(const std::string &path, const Header &expected,
                     size_t bytes) {
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;
    struct stat status {};
    Header observed;
    const bool valid = fstat(fd, &status) == 0 &&
        static_cast<uint64_t>(status.st_size) == bytes &&
        pread(fd, &observed, sizeof(observed), 0) == sizeof(observed) &&
        std::memcmp(&observed, &expected, sizeof(expected)) == 0;
    if (!valid) { close(fd); return false; }
    void *mapping = mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (mapping == MAP_FAILED) return false;
    mapping_ = static_cast<const uint8_t *>(mapping);
    payload_ = mapping_ + sizeof(Header);
    bytes_ = bytes;
    stride_ = expected.stride;
    return true;
  }
  void Close() {
    if (mapping_ != nullptr) munmap(const_cast<uint8_t *>(mapping_), bytes_);
    mapping_ = payload_ = nullptr;
    bytes_ = stride_ = 0;
    built_ = false;
  }
  const uint8_t *mapping_{nullptr}, *payload_{nullptr};
  size_t bytes_{0}, stride_{0};
  bool built_{false};
};

#endif
