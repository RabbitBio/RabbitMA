#ifndef MEGAHIT_SEQUENCE_IO_READ_CHUNK_INDEX_H_
#define MEGAHIT_SEQUENCE_IO_READ_CHUNK_INDEX_H_

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// A tiny persistent index over the variable-length [length][packed bases]
// records in reads.lib.bin.  The entries contain only exact record boundaries
// and additive metadata; they do not depend on k, thread topology, or an
// assembly policy.  Buildlib can therefore publish them for every later
// count/index pass without changing biological semantics.
struct PackedReadChunk {
  PackedReadChunk() = default;
  PackedReadChunk(uint64_t word_begin_arg, uint64_t word_end_arg,
                  uint64_t read_begin_arg, uint64_t read_end_arg,
                  uint64_t num_bases_arg = 0,
                  unsigned max_read_len_arg = 0)
      : word_begin(word_begin_arg),
        word_end(word_end_arg),
        read_begin(read_begin_arg),
        read_end(read_end_arg),
        num_bases(num_bases_arg),
        max_read_len(max_read_len_arg) {}

  uint64_t word_begin{0};
  uint64_t word_end{0};
  uint64_t read_begin{0};
  uint64_t read_end{0};
  uint64_t num_bases{0};
  unsigned max_read_len{0};
};

struct PackedReadChunkIndex {
  uint64_t num_reads{0};
  uint64_t num_bases{0};
  unsigned max_read_len{0};
  std::vector<PackedReadChunk> chunks;
};

inline std::string PackedReadChunkIndexPath(const std::string &read_path) {
  return read_path + ".ridx.chunks";
}

inline int64_t PackedReadMtimeSeconds(const struct stat &status) {
#if defined(__APPLE__)
  return static_cast<int64_t>(status.st_mtimespec.tv_sec);
#else
  return static_cast<int64_t>(status.st_mtim.tv_sec);
#endif
}

inline int64_t PackedReadMtimeNanoseconds(const struct stat &status) {
#if defined(__APPLE__)
  return static_cast<int64_t>(status.st_mtimespec.tv_nsec);
#else
  return static_cast<int64_t>(status.st_mtim.tv_nsec);
#endif
}

inline bool LoadPackedReadChunkIndex(const std::string &read_path,
                                     PackedReadChunkIndex *index) {
  struct stat status;
  if (stat(read_path.c_str(), &status) != 0 || status.st_size <= 0 ||
      status.st_size % static_cast<off_t>(sizeof(uint32_t)) != 0) {
    return false;
  }
  std::ifstream input(PackedReadChunkIndexPath(read_path).c_str());
  std::string magic;
  uint64_t stored_bytes = 0;
  int64_t stored_mtime_seconds = 0;
  int64_t stored_mtime_nanoseconds = 0;
  uint64_t chunk_count = 0;
  PackedReadChunkIndex loaded;
  if (!(input >> magic >> stored_bytes >> stored_mtime_seconds >>
        stored_mtime_nanoseconds >> loaded.num_reads >> loaded.num_bases >>
        loaded.max_read_len >> chunk_count) ||
      magic != "MEGAHIT_READ_CHUNKS_V1" ||
      stored_bytes != static_cast<uint64_t>(status.st_size) ||
      stored_mtime_seconds != PackedReadMtimeSeconds(status) ||
      stored_mtime_nanoseconds != PackedReadMtimeNanoseconds(status) ||
      chunk_count == 0 ||
      chunk_count > static_cast<uint64_t>(
                        std::numeric_limits<size_t>::max())) {
    return false;
  }
  loaded.chunks.resize(static_cast<size_t>(chunk_count));
  uint64_t observed_bases = 0;
  unsigned observed_max = 0;
  for (size_t i = 0; i < loaded.chunks.size(); ++i) {
    PackedReadChunk &chunk = loaded.chunks[i];
    if (!(input >> chunk.word_begin >> chunk.word_end >> chunk.read_begin >>
          chunk.read_end >> chunk.num_bases >> chunk.max_read_len) ||
        chunk.word_begin >= chunk.word_end ||
        chunk.read_begin >= chunk.read_end ||
        (i == 0 && (chunk.word_begin != 0 || chunk.read_begin != 0)) ||
        (i != 0 &&
         (chunk.word_begin != loaded.chunks[i - 1u].word_end ||
          chunk.read_begin != loaded.chunks[i - 1u].read_end))) {
      return false;
    }
    if (chunk.num_bases >
        std::numeric_limits<uint64_t>::max() - observed_bases) {
      return false;
    }
    observed_bases += chunk.num_bases;
    observed_max = std::max(observed_max, chunk.max_read_len);
  }
  const PackedReadChunk &last = loaded.chunks.back();
  if (last.word_end !=
          static_cast<uint64_t>(status.st_size) / sizeof(uint32_t) ||
      last.read_end != loaded.num_reads ||
      observed_bases != loaded.num_bases ||
      observed_max != loaded.max_read_len) {
    return false;
  }
  *index = std::move(loaded);
  return true;
}

inline bool PublishPackedReadChunkIndex(
    const std::string &read_path, uint64_t num_reads, uint64_t num_bases,
    unsigned max_read_len, const std::vector<PackedReadChunk> &chunks) {
  if (chunks.empty()) return false;
  struct stat status;
  if (stat(read_path.c_str(), &status) != 0 || status.st_size <= 0) {
    return false;
  }
  const std::string path = PackedReadChunkIndexPath(read_path);
  const std::string temporary_path =
      path + ".tmp." +
      std::to_string(static_cast<unsigned long long>(getpid()));
  std::ofstream output(temporary_path.c_str(),
                       std::ofstream::out | std::ofstream::trunc);
  if (!output) return false;
  output << "MEGAHIT_READ_CHUNKS_V1 "
         << static_cast<uint64_t>(status.st_size) << ' '
         << PackedReadMtimeSeconds(status) << ' '
         << PackedReadMtimeNanoseconds(status) << ' ' << num_reads << ' '
         << num_bases << ' ' << max_read_len << ' ' << chunks.size() << '\n';
  for (const PackedReadChunk &chunk : chunks) {
    output << chunk.word_begin << ' ' << chunk.word_end << ' '
           << chunk.read_begin << ' ' << chunk.read_end << ' '
           << chunk.num_bases << ' ' << chunk.max_read_len << '\n';
  }
  output.close();
  if (!output || rename(temporary_path.c_str(), path.c_str()) != 0) {
    std::remove(temporary_path.c_str());
    return false;
  }
  return true;
}

#endif  // MEGAHIT_SEQUENCE_IO_READ_CHUNK_INDEX_H_
