#ifndef MEGAHIT_UTILS_COALESCED_READS_H
#define MEGAHIT_UTILS_COALESCED_READS_H

#include <sys/uio.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <limits>
#include <vector>

struct PositionedRead {
  uint64_t offset;
  size_t bytes;
  uint8_t *destination;
};

struct CoalescedReadStats {
  uint64_t calls{0};
  uint64_t bytes{0};
};

// Scatter contiguous file ranges directly into their owning bucket buffers.
// No intermediate copy and no reads of gaps between requested extents. The
// destination of each request remains fixed even though requests are sorted.
inline bool ReadCoalescedExtents(int fd, std::vector<PositionedRead> *requests,
                                 CoalescedReadStats *stats) {
  std::sort(requests->begin(), requests->end(),
            [](const PositionedRead &a, const PositionedRead &b) {
              return a.offset < b.offset;
            });
  long limit = sysconf(_SC_IOV_MAX);
  if (limit <= 0) limit = 1024;
  std::vector<struct iovec> vectors;
  vectors.reserve(std::min<size_t>(requests->size(), limit));
  for (size_t begin = 0; begin < requests->size();) {
    vectors.clear();
    uint64_t offset = (*requests)[begin].offset;
    uint64_t end_offset = offset;
    size_t end = begin;
    while (end < requests->size() && vectors.size() < size_t(limit)) {
      const PositionedRead &request = (*requests)[end];
      if (request.offset != end_offset) break;
      if (request.bytes > uint64_t(std::numeric_limits<off_t>::max()) ||
          request.offset > uint64_t(std::numeric_limits<off_t>::max()) -
                               request.bytes) {
        return false;
      }
      if (request.bytes != 0) {
        vectors.push_back(iovec{request.destination, request.bytes});
        end_offset += request.bytes;
      }
      ++end;
    }
    size_t first = 0;
    while (first < vectors.size()) {
      ++stats->calls;
      const ssize_t amount = preadv(fd, vectors.data() + first,
                                    int(vectors.size() - first), off_t(offset));
      if (amount < 0 && errno == EINTR) continue;
      if (amount <= 0) return false;
      stats->bytes += uint64_t(amount);
      offset += uint64_t(amount);
      size_t consumed = size_t(amount);
      while (first < vectors.size() && consumed >= vectors[first].iov_len) {
        consumed -= vectors[first].iov_len;
        ++first;
      }
      if (consumed != 0) {
        vectors[first].iov_base =
            static_cast<uint8_t *>(vectors[first].iov_base) + consumed;
        vectors[first].iov_len -= consumed;
      }
    }
    begin = end;
  }
  return true;
}

#endif
