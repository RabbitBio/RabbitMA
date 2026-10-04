#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "utils/coalesced_reads.h"

void Require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  char path[] = "/tmp/rabbitma-coalesced-read-XXXXXX";
  const int fd = mkstemp(path);
  Require(fd >= 0, "cannot create test file");
  unlink(path);
  std::vector<uint8_t> data(1u << 20u);
  for (size_t i = 0; i < data.size(); ++i) data[i] = (i * 71 + i / 13) & 255;
  Require(write(fd, data.data(), data.size()) == ssize_t(data.size()),
          "cannot populate test file");
  for (bool gaps : {false, true}) {
    const size_t count = 4097;
    const size_t width = 17;
    std::vector<uint8_t> output(count * (width + 2), 0xA5);
    std::vector<PositionedRead> requests;
    for (size_t i = 0; i < count; ++i) {
      const size_t offset = i * width + (gaps ? i / 77 * 5 : 0);
      requests.push_back(PositionedRead{
          offset, width, output.data() + i * (width + 2) + 1});
    }
    std::reverse(requests.begin(), requests.end());
    CoalescedReadStats stats;
    Require(ReadCoalescedExtents(fd, &requests, &stats), "scatter read failed");
    Require(stats.bytes == count * width, "gap bytes were read or bytes lost");
    Require(stats.calls < count / 10, "contiguous requests were not coalesced");
    for (size_t i = 0; i < count; ++i) {
      const size_t offset = i * width + (gaps ? i / 77 * 5 : 0);
      Require(std::equal(data.begin() + offset, data.begin() + offset + width,
                         output.begin() + i * (width + 2) + 1),
              "bytes reached the wrong destination");
      Require(output[i * (width + 2)] == 0xA5 &&
                  output[i * (width + 2) + width + 1] == 0xA5,
              "scatter overwrote a destination boundary");
    }
  }

  // preadv first returns a short read inside an iovec, then EOF. Only bytes
  // actually present in the file may be written into the destination.
  std::vector<uint8_t> tail(32, 0xA5);
  std::vector<PositionedRead> truncated{
      {data.size() - 8, 17, tail.data()},
      {data.size() + 9, 15, tail.data() + 17}};
  CoalescedReadStats short_stats;
  Require(!ReadCoalescedExtents(fd, &truncated, &short_stats),
          "truncated file must fail");
  Require(short_stats.bytes == 8 && short_stats.calls == 2,
          "short-read progress is wrong");
  Require(std::equal(data.end() - 8, data.end(), tail.begin()),
          "short-read bytes differ");
  Require(std::all_of(tail.begin() + 8, tail.end(),
                       [](uint8_t c) { return c == 0xA5; }),
          "short read wrote nonexistent bytes");

  std::vector<PositionedRead> empty{{0, 0, nullptr}};
  CoalescedReadStats empty_stats;
  Require(ReadCoalescedExtents(fd, &empty, &empty_stats) && empty_stats.calls == 0,
          "zero-length extent must not read");
  std::vector<PositionedRead> overflow{
      {std::numeric_limits<uint64_t>::max(), 1, tail.data()}};
  Require(!ReadCoalescedExtents(fd, &overflow, &empty_stats),
          "unrepresentable file offset must fail");
  close(fd);
  std::cout << "Coalesced reads preserve exact destinations, gaps, iovec limits, "
               "short reads and guards\n";
}
