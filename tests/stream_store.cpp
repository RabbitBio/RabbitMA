#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "utils/stream_store.h"

int main() {
  std::vector<uint8_t> source(2048), expected(2048), actual(2048);
  for (size_t i = 0; i < source.size(); ++i) source[i] = (i * 71u + i / 13u) & 255u;
  uint64_t comparisons = 0;
  for (size_t from = 0; from < 64; ++from) {
    for (size_t to = 0; to < 64; ++to) {
      for (size_t length : {0u, 1u, 3u, 15u, 16u, 31u, 32u, 63u, 64u, 65u,
                            127u, 128u, 129u, 255u, 256u, 320u, 384u, 1025u}) {
        std::fill(expected.begin(), expected.end(), 0xA5u);
        actual = expected;
        std::memcpy(expected.data() + to, source.data() + from, length);
        StreamStoreCopy(actual.data() + to, source.data() + from, length);
        StreamStoreFence();
        if (actual != expected) throw std::runtime_error("bytes or boundary guard changed");
        ++comparisons;
      }
    }
  }
  // Adjacent ranges can share their boundary cache line. Every interior
  // streaming store must remain strictly inside its owning range.
  for (size_t length : {17u, 65u, 129u, 256u, 320u, 384u}) {
    std::vector<uint8_t> output(length * 8u + 128u, 0xA5u);
    for (unsigned repeat = 0; repeat < 100; ++repeat) {
#pragma omp parallel for num_threads(8)
      for (unsigned thread = 0; thread < 8; ++thread) {
        StreamStoreCopy(output.data() + 3u + thread * length,
                         source.data() + thread, length);
        StreamStoreFence();
      }
      for (unsigned thread = 0; thread < 8; ++thread)
        if (std::memcmp(output.data() + 3u + thread * length,
                        source.data() + thread, length) != 0)
          throw std::runtime_error("adjacent worker range changed");
      if (output[2] != 0xA5u || output[3u + 8u * length] != 0xA5u)
        throw std::runtime_error("worker boundary guard changed");
      ++comparisons;
    }
  }
  std::cout << comparisons << " exact byte and boundary comparisons\n";
}
