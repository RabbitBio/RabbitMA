#ifndef MEGAHIT_EDGE_BUCKET_HISTOGRAM_H
#define MEGAHIT_EDGE_BUCKET_HISTOGRAM_H

#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// A count->SeqToSdbg handoff for the first radix histogram.  Rows correspond
// exactly to the sharded edge files and columns to the 8-base radix prefix.
// The data is deliberately a small, versioned sidecar: edge files remain the
// source of truth and older/newer executables simply fall back to scanning.
struct EdgeBucketHistogram {
  static constexpr uint32_t kVersion = 1;
  static constexpr uint32_t kNumBuckets = 1u << 16u;

  uint32_t kmer_size{0};
  uint32_t num_files{0};
  uint64_t num_edges{0};
  std::vector<uint64_t> counts;

  bool Write(const std::string &prefix) const {
    if (num_files == 0 || counts.size() !=
                              static_cast<size_t>(num_files) * kNumBuckets) {
      return false;
    }
    std::ofstream output(prefix + ".edges.lv0",
                         std::ios::binary | std::ios::trunc);
    if (!output) {
      return false;
    }
    static const char kMagic[8] = {'M', 'H', 'L', 'V', '0', 'H', '1', '\0'};
    const uint32_t version = kVersion;
    const uint32_t num_buckets = kNumBuckets;
    output.write(kMagic, sizeof(kMagic));
    output.write(reinterpret_cast<const char *>(&version), sizeof(version));
    output.write(reinterpret_cast<const char *>(&kmer_size),
                 sizeof(kmer_size));
    output.write(reinterpret_cast<const char *>(&num_files),
                 sizeof(num_files));
    output.write(reinterpret_cast<const char *>(&num_buckets),
                 sizeof(num_buckets));
    output.write(reinterpret_cast<const char *>(&num_edges),
                 sizeof(num_edges));
    output.write(reinterpret_cast<const char *>(counts.data()),
                 counts.size() * sizeof(counts[0]));
    return output.good();
  }

  bool Read(const std::string &prefix) {
    std::ifstream input(prefix + ".edges.lv0", std::ios::binary);
    if (!input) {
      return false;
    }
    char magic[8];
    uint32_t version = 0;
    uint32_t num_buckets = 0;
    input.read(magic, sizeof(magic));
    input.read(reinterpret_cast<char *>(&version), sizeof(version));
    input.read(reinterpret_cast<char *>(&kmer_size), sizeof(kmer_size));
    input.read(reinterpret_cast<char *>(&num_files), sizeof(num_files));
    input.read(reinterpret_cast<char *>(&num_buckets), sizeof(num_buckets));
    input.read(reinterpret_cast<char *>(&num_edges), sizeof(num_edges));
    static const char kMagic[8] = {'M', 'H', 'L', 'V', '0', 'H', '1', '\0'};
    if (!input || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0 ||
        version != kVersion || num_buckets != kNumBuckets || num_files == 0 ||
        num_files > std::numeric_limits<size_t>::max() / kNumBuckets) {
      counts.clear();
      return false;
    }
    counts.resize(static_cast<size_t>(num_files) * kNumBuckets);
    input.read(reinterpret_cast<char *>(counts.data()),
               counts.size() * sizeof(counts[0]));
    if (!input) {
      counts.clear();
      return false;
    }
    // Reject trailing/truncated format variants rather than silently using a
    // histogram whose row interpretation may differ.
    char trailing;
    if (input.read(&trailing, 1)) {
      counts.clear();
      return false;
    }
    return input.eof();
  }
};

#endif  // MEGAHIT_EDGE_BUCKET_HISTOGRAM_H
