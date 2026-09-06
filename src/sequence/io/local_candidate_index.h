#ifndef MEGAHIT_LOCAL_CANDIDATE_INDEX_H
#define MEGAHIT_LOCAL_CANDIDATE_INDEX_H

#include <cstdint>
#include <cstring>
#include <type_traits>

// A local-assembly candidate is the uint32-word offset of a packed read
// record's length header in reads.lib.bin.  Offsets are stored sorted and
// unique, so local mapping can merge them with the normal sequential record
// walk without one hash lookup per read.
struct LocalCandidateFileHeader {
  char magic[8];
  uint32_t version;
  uint32_t header_bytes;
  uint64_t source_words;
  uint64_t source_reads;
  uint64_t candidate_count;
};

static_assert(std::is_trivially_copyable<LocalCandidateFileHeader>::value,
              "candidate header must be directly serializable");
static_assert(sizeof(LocalCandidateFileHeader) == 40u,
              "candidate header layout changed");

inline LocalCandidateFileHeader MakeLocalCandidateFileHeader(
    uint64_t source_words, uint64_t source_reads,
    uint64_t candidate_count) {
  LocalCandidateFileHeader header{};
  const char magic[8] = {'M', 'H', 'L', 'C', 'A', 'N', 'D', '1'};
  std::memcpy(header.magic, magic, sizeof(magic));
  header.version = 1u;
  header.header_bytes = sizeof(LocalCandidateFileHeader);
  header.source_words = source_words;
  header.source_reads = source_reads;
  header.candidate_count = candidate_count;
  return header;
}

inline bool IsValidLocalCandidateFileHeader(
    const LocalCandidateFileHeader &header) {
  const char magic[8] = {'M', 'H', 'L', 'C', 'A', 'N', 'D', '1'};
  return std::memcmp(header.magic, magic, sizeof(magic)) == 0 &&
         header.version == 1u &&
         header.header_bytes == sizeof(LocalCandidateFileHeader);
}

#endif  // MEGAHIT_LOCAL_CANDIDATE_INDEX_H
