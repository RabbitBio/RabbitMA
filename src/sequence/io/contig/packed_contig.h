#ifndef MEGAHIT_PACKED_CONTIG_H
#define MEGAHIT_PACKED_CONTIG_H

#include <cstdint>
#include <cstring>
#include <string>

namespace packed_contig {

static const char kMagic[8] = {'M', 'G', 'C', 'T', 'G', '0', '1', '\0'};
static constexpr uint32_t kVersion = 1u;

struct FileHeader {
  char magic[8];
  uint32_t version;
  uint32_t reserved;
};

struct RecordHeader {
  uint32_t length;
  uint32_t kmer_size;
  int32_t flag;
  float multiplicity;
  int64_t id;
};

static_assert(sizeof(FileHeader) == 16u, "packed contig header layout");
static_assert(sizeof(RecordHeader) == 24u, "packed contig record layout");

inline FileHeader MakeFileHeader() {
  FileHeader header{};
  std::memcpy(header.magic, kMagic, sizeof(kMagic));
  header.version = kVersion;
  return header;
}

inline bool IsValid(const FileHeader &header) {
  return std::memcmp(header.magic, kMagic, sizeof(kMagic)) == 0 &&
         header.version == kVersion;
}

inline std::string Path(const std::string &fasta_path) {
  return fasta_path + ".mgb";
}

}  // namespace packed_contig

#endif  // MEGAHIT_PACKED_CONTIG_H
