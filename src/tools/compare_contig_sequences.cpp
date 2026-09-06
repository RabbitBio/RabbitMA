// Compare two FASTA or MEGAHIT packed-contig files after canonicalizing
// strand and circular origin.  This is intentionally a standalone utility:
//
//   c++ -O3 -std=c++11 compare_contig_sequences.cpp -o compare-contigs
//   compare-contigs left.fa[.mgb] right.fa[.mgb] [examples]

// The comparison is exact (canonical strings are sorted, not fingerprinted)
// and retains duplicate multiplicity.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct PackedFileHeader {
  char magic[8];
  uint32_t version;
  uint32_t reserved;
};

struct PackedRecordHeader {
  uint32_t length;
  uint32_t kmer_size;
  int32_t flag;
  float multiplicity;
  int64_t id;
};

static_assert(sizeof(PackedFileHeader) == 16, "packed header layout");
static_assert(sizeof(PackedRecordHeader) == 24, "packed record layout");

char Complement(char base) {
  switch (base) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    default: throw std::runtime_error("invalid DNA base");
  }
}

std::string ReverseComplement(const std::string &sequence) {
  std::string result(sequence.size(), 'N');
  for (size_t i = 0; i < sequence.size(); ++i) {
    result[sequence.size() - 1u - i] = Complement(sequence[i]);
  }
  return result;
}

size_t MinimalRotation(const std::string &sequence) {
  if (sequence.empty()) return 0;
  const size_t size = sequence.size();
  size_t first = 0;
  size_t second = 1;
  size_t matched = 0;
  while (first < size && second < size && matched < size) {
    const char lhs = sequence[(first + matched) % size];
    const char rhs = sequence[(second + matched) % size];
    if (lhs == rhs) {
      ++matched;
      continue;
    }
    if (lhs > rhs) {
      first += matched + 1u;
      if (first <= second) first = second + 1u;
    } else {
      second += matched + 1u;
      if (second <= first) second = first + 1u;
    }
    matched = 0;
  }
  return std::min(first, second);
}

std::string Rotate(const std::string &sequence, size_t offset) {
  if (sequence.empty() || offset == 0) return sequence;
  return sequence.substr(offset) + sequence.substr(0, offset);
}

std::string CanonicalLinear(std::string sequence) {
  std::string reverse = ReverseComplement(sequence);
  if (reverse < sequence) sequence.swap(reverse);
  return sequence;
}

std::string CanonicalCircular(std::string sequence, uint32_t kmer_size) {
  // Packed loop records repeat exactly k bases at the end.  The biological
  // circular sequence is the non-overlapping core.
  if (sequence.size() > kmer_size) sequence.resize(sequence.size() - kmer_size);
  std::string reverse = ReverseComplement(sequence);
  sequence = Rotate(sequence, MinimalRotation(sequence));
  reverse = Rotate(reverse, MinimalRotation(reverse));
  if (reverse < sequence) sequence.swap(reverse);
  return sequence;
}

bool LooksPacked(const std::string &path) {
  std::ifstream input(path.c_str(), std::ios::binary);
  char magic[8]{};
  input.read(magic, sizeof(magic));
  static const char expected[8] = {'M', 'G', 'C', 'T', 'G', '0', '1', '\0'};
  return input.gcount() == static_cast<std::streamsize>(sizeof(magic)) &&
         std::memcmp(magic, expected, sizeof(magic)) == 0;
}

void ReadPacked(const std::string &path, std::vector<std::string> *output) {
  std::ifstream input(path.c_str(), std::ios::binary);
  PackedFileHeader file_header{};
  input.read(reinterpret_cast<char *>(&file_header), sizeof(file_header));
  if (!input || file_header.version != 1u) {
    throw std::runtime_error("invalid packed contig file: " + path);
  }
  std::vector<uint32_t> words;
  while (true) {
    PackedRecordHeader header{};
    input.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (input.gcount() == 0 && input.eof()) break;
    if (!input) throw std::runtime_error("truncated packed header: " + path);
    words.resize((static_cast<size_t>(header.length) + 15u) / 16u);
    input.read(reinterpret_cast<char *>(words.data()),
               words.size() * sizeof(uint32_t));
    if (!input) throw std::runtime_error("truncated packed sequence: " + path);
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string sequence(header.length, 'A');
    for (uint32_t i = 0; i < header.length; ++i) {
      sequence[i] = bases[(words[i >> 4u] >>
                           ((15u - (i & 15u)) << 1u)) & 3u];
    }
    output->push_back((header.flag & 2) != 0
                          ? CanonicalCircular(std::move(sequence),
                                              header.kmer_size)
                          : CanonicalLinear(std::move(sequence)));
  }
}

void AppendFastaRecord(std::string *sequence,
                       std::vector<std::string> *output) {
  if (sequence->empty()) return;
  // A general FASTA has no machine-readable loop flag.  Strand is still
  // normalized exactly; packed MEGAHIT output should be used when circular
  // origin equivalence is required.
  output->push_back(CanonicalLinear(std::move(*sequence)));
  sequence->clear();
}

void ReadFasta(const std::string &path, std::vector<std::string> *output) {
  std::ifstream input(path.c_str());
  if (!input) throw std::runtime_error("cannot open FASTA: " + path);
  std::string line;
  std::string sequence;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line.front() == '>') {
      AppendFastaRecord(&sequence, output);
    } else {
      for (char base : line) {
        if (base >= 'a' && base <= 'z') base -= 'a' - 'A';
        if (base != ' ' && base != '\t') sequence.push_back(base);
      }
    }
  }
  AppendFastaRecord(&sequence, output);
}

std::vector<std::string> ReadCanonical(const std::string &path) {
  std::vector<std::string> result;
  if (LooksPacked(path)) ReadPacked(path, &result);
  else ReadFasta(path, &result);
  std::sort(result.begin(), result.end());
  return result;
}

std::string Abbreviate(const std::string &sequence) {
  const size_t limit = 160u;
  if (sequence.size() <= limit) return sequence;
  return sequence.substr(0, limit / 2u) + "..." +
         sequence.substr(sequence.size() - limit / 2u);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr << "usage: " << argv[0]
              << " LEFT.fa[.mgb] RIGHT.fa[.mgb] [examples]\n";
    return 2;
  }
  try {
    const size_t example_limit =
        argc == 4 ? static_cast<size_t>(std::stoull(argv[3])) : 5u;
    std::vector<std::string> left = ReadCanonical(argv[1]);
    std::vector<std::string> right = ReadCanonical(argv[2]);
    size_t i = 0;
    size_t j = 0;
    uint64_t common = 0;
    uint64_t only_left = 0;
    uint64_t only_right = 0;
    std::vector<std::string> left_examples;
    std::vector<std::string> right_examples;
    while (i < left.size() || j < right.size()) {
      if (j == right.size() ||
          (i < left.size() && left[i] < right[j])) {
        if (left_examples.size() < example_limit) left_examples.push_back(left[i]);
        ++only_left;
        ++i;
      } else if (i == left.size() || right[j] < left[i]) {
        if (right_examples.size() < example_limit) right_examples.push_back(right[j]);
        ++only_right;
        ++j;
      } else {
        ++common;
        ++i;
        ++j;
      }
    }
    const uint64_t denominator = std::max(left.size(), right.size());
    const double mismatch_percent = denominator == 0
        ? 0.0
        : 100.0 * static_cast<double>(std::max(only_left, only_right)) /
              static_cast<double>(denominator);
    std::cout << "left=" << left.size() << " right=" << right.size()
              << " common=" << common << " only_left=" << only_left
              << " only_right=" << only_right << " mismatch="
              << std::fixed << std::setprecision(9) << mismatch_percent
              << "%\n";
    for (const std::string &sequence : left_examples) {
      std::cout << "only_left len=" << sequence.size() << " "
                << Abbreviate(sequence) << '\n';
    }
    for (const std::string &sequence : right_examples) {
      std::cout << "only_right len=" << sequence.size() << " "
                << Abbreviate(sequence) << '\n';
    }
    return only_left == 0 && only_right == 0 ? 0 : 1;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
