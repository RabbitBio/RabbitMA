#ifndef MEGAHIT_SEQUENCE_IO_READ_ANCHOR_POSITIONS_H_
#define MEGAHIT_SEQUENCE_IO_READ_ANCHOR_POSITIONS_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// A buildlib side stream containing only the minimizer positions selected in
// each packed read.  Keys are deliberately not materialized: read-index can
// recover an a-mer from the packed library in at most three word operations.
// Keeping this representation independent of locator and bucket widths lets
// the later index builder retain all of its resource-adaptive choices.
constexpr uint64_t kReadAnchorPositionMagic =
    UINT64_C(0x5241504f53563131);  // "RAPOSV11"
constexpr uint32_t kReadAnchorPositionVersion = 1u;

struct ReadAnchorPositionHeader {
  uint64_t magic{kReadAnchorPositionMagic};
  uint32_t version{kReadAnchorPositionVersion};
  uint32_t header_bytes{sizeof(ReadAnchorPositionHeader)};
  uint32_t anchor_len{0};
  uint32_t window_len{0};
  uint64_t source_bytes{0};
  uint64_t num_reads{0};
  uint64_t num_bases{0};
  uint64_t payload_bytes{0};
  uint64_t chunk_index_offset{0};
  uint64_t num_chunks{0};
};

struct ReadAnchorPositionChunk {
  ReadAnchorPositionChunk() = default;
  ReadAnchorPositionChunk(uint64_t word_begin_arg, uint64_t word_end_arg,
                          uint64_t read_begin_arg, uint64_t read_end_arg,
                          uint64_t payload_begin_arg,
                          uint64_t payload_end_arg)
      : word_begin(word_begin_arg),
        word_end(word_end_arg),
        read_begin(read_begin_arg),
        read_end(read_end_arg),
        payload_begin(payload_begin_arg),
        payload_end(payload_end_arg) {}
  uint64_t word_begin{0};
  uint64_t word_end{0};
  uint64_t read_begin{0};
  uint64_t read_end{0};
  uint64_t payload_begin{0};
  uint64_t payload_end{0};
};

inline std::string ReadAnchorPositionPath(const std::string &read_path) {
  return read_path + ".ridx.positions";
}

inline uint64_t ReadAnchorMix64(uint64_t value) {
  value ^= value >> 30u;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27u;
  value *= UINT64_C(0x94d049bb133111eb);
  value ^= value >> 31u;
  return value;
}

inline uint64_t ReadAnchorRankHash(uint64_t key) {
  return ReadAnchorMix64(key ^ UINT64_C(0x9e3779b97f4a7c15));
}

inline uint8_t PackedReadBase(const uint32_t *sequence, uint32_t position) {
  const unsigned lane = position & 15u;
  return static_cast<uint8_t>(
      (sequence[position >> 4u] >> ((15u - lane) << 1u)) & 3u);
}

// Extract the exact forward a-mer representation produced by the rolling
// selector.  At most three packed uint32 words are touched for a <= 31.
inline uint64_t ExtractPackedReadKey(const uint32_t *sequence,
                                     uint32_t start,
                                     unsigned anchor_len) {
  uint64_t key = 0;
  unsigned remaining = anchor_len;
  uint32_t word = start >> 4u;
  unsigned lane = start & 15u;
  while (remaining != 0u) {
    const unsigned take = std::min<unsigned>(remaining, 16u - lane);
    const unsigned used_bits = take << 1u;
    const uint32_t aligned = sequence[word] << (lane << 1u);
    const uint32_t segment =
        used_bits == 32u ? aligned : aligned >> (32u - used_bits);
    key = (key << used_bits) | segment;
    remaining -= take;
    ++word;
    lane = 0u;
  }
  return key;
}

template <class Visitor>
inline void ForEachPackedReadAnchor(
    const uint32_t *sequence, unsigned length, unsigned anchor_len,
    unsigned window_len, std::vector<uint32_t> *queue_pos,
    std::vector<uint64_t> *queue_hash, std::vector<uint64_t> *queue_key,
    const Visitor &visitor) {
  if (length < window_len) return;
  const unsigned num_anchor_positions = length - anchor_len + 1u;
  const unsigned anchors_per_window = window_len - anchor_len + 1u;
  queue_pos->resize(num_anchor_positions);
  queue_hash->resize(num_anchor_positions);
  queue_key->resize(num_anchor_positions);

  const uint64_t key_mask =
      anchor_len == 32u
          ? std::numeric_limits<uint64_t>::max()
          : (uint64_t{1} << (2u * anchor_len)) - 1u;
  uint64_t key = 0;
  size_t queue_head = 0;
  size_t queue_tail = 0;
  uint32_t last_emitted = std::numeric_limits<uint32_t>::max();
  for (unsigned base_pos = 0; base_pos < length; ++base_pos) {
    key = ((key << 2u) | PackedReadBase(sequence, base_pos)) & key_mask;
    if (base_pos + 1u < anchor_len) continue;
    const uint32_t anchor_pos = base_pos + 1u - anchor_len;
    const uint64_t hash = ReadAnchorRankHash(key);
    while (queue_head < queue_tail &&
           (*queue_pos)[queue_head] + anchors_per_window <= anchor_pos) {
      ++queue_head;
    }
    // Right-most tie breaking is part of the exact on-disk index semantics.
    while (queue_head < queue_tail &&
           hash <= (*queue_hash)[queue_tail - 1u]) {
      --queue_tail;
    }
    (*queue_pos)[queue_tail] = anchor_pos;
    (*queue_hash)[queue_tail] = hash;
    (*queue_key)[queue_tail] = key;
    ++queue_tail;
    if (anchor_pos + 1u < anchors_per_window) continue;
    const uint32_t selected_pos = (*queue_pos)[queue_head];
    if (selected_pos == last_emitted) continue;
    last_emitted = selected_pos;
    visitor((*queue_key)[queue_head], selected_pos);
  }
}

inline void AppendReadAnchorVarint(uint32_t value,
                                   std::vector<uint8_t> *output) {
  while (value >= 0x80u) {
    output->push_back(static_cast<uint8_t>(value | 0x80u));
    value >>= 7u;
  }
  output->push_back(static_cast<uint8_t>(value));
}

inline bool ReadReadAnchorVarint(const uint8_t **cursor, const uint8_t *end,
                                 uint32_t *value) {
  uint32_t result = 0;
  unsigned shift = 0;
  while (*cursor != end && shift < 35u) {
    const uint8_t byte = *(*cursor)++;
    result |= static_cast<uint32_t>(byte & 0x7fu) << shift;
    if ((byte & 0x80u) == 0u) {
      *value = result;
      return true;
    }
    shift += 7u;
  }
  return false;
}

#endif  // MEGAHIT_SEQUENCE_IO_READ_ANCHOR_POSITIONS_H_
