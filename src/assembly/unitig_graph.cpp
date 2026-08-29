//
// Created by vout on 11/19/18.
//

#include "unitig_graph.h"
#include <omp.h>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif

#include "kmlib/kmbitvector.h"
#include "utils/mutex.h"
#include "utils/utils.h"

namespace {

#ifdef __linux__
bool ReadMemAvailable(uint64_t *available_bytes) {
  std::ifstream input("/proc/meminfo");
  std::string line;
  while (std::getline(input, line)) {
    if (line.compare(0, 13, "MemAvailable:") != 0) {
      continue;
    }

    std::istringstream fields(line);
    std::string name;
    std::string unit;
    uint64_t kibibytes = 0;
    if (!(fields >> name >> kibibytes >> unit) || name != "MemAvailable:" ||
        unit != "kB" ||
        kibibytes > std::numeric_limits<uint64_t>::max() / 1024) {
      return false;
    }
    *available_bytes = kibibytes * 1024;
    return true;
  }
  return false;
}

bool ReadCgroupValue(const std::string &path, uint64_t *value) {
  std::ifstream input(path);
  std::string text;
  if (!(input >> text) || text == "max") {
    return false;
  }

  char *end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno != 0 || end == text.c_str() || *end != '\0') {
    return false;
  }
  *value = static_cast<uint64_t>(parsed);
  return true;
}

std::string CurrentMemoryCgroup(bool unified) {
  std::ifstream input("/proc/self/cgroup");
  std::string line;
  while (std::getline(input, line)) {
    const size_t first_colon = line.find(':');
    const size_t second_colon =
        first_colon == std::string::npos
            ? std::string::npos
            : line.find(':', first_colon + 1);
    if (second_colon == std::string::npos) {
      continue;
    }
    const std::string controllers =
        line.substr(first_colon + 1, second_colon - first_colon - 1);
    bool is_memory_controller = false;
    size_t begin = 0;
    while (begin <= controllers.size()) {
      const size_t comma = controllers.find(',', begin);
      const size_t end =
          comma == std::string::npos ? controllers.size() : comma;
      if (controllers.compare(begin, end - begin, "memory") == 0) {
        is_memory_controller = true;
        break;
      }
      if (comma == std::string::npos) {
        break;
      }
      begin = comma + 1;
    }
    if ((unified && controllers.empty()) ||
        (!unified && is_memory_controller)) {
      return line.substr(second_colon + 1);
    }
  }
  return {};
}

void TightenAvailableMemoryFromCgroupHierarchy(
    const std::string &root, const std::string &relative_path,
    const char *limit_name, const char *usage_name, uint64_t *available_bytes,
    bool *found_available_memory) {
  std::string directory = root;
  if (!relative_path.empty() && relative_path != "/") {
    if (relative_path.front() != '/') {
      directory.push_back('/');
    }
    directory += relative_path;
    while (directory.size() > root.size() && directory.back() == '/') {
      directory.pop_back();
    }
  }

  while (directory.size() >= root.size() &&
         directory.compare(0, root.size(), root) == 0) {
    uint64_t limit = 0;
    uint64_t usage = 0;
    if (ReadCgroupValue(directory + "/" + limit_name, &limit) &&
        ReadCgroupValue(directory + "/" + usage_name, &usage) &&
        // cgroup v1 represents an unlimited limit as a value close to
        // INT64_MAX.  No physical host can provide one exbibyte here, so such
        // values are not finite constraints.
        limit < (uint64_t{1} << 60u)) {
      const uint64_t remaining = limit > usage ? limit - usage : 0;
      if (!*found_available_memory || remaining < *available_bytes) {
        *available_bytes = remaining;
        *found_available_memory = true;
      }
    }

    if (directory == root) {
      break;
    }
    const size_t slash = directory.find_last_of('/');
    if (slash == std::string::npos || slash < root.size()) {
      directory = root;
    } else {
      directory.resize(slash);
    }
  }
}

void TightenAvailableMemoryFromCgroup(uint64_t *available_bytes,
                                      bool *found_available_memory) {
  const std::string unified_path = CurrentMemoryCgroup(true);
  if (!unified_path.empty()) {
    TightenAvailableMemoryFromCgroupHierarchy(
        "/sys/fs/cgroup", unified_path, "memory.max", "memory.current",
        available_bytes, found_available_memory);
  }

  const std::string legacy_path = CurrentMemoryCgroup(false);
  if (!legacy_path.empty()) {
    TightenAvailableMemoryFromCgroupHierarchy(
        "/sys/fs/cgroup/memory", legacy_path, "memory.limit_in_bytes",
        "memory.usage_in_bytes", available_bytes, found_available_memory);
  }
}
#endif

constexpr uint32_t kNoSimpleNeighbor = std::numeric_limits<uint32_t>::max();
constexpr uint32_t kNoCompactSimpleNeighbor =
    (static_cast<uint32_t>(1) << 28u) - 1u;
constexpr uint64_t kCompactSimpleNeighborBytesPerEdge = 7;
constexpr uint64_t kBlockSimpleNeighborEdges = 256;

struct WideDoubleNeighborAccessor {
  uint32_t *next;
  uint32_t *prev;

  void StoreNext(uint64_t edge, uint32_t value) const { next[edge] = value; }
  void StorePrev(uint64_t edge, uint32_t value) const { prev[edge] = value; }
  uint32_t LoadNext(uint64_t edge) const { return next[edge]; }
  uint32_t LoadPrev(uint64_t edge) const { return prev[edge]; }
  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadNext(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
  uint64_t Prev(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadPrev(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

struct CompactDoubleNeighborAccessor {
  uint8_t *records;

  void StoreNext(uint64_t edge, uint32_t value) const {
    assert(value == kNoSimpleNeighbor || value < kNoCompactSimpleNeighbor);
    const uint32_t packed = value == kNoSimpleNeighbor
                                ? kNoCompactSimpleNeighbor
                                : value;
    uint8_t *record =
        records + edge * kCompactSimpleNeighborBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    std::memcpy(record, &packed, sizeof(packed));
#else
    record[0] = static_cast<uint8_t>(packed);
    record[1] = static_cast<uint8_t>(packed >> 8u);
    record[2] = static_cast<uint8_t>(packed >> 16u);
    record[3] = static_cast<uint8_t>((packed >> 24u) & 0x0Fu);
#endif
  }

  void StorePrev(uint64_t edge, uint32_t value) const {
    assert(value == kNoSimpleNeighbor || value < kNoCompactSimpleNeighbor);
    const uint32_t packed = value == kNoSimpleNeighbor
                                ? kNoCompactSimpleNeighbor
                                : value;
    uint8_t *record =
        records + edge * kCompactSimpleNeighborBytesPerEdge;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    const uint32_t encoded =
        (packed << 4u) | static_cast<uint32_t>(record[3] & 0x0Fu);
    std::memcpy(record + 3, &encoded, sizeof(encoded));
#else
    record[3] = static_cast<uint8_t>((record[3] & 0x0Fu) |
                                     ((packed & 0x0Fu) << 4u));
    record[4] = static_cast<uint8_t>(packed >> 4u);
    record[5] = static_cast<uint8_t>(packed >> 12u);
    record[6] = static_cast<uint8_t>(packed >> 20u);
#endif
  }

  uint32_t LoadNext(uint64_t edge) const {
    const uint8_t *record =
        records + edge * kCompactSimpleNeighborBytesPerEdge;
    uint32_t packed;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    std::memcpy(&packed, record, sizeof(packed));
    packed &= kNoCompactSimpleNeighbor;
#else
    packed = static_cast<uint32_t>(record[0]) |
             (static_cast<uint32_t>(record[1]) << 8u) |
             (static_cast<uint32_t>(record[2]) << 16u) |
             (static_cast<uint32_t>(record[3] & 0x0Fu) << 24u);
#endif
    return packed == kNoCompactSimpleNeighbor ? kNoSimpleNeighbor : packed;
  }

  uint32_t LoadPrev(uint64_t edge) const {
    const uint8_t *record =
        records + edge * kCompactSimpleNeighborBytesPerEdge;
    uint32_t packed;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    std::memcpy(&packed, record + 3, sizeof(packed));
    packed >>= 4u;
#else
    packed = static_cast<uint32_t>(record[3] >> 4u) |
             (static_cast<uint32_t>(record[4]) << 4u) |
             (static_cast<uint32_t>(record[5]) << 12u) |
             (static_cast<uint32_t>(record[6]) << 20u);
#endif
    return packed == kNoCompactSimpleNeighbor ? kNoSimpleNeighbor : packed;
  }

  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadNext(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
  uint64_t Prev(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadPrev(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

struct XorNeighborAccessor {
  const uint32_t *link;

  uint64_t Next(uint64_t edge, uint64_t came_from) const {
    return Other(edge, came_from);
  }
  uint64_t Prev(uint64_t edge, uint64_t came_from) const {
    return Other(edge, came_from);
  }

 private:
  uint64_t Other(uint64_t edge, uint64_t came_from) const {
    const uint32_t encoded_came_from =
        came_from == SDBG::kNullID ? kNoSimpleNeighbor
                                   : static_cast<uint32_t>(came_from);
    const uint32_t value = link[edge] ^ encoded_came_from;
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

struct ForwardNeighborAccessor {
  const uint32_t *next;

  uint32_t LoadNext32(uint64_t edge) const { return next[edge]; }
  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadNext32(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

// A simple outgoing edge is located within the small BOSS edge group ending
// at Forward(edge).  Zero means no simple successor; values 1..15 encode the
// offset from that group end.  Packing two codes per byte reduces the
// temporary unitig snapshot from four bytes to half a byte per SDBG edge.
struct OffsetNeighborAccessor {
  SDBG *sdbg;
  const uint8_t *codes;

  uint8_t LoadCode(uint64_t edge) const {
    const uint8_t packed = codes[edge >> 1u];
    return (edge & 1u) != 0u ? packed >> 4u : packed & 0x0fu;
  }

  uint32_t LoadNext32(uint64_t edge) const {
    const uint8_t code = LoadCode(edge);
    if (code == 0) {
      return kNoSimpleNeighbor;
    }
    const uint64_t first_outgoing = sdbg->Forward(edge);
    assert(first_outgoing != SDBG::kNullID);
    assert(first_outgoing >= static_cast<uint64_t>(code - 1u));
    const uint64_t next = first_outgoing - (code - 1u);
    assert(next < static_cast<uint64_t>(kNoSimpleNeighbor));
    return static_cast<uint32_t>(next);
  }

  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadNext32(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

// Direct simple-next lookup with one byte per edge.  Absolute edge IDs are
// monotone for each normalized W symbol, so a small logical edge block shares
// one 32-bit base per alphabet symbol.  Code zero is null, 1..254 are exact
// offsets and 255 retains an exact (rare) succinct-topology fallback.
struct BlockNeighborAccessor {
  SDBG *sdbg;
  const uint8_t *codes;
  const uint32_t *bases;

  uint32_t LoadNext32(uint64_t edge) const {
    const uint8_t code = codes[edge];
    if (code == 0) {
      return kNoSimpleNeighbor;
    }
    if (code == 255u) {
      const uint64_t next = sdbg->NextSimplePathEdge(edge);
      return next == SDBG::kNullID ? kNoSimpleNeighbor
                                   : static_cast<uint32_t>(next);
    }
    uint8_t c = sdbg->GetW(edge);
    if (c > kAlphabetSize) {
      c -= kAlphabetSize;
    }
    const uint32_t base =
        bases[(edge / kBlockSimpleNeighborEdges) * (kAlphabetSize + 1u) + c];
    assert(base != kNoSimpleNeighbor);
    return base + static_cast<uint32_t>(code - 1u);
  }

  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t value = LoadNext32(edge);
    return value == kNoSimpleNeighbor ? SDBG::kNullID
                                      : static_cast<uint64_t>(value);
  }
};

// A snapshot built before SDBG tip pruning contains exact non-zero links that
// survived pruning.  Missing/stale entries are resolved from the live graph,
// so links newly exposed by pruning are represented without rebuilding the
// edge-sized cache.
struct SharedSnapshotNeighborAccessor {
  SDBG *sdbg;

  uint32_t LoadNext32(uint64_t edge) const {
    const uint64_t next = sdbg->CachedOrLiveNextSimplePathEdge(edge);
    return next == SDBG::kNullID ? kNoSimpleNeighbor
                                 : static_cast<uint32_t>(next);
  }

  uint64_t Next(uint64_t edge, uint64_t) const {
    const uint32_t next = LoadNext32(edge);
    return next == kNoSimpleNeighbor ? SDBG::kNullID
                                     : static_cast<uint64_t>(next);
  }
};

struct UncachedNeighborAccessor {
  SDBG *sdbg;
  const AtomicBitVector *has_simple_next;

  uint64_t Next(uint64_t edge, uint64_t) const {
    return has_simple_next->at(edge) ? sdbg->UniqueNextEdge(edge)
                                     : SDBG::kNullID;
  }
  uint64_t Prev(uint64_t edge, uint64_t) const {
    const uint64_t prev = sdbg->UniquePrevEdge(edge);
    return prev != SDBG::kNullID && has_simple_next->at(prev)
               ? prev
               : SDBG::kNullID;
  }
};

template <class NeighborAccessor>
void BuildDoubleSimpleNeighbors(
    SDBG *sdbg, uint64_t num_edges, int num_threads,
    const NeighborAccessor &neighbors,
    std::vector<std::vector<uint64_t>> *thread_terminals) {
#pragma omp parallel num_threads(num_threads)
  {
    const int tid = omp_get_thread_num();
    auto &local_terminals = (*thread_terminals)[tid];
    local_terminals.reserve(
        static_cast<size_t>(num_edges / std::max(1, num_threads) / 64 + 1));

#pragma omp for schedule(static)
    for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
      uint32_t next_id = kNoSimpleNeighbor;
      if (sdbg->IsValidEdge(edge_idx)) {
        const uint64_t next = sdbg->NextSimplePathEdge(edge_idx);
        if (next == SDBG::kNullID) {
          local_terminals.push_back(edge_idx);
        } else {
          assert(next < static_cast<uint64_t>(kNoSimpleNeighbor));
          next_id = static_cast<uint32_t>(next);
        }
      }
      // Store next before prev: compact records share byte 3, and the prev
      // setter deliberately preserves next's low nibble.
      neighbors.StoreNext(edge_idx, next_id);
      neighbors.StorePrev(edge_idx, kNoSimpleNeighbor);
    }
  }

  // A simple-next target has exactly one simple predecessor, so each slot has
  // one writer even though the inverse fill is parallel.
#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
    const uint32_t next = neighbors.LoadNext(edge_idx);
    if (next != kNoSimpleNeighbor) {
      assert(neighbors.LoadPrev(next) == kNoSimpleNeighbor);
      neighbors.StorePrev(next, static_cast<uint32_t>(edge_idx));
    }
  }
}

template <class NeighborAccessor>
size_t AssembleNonLoopPath(
    uint64_t edge_idx, int tid, SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const NeighborAccessor &neighbors) {
  if (!locks->try_lock(edge_idx)) {
    return 0;
  }

  bool will_be_added = true;
  uint64_t cur_edge = edge_idx;
  uint64_t prev_edge;
  uint64_t came_from = SDBG::kNullID;
  int64_t depth = sdbg->EdgeMultiplicity(edge_idx);
  uint32_t length = 1;

  while ((prev_edge = neighbors.Prev(cur_edge, came_from)) != SDBG::kNullID) {
    came_from = cur_edge;
    cur_edge = prev_edge;
    if (!locks->try_lock(cur_edge)) {
      will_be_added = false;
      break;
    }
    depth += sdbg->EdgeMultiplicity(cur_edge);
    ++length;
  }

  if (!will_be_added) {
    return 0;
  }

  const uint64_t rc_start = sdbg->EdgeReverseComplement(edge_idx);
  uint64_t rc_end;
  assert(rc_start != SDBG::kNullID);

  if (!locks->try_lock(rc_start)) {
    rc_end = sdbg->EdgeReverseComplement(cur_edge);
    if (std::max(edge_idx, cur_edge) < std::max(rc_start, rc_end)) {
      will_be_added = false;
    }
  } else {
    uint64_t rc_cur_edge = rc_start;
    uint64_t rc_came_from = SDBG::kNullID;
    rc_end = rc_cur_edge;
    bool extend_full = true;
    while (true) {
      const uint64_t rc_next_edge =
          neighbors.Next(rc_cur_edge, rc_came_from);
      if (rc_next_edge == SDBG::kNullID) {
        break;
      }
      rc_came_from = rc_cur_edge;
      rc_cur_edge = rc_next_edge;
      rc_end = rc_next_edge;
      if (!locks->try_lock(rc_next_edge)) {
        extend_full = false;
        break;
      }
    }
    if (!extend_full) {
      rc_end = sdbg->EdgeReverseComplement(cur_edge);
      assert(rc_end != SDBG::kNullID);
    }
  }

  if (will_be_added) {
    (*thread_vertices)[tid].emplace_back(cur_edge, edge_idx, rc_start, rc_end,
                                         depth, length);
    return cur_edge == rc_start;
  }
  return 0;
}

template <class NeighborAccessor>
size_t AssembleForwardNonLoopPath(
    uint64_t start_edge, int tid, SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const NeighborAccessor &neighbors) {
  if (!locks->try_lock(start_edge)) {
    return 0;
  }

  bool will_be_added = true;
  uint64_t cur_edge = start_edge;
  uint64_t came_from = SDBG::kNullID;
  int64_t depth = sdbg->EdgeMultiplicity(start_edge);
  uint32_t length = 1;

  while (true) {
    const uint64_t next_edge = neighbors.Next(cur_edge, came_from);
    if (next_edge == SDBG::kNullID) {
      break;
    }
    came_from = cur_edge;
    cur_edge = next_edge;
    if (!locks->try_lock(cur_edge)) {
      will_be_added = false;
      break;
    }
    depth += sdbg->EdgeMultiplicity(cur_edge);
    ++length;
  }

  if (!will_be_added) {
    return 0;
  }

  const uint64_t end_edge = cur_edge;
  const uint64_t rc_start = sdbg->EdgeReverseComplement(end_edge);
  uint64_t rc_end;
  assert(rc_start != SDBG::kNullID);

  if (!locks->try_lock(rc_start)) {
    rc_end = sdbg->EdgeReverseComplement(start_edge);
    if (std::max(end_edge, start_edge) < std::max(rc_start, rc_end)) {
      will_be_added = false;
    }
  } else {
    uint64_t rc_cur_edge = rc_start;
    uint64_t rc_came_from = SDBG::kNullID;
    rc_end = rc_cur_edge;
    bool extend_full = true;
    while (true) {
      const uint64_t rc_next_edge =
          neighbors.Next(rc_cur_edge, rc_came_from);
      if (rc_next_edge == SDBG::kNullID) {
        break;
      }
      rc_came_from = rc_cur_edge;
      rc_cur_edge = rc_next_edge;
      rc_end = rc_next_edge;
      if (!locks->try_lock(rc_next_edge)) {
        extend_full = false;
        break;
      }
    }
    if (!extend_full) {
      rc_end = sdbg->EdgeReverseComplement(start_edge);
      assert(rc_end != SDBG::kNullID);
    }
  }

  if (will_be_added) {
    (*thread_vertices)[tid].emplace_back(start_edge, end_edge, rc_start,
                                         rc_end, depth, length);
    return start_edge == rc_start;
  }
  return 0;
}

template <class NeighborAccessor>
size_t AssembleCachedNonLoopPaths(
    SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const std::vector<uint64_t> &terminal_ids, int num_threads,
    const NeighborAccessor &neighbors) {
  size_t count_palindrome = 0;
#pragma omp parallel num_threads(num_threads) reduction(+ : count_palindrome)
  {
    const int tid = omp_get_thread_num();
    const uint64_t team_size = omp_get_num_threads();
    const uint64_t block_size = sdbg->size() / team_size;
    const uint64_t remainder = sdbg->size() % team_size;
    const uint64_t edge_begin =
        static_cast<uint64_t>(tid) * block_size +
        std::min<uint64_t>(tid, remainder);
    const uint64_t edge_end =
        edge_begin + block_size + (static_cast<uint64_t>(tid) < remainder);
    auto terminal_begin =
        std::lower_bound(terminal_ids.begin(), terminal_ids.end(), edge_begin);
    auto terminal_end =
        std::lower_bound(terminal_begin, terminal_ids.end(), edge_end);
    for (auto it = terminal_begin; it != terminal_end; ++it) {
      count_palindrome += AssembleNonLoopPath(
          *it, tid, sdbg, locks, thread_vertices, neighbors);
    }
  }
  return count_palindrome;
}

template <class NeighborAccessor>
size_t AssembleForwardNonLoopPaths(
    SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const std::vector<uint64_t> &start_ids, int num_threads,
    const NeighborAccessor &neighbors) {
  size_t count_palindrome = 0;
#pragma omp parallel num_threads(num_threads) reduction(+ : count_palindrome)
  {
    const int tid = omp_get_thread_num();
    const uint64_t team_size = omp_get_num_threads();
    const uint64_t block_size = sdbg->size() / team_size;
    const uint64_t remainder = sdbg->size() % team_size;
    const uint64_t edge_begin =
        static_cast<uint64_t>(tid) * block_size +
        std::min<uint64_t>(tid, remainder);
    const uint64_t edge_end =
        edge_begin + block_size + (static_cast<uint64_t>(tid) < remainder);
    auto start_begin =
        std::lower_bound(start_ids.begin(), start_ids.end(), edge_begin);
    auto start_end = std::lower_bound(start_begin, start_ids.end(), edge_end);
    for (auto it = start_begin; it != start_end; ++it) {
      count_palindrome += AssembleForwardNonLoopPath(
          *it, tid, sdbg, locks, thread_vertices, neighbors);
    }
  }
  return count_palindrome;
}

// Interleave independent simple-path walks within each worker.  A single path
// is a dependency chain (next[cur] -> next[next[cur]]), but different starts
// have no such dependency.  Advancing several states round-robin exposes
// memory-level parallelism for the dense next array, SDBG multiplicities and
// lock words without changing graph semantics.
template <class ForwardAccessor>
size_t AssembleForwardNonLoopPathsMultistream(
    SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const std::vector<uint64_t> &start_ids, int num_threads,
    const ForwardAccessor &neighbors, unsigned stream_width) {
  struct PathState {
    bool active{false};
    bool reverse_phase{false};
    uint64_t start{SDBG::kNullID};
    uint64_t cur{SDBG::kNullID};
    uint64_t end{SDBG::kNullID};
    uint64_t rc_start{SDBG::kNullID};
    uint64_t rc_cur{SDBG::kNullID};
    uint64_t rc_end{SDBG::kNullID};
    uint32_t pending_next{kNoSimpleNeighbor};
    int64_t depth{0};
    uint32_t length{0};
  };

  stream_width = std::max(1u, stream_width);
  const size_t chunks_per_worker = 64;
  const size_t chunk_size = std::max<size_t>(
      64, DivCeiling(start_ids.size(),
                     static_cast<size_t>(std::max(1, num_threads)) *
                         chunks_per_worker));
  std::atomic<size_t> next_chunk{0};
  size_t count_palindrome = 0;

#pragma omp parallel num_threads(num_threads) reduction(+ : count_palindrome)
  {
    const int tid = omp_get_thread_num();
    std::vector<PathState> states(stream_width);
    size_t local_next = 0;
    size_t local_end = 0;
    bool input_exhausted = false;

    auto take_start = [&](uint64_t *start) {
      while (local_next == local_end) {
        const size_t begin =
            next_chunk.fetch_add(chunk_size, std::memory_order_relaxed);
        if (begin >= start_ids.size()) {
          return false;
        }
        local_next = begin;
        local_end = std::min(start_ids.size(), begin + chunk_size);
      }
      *start = start_ids[local_next++];
      return true;
    };

    auto emit = [&](PathState *state) {
      (*thread_vertices)[tid].emplace_back(
          state->start, state->end, state->rc_start, state->rc_end,
          state->depth, state->length);
      count_palindrome += state->start == state->rc_start;
      state->active = false;
    };

    for (;;) {
      size_t active_states = 0;
      for (PathState &state : states) {
        while (!state.active && !input_exhausted) {
          uint64_t start;
          if (!take_start(&start)) {
            input_exhausted = true;
            break;
          }
          if (!locks->try_lock(start)) {
            continue;
          }
          state.active = true;
          state.reverse_phase = false;
          state.start = start;
          state.cur = start;
          state.end = SDBG::kNullID;
          state.rc_start = SDBG::kNullID;
          state.rc_cur = SDBG::kNullID;
          state.rc_end = SDBG::kNullID;
          state.depth = sdbg->EdgeMultiplicity(start);
          state.length = 1;
        }
        active_states += state.active;
      }
      if (active_states == 0) {
        break;
      }

      // Materialize every independent next-link load in a separate stage.
      // Unlike a prefetch immediately followed by a dependent load, this
      // gives the machine a full stream-width window of real misses to overlap.
      for (PathState &state : states) {
        if (state.active) {
          const uint64_t edge =
              state.reverse_phase ? state.rc_cur : state.cur;
          state.pending_next = neighbors.LoadNext32(edge);
        }
      }

      // The decoded next IDs make the lock and multiplicity addresses known.
      // Start those dependent requests before the mutation/accumulation pass.
      for (PathState &state : states) {
        if (state.active && state.pending_next != kNoSimpleNeighbor) {
          locks->prefetch_for_write(state.pending_next);
          if (!state.reverse_phase) {
            sdbg->PrefetchMultiplicity(state.pending_next);
          }
        }
      }

      for (PathState &state : states) {
        if (!state.active) {
          continue;
        }
        if (!state.reverse_phase) {
          const uint64_t next =
              state.pending_next == kNoSimpleNeighbor
                  ? SDBG::kNullID
                  : static_cast<uint64_t>(state.pending_next);
          if (next != SDBG::kNullID) {
            state.cur = next;
            if (!locks->try_lock(next)) {
              state.active = false;
              continue;
            }
            state.depth += sdbg->EdgeMultiplicity(next);
            ++state.length;
            continue;
          }

          state.end = state.cur;
          state.rc_start = sdbg->EdgeReverseComplement(state.end);
          assert(state.rc_start != SDBG::kNullID);
          if (!locks->try_lock(state.rc_start)) {
            state.rc_end = sdbg->EdgeReverseComplement(state.start);
            if (std::max(state.end, state.start) <
                std::max(state.rc_start, state.rc_end)) {
              state.active = false;
            } else {
              emit(&state);
            }
            continue;
          }
          state.reverse_phase = true;
          state.rc_cur = state.rc_start;
          state.rc_end = state.rc_start;
          continue;
        }

        const uint64_t rc_next =
            state.pending_next == kNoSimpleNeighbor
                ? SDBG::kNullID
                : static_cast<uint64_t>(state.pending_next);
        if (rc_next == SDBG::kNullID) {
          emit(&state);
          continue;
        }
        state.rc_cur = rc_next;
        state.rc_end = rc_next;
        if (!locks->try_lock(rc_next)) {
          state.rc_end = sdbg->EdgeReverseComplement(state.start);
          emit(&state);
        }
      }
    }
  }
  return count_palindrome;
}

template <class NeighborAccessor>
size_t AssembleUncachedNonLoopPaths(
    SDBG *sdbg, AtomicBitVector *locks,
    std::vector<std::vector<UnitigGraphVertex>> *thread_vertices,
    const AtomicBitVector &has_simple_next,
    const NeighborAccessor &neighbors) {
  size_t count_palindrome = 0;
#pragma omp parallel for reduction(+ : count_palindrome)
  for (uint64_t edge_idx = 0; edge_idx < sdbg->size(); ++edge_idx) {
    if (sdbg->IsValidEdge(edge_idx) && !has_simple_next.at(edge_idx)) {
      count_palindrome += AssembleNonLoopPath(
          edge_idx, omp_get_thread_num(), sdbg, locks, thread_vertices,
          neighbors);
    }
  }
  return count_palindrome;
}

}  // namespace

UnitigGraph::UnitigGraph(SDBG *sdbg)
    : sdbg_(sdbg), adapter_impl_(this), sudo_adapter_impl_(this) {
  id_map_.clear();
  dense_id_map_.clear();
  vertices_.clear();
  active_ids_.clear();
  AtomicBitVector locks(sdbg_->size());
  AtomicBitVector has_simple_next;
  AtomicBitVector has_simple_prev;
  std::unique_ptr<uint32_t[]> owned_simple_next;
  uint32_t *simple_next = nullptr;
  bool reusing_forward_topology = false;
  std::unique_ptr<uint32_t[]> simple_prev;
  std::unique_ptr<uint8_t[]> compact_simple_neighbors;
  std::unique_ptr<uint8_t[]> simple_offsets;
  std::unique_ptr<uint8_t[]> block_simple_codes;
  std::unique_ptr<uint32_t[]> block_simple_bases;
  std::unique_ptr<uint32_t[]> simple_link;
  const uint64_t num_edges = sdbg_->size();
  const bool compact_edge_ids =
      num_edges < static_cast<uint64_t>(kNoSimpleNeighbor);
  const bool use_compact_simple_neighbors =
      num_edges < static_cast<uint64_t>(kNoCompactSimpleNeighbor) &&
      std::getenv("MEGAHIT_ENABLE_COMPACT_SIMPLE_NEIGHBORS") != nullptr;
  const uint64_t bitvector_word_bits = AtomicBitVector::bits_per_word();
  const uint64_t bitvector_bytes =
      (num_edges / bitvector_word_bits +
       static_cast<uint64_t>(num_edges % bitvector_word_bits != 0)) *
      sizeof(AtomicBitVector::word_type);
  const uint64_t double_neighbor_bytes =
      compact_edge_ids
          ? num_edges * (use_compact_simple_neighbors
                             ? kCompactSimpleNeighborBytesPerEdge
                             : 2 * sizeof(uint32_t))
                       : std::numeric_limits<uint64_t>::max();
  const uint64_t xor_neighbor_bytes =
      compact_edge_ids ? num_edges * sizeof(uint32_t) + bitvector_bytes
                       : std::numeric_limits<uint64_t>::max();
  const uint64_t forward_neighbor_bytes = xor_neighbor_bytes;
  const uint64_t offset_neighbor_bytes =
      (num_edges + 1u) / 2u + bitvector_bytes;
  const uint64_t num_simple_blocks =
      (num_edges + kBlockSimpleNeighborEdges - 1u) /
      kBlockSimpleNeighborEdges;
  const uint64_t block_neighbor_bytes =
      num_edges + num_simple_blocks * (kAlphabetSize + 1u) *
                      sizeof(uint32_t) +
      bitvector_bytes;

  enum class SimpleNeighborMode {
    kOff,
    kOffset,
    kBlock,
    kXor,
    kForward,
    kDouble
  };
  SimpleNeighborMode simple_neighbor_mode = SimpleNeighborMode::kOff;
  bool auto_select_mode = true;
  const char *requested_mode =
      std::getenv("MEGAHIT_SIMPLE_NEIGHBOR_MODE");
  if (std::getenv("MEGAHIT_DISABLE_SIMPLE_NEIGHBOR_CACHE") != nullptr) {
    auto_select_mode = false;
  } else if (requested_mode != nullptr) {
    auto_select_mode = false;
    if (std::strcmp(requested_mode, "double") == 0) {
      simple_neighbor_mode = SimpleNeighborMode::kDouble;
    } else if (std::strcmp(requested_mode, "forward") == 0) {
      simple_neighbor_mode = SimpleNeighborMode::kForward;
    } else if (std::strcmp(requested_mode, "offset") == 0) {
      simple_neighbor_mode = SimpleNeighborMode::kOffset;
    } else if (std::strcmp(requested_mode, "block") == 0) {
      simple_neighbor_mode = SimpleNeighborMode::kBlock;
    } else if (std::strcmp(requested_mode, "xor") == 0) {
      simple_neighbor_mode = SimpleNeighborMode::kXor;
    } else if (std::strcmp(requested_mode, "off") != 0) {
      xwarn("Unknown MEGAHIT_SIMPLE_NEIGHBOR_MODE '{s}'; using auto mode\n",
            requested_mode);
      auto_select_mode = true;
    }
  }

  if (auto_select_mode) {
    uint64_t available_bytes = 0;
    bool found_available_memory = false;
#ifdef __linux__
    // MemAvailable includes reclaimable page cache and is therefore stable
    // across cold and warm input-cache states.  Older kernels without it fall
    // back to the former sysinfo estimate.  A finite cgroup limit always
    // tightens the host-wide value afterwards.
    found_available_memory = ReadMemAvailable(&available_bytes);
    struct sysinfo memory_info;
    if (!found_available_memory && ::sysinfo(&memory_info) == 0) {
      available_bytes =
          (static_cast<uint64_t>(memory_info.freeram) +
           static_cast<uint64_t>(memory_info.bufferram)) *
          static_cast<uint64_t>(memory_info.mem_unit);
      found_available_memory = true;
    }
    TightenAvailableMemoryFromCgroup(&available_bytes,
                                     &found_available_memory);
#endif
    // --host_mem is the assembly-wide contract, while MemAvailable/cgroup is
    // the amount the process can still obtain right now.  Honor the tighter
    // constraint.  This keeps the choice portable across machines instead of
    // keying it to a particular socket count or physical-RAM size.
    const uint64_t host_budget = sdbg_->HostMemoryBudget();
    if (host_budget != 0 &&
        (!found_available_memory || host_budget < available_bytes)) {
      available_bytes = host_budget;
      found_available_memory = true;
    }

    // A caller that explicitly materialized the diagnostic wide topology
    // cache has already paid for a uint32 Forward array; reuse that storage in
    // place.  The normal/default path never allocates a new four-byte-per-edge
    // table merely because a large-memory host happens to be available.
    if (sdbg_->HasSimplePathSnapshot()) {
      simple_neighbor_mode = SimpleNeighborMode::kBlock;
    } else if (sdbg_->HasReusableForwardTopologyCache()) {
      simple_neighbor_mode = SimpleNeighborMode::kForward;
    } else if (compact_edge_ids && found_available_memory &&
               block_neighbor_bytes <= available_bytes / 4u) {
      // One byte per edge plus a small block directory removes the dependent
      // rank/select walk from unitig traversal and is the default balanced
      // representation.
      simple_neighbor_mode = SimpleNeighborMode::kBlock;
    } else if (compact_edge_ids &&
               (!found_available_memory ||
                offset_neighbor_bytes <= available_bytes / 4u)) {
      // Half a byte per edge is the bounded-memory fallback.  It preserves the
      // exact graph topology and pays an extra Forward lookup only while
      // walking a chain.
      simple_neighbor_mode = SimpleNeighborMode::kOffset;
    }
  }

  if (!compact_edge_ids && simple_neighbor_mode != SimpleNeighborMode::kOff) {
    xwarn("Simple-neighbor cache requires fewer than {} edges; using off mode\n",
          kNoSimpleNeighbor);
    simple_neighbor_mode = SimpleNeighborMode::kOff;
  }

  const bool reusing_simple_path_snapshot =
      simple_neighbor_mode == SimpleNeighborMode::kBlock &&
      sdbg_->HasSimplePathSnapshot();

  if (simple_neighbor_mode == SimpleNeighborMode::kOffset) {
    try {
      simple_offsets.reset(new uint8_t[(num_edges + 1u) / 2u]);
      has_simple_prev = AtomicBitVector(num_edges);
    } catch (const std::bad_alloc &) {
      simple_offsets.reset();
      has_simple_prev = AtomicBitVector();
      if (!auto_select_mode) {
        throw;
      }
      simple_neighbor_mode = SimpleNeighborMode::kOff;
      xwarn("Offset simple-neighbor allocation failed; using off mode\n");
    }
  }
  if (simple_neighbor_mode == SimpleNeighborMode::kBlock) {
    try {
      if (!reusing_simple_path_snapshot) {
        block_simple_codes.reset(new uint8_t[num_edges]);
        block_simple_bases.reset(
            new uint32_t[num_simple_blocks * (kAlphabetSize + 1u)]);
      }
      has_simple_prev = AtomicBitVector(num_edges);
    } catch (const std::bad_alloc &) {
      block_simple_codes.reset();
      block_simple_bases.reset();
      has_simple_prev = AtomicBitVector();
      if (!auto_select_mode) {
        throw;
      }
      simple_neighbor_mode = SimpleNeighborMode::kOff;
      xwarn("Block simple-neighbor allocation failed; using off mode\n");
    }
  }

  // Auto mode is opportunistic: if the available-memory estimate becomes
  // stale between probing and allocation, retain the exact off-mode path
  // instead of turning an optional cache into a fatal allocation failure.
  // Explicitly requested modes keep their fail-fast semantics.
  if (auto_select_mode && simple_neighbor_mode == SimpleNeighborMode::kDouble) {
    try {
      if (use_compact_simple_neighbors) {
        compact_simple_neighbors.reset(
            new uint8_t[num_edges * kCompactSimpleNeighborBytesPerEdge]);
      } else {
        owned_simple_next.reset(new uint32_t[num_edges]);
        simple_next = owned_simple_next.get();
        simple_prev.reset(new uint32_t[num_edges]);
      }
    } catch (const std::bad_alloc &) {
      compact_simple_neighbors.reset();
      owned_simple_next.reset();
      simple_next = nullptr;
      simple_prev.reset();
      simple_neighbor_mode = SimpleNeighborMode::kXor;
      xwarn("Double simple-neighbor allocation failed; trying xor mode\n");
    }
  }
  if (auto_select_mode && simple_neighbor_mode == SimpleNeighborMode::kXor) {
    try {
      simple_link.reset(new uint32_t[num_edges]);
      has_simple_prev = AtomicBitVector(num_edges);
    } catch (const std::bad_alloc &) {
      simple_link.reset();
      has_simple_prev = AtomicBitVector();
      simple_neighbor_mode = SimpleNeighborMode::kOff;
      xwarn("Xor simple-neighbor allocation failed; using off mode\n");
    }
  }

  if (auto_select_mode &&
      simple_neighbor_mode == SimpleNeighborMode::kForward) {
    simple_next = sdbg_->BorrowForwardCacheForSimpleNeighbors();
    reusing_forward_topology = simple_next != nullptr;
    if (!simple_next) {
      try {
        owned_simple_next.reset(new uint32_t[num_edges]);
        simple_next = owned_simple_next.get();
      } catch (const std::bad_alloc &) {
        owned_simple_next.reset();
        simple_next = nullptr;
        simple_neighbor_mode = SimpleNeighborMode::kOff;
        xwarn("Forward simple-neighbor allocation failed; using off mode\n");
      }
    }
    if (simple_neighbor_mode == SimpleNeighborMode::kForward) {
      try {
        has_simple_prev = AtomicBitVector(num_edges);
      } catch (const std::bad_alloc &) {
        if (reusing_forward_topology) {
          sdbg_->RestoreBorrowedForwardTopologyCache();
          reusing_forward_topology = false;
        }
        owned_simple_next.reset();
        simple_next = nullptr;
        has_simple_prev = AtomicBitVector();
        simple_neighbor_mode = SimpleNeighborMode::kOff;
        xwarn("Forward simple-neighbor allocation failed; using off mode\n");
      }
    }
  }
  if (simple_neighbor_mode == SimpleNeighborMode::kForward && !simple_next) {
    simple_next = sdbg_->BorrowForwardCacheForSimpleNeighbors();
    reusing_forward_topology = simple_next != nullptr;
    if (!simple_next) {
      owned_simple_next.reset(new uint32_t[num_edges]);
      simple_next = owned_simple_next.get();
    }
    has_simple_prev = AtomicBitVector(num_edges);
  }

  const char *simple_neighbor_mode_name =
      simple_neighbor_mode == SimpleNeighborMode::kDouble
          ? (use_compact_simple_neighbors ? "compact-double" : "double")
          : (simple_neighbor_mode == SimpleNeighborMode::kForward
                 ? (reusing_forward_topology ? "forward-reused" : "forward")
                 : (simple_neighbor_mode == SimpleNeighborMode::kOffset
                        ? "offset"
                 : (simple_neighbor_mode == SimpleNeighborMode::kBlock
                        ? "block"
                 : (simple_neighbor_mode == SimpleNeighborMode::kXor
                        ? "xor"
                        : "off"))));
  const uint64_t simple_neighbor_temporary_bytes =
      simple_neighbor_mode == SimpleNeighborMode::kDouble
          ? double_neighbor_bytes
          : (simple_neighbor_mode == SimpleNeighborMode::kForward
                 ? (reusing_forward_topology ? bitvector_bytes
                                             : forward_neighbor_bytes)
                 : (simple_neighbor_mode == SimpleNeighborMode::kOffset
                        ? offset_neighbor_bytes
                 : (simple_neighbor_mode == SimpleNeighborMode::kBlock
                        ? (reusing_simple_path_snapshot ? bitvector_bytes
                                                        : block_neighbor_bytes)
                 : (simple_neighbor_mode == SimpleNeighborMode::kXor
                        ? xor_neighbor_bytes
                        : bitvector_bytes))));
  xinfo("Simple-neighbor mode: {s}; additional temporary bytes: {}\n",
        simple_neighbor_mode_name, simple_neighbor_temporary_bytes);

  size_t count_palindrome = 0;
  const int num_threads = omp_get_max_threads();
  std::vector<std::vector<UnitigGraphVertex>> thread_vertices(num_threads);
  std::vector<uint64_t> terminal_ids;
  std::vector<size_t> terminal_thread_offsets;
  std::vector<std::vector<uint64_t>> thread_terminals;
  if (simple_neighbor_mode != SimpleNeighborMode::kOff) {
    thread_terminals.resize(num_threads);
  }
  SimpleTimer unitig_stage_timer;
  unitig_stage_timer.start();

  // Both cache modes materialize NextSimplePathEdge once and preserve the
  // former full scan's ascending terminal order.  The double mode is fastest;
  // the XOR-linked mode trades one extra chain conversion for half the dense
  // storage.  Arrays are allocated without value-initialization so their
  // parallel fills first-touch pages on the worker owning each edge range.
  if (simple_neighbor_mode == SimpleNeighborMode::kForward) {
    std::vector<std::vector<uint64_t>> thread_starts(num_threads);

#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      auto &local_starts = thread_starts[tid];
      local_starts.reserve(
          static_cast<size_t>(num_edges / std::max(1, num_threads) / 64 + 1));

#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        uint32_t next_id = kNoSimpleNeighbor;
        if (sdbg_->IsValidEdge(edge_idx)) {
          const uint64_t next =
              reusing_forward_topology
                  ? sdbg_->ConvertBorrowedForwardEntryToSimpleNeighbor(edge_idx)
                  : sdbg_->NextSimplePathEdge(edge_idx);
          if (next != SDBG::kNullID) {
            assert(next < static_cast<uint64_t>(kNoSimpleNeighbor));
            next_id = static_cast<uint32_t>(next);
            has_simple_prev.set(next);
          }
        }
        simple_next[edge_idx] = next_id;
      }

      // The implicit barrier above makes the predecessor bitmap complete.
#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) &&
            !has_simple_prev.at(edge_idx)) {
          local_starts.push_back(edge_idx);
        }
      }
    }
    thread_terminals.swap(thread_starts);
  } else if (simple_neighbor_mode == SimpleNeighborMode::kOffset) {
    std::vector<std::vector<uint64_t>> thread_starts(num_threads);
    const uint64_t num_code_bytes = (num_edges + 1u) / 2u;

#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      auto &local_starts = thread_starts[tid];
      local_starts.reserve(
          static_cast<size_t>(num_edges / std::max(1, num_threads) / 64 + 1));

#pragma omp for schedule(static)
      for (uint64_t byte_idx = 0; byte_idx < num_code_bytes; ++byte_idx) {
        uint8_t packed = 0;
        const uint64_t first_edge = byte_idx * 2u;
        for (unsigned lane = 0; lane < 2u; ++lane) {
          const uint64_t edge_idx = first_edge + lane;
          if (edge_idx >= num_edges || !sdbg_->IsValidEdge(edge_idx)) {
            continue;
          }
          uint8_t offset = 0;
          const uint64_t next =
              sdbg_->NextSimplePathEdgeWithOffset(edge_idx, &offset);
          if (next != SDBG::kNullID) {
            const uint8_t code = static_cast<uint8_t>(offset + 1u);
            packed |= static_cast<uint8_t>(code << (lane * 4u));
            has_simple_prev.set(next);
          }
        }
        simple_offsets[byte_idx] = packed;
      }

#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) &&
            !has_simple_prev.at(edge_idx)) {
          local_starts.push_back(edge_idx);
        }
      }
    }
    thread_terminals.swap(thread_starts);
  } else if (simple_neighbor_mode == SimpleNeighborMode::kBlock) {
    std::vector<std::vector<uint64_t>> thread_starts(num_threads);
    uint64_t overflow_count = 0;

    if (reusing_simple_path_snapshot) {
#pragma omp parallel num_threads(num_threads)
      {
        const int tid = omp_get_thread_num();
        auto &local_terminals = thread_terminals[tid];
        auto &local_starts = thread_starts[tid];
        local_terminals.reserve(static_cast<size_t>(
            num_edges / std::max(1, num_threads) / 64 + 1));
        local_starts.reserve(local_terminals.capacity());

#pragma omp for schedule(static)
        for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
          if (!sdbg_->IsValidEdge(edge_idx)) continue;
          const uint64_t next =
              sdbg_->CachedOrLiveNextSimplePathEdge(edge_idx);
          if (next == SDBG::kNullID) {
            local_terminals.push_back(edge_idx);
          } else {
            has_simple_prev.set(next);
          }
        }

#pragma omp for schedule(static)
        for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
          if (sdbg_->IsValidEdge(edge_idx) &&
              !has_simple_prev.at(edge_idx)) {
            local_starts.push_back(edge_idx);
          }
        }
      }
      xinfo("Reused pre-pruning simple-path snapshot; {} original overflow "
            "edges use live topology\n",
            sdbg_->SimplePathSnapshotOverflowCount());
    } else {
#pragma omp parallel reduction(+ : overflow_count) num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      auto &local_starts = thread_starts[tid];
      local_starts.reserve(
          static_cast<size_t>(num_edges / std::max(1, num_threads) / 64 + 1));
      std::array<uint32_t, kBlockSimpleNeighborEdges> local_next;

#pragma omp for schedule(static)
      for (uint64_t block_id = 0; block_id < num_simple_blocks; ++block_id) {
        uint32_t *block_bases =
            block_simple_bases.get() + block_id * (kAlphabetSize + 1u);
        std::fill(block_bases, block_bases + kAlphabetSize + 1,
                  kNoSimpleNeighbor);
        const uint64_t edge_begin = block_id * kBlockSimpleNeighborEdges;
        const uint64_t edge_end = std::min<uint64_t>(
            num_edges, edge_begin + kBlockSimpleNeighborEdges);

        for (uint64_t edge_idx = edge_begin; edge_idx < edge_end;
             ++edge_idx) {
          const uint64_t lane = edge_idx - edge_begin;
          uint32_t next_id = kNoSimpleNeighbor;
          if (sdbg_->IsValidEdge(edge_idx)) {
            const uint64_t next = sdbg_->NextSimplePathEdge(edge_idx);
            if (next != SDBG::kNullID) {
              next_id = static_cast<uint32_t>(next);
              uint8_t c = sdbg_->GetW(edge_idx);
              if (c > kAlphabetSize) {
                c -= kAlphabetSize;
              }
              block_bases[c] = std::min(block_bases[c], next_id);
              has_simple_prev.set(next);
            }
          }
          local_next[lane] = next_id;
        }

        for (uint64_t edge_idx = edge_begin; edge_idx < edge_end;
             ++edge_idx) {
          const uint32_t next_id = local_next[edge_idx - edge_begin];
          if (next_id == kNoSimpleNeighbor) {
            block_simple_codes[edge_idx] = 0;
            continue;
          }
          uint8_t c = sdbg_->GetW(edge_idx);
          if (c > kAlphabetSize) {
            c -= kAlphabetSize;
          }
          const uint32_t delta = next_id - block_bases[c];
          if (delta < 254u) {
            block_simple_codes[edge_idx] =
                static_cast<uint8_t>(delta + 1u);
          } else {
            block_simple_codes[edge_idx] = 255u;
            ++overflow_count;
          }
        }
      }

#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) &&
            !has_simple_prev.at(edge_idx)) {
          local_starts.push_back(edge_idx);
        }
      }
    }
      xinfo("Block simple-neighbor overflow edges: {} / {}\n",
            overflow_count, num_edges);
    }
    thread_terminals.swap(thread_starts);
  } else if (simple_neighbor_mode == SimpleNeighborMode::kDouble) {
    if (!compact_simple_neighbors && !simple_next) {
      if (use_compact_simple_neighbors) {
        compact_simple_neighbors.reset(
            new uint8_t[num_edges * kCompactSimpleNeighborBytesPerEdge]);
      } else {
        owned_simple_next.reset(new uint32_t[num_edges]);
        simple_next = owned_simple_next.get();
        simple_prev.reset(new uint32_t[num_edges]);
      }
    }

    if (compact_simple_neighbors) {
      const CompactDoubleNeighborAccessor neighbors{
          compact_simple_neighbors.get()};
      BuildDoubleSimpleNeighbors(sdbg_, num_edges, num_threads, neighbors,
                                 &thread_terminals);
    } else {
      const WideDoubleNeighborAccessor neighbors{
          simple_next, simple_prev.get()};
      BuildDoubleSimpleNeighbors(sdbg_, num_edges, num_threads, neighbors,
                                 &thread_terminals);
    }
  } else if (simple_neighbor_mode == SimpleNeighborMode::kXor) {
    if (!simple_link) {
      simple_link.reset(new uint32_t[num_edges]);
      has_simple_prev = AtomicBitVector(num_edges);
    }
    std::vector<std::vector<uint64_t>> thread_starts(num_threads);

#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      auto &local_terminals = thread_terminals[tid];
      auto &local_starts = thread_starts[tid];
      local_terminals.reserve(
          static_cast<size_t>(num_edges / std::max(1, num_threads) / 64 + 1));
      local_starts.reserve(local_terminals.capacity());

#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        simple_link[edge_idx] = kNoSimpleNeighbor;
        if (!sdbg_->IsValidEdge(edge_idx)) {
          continue;
        }

        const uint64_t next = sdbg_->NextSimplePathEdge(edge_idx);
        if (next == SDBG::kNullID) {
          local_terminals.push_back(edge_idx);
        } else {
          assert(next < static_cast<uint64_t>(kNoSimpleNeighbor));
          simple_link[edge_idx] = static_cast<uint32_t>(next);
          has_simple_prev.set(next);
        }
      }

      // The implicit barrier above makes the incoming-edge bitmap complete.
      // A valid edge without a simple predecessor starts one non-loop chain.
#pragma omp for schedule(static)
      for (uint64_t edge_idx = 0; edge_idx < num_edges; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) && !has_simple_prev.at(edge_idx)) {
          local_starts.push_back(edge_idx);
        }
      }

      // Distinct starts lead through disjoint chains: NextSimplePathEdge is
      // injective because every target has a unique incoming edge.  Therefore
      // each link has one non-atomic writer.  Pure cycles have no start and are
      // intentionally left in next-ID form; the loop phase uses the SDBG.
      for (uint64_t start : local_starts) {
        uint32_t prev = kNoSimpleNeighbor;
        uint32_t cur = static_cast<uint32_t>(start);
        while (true) {
          const uint32_t next = simple_link[cur];
          simple_link[cur] = prev ^ next;
          if (next == kNoSimpleNeighbor) {
            break;
          }
          prev = cur;
          cur = next;
        }
      }
    }
  } else {
    has_simple_next = AtomicBitVector(num_edges);
#pragma omp parallel for schedule(static)
    for (uint64_t word_idx = 0; word_idx < has_simple_next.word_count();
         ++word_idx) {
      unsigned long word = 0;
      const uint64_t edge_begin = word_idx * has_simple_next.bits_per_word();
      const uint64_t edge_end = std::min<uint64_t>(
          edge_begin + has_simple_next.bits_per_word(), num_edges);
      for (uint64_t edge_idx = edge_begin; edge_idx < edge_end; ++edge_idx) {
        if (sdbg_->IsValidEdge(edge_idx) &&
            sdbg_->NextSimplePathEdge(edge_idx) != SDBG::kNullID) {
          word |= 1ul << (edge_idx - edge_begin);
        }
      }
      has_simple_next.store_word(word_idx, word);
    }
  }

  if (simple_neighbor_mode != SimpleNeighborMode::kOff) {
    // Static scheduling gives every thread a contiguous edge-ID interval.
    // Concatenating those intervals in thread order preserves global order.
    terminal_thread_offsets.resize(num_threads + 1, 0);
    for (int tid = 0; tid < num_threads; ++tid) {
      terminal_thread_offsets[tid + 1] =
          terminal_thread_offsets[tid] + thread_terminals[tid].size();
    }
    terminal_ids.resize(terminal_thread_offsets.back());
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int tid = 0; tid < num_threads; ++tid) {
      std::copy(thread_terminals[tid].begin(), thread_terminals[tid].end(),
                terminal_ids.begin() + terminal_thread_offsets[tid]);
    }
    std::vector<std::vector<uint64_t>>().swap(thread_terminals);
  }
  unitig_stage_timer.stop();
  xinfo("Unitig neighbor discovery time: {.4}\n",
        unitig_stage_timer.elapsed());
  unitig_stage_timer.reset();
  unitig_stage_timer.start();

  // Assemble non-loop paths.  Cached modes process only the collected
  // terminal IDs while retaining each static edge range's original owner.
  // Dispatch once here so inner chain walks are fully specialized for their
  // representation and contain no per-edge mode branch.
  if (simple_neighbor_mode == SimpleNeighborMode::kForward) {
    const ForwardNeighborAccessor neighbors{simple_next};
    if (std::getenv("MEGAHIT_DISABLE_MULTISTREAM_UNITIG") != nullptr) {
      count_palindrome = AssembleForwardNonLoopPaths(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads,
          neighbors);
    } else {
      // Pointer-chasing paths expose one outstanding dependent read each.
      // Keep a modest architecture-neutral pool of independent paths in
      // flight; the state is only a few KiB per worker and the width remains
      // bounded when millions of starts are available.  The environment
      // override remains useful for unusual memory systems.
      unsigned stream_width = static_cast<unsigned>(std::min<size_t>(
          32u, std::max<size_t>(1u, DivCeiling(
                                    terminal_ids.size(),
                                    static_cast<size_t>(num_threads)))));
      if (const char *value = std::getenv("MEGAHIT_UNITIG_MLP_WIDTH")) {
        char *end = nullptr;
        const unsigned long parsed = std::strtoul(value, &end, 10);
        if (end != value && *end == '\0' && parsed >= 1 && parsed <= 64) {
          stream_width = static_cast<unsigned>(parsed);
        }
      }
      xinfo("Unitig chain traversal: {} interleaved paths per worker\n",
            stream_width);
      count_palindrome = AssembleForwardNonLoopPathsMultistream(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads, neighbors,
          stream_width);
    }
  } else if (simple_neighbor_mode == SimpleNeighborMode::kOffset) {
    const OffsetNeighborAccessor neighbors{sdbg_, simple_offsets.get()};
    unsigned stream_width = static_cast<unsigned>(std::min<size_t>(
        32u, std::max<size_t>(1u, DivCeiling(
                                  terminal_ids.size(),
                                  static_cast<size_t>(num_threads)))));
    if (const char *value = std::getenv("MEGAHIT_UNITIG_MLP_WIDTH")) {
      char *end = nullptr;
      const unsigned long parsed = std::strtoul(value, &end, 10);
      if (end != value && *end == '\0' && parsed >= 1 && parsed <= 64) {
        stream_width = static_cast<unsigned>(parsed);
      }
    }
    xinfo("Unitig chain traversal: {} interleaved paths per worker\n",
          stream_width);
    count_palindrome = AssembleForwardNonLoopPathsMultistream(
        sdbg_, &locks, &thread_vertices, terminal_ids, num_threads, neighbors,
        stream_width);
  } else if (simple_neighbor_mode == SimpleNeighborMode::kBlock) {
    unsigned stream_width = static_cast<unsigned>(std::min<size_t>(
        32u, std::max<size_t>(1u, DivCeiling(
                                  terminal_ids.size(),
                                  static_cast<size_t>(num_threads)))));
    xinfo("Unitig chain traversal: {} interleaved paths per worker\n",
          stream_width);
    if (reusing_simple_path_snapshot) {
      const SharedSnapshotNeighborAccessor neighbors{sdbg_};
      count_palindrome = AssembleForwardNonLoopPathsMultistream(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads,
          neighbors, stream_width);
    } else {
      const BlockNeighborAccessor neighbors{
          sdbg_, block_simple_codes.get(), block_simple_bases.get()};
      count_palindrome = AssembleForwardNonLoopPathsMultistream(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads,
          neighbors, stream_width);
    }
  } else if (simple_neighbor_mode == SimpleNeighborMode::kDouble) {
    if (compact_simple_neighbors) {
      const CompactDoubleNeighborAccessor neighbors{
          compact_simple_neighbors.get()};
      count_palindrome = AssembleCachedNonLoopPaths(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads,
          neighbors);
    } else {
      const WideDoubleNeighborAccessor neighbors{
          simple_next, simple_prev.get()};
      count_palindrome = AssembleCachedNonLoopPaths(
          sdbg_, &locks, &thread_vertices, terminal_ids, num_threads,
          neighbors);
    }
  } else if (simple_neighbor_mode == SimpleNeighborMode::kXor) {
    const XorNeighborAccessor neighbors{simple_link.get()};
    count_palindrome = AssembleCachedNonLoopPaths(
        sdbg_, &locks, &thread_vertices, terminal_ids, num_threads, neighbors);
  } else {
    const UncachedNeighborAccessor neighbors{sdbg_, &has_simple_next};
    count_palindrome = AssembleUncachedNonLoopPaths(
        sdbg_, &locks, &thread_vertices, has_simple_next, neighbors);
  }

  size_t num_non_loop_vertices = 0;
  for (const auto &local : thread_vertices) {
    num_non_loop_vertices += local.size();
  }
  vertices_.reserve(num_non_loop_vertices);
  for (auto &local : thread_vertices) {
    vertices_.insert(vertices_.end(), std::make_move_iterator(local.begin()),
                     std::make_move_iterator(local.end()));
    std::vector<UnitigGraphVertex>().swap(local);
  }
  unitig_stage_timer.stop();
  xinfo("Unitig non-loop chain assembly time: {.4}\n",
        unitig_stage_timer.elapsed());
  if (simple_neighbor_mode == SimpleNeighborMode::kForward &&
      num_threads == 1) {
    // Forward-only links discover chains by their starts.  Restore the
    // historical terminal-edge order so single-thread output remains byte
    // identical.  Multi-thread output ordering was never stable, so avoid an
    // O(V log V) serial sort on the performance-critical parallel path.
    std::sort(vertices_.begin(), vertices_.end(),
              [](const UnitigGraphVertex &lhs,
                 const UnitigGraphVertex &rhs) {
                return lhs.TerminalOrderKey() < rhs.TerminalOrderKey();
              });
  }
  xinfo("Graph size without loops: {}, palindrome: {}\n", vertices_.size(),
        count_palindrome);

  // The simple-neighbor snapshot is dead before loop discovery.  In the
  // common wide-cache path its storage is the temporarily repurposed Forward
  // table, so restore that table in place now.  This keeps Forward, Backward
  // and graph-cleaning performance while eliminating their edge-sized peak
  // overlap with a separate simple-next allocation.
  if (reusing_forward_topology) {
    sdbg_->RestoreBorrowedForwardTopologyCache();
    reusing_forward_topology = false;
  }
  owned_simple_next.reset();
  simple_next = nullptr;
  simple_prev.reset();
  compact_simple_neighbors.reset();
  simple_offsets.reset();
  if (simple_neighbor_mode == SimpleNeighborMode::kBlock &&
      !reusing_simple_path_snapshot &&
      std::getenv("MEGAHIT_DISABLE_RETAINED_SIMPLE_NEXT") == nullptr) {
    materialization_simple_codes_ = std::move(block_simple_codes);
    materialization_simple_bases_ = std::move(block_simple_bases);
    materialization_simple_edge_count_ = num_edges;
    const uint64_t retained_bytes =
        num_edges + num_simple_blocks * (kAlphabetSize + 1u) *
                        sizeof(uint32_t);
    xinfo("Retained simple-successor snapshot: {.2} MiB\n",
          static_cast<double>(retained_bytes) / (1u << 20u));
  } else {
    block_simple_codes.reset();
    block_simple_bases.reset();
    if (reusing_simple_path_snapshot) {
      xinfo("Retained shared simple-successor snapshot: {.2} MiB\n",
            static_cast<double>(sdbg_->SimplePathSnapshotBytes()) /
                (1u << 20u));
    }
  }
  simple_link.reset();
  has_simple_prev = AtomicBitVector();
  has_simple_next = AtomicBitVector();
  std::vector<uint64_t>().swap(terminal_ids);
  std::vector<size_t>().swap(terminal_thread_offsets);

  // assemble looped paths
  unitig_stage_timer.reset();
  unitig_stage_timer.start();
  // The historical loop walk is serialized by a global mutex, but the
  // winning seed used to depend on which OpenMP worker reached that mutex
  // first.  Besides rotating the emitted circular sequence, that also changes
  // its historical seed-weighted depth (the seed edge is intentionally
  // counted exactly as in the legacy implementation).  Scan candidates in
  // parallel, then claim components in ascending edge order: this is the
  // legacy one-thread order without turning the multi-billion-edge scan into
  // a serial pass.
  const int loop_scan_threads = std::max(1, omp_get_max_threads());
  std::vector<std::vector<uint64_t>> local_loop_candidates(loop_scan_threads);
#pragma omp parallel num_threads(loop_scan_threads)
  {
    auto &local = local_loop_candidates[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t edge_idx = 0; edge_idx < sdbg_->size(); ++edge_idx) {
      if (!locks.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
        local.push_back(edge_idx);
      }
    }
  }
  size_t num_loop_candidates = 0;
  for (const auto &local : local_loop_candidates) {
    num_loop_candidates += local.size();
  }
  std::vector<uint64_t> loop_seed_candidates;
  loop_seed_candidates.reserve(num_loop_candidates);
  for (auto &local : local_loop_candidates) {
    loop_seed_candidates.insert(loop_seed_candidates.end(), local.begin(),
                                local.end());
    std::vector<uint64_t>().swap(local);
  }
  std::sort(loop_seed_candidates.begin(), loop_seed_candidates.end());

  size_t count_loop = 0;
  for (const uint64_t edge_idx : loop_seed_candidates) {
    if (!locks.at(edge_idx) && sdbg_->IsValidEdge(edge_idx)) {
      if (!locks.at(edge_idx)) {
        uint64_t cur_edge = edge_idx;
        uint64_t rc_edge = sdbg_->EdgeReverseComplement(edge_idx);
        uint64_t depth = sdbg_->EdgeMultiplicity(edge_idx);
        uint32_t length = 0;
        // whether it is marked before entering the loop
        bool rc_marked = locks.at(rc_edge);

        while (!locks.at(cur_edge)) {
          locks.set(cur_edge);
          depth += sdbg_->EdgeMultiplicity(cur_edge);
          ++length;
          // Keep the loop assembly path on the original SDBG queries.  Loop
          // seed selection is intentionally unchanged in this low-risk pass;
          // the dense snapshot only accelerates non-loop construction.
          cur_edge = sdbg_->PrevSimplePathEdge(cur_edge);
          assert(cur_edge != SDBG::kNullID);
        }
        assert(cur_edge == edge_idx);

        if (!rc_marked) {
          uint64_t start = sdbg_->NextSimplePathEdge(edge_idx);
          uint64_t end = edge_idx;
          vertices_.emplace_back(start, end, sdbg_->EdgeReverseComplement(end),
                                 sdbg_->EdgeReverseComplement(start), depth,
                                 length, true);
          count_loop += 1;
        }
      }
    }
  }

  if (vertices_.size() >= kMaxNumVertices) {
    xfatal(
        "Too many vertices in the unitig graph ({} >= {}), "
        "you may increase the kmer size to remove tons of erroneous kmers.\n",
        vertices_.size(), kMaxNumVertices);
  }
  unitig_stage_timer.stop();
  xinfo("Unitig loop discovery time: {.4}\n", unitig_stage_timer.elapsed());

  unitig_stage_timer.reset();
  unitig_stage_timer.start();
  const bool profile_unitig_build =
      std::getenv("MEGAHIT_PROFILE_UNITIG_BUILD") != nullptr;
  double endpoint_stage_begin = omp_get_wtime();
  double endpoint_release_time = 0;
  double endpoint_allocate_time = 0;
  double endpoint_insert_time = 0;
  sdbg_->FreeMultiplicity();
  endpoint_release_time = omp_get_wtime() - endpoint_stage_begin;
  const uint64_t num_endpoints = vertices_.size() * 2 - count_palindrome;

  // RabbitTClust's dense-ID pattern is substantially cheaper than a hash
  // probe on every graph traversal.  Use it when its storage is no larger
  // than a conservative estimate for the flat hash table; retain the sparse
  // map for long-unitig graphs where a dense edge-sized array would waste
  // memory.
  constexpr uint64_t kMaxDenseEdgesPerEndpoint = 4;
  use_dense_id_map_ =
      num_endpoints != 0 &&
      sdbg_->size() <= num_endpoints * kMaxDenseEdgesPerEndpoint;
  if (use_dense_id_map_) {
    dense_id_map_.assign(sdbg_->size(), static_cast<size_type>(-1));
    xinfo("Using dense edge-to-unitig map: {} bytes\n",
          dense_id_map_.size() * sizeof(dense_id_map_[0]));
  } else {
    id_map_.reserve(num_endpoints, num_threads, sdbg_->size());
    xinfo("Using sharded sparse edge-to-unitig map: {} endpoints, {}-bit "
          "keys, {} shards, batch {}\n",
          num_endpoints, id_map_.key_bits(), id_map_.num_shards(),
          id_map_.batch_size());
  }
  endpoint_allocate_time = omp_get_wtime() - endpoint_stage_begin -
                           endpoint_release_time;

  // Endpoint ownership is immutable during construction.  The sparse map is
  // already split into independently locked phmap shards, so build it with the
  // full OpenMP team instead of serially inserting tens of millions of keys.
  // Dense mode has distinct endpoint writers and needs no lock.
#pragma omp parallel num_threads(num_threads)
  {
    EndpointMap::BufferedInserter buffered(&id_map_);
#pragma omp for schedule(static)
    for (size_type i = 0; i < vertices_.size(); ++i) {
      VertexAdapter adapter(vertices_[i]);
      if (use_dense_id_map_) {
        dense_id_map_[adapter.b()] = i;
        dense_id_map_[adapter.rb()] = i;
      } else {
        buffered.Insert(adapter.b(), i);
        buffered.Insert(adapter.rb(), i);
      }
    }
    if (!use_dense_id_map_) {
      buffered.FlushAll();
    }
  }
  endpoint_insert_time = omp_get_wtime() - endpoint_stage_begin -
                         endpoint_release_time - endpoint_allocate_time;
  active_ids_.resize(vertices_.size());
  std::iota(active_ids_.begin(), active_ids_.end(), size_type{0});

  // Preserve the deterministic order of the historical one-thread builder
  // separately from the multistream builder's schedule-dependent slots.  A
  // Refresh mutates terminal edges, hence this value must be captured once
  // here rather than recomputed from the current vertex later.  Initial loop
  // components followed all non-loop paths in the legacy constructor.
  legacy_order_keys_.resize(vertices_.size());
#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (size_type i = 0; i < vertices_.size(); ++i) {
    const uint64_t terminal_key = vertices_[i].TerminalOrderKey();
    assert(terminal_key < (uint64_t{1} << 63u));
    legacy_order_keys_[i] =
        terminal_key |
        (VertexAdapter(vertices_[i]).IsLoop() ? (uint64_t{1} << 63u) : 0u);
  }

  // Direct handles require one orientation bit in addition to the stable
  // 32-bit slot ID.  For larger theoretical graphs retain the exact SDBG +
  // endpoint-map traversal without truncating IDs.
  if (vertices_.size() < (size_t{1} << 31u) &&
      std::getenv("MEGAHIT_DISABLE_DIRECT_UNITIG_ADJACENCY") == nullptr) {
    try {
      direct_adjacency_.reset(new DirectAdjacency[vertices_.size()]);
      direct_adjacency_epoch_.resize(vertices_.size());
#pragma omp parallel for schedule(static)
      for (size_type i = 0; i < vertices_.size(); ++i) {
        direct_adjacency_epoch_[i].v.store(0, std::memory_order_relaxed);
      }
      use_direct_adjacency_ = true;
      xinfo("Using lazy direct unitig adjacency: {} bytes for {} stable "
            "slots\n",
            vertices_.size() *
                (sizeof(DirectAdjacency) +
                 sizeof(direct_adjacency_epoch_[0])),
            vertices_.size());
    } catch (const std::bad_alloc &) {
      direct_adjacency_.reset();
      direct_adjacency_epoch_.clear();
      use_direct_adjacency_ = false;
      xwarn("Direct unitig adjacency allocation failed; using SDBG "
            "traversal\n");
    }
  }
  if (!use_dense_id_map_) {
    assert(num_endpoints >= id_map_.size());
  }
  if (profile_unitig_build) {
    const double endpoint_total = omp_get_wtime() - endpoint_stage_begin;
    xinfo("Endpoint build profile: release={.4}, reserve={.4}, insert={.4}, "
          "active+adjacency={.4}\n",
          endpoint_release_time, endpoint_allocate_time, endpoint_insert_time,
          endpoint_total - endpoint_release_time - endpoint_allocate_time -
              endpoint_insert_time);
  }
  unitig_stage_timer.stop();
  xinfo("Unitig endpoint-map construction time: {.4}\n",
        unitig_stage_timer.elapsed());
}

void UnitigGraph::BuildDirectAdjacency(size_type id) {
  DirectAdjacency &adjacency = direct_adjacency_[id];
  for (unsigned strand = 0; strand < 2; ++strand) {
    VertexAdapter source(vertices_[id], strand, id);
    uint64_t next_starts[4];
    const int degree = sdbg_->OutgoingEdges(source.e(), next_starts);
    if (degree < 0 || degree > 4) {
      xfatal("Invalid SDBG degree {} while caching unitig {}:{}\n", degree,
             id, strand);
    }
    uint32_t *targets = adjacency.targets + strand * 4u;
    for (int i = 0; i < degree; ++i) {
      const size_type target_id = VertexIdForSdbgId(next_starts[i]);
      if (target_id >= vertices_.size()) {
        xfatal("Endpoint {} mapped to invalid unitig {} while caching {}:{}\n",
               next_starts[i], target_id, id, strand);
      }
      VertexAdapter target(vertices_[target_id], 0, target_id);
      const unsigned target_strand = target.b() == next_starts[i] ? 0u : 1u;
      targets[i] = (target_id << 1u) | target_strand;
    }
    if (degree < 4) {
      targets[degree] = std::numeric_limits<uint32_t>::max();
    }
  }
}

void UnitigGraph::EnsureDirectAdjacency(size_type id) {
  assert(use_direct_adjacency_);
  static constexpr uint32_t kBuilding =
      std::numeric_limits<uint32_t>::max();
  auto &epoch = direct_adjacency_epoch_[id].v;
  uint32_t observed = epoch.load(std::memory_order_acquire);
  while (observed != direct_adjacency_generation_) {
    if (observed != kBuilding &&
        epoch.compare_exchange_weak(observed, kBuilding,
                                    std::memory_order_acquire,
                                    std::memory_order_relaxed)) {
      BuildDirectAdjacency(id);
      epoch.store(direct_adjacency_generation_, std::memory_order_release);
      return;
    }
    if (observed == kBuilding) {
      std::this_thread::yield();
    }
    observed = epoch.load(std::memory_order_acquire);
  }
}

void UnitigGraph::InvalidateDirectAdjacency() {
  if (!use_direct_adjacency_) {
    return;
  }
  if (direct_adjacency_generation_ + 1u ==
      std::numeric_limits<uint32_t>::max()) {
#pragma omp parallel for schedule(static)
    for (size_type i = 0; i < vertices_.size(); ++i) {
      direct_adjacency_epoch_[i].v.store(0, std::memory_order_relaxed);
    }
    direct_adjacency_generation_ = 1;
  } else {
    ++direct_adjacency_generation_;
  }
}

std::vector<UnitigGraph::size_type> UnitigGraph::CollectRefreshAffected() {
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<size_type>> local_affected(max_threads);
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const size_t begin = active_ids_.size() * static_cast<size_t>(tid) /
                         static_cast<size_t>(thread_count);
    const size_t end = active_ids_.size() * static_cast<size_t>(tid + 1) /
                       static_cast<size_t>(thread_count);
    auto &affected = local_affected[tid];
    for (size_t active_index = begin; active_index < end; ++active_index) {
      const size_type id = active_ids_[active_index];
      SudoVertexAdapter adapter(vertices_[id], 0, id);
      bool flagged = adapter.IsToDelete() || adapter.IsToDisconnect();
      adapter.ReverseComplement();
      flagged = flagged || adapter.IsToDisconnect();
      if (flagged) {
        affected.push_back(id);
      }
    }
  }

  size_t affected_size = 0;
  for (const auto &ids : local_affected) {
    affected_size += ids.size();
  }
  std::vector<size_type> affected;
  affected.reserve(affected_size);
  // active_ids_ is stable-ordered and each worker owns one contiguous range,
  // so concatenation retains the same deterministic order as a full scan.
  for (auto &ids : local_affected) {
    affected.insert(affected.end(), ids.begin(), ids.end());
  }
  return affected;
}

std::vector<UnitigGraph::size_type> UnitigGraph::CollectRefreshSeeds(
    const std::vector<size_type> &affected) {
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<size_type>> local_seeds(max_threads);
#pragma omp parallel
  {
    auto &seeds = local_seeds[omp_get_thread_num()];
    seeds.reserve(affected.size() / std::max(1, omp_get_num_threads()) * 4);
#pragma omp for schedule(static)
    for (size_t affected_index = 0; affected_index < affected.size();
         ++affected_index) {
      const size_type id = affected[affected_index];
      seeds.push_back(id);
      auto adapter = MakeVertexAdapter(id);
      for (int strand = 0; strand < 2;
           ++strand, adapter.ReverseComplement()) {
        VertexAdapter neighbors[4];
        const int degree = GetNextAdapters(adapter, neighbors);
        for (int i = 0; i < degree; ++i) {
          seeds.push_back(neighbors[i].UnitigId());
        }
      }
    }
  }

  size_t seed_count = 0;
  for (const auto &ids : local_seeds) {
    seed_count += ids.size();
  }
  std::vector<size_type> seeds;
  seeds.reserve(seed_count);
  for (auto &ids : local_seeds) {
    seeds.insert(seeds.end(), ids.begin(), ids.end());
  }
  std::sort(seeds.begin(), seeds.end());
  seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
  return seeds;
}

void UnitigGraph::RefreshDisconnected(
    const std::vector<size_type> &affected) {
#pragma omp parallel for schedule(static)
  for (size_t affected_index = 0; affected_index < affected.size();
       ++affected_index) {
    const size_type i = affected[affected_index];
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsToDelete() || adapter.IsPalindrome() || adapter.IsLoop()) {
      continue;
    }

    uint8_t to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();
    uint8_t rc_to_disconnect = adapter.IsToDisconnect();
    adapter.ReverseComplement();

    if (!to_disconnect && !rc_to_disconnect) {
      continue;
    }

    if (adapter.GetLength() <= to_disconnect + rc_to_disconnect) {
      adapter.SetToDelete();
      continue;
    }

    auto old_start = adapter.b();
    auto old_end = adapter.e();
    auto old_rc_start = adapter.rb();
    auto old_rc_end = adapter.re();
    uint64_t new_start, new_end, new_rc_start, new_rc_end;

    if (to_disconnect) {
      new_start = sdbg_->NextSimplePathEdge(old_start);
      new_rc_end = sdbg_->PrevSimplePathEdge(old_rc_end);
      assert(new_start != SDBG::kNullID && new_rc_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_start);
      sdbg_->SetInvalidEdge(old_rc_end);
    } else {
      new_start = old_start;
      new_rc_end = old_rc_end;
    }

    if (rc_to_disconnect) {
      new_rc_start = sdbg_->NextSimplePathEdge(old_rc_start);
      new_end = sdbg_->PrevSimplePathEdge(old_end);
      assert(new_rc_start != SDBG::kNullID && new_end != SDBG::kNullID);
      sdbg_->SetInvalidEdge(old_rc_start);
      sdbg_->SetInvalidEdge(old_end);
    } else {
      new_rc_start = old_rc_start;
      new_end = old_end;
    }

    uint32_t new_length =
        adapter.GetLength() - to_disconnect - rc_to_disconnect;
    uint64_t new_total_depth = lround(adapter.GetAvgDepth() * new_length);
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(new_length);
    adapter.SetTotalDepth(new_total_depth);

    if (to_disconnect) {
      ReplaceVertexIdForSdbgId(old_start, new_start, i);
    }
    if (rc_to_disconnect) {
      ReplaceVertexIdForSdbgId(old_rc_start, new_rc_start, i);
    }
  }
}

void UnitigGraph::RefreshDelta(bool set_changed,
                               const std::vector<size_type> &affected,
                               const std::vector<size_type> &seeds) {
  static const uint8_t kDeleted = 0x1;
  static const uint8_t kVisited = 0x2;
  static const bool profile_refresh =
      std::getenv("MEGAHIT_PROFILE_GRAPH_REFRESH") != nullptr;
  SimpleTimer timer;
  double stages[6] = {};

  // Delta refresh operates against the exact SDBG topology while mutating
  // endpoints.  The direct cache remains available before and after this
  // interval; only records in the changed neighborhood are invalidated.
  direct_adjacency_suspended_ = true;
  timer.start();
  RefreshDisconnected(affected);
  timer.stop();
  stages[0] = timer.elapsed();

  timer.reset();
  timer.start();
#pragma omp parallel for schedule(static)
  for (size_t affected_index = 0; affected_index < affected.size();
       ++affected_index) {
    const size_type id = affected[affected_index];
    auto adapter = MakeSudoAdapter(id);
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2;
         ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.e();
      for (size_t j = 1; j < adapter.GetLength(); ++j) {
        const uint64_t prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.b());
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.IsPalindrome()) {
        break;
      }
    }
  }
  timer.stop();
  stages[1] = timer.elapsed();

  timer.reset();
  timer.start();
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint64_t>> local_starts(max_threads);
  std::vector<std::vector<size_type>> local_loops(max_threads);
#pragma omp parallel
  {
    auto &starts = local_starts[omp_get_thread_num()];
    auto &loops = local_loops[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
      const size_type seed_id = seeds[seed_index];
      auto seed_adapter = MakeSudoAdapter(seed_id);
      if (seed_adapter.IsToDelete() ||
          (seed_adapter.GetFlag() & kDeleted) ||
          seed_adapter.IsStandalone()) {
        continue;
      }
      for (unsigned strand = 0; strand < 2; ++strand) {
        auto origin = MakeSudoAdapter(seed_id, strand);
        auto cur = origin;
        size_type loop_owner = origin.UnitigId();
        bool classified = false;
        for (size_t hop = 0; hop <= active_ids_.size(); ++hop) {
          auto prev = PrevSimplePathAdapter(cur);
          if (!prev.IsValid()) {
            starts.push_back((uint64_t{cur.UnitigId()} << 1u) |
                             cur.strand());
            classified = true;
            break;
          }
          if (LegacyOrderLess(prev.UnitigId(), loop_owner)) {
            loop_owner = prev.UnitigId();
          }
          if (prev.UnitigId() == origin.UnitigId()) {
            // Delta refresh may be entered from only a subset of a newly
            // formed circular component.  Choosing the seed itself would
            // make the linearized cut depend on which dirty vertex happened
            // to reach this pass.  A full historical Refresh visits every
            // member and lets the earliest legacy vertex own the loop, so
            // recover that same owner while this predecessor walk is already
            // traversing the complete component.
            loops.push_back(loop_owner);
            classified = true;
            break;
          }
          cur = prev;
        }
        if (!classified) {
          xfatal("Delta refresh predecessor walk did not converge from {}:{}\n",
                 seed_id, strand);
        }
      }
    }
  }

  size_t start_count = 0;
  size_t loop_count = 0;
  for (int tid = 0; tid < max_threads; ++tid) {
    start_count += local_starts[tid].size();
    loop_count += local_loops[tid].size();
  }
  std::vector<uint64_t> starts;
  std::vector<size_type> loop_candidates;
  starts.reserve(start_count);
  loop_candidates.reserve(loop_count);
  for (int tid = 0; tid < max_threads; ++tid) {
    starts.insert(starts.end(), local_starts[tid].begin(),
                  local_starts[tid].end());
    loop_candidates.insert(loop_candidates.end(), local_loops[tid].begin(),
                           local_loops[tid].end());
  }
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  const bool experimental_canonical_loop_order =
      std::getenv("MEGAHIT_EXPERIMENTAL_CANONICAL_REFRESH_LOOP_ORDER") !=
      nullptr;
  if (!experimental_canonical_loop_order) {
    // Use the immutable historical construction order.  Recomputing a
    // terminal key here is insufficient because earlier Refresh calls mutate
    // the endpoints of the surviving owner vertex.
    std::sort(loop_candidates.begin(), loop_candidates.end(),
              [this](size_type lhs, size_type rhs) {
                return LegacyOrderLess(lhs, rhs);
              });
  } else {
    // Keep canonical ownership available only for experiments.  It changes
    // the decomposition of rare self-complementary cycle components and is
    // therefore not equivalent to the historical Refresh order.
    std::sort(loop_candidates.begin(), loop_candidates.end(),
              [this](size_type lhs, size_type rhs) {
                auto lhs_adapter = MakeSudoAdapter(lhs);
                auto rhs_adapter = MakeSudoAdapter(rhs);
                const uint64_t lhs_key = lhs_adapter.canonical_id();
                const uint64_t rhs_key = rhs_adapter.canonical_id();
                return lhs_key != rhs_key ? lhs_key < rhs_key : lhs < rhs;
              });
  }
  loop_candidates.erase(
      std::unique(loop_candidates.begin(), loop_candidates.end()),
      loop_candidates.end());

  AtomicBitVector locks(vertices_.size());
#pragma omp parallel for schedule(dynamic, 16)
  for (size_t start_index = 0; start_index < starts.size(); ++start_index) {
    const uint64_t start_handle = starts[start_index];
    const size_type id = static_cast<size_type>(start_handle >> 1u);
    auto adapter = MakeSudoAdapter(id, start_handle & 1u);
    if (adapter.IsToDelete() || (adapter.GetFlag() & kDeleted) ||
        adapter.IsStandalone() ||
        PrevSimplePathAdapter(adapter).IsValid()) {
      continue;
    }
    if (!locks.try_lock(id)) {
      continue;
    }

    std::vector<SudoVertexAdapter> linear_path;
    for (auto cur = NextSimplePathAdapter(adapter); cur.IsValid();
         cur = NextSimplePathAdapter(cur)) {
      linear_path.emplace_back(cur);
    }

    if (linear_path.empty()) {
      adapter.SetFlag(kVisited);
      continue;
    }

    const size_type back_id = linear_path.back().UnitigId();
    if (back_id != id && !locks.try_lock(back_id)) {
      if (LegacyOrderLess(id, back_id)) {
        locks.unlock(id);
        continue;
      }
      locks.lock(back_id);
    }

    uint32_t new_length = adapter.GetLength();
    uint64_t new_total_depth = adapter.GetTotalDepth();
    adapter.SetFlag(kVisited);
    for (auto &vertex : linear_path) {
      new_length += vertex.GetLength();
      new_total_depth += vertex.GetTotalDepth();
      if (vertex.canonical_id() != adapter.canonical_id()) {
        vertex.SetFlag(kDeleted);
      }
    }

    const uint64_t new_start = adapter.b();
    const uint64_t new_rc_end = adapter.re();
    const uint64_t new_rc_start = linear_path.back().rb();
    const uint64_t new_end = linear_path.back().e();
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(new_length);
    adapter.SetTotalDepth(new_total_depth);
    if (set_changed) {
      adapter.SetChanged();
    }
  }
  timer.stop();
  stages[2] = timer.elapsed();

  timer.reset();
  timer.start();
  for (size_t loop_index = 0; loop_index < loop_candidates.size();
       ++loop_index) {
    const size_type id = loop_candidates[loop_index];
    auto adapter = MakeSudoAdapter(id);
    if (adapter.IsStandalone() || adapter.GetFlag()) {
      continue;
    }
    if (experimental_canonical_loop_order) adapter.ToUniqueFormat();
    uint32_t length = adapter.GetLength();
    uint64_t total_depth = adapter.GetTotalDepth();
    SudoVertexAdapter next_adapter = adapter;
    while (true) {
      next_adapter = NextSimplePathAdapter(next_adapter);
      assert(next_adapter.IsValid());
      if (next_adapter.b() == adapter.b()) {
        break;
      }
      // A self-reverse-complement circular component can revisit the owner
      // through its opposite orientation before returning to the original
      // oriented endpoint.  Both adapters refer to the same logical vertex;
      // deleting that reverse view also deletes the owner and silently drops
      // the whole circle during compaction.  Other vertices in the component
      // still have to be retired exactly as before.
      if (next_adapter.canonical_id() != adapter.canonical_id()) {
        next_adapter.SetFlag(kDeleted);
      }
      length += next_adapter.GetLength();
      total_depth += next_adapter.GetTotalDepth();
    }

    const uint64_t new_start = adapter.b();
    const uint64_t new_end = sdbg_->PrevSimplePathEdge(new_start);
    const uint64_t new_rc_end = adapter.re();
    const uint64_t new_rc_start = sdbg_->NextSimplePathEdge(new_rc_end);
    assert(new_start == sdbg_->EdgeReverseComplement(new_rc_end));
    assert(new_end == sdbg_->EdgeReverseComplement(new_rc_start));
    adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
    adapter.SetLength(length);
    adapter.SetTotalDepth(total_depth);
    adapter.SetLooped();
    if (set_changed) {
      adapter.SetChanged();
    }
  }
  timer.stop();
  stages[3] = timer.elapsed();

  timer.reset();
  timer.start();
  const size_t old_vertex_count = active_ids_.size();
  const int compact_threads = omp_get_max_threads();
  if (compact_threads > 1 &&
      old_vertex_count / static_cast<size_t>(compact_threads) >= 4096) {
    std::vector<size_t> range_begin(compact_threads + 1, 0);
    std::vector<size_t> kept_offsets(compact_threads + 1, 0);
    for (int tid = 0; tid <= compact_threads; ++tid) {
      range_begin[tid] = old_vertex_count * static_cast<size_t>(tid) /
                         static_cast<size_t>(compact_threads);
    }
#pragma omp parallel for schedule(static) num_threads(compact_threads)
    for (int tid = 0; tid < compact_threads; ++tid) {
      size_t kept = 0;
      for (size_t i = range_begin[tid]; i < range_begin[tid + 1]; ++i) {
        const size_type stable_id = active_ids_[i];
        kept +=
            !(SudoVertexAdapter(vertices_[stable_id]).GetFlag() & kDeleted);
      }
      kept_offsets[tid + 1] = kept;
    }
    for (int tid = 0; tid < compact_threads; ++tid) {
      kept_offsets[tid + 1] += kept_offsets[tid];
    }
    std::vector<size_type> compacted_ids(kept_offsets.back());
#pragma omp parallel for schedule(static) num_threads(compact_threads)
    for (int tid = 0; tid < compact_threads; ++tid) {
      size_t output = kept_offsets[tid];
      for (size_t i = range_begin[tid]; i < range_begin[tid + 1]; ++i) {
        const size_type stable_id = active_ids_[i];
        if (!(SudoVertexAdapter(vertices_[stable_id]).GetFlag() & kDeleted)) {
          compacted_ids[output++] = stable_id;
        }
      }
      assert(output == kept_offsets[tid + 1]);
    }
    active_ids_.swap(compacted_ids);
  } else {
    size_t output = 0;
    for (size_type stable_id : active_ids_) {
      if (!(SudoVertexAdapter(vertices_[stable_id]).GetFlag() & kDeleted)) {
        active_ids_[output++] = stable_id;
      }
    }
    active_ids_.resize(output);
  }
  timer.stop();
  stages[4] = timer.elapsed();

  timer.reset();
  timer.start();
  std::vector<size_type> remap_ids = seeds;
  remap_ids.reserve(remap_ids.size() + starts.size() +
                    loop_candidates.size());
  for (uint64_t handle : starts) {
    remap_ids.push_back(static_cast<size_type>(handle >> 1u));
  }
  remap_ids.insert(remap_ids.end(), loop_candidates.begin(),
                   loop_candidates.end());
  std::sort(remap_ids.begin(), remap_ids.end());
  remap_ids.erase(std::unique(remap_ids.begin(), remap_ids.end()),
                  remap_ids.end());

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < remap_ids.size(); ++i) {
    const size_type id = remap_ids[i];
    auto adapter = MakeSudoAdapter(id);
    if (adapter.GetFlag() & kDeleted) {
      continue;
    }
    UpdateVertexIdForSdbgIdConcurrent(adapter.b(), id);
    UpdateVertexIdForSdbgIdConcurrent(adapter.rb(), id);
    adapter.SetFlag(0);
  }

  std::vector<std::vector<size_type>> local_post_neighbors(max_threads);
#pragma omp parallel
  {
    auto &neighbors_out = local_post_neighbors[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t i = 0; i < remap_ids.size(); ++i) {
      const size_type id = remap_ids[i];
      auto adapter = MakeSudoAdapter(id);
      if (adapter.GetFlag() & kDeleted) {
        continue;
      }
      for (int strand = 0; strand < 2;
           ++strand, adapter.ReverseComplement()) {
        SudoVertexAdapter neighbors[4];
        const int degree = GetNextAdapters(adapter, neighbors);
        for (int j = 0; j < degree; ++j) {
          neighbors_out.push_back(neighbors[j].UnitigId());
        }
      }
    }
  }
  std::vector<size_type> adjacency_dirty = remap_ids;
  for (auto &ids : local_post_neighbors) {
    adjacency_dirty.insert(adjacency_dirty.end(), ids.begin(), ids.end());
  }
  std::sort(adjacency_dirty.begin(), adjacency_dirty.end());
  adjacency_dirty.erase(
      std::unique(adjacency_dirty.begin(), adjacency_dirty.end()),
      adjacency_dirty.end());
  if (use_direct_adjacency_) {
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < adjacency_dirty.size(); ++i) {
      direct_adjacency_epoch_[adjacency_dirty[i]].v.store(
          0, std::memory_order_relaxed);
    }
  }
  direct_adjacency_suspended_ = false;
  timer.stop();
  stages[5] = timer.elapsed();

  if (profile_refresh) {
    xinfo(
        "DeltaRefresh profile: affected={}, seeds={}, starts={}, "
        "disconnect={.4}, invalidate={.4}, merge={.4}, loops={.4}, "
        "compact={.4}, remap={.4}\n",
        affected.size(), seeds.size(), starts.size(), stages[0], stages[1],
        stages[2], stages[3], stages[4], stages[5]);
  }
}

void UnitigGraph::Refresh(bool set_changed) {
  static const uint8_t kDeleted = 0x1;
  static const uint8_t kVisited = 0x2;
  static const bool profile_refresh =
      std::getenv("MEGAHIT_PROFILE_GRAPH_REFRESH") != nullptr;
  SimpleTimer refresh_timer;
  double refresh_stages[6] = {};
  static const bool delta_refresh_enabled =
      std::getenv("MEGAHIT_DISABLE_DELTA_REFRESH") == nullptr;
  static const bool force_delta_refresh =
      std::getenv("MEGAHIT_FORCE_DELTA_REFRESH") != nullptr;
  refresh_timer.start();
  const std::vector<size_type> refresh_affected = CollectRefreshAffected();

  // One changed unitig has at most eight incident oriented neighbors.  Avoid
  // even constructing a local region when that degree-derived upper bound is
  // already comparable to a dense pass.  After deduplication, select delta
  // only when four local streams are still cheaper than one graph-wide
  // topology/remap sweep.  A local seed can issue two oriented predecessor
  // walks, two successor walks, two endpoint remaps, and up to eight cache
  // invalidations, so use sixteen random graph operations as the conservative
  // cost of one seed versus one contiguous active-slot visit.  These are
  // graph-complexity decisions and do not
  // encode a NUMA-node or thread-count assumption.
  static const size_t kMaxIncidentNeighbors = 8;
  std::vector<size_type> refresh_seeds;
  const bool delta_region_can_be_sparse =
      !refresh_affected.empty() &&
      refresh_affected.size() <=
          active_ids_.size() / (kMaxIncidentNeighbors + 1);
  if (delta_refresh_enabled &&
      (force_delta_refresh || delta_region_can_be_sparse)) {
    refresh_seeds = CollectRefreshSeeds(refresh_affected);
    if (force_delta_refresh ||
        refresh_seeds.size() <= active_ids_.size() / 16) {
      RefreshDelta(set_changed, refresh_affected, refresh_seeds);
      return;
    }
  }

  // Endpoint fields are rewritten while chains are coalesced.  A cache entry
  // materialized midway through that mutation would be internally consistent
  // but stale by the end of the same generation.  Keep Refresh on the exact
  // SDBG path; direct handles resume after the mutation barrier and endpoint
  // remap below.
  direct_adjacency_suspended_ = true;
  RefreshDisconnected(refresh_affected);
  refresh_timer.stop();
  refresh_stages[0] = refresh_timer.elapsed();
  refresh_timer.reset();
  refresh_timer.start();
#pragma omp parallel for schedule(static)
  for (size_t affected_index = 0;
       affected_index < refresh_affected.size(); ++affected_index) {
    const size_type i = refresh_affected[affected_index];
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsToDelete()) {
      continue;
    }
    adapter.SetFlag(kDeleted);
    if (adapter.IsStandalone()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      uint64_t cur_edge = adapter.e();
      for (size_t j = 1; j < adapter.GetLength(); ++j) {
        auto prev = sdbg_->UniquePrevEdge(cur_edge);
        sdbg_->SetInvalidEdge(cur_edge);
        cur_edge = prev;
        assert(cur_edge != SDBG::kNullID);
      }
      assert(cur_edge == adapter.b());
      sdbg_->SetInvalidEdge(cur_edge);
      if (adapter.IsPalindrome()) {
        break;
      }
    }
  }
  refresh_timer.stop();
  refresh_stages[1] = refresh_timer.elapsed();
  // Disconnect/delete has now published the new SDBG topology.  Any direct
  // record touched by the merge pass is rebuilt once for this generation.
  InvalidateDirectAdjacency();
  refresh_timer.reset();
  refresh_timer.start();

  AtomicBitVector locks(vertices_.size());
#pragma omp parallel for
  for (size_type active_index = 0; active_index < active_ids_.size();
       ++active_index) {
    const size_type i = active_ids_[active_index];
    auto adapter = MakeSudoAdapter(i);
    if (adapter.IsStandalone() || (adapter.GetFlag() & kDeleted)) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      if (PrevSimplePathAdapter(adapter).IsValid()) {
        continue;
      }
      if (!locks.try_lock(i)) {
        break;
      }
      std::vector<SudoVertexAdapter> linear_path;
      for (auto cur = NextSimplePathAdapter(adapter); cur.IsValid();
           cur = NextSimplePathAdapter(cur)) {
        linear_path.emplace_back(cur);
      }

      if (linear_path.empty()) {
        adapter.SetFlag(kVisited);
        break;
      }

      size_type back_id = linear_path.back().UnitigId();
      if (back_id != i && !locks.try_lock(back_id)) {
        if (LegacyOrderLess(i, back_id)) {
          locks.unlock(i);
          break;
        } else {
          locks.lock(back_id);
        }
      }

      auto new_length = adapter.GetLength();
      auto new_total_depth = adapter.GetTotalDepth();
      adapter.SetFlag(kVisited);

      for (auto &v : linear_path) {
        new_length += v.GetLength();
        new_total_depth += v.GetTotalDepth();
        if (v.canonical_id() != adapter.canonical_id()) v.SetFlag(kDeleted);
      }

      auto new_start = adapter.b();
      auto new_rc_end = adapter.re();
      auto new_rc_start = linear_path.back().rb();
      auto new_end = linear_path.back().e();

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(new_length);
      adapter.SetTotalDepth(new_total_depth);
      if (set_changed) adapter.SetChanged();
      break;
    }
  }
  refresh_timer.stop();
  refresh_stages[2] = refresh_timer.elapsed();
  refresh_timer.reset();
  refresh_timer.start();

  // looped path
  const int loop_threads = std::max(1, omp_get_max_threads());
  std::vector<std::vector<size_type>> local_loop_ids(loop_threads);
#pragma omp parallel num_threads(loop_threads)
  {
    auto &local = local_loop_ids[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_type active_index = 0; active_index < active_ids_.size();
         ++active_index) {
      const size_type i = active_ids_[active_index];
      auto adapter = MakeSudoAdapter(i);
      if (!adapter.IsStandalone() && !adapter.GetFlag()) {
        local.push_back(i);
      }
    }
  }
  size_t num_loop_ids = 0;
  for (const auto &local : local_loop_ids) num_loop_ids += local.size();
  std::vector<size_type> loop_ids;
  loop_ids.reserve(num_loop_ids);
  for (auto &local : local_loop_ids) {
    loop_ids.insert(loop_ids.end(), local.begin(), local.end());
  }
  const bool experimental_canonical_loop_order =
      std::getenv("MEGAHIT_EXPERIMENTAL_CANONICAL_REFRESH_LOOP_ORDER") !=
      nullptr;
  if (!experimental_canonical_loop_order) {
    std::sort(loop_ids.begin(), loop_ids.end(),
              [this](size_type lhs, size_type rhs) {
                return LegacyOrderLess(lhs, rhs);
              });
  } else {
    std::sort(loop_ids.begin(), loop_ids.end(),
              [this](size_type lhs, size_type rhs) {
                auto lhs_adapter = MakeSudoAdapter(lhs);
                auto rhs_adapter = MakeSudoAdapter(rhs);
                const uint64_t lhs_key = lhs_adapter.canonical_id();
                const uint64_t rhs_key = rhs_adapter.canonical_id();
                return lhs_key != rhs_key ? lhs_key < rhs_key : lhs < rhs;
              });
  }
  loop_ids.erase(std::unique(loop_ids.begin(), loop_ids.end()), loop_ids.end());

  for (const size_type i : loop_ids) {
    auto adapter = MakeSudoAdapter(i);
    if (!adapter.IsStandalone() && !adapter.GetFlag()) {
      if (experimental_canonical_loop_order) adapter.ToUniqueFormat();
      uint32_t length = adapter.GetLength();
      uint64_t total_depth = adapter.GetTotalDepth();
      SudoVertexAdapter next_adapter = adapter;
      while (true) {
        next_adapter = NextSimplePathAdapter(next_adapter);
        assert(next_adapter.IsValid());
        if (next_adapter.b() == adapter.b()) {
          break;
        }
        // Do not delete the loop owner when the walk reaches its reverse
        // orientation in a self-reverse-complement circular component.
        if (next_adapter.canonical_id() != adapter.canonical_id()) {
          next_adapter.SetFlag(kDeleted);
        }
        length += next_adapter.GetLength();
        total_depth += next_adapter.GetTotalDepth();
      }

      auto new_start = adapter.b();
      auto new_end = sdbg_->PrevSimplePathEdge(new_start);
      auto new_rc_end = adapter.re();
      auto new_rc_start = sdbg_->NextSimplePathEdge(new_rc_end);
      assert(new_start == sdbg_->EdgeReverseComplement(new_rc_end));
      assert(new_end == sdbg_->EdgeReverseComplement(new_rc_start));

      adapter.SetBeginEnd(new_start, new_end, new_rc_start, new_rc_end);
      adapter.SetLength(length);
      adapter.SetTotalDepth(total_depth);
      adapter.SetLooped();
      if (set_changed) adapter.SetChanged();
    }
  }
  refresh_timer.stop();
  refresh_stages[3] = refresh_timer.elapsed();
  refresh_timer.reset();
  refresh_timer.start();

  const size_t old_vertex_count = active_ids_.size();
  const int compact_threads = omp_get_max_threads();
  static const bool parallel_stable_compaction =
      std::getenv("MEGAHIT_DISABLE_PARALLEL_VERTEX_COMPACTION") == nullptr;
  bool remapped_during_compaction = false;
  if (parallel_stable_compaction && compact_threads > 1 &&
      old_vertex_count >= (size_t{1} << 18u)) {
    std::vector<size_t> range_begin(compact_threads + 1, 0);
    std::vector<size_t> kept_offsets(compact_threads + 1, 0);
    for (int tid = 0; tid <= compact_threads; ++tid) {
      range_begin[tid] = old_vertex_count * static_cast<size_t>(tid) /
                         static_cast<size_t>(compact_threads);
    }
#pragma omp parallel for schedule(static) num_threads(compact_threads)
    for (int tid = 0; tid < compact_threads; ++tid) {
      size_t kept = 0;
      for (size_t i = range_begin[tid]; i < range_begin[tid + 1]; ++i) {
        const size_type stable_id = active_ids_[i];
        kept +=
            !(SudoVertexAdapter(vertices_[stable_id]).GetFlag() & kDeleted);
      }
      kept_offsets[tid + 1] = kept;
    }
    for (int tid = 0; tid < compact_threads; ++tid) {
      kept_offsets[tid + 1] += kept_offsets[tid];
    }

    std::vector<size_type> compacted_ids(kept_offsets.back());
#pragma omp parallel for schedule(static) num_threads(compact_threads)
    for (int tid = 0; tid < compact_threads; ++tid) {
      size_t output = kept_offsets[tid];
      for (size_t i = range_begin[tid]; i < range_begin[tid + 1]; ++i) {
        const size_type stable_id = active_ids_[i];
        auto adapter = MakeSudoAdapter(stable_id);
        if (adapter.GetFlag() & kDeleted) {
          continue;
        }
        assert(adapter.IsStandalone() || adapter.GetFlag());
        UpdateVertexIdForSdbgIdConcurrent(adapter.b(), stable_id);
        UpdateVertexIdForSdbgIdConcurrent(adapter.rb(), stable_id);
        adapter.SetFlag(0);
        compacted_ids[output++] = stable_id;
      }
      assert(output == kept_offsets[tid + 1]);
    }
    active_ids_.swap(compacted_ids);
    remapped_during_compaction = true;
  } else {
    size_t output = 0;
    for (size_type stable_id : active_ids_) {
      auto adapter = MakeSudoAdapter(stable_id);
      if (adapter.GetFlag() & kDeleted) {
        continue;
      }
      active_ids_[output++] = stable_id;
    }
    active_ids_.resize(output);
  }
  refresh_timer.stop();
  refresh_stages[4] = refresh_timer.elapsed();
  refresh_timer.reset();
  refresh_timer.start();

  if (!remapped_during_compaction) {
#pragma omp parallel for
    for (size_type active_index = 0; active_index < active_ids_.size();
         ++active_index) {
      const size_type stable_id = active_ids_[active_index];
      auto adapter = MakeSudoAdapter(stable_id);
      assert(adapter.IsStandalone() || adapter.GetFlag());
      UpdateVertexIdForSdbgIdConcurrent(adapter.b(), stable_id);
      UpdateVertexIdForSdbgIdConcurrent(adapter.rb(), stable_id);
      adapter.SetFlag(0);
    }
  }
  // Merge/loop construction changed logical endpoint ownership.  Endpoint
  // remapping above is complete, so future users may rebuild direct handles
  // against the new stable IDs.
  InvalidateDirectAdjacency();
  direct_adjacency_suspended_ = false;
  refresh_timer.stop();
  refresh_stages[5] = refresh_timer.elapsed();
  if (profile_refresh) {
    xinfo(
        "Refresh profile: disconnect={.4}, invalidate={.4}, merge={.4}, "
        "loops={.4}, compact={.4}, remap={.4}\n",
        refresh_stages[0], refresh_stages[1], refresh_stages[2],
        refresh_stages[3], refresh_stages[4], refresh_stages[5]);
  }
}

std::string UnitigGraph::VertexToDNAString(VertexAdapter v) {
  v.ToUniqueFormat();
  const size_t vertex_length = v.GetLength();
  std::string label(k() + vertex_length, '\0');
  uint8_t seq[kMaxK];
  sdbg_->GetLabel(v.b(), seq);
  for (unsigned i = 0; i < sdbg_->k(); ++i) {
    assert(seq[i] >= 1 && seq[i] <= 4);
    label[i] = "ACGT"[seq[i] - 1];
  }

  // A UnitigGraph vertex is already known to be a simple path. Materialize it
  // in forward order so each edge only needs its unique successor; the former
  // reverse walk redundantly re-checked both degrees at every base and then
  // reversed the complete string.
  uint64_t cur_edge = v.b();
  for (size_t i = 0; i < vertex_length; ++i) {
    int8_t cur_char = sdbg_->GetW(cur_edge);
    label[sdbg_->k() + i] =
        "ACGT"[cur_char > 4 ? (cur_char - 5) : (cur_char - 1)];

    if (i + 1 < vertex_length) {
      cur_edge = SimpleNextForMaterialization(cur_edge);
      if (cur_edge == SDBG::kNullID) {
        xfatal("{}, {}, {}, {}, ({}, {}), {}, {}\n", v.b(), v.e(), v.rb(),
               v.re(), sdbg_->EdgeReverseComplement(v.e()),
               sdbg_->EdgeReverseComplement(v.b()), v.GetLength(), i + 1);
      }
    }
  }

  if (cur_edge != v.e()) {
    xfatal("fwd: {}, {}, rev: {}, {}, ({}, {}) length: {}\n", v.b(), v.e(),
           v.rb(), v.re(), sdbg_->EdgeReverseComplement(v.e()),
           sdbg_->EdgeReverseComplement(v.b()), v.GetLength());
  }

  return label;
}

uint64_t UnitigGraph::SimpleNextForMaterialization(uint64_t edge) const {
  uint64_t shared_next = SDBG::kNullID;
  if (sdbg_->TryCachedNextSimplePathEdge(edge, &shared_next)) {
    return shared_next;
  }
  if (materialization_simple_codes_ &&
      edge < materialization_simple_edge_count_) {
    const uint8_t code = materialization_simple_codes_[edge];
    if (code != 0 && code != 255u) {
      uint8_t c = sdbg_->GetW(edge);
      if (c > kAlphabetSize) {
        c -= kAlphabetSize;
      }
      const uint32_t base = materialization_simple_bases_[
          (edge / kBlockSimpleNeighborEdges) * (kAlphabetSize + 1u) + c];
      if (base != kNoSimpleNeighbor) {
        const uint64_t next =
            static_cast<uint64_t>(base) + static_cast<uint64_t>(code - 1u);
        if (next < materialization_simple_edge_count_ &&
            sdbg_->IsValidEdge(next)) {
          return next;
        }
      }
    } else if (code == 255u) {
      const uint64_t next = sdbg_->NextSimplePathEdge(edge);
      if (next != SDBG::kNullID && sdbg_->IsValidEdge(next)) {
        return next;
      }
    }
  }
  return sdbg_->UniqueNextEdge(edge);
}
