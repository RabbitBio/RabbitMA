#include "startup_affinity.h"

#include <omp.h>
#include <pthread.h>
#include <sched.h>

#include <array>
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef __linux__
#include <linux/mempolicy.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

namespace {

// libgomp's own ELF constructor parses OMP_PROC_BIND/OMP_PLACES and binds the
// initial thread before any constructor of this binary runs, so capturing
// sched_getaffinity at startup would only see place 0.  The union of all
// OpenMP places is computed from the full startup mask, so it reconstructs
// the set of CPUs this process may use (and still respects taskset/cgroups,
// within which libgomp builds its place list).  Without binding, places are
// absent and the inherited mask is already correct.
cpu_set_t ComputeWideMask() {
  cpu_set_t mask;
  CPU_ZERO(&mask);
  bool any = false;
  const int num_places = omp_get_num_places();
  for (int p = 0; p < num_places; ++p) {
    const int n = omp_get_place_num_procs(p);
    if (n <= 0) {
      continue;
    }
    std::vector<int> ids(n);
    omp_get_place_proc_ids(p, ids.data());
    for (int id : ids) {
      if (id >= 0 && id < CPU_SETSIZE) {
        CPU_SET(id, &mask);
        any = true;
      }
    }
  }
  if (!any) {
    sched_getaffinity(0, sizeof(mask), &mask);
  }
  return mask;
}

bool ParseIndexList(const std::string &text, std::vector<int> *indices) {
  indices->clear();
  const char *cursor = text.c_str();
  while (*cursor != '\0' && *cursor != '\n') {
    char *end = nullptr;
    errno = 0;
    const long first = std::strtol(cursor, &end, 10);
    if (errno != 0 || end == cursor || first < 0 ||
        first > std::numeric_limits<int>::max()) {
      return false;
    }
    long last = first;
    cursor = end;
    if (*cursor == '-') {
      ++cursor;
      errno = 0;
      last = std::strtol(cursor, &end, 10);
      if (errno != 0 || end == cursor || last < first ||
          last > std::numeric_limits<int>::max()) {
        return false;
      }
      cursor = end;
    }
    for (long value = first; value <= last; ++value) {
      indices->push_back(static_cast<int>(value));
    }
    if (*cursor == ',') {
      ++cursor;
    } else if (*cursor != '\0' && *cursor != '\n') {
      return false;
    }
  }
  return !indices->empty();
}

bool ReadTextFile(const std::string &path, std::string *text) {
  std::ifstream input(path);
  return input.good() && static_cast<bool>(std::getline(input, *text));
}

uint64_t ParseCacheSize(const std::string &text) {
  if (text.empty()) {
    return 0;
  }
  char *end = nullptr;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (errno != 0 || end == text.c_str()) {
    return 0;
  }
  uint64_t multiplier = 1;
  if (*end == 'K' || *end == 'k') {
    multiplier = uint64_t{1} << 10u;
  } else if (*end == 'M' || *end == 'm') {
    multiplier = uint64_t{1} << 20u;
  } else if (*end == 'G' || *end == 'g') {
    multiplier = uint64_t{1} << 30u;
  } else if (*end != '\0') {
    return 0;
  }
  return value > std::numeric_limits<uint64_t>::max() / multiplier
             ? 0
             : static_cast<uint64_t>(value) * multiplier;
}

uint64_t LastLevelCacheForCpu(int cpu) {
#ifdef __linux__
  uint64_t largest = 0;
  // Cache indices are a dense, very small sysfs namespace.  Stop at the first
  // missing entry instead of assuming a particular number of cache levels.
  for (unsigned index = 0;; ++index) {
    const std::string prefix = "/sys/devices/system/cpu/cpu" +
                               std::to_string(cpu) + "/cache/index" +
                               std::to_string(index) + "/";
    std::string level;
    std::string type;
    std::string size;
    if (!ReadTextFile(prefix + "level", &level)) {
      break;
    }
    if (!ReadTextFile(prefix + "type", &type) ||
        !ReadTextFile(prefix + "size", &size)) {
      continue;
    }
    if (type == "Unified" || type == "Data") {
      largest = std::max(largest, ParseCacheSize(size));
    }
  }
  return largest;
#else
  (void)cpu;
  return 0;
#endif
}

unsigned PhysicalCoreCountForCpus(const std::vector<int> &cpus) {
#ifdef __linux__
  std::set<std::pair<int, int>> cores;
  for (int cpu : cpus) {
    const std::string prefix = "/sys/devices/system/cpu/cpu" +
                               std::to_string(cpu) + "/topology/";
    std::string package_text;
    std::string core_text;
    if (!ReadTextFile(prefix + "physical_package_id", &package_text) ||
        !ReadTextFile(prefix + "core_id", &core_text)) {
      return cpus.size();
    }
    char *package_end = nullptr;
    char *core_end = nullptr;
    const long package = std::strtol(package_text.c_str(), &package_end, 10);
    const long core = std::strtol(core_text.c_str(), &core_end, 10);
    if (package_end == package_text.c_str() || *package_end != '\0' ||
        core_end == core_text.c_str() || *core_end != '\0') {
      return cpus.size();
    }
    cores.emplace(static_cast<int>(package), static_cast<int>(core));
  }
  return cores.empty() ? cpus.size() : cores.size();
#else
  return cpus.size();
#endif
}

#ifdef __linux__
bool IsMemoryNodeAllowed(int node, int max_node) {
#if defined(SYS_get_mempolicy)
  constexpr unsigned kBitsPerWord = sizeof(unsigned long) * 8u;
  const unsigned num_words =
      (static_cast<unsigned>(max_node) + kBitsPerWord) / kBitsPerWord;
  std::vector<unsigned long> allowed(num_words, 0);
  int current_mode = 0;
  const unsigned maxnode_bits = num_words * kBitsPerWord;
  if (::syscall(SYS_get_mempolicy, &current_mode, allowed.data(), maxnode_bits,
                nullptr, MPOL_F_MEMS_ALLOWED) != 0) {
    return true;  // topology is still useful if the policy query is blocked
  }
  return (allowed[static_cast<unsigned>(node) / kBitsPerWord] >>
          (static_cast<unsigned>(node) % kBitsPerWord)) & 1ul;
#else
  (void)node;
  (void)max_node;
  return true;
#endif
}
#endif

NumaTopology DiscoverNumaTopology() {
  NumaTopology topology;
  topology.cpu_to_domain.assign(CPU_SETSIZE, -1);
  const cpu_set_t allowed_cpus = ComputeWideMask();

#ifdef __linux__
  std::string online_text;
  std::vector<int> online_nodes;
  if (ReadTextFile("/sys/devices/system/node/online", &online_text) &&
      ParseIndexList(online_text, &online_nodes)) {
    const int max_node = *std::max_element(online_nodes.begin(),
                                           online_nodes.end());
    for (int node : online_nodes) {
      if (!IsMemoryNodeAllowed(node, max_node)) {
        continue;
      }
      std::string cpu_text;
      std::vector<int> node_cpus;
      if (!ReadTextFile("/sys/devices/system/node/node" +
                            std::to_string(node) + "/cpulist",
                        &cpu_text) ||
          !ParseIndexList(cpu_text, &node_cpus)) {
        continue;
      }
      std::vector<int> usable_cpus;
      for (int cpu : node_cpus) {
        if (cpu >= 0 && cpu < CPU_SETSIZE && CPU_ISSET(cpu, &allowed_cpus)) {
          usable_cpus.push_back(cpu);
        }
      }
      if (usable_cpus.empty()) {
        continue;
      }
      const unsigned domain = topology.node_ids.size();
      topology.node_ids.push_back(node);
      topology.cpu_counts.push_back(usable_cpus.size());
      topology.physical_core_counts.push_back(
          PhysicalCoreCountForCpus(usable_cpus));
      topology.last_level_cache_bytes.push_back(
          LastLevelCacheForCpu(usable_cpus.front()));
      for (int cpu : usable_cpus) {
        topology.cpu_to_domain[cpu] = static_cast<int>(domain);
      }
    }
  }
#endif

  if (topology.node_ids.empty()) {
    topology.node_ids.push_back(-1);
    unsigned cpu_count = 0;
    int first_cpu = -1;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (CPU_ISSET(cpu, &allowed_cpus)) {
        topology.cpu_to_domain[cpu] = 0;
        ++cpu_count;
        if (first_cpu < 0) {
          first_cpu = cpu;
        }
      }
    }
    topology.cpu_counts.push_back(std::max(1u, cpu_count));
    std::vector<int> usable_cpus;
    usable_cpus.reserve(cpu_count);
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (CPU_ISSET(cpu, &allowed_cpus)) {
        usable_cpus.push_back(cpu);
      }
    }
    topology.physical_core_counts.push_back(
        std::max(1u, PhysicalCoreCountForCpus(usable_cpus)));
    topology.last_level_cache_bytes.push_back(
        first_cpu < 0 ? 0 : LastLevelCacheForCpu(first_cpu));
  }
  return topology;
}

#ifdef __linux__
bool PageInterior(void *address, size_t bytes, void **aligned_address,
                  size_t *aligned_bytes) {
  if (address == nullptr || bytes == 0) {
    return false;
  }
  const long page_size_long = ::sysconf(_SC_PAGESIZE);
  if (page_size_long <= 0) {
    return false;
  }
  const uintptr_t page_size = static_cast<uintptr_t>(page_size_long);
  const uintptr_t raw_begin = reinterpret_cast<uintptr_t>(address);
  if (bytes > std::numeric_limits<uintptr_t>::max() - raw_begin) {
    return false;
  }
  const uintptr_t raw_end = raw_begin + bytes;
  const uintptr_t aligned_begin =
      raw_begin + ((page_size - raw_begin % page_size) % page_size);
  const uintptr_t aligned_end = raw_end - raw_end % page_size;
  if (aligned_begin >= aligned_end) {
    return false;
  }
  *aligned_address = reinterpret_cast<void *>(aligned_begin);
  *aligned_bytes = aligned_end - aligned_begin;
  return true;
}
#endif

}  // namespace

uint64_t NumaTopology::total_last_level_cache_bytes() const {
  uint64_t total = 0;
  for (uint64_t bytes : last_level_cache_bytes) {
    if (bytes > std::numeric_limits<uint64_t>::max() - total) {
      return std::numeric_limits<uint64_t>::max();
    }
    total += bytes;
  }
  return total;
}

unsigned NumaTopology::total_physical_core_count() const {
  uint64_t total = 0;
  for (unsigned count : physical_core_counts) {
    total += count;
  }
  return static_cast<unsigned>(std::min<uint64_t>(
      total, std::numeric_limits<unsigned>::max()));
}

void ResetThreadAffinityToStartupMask() {
  static const cpu_set_t mask = ComputeWideMask();
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask);
}

const NumaTopology &GetNumaTopology() {
  static const NumaTopology topology = DiscoverNumaTopology();
  return topology;
}

unsigned CurrentNumaDomain() {
#ifdef __linux__
  const int cpu = sched_getcpu();
  const auto &mapping = GetNumaTopology().cpu_to_domain;
  if (cpu >= 0 && static_cast<size_t>(cpu) < mapping.size() &&
      mapping[cpu] >= 0) {
    return static_cast<unsigned>(mapping[cpu]);
  }
#endif
  return 0;
}

bool InterleaveMemoryPages(void *address, size_t bytes) {
#if defined(__linux__) && defined(SYS_get_mempolicy) && defined(SYS_mbind)
  void *aligned_address = nullptr;
  size_t aligned_bytes = 0;
  if (!PageInterior(address, bytes, &aligned_address, &aligned_bytes)) {
    return false;
  }

  constexpr unsigned kMaxNodes = 1024;
  constexpr unsigned kBitsPerWord = sizeof(unsigned long) * 8u;
  std::array<unsigned long, kMaxNodes / kBitsPerWord> allowed{};
  int current_mode = 0;
  if (::syscall(SYS_get_mempolicy, &current_mode, allowed.data(), kMaxNodes,
                nullptr, MPOL_F_MEMS_ALLOWED) != 0) {
    return false;
  }
  unsigned allowed_nodes = 0;
  for (unsigned long word : allowed) {
    allowed_nodes += static_cast<unsigned>(__builtin_popcountl(word));
  }
  if (allowed_nodes < 2) {
    return false;
  }
  return ::syscall(SYS_mbind, aligned_address, aligned_bytes, MPOL_INTERLEAVE,
                   allowed.data(), kMaxNodes, 0) == 0;
#else
  (void)address;
  (void)bytes;
  return false;
#endif
}

bool BindMemoryPagesToNumaDomain(void *address, size_t bytes,
                                 unsigned domain) {
#if defined(__linux__) && defined(SYS_mbind)
  const NumaTopology &topology = GetNumaTopology();
  if (domain >= topology.node_ids.size() || topology.node_ids[domain] < 0) {
    return false;
  }
  void *aligned_address = nullptr;
  size_t aligned_bytes = 0;
  if (!PageInterior(address, bytes, &aligned_address, &aligned_bytes)) {
    return false;
  }
  const unsigned node = static_cast<unsigned>(topology.node_ids[domain]);
  constexpr unsigned kBitsPerWord = sizeof(unsigned long) * 8u;
  const unsigned num_words = node / kBitsPerWord + 1u;
  std::vector<unsigned long> mask(num_words, 0);
  mask[node / kBitsPerWord] |= 1ul << (node % kBitsPerWord);
  return ::syscall(SYS_mbind, aligned_address, aligned_bytes, MPOL_BIND,
                   mask.data(), num_words * kBitsPerWord, 0) == 0;
#else
  (void)address;
  (void)bytes;
  (void)domain;
  return false;
#endif
}

bool AdviseHugePages(void *address, size_t bytes) {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  void *aligned_address = nullptr;
  size_t aligned_bytes = 0;
  if (!PageInterior(address, bytes, &aligned_address, &aligned_bytes)) {
    return false;
  }
  return ::madvise(aligned_address, aligned_bytes, MADV_HUGEPAGE) == 0;
#else
  (void)address;
  (void)bytes;
  return false;
#endif
}

bool DiscardMemoryPages(void *address, size_t bytes) {
#if defined(__linux__) && defined(MADV_DONTNEED)
  void *aligned_address = nullptr;
  size_t aligned_bytes = 0;
  if (!PageInterior(address, bytes, &aligned_address, &aligned_bytes)) {
    return false;
  }
  return ::madvise(aligned_address, aligned_bytes, MADV_DONTNEED) == 0;
#else
  (void)address;
  (void)bytes;
  return false;
#endif
}
