#include "startup_affinity.h"

#include <omp.h>
#include <pthread.h>
#include <sched.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
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

struct CacheDescription {
  unsigned level{0};
  uint64_t bytes{0};
  std::string identity;
};

CacheDescription LastLevelCacheForCpu(int cpu) {
  CacheDescription result;
#ifdef __linux__
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
    if (type != "Unified" && type != "Data") {
      continue;
    }
    char *level_end = nullptr;
    const unsigned long parsed_level =
        std::strtoul(level.c_str(), &level_end, 10);
    const uint64_t parsed_size = ParseCacheSize(size);
    if (level_end == level.c_str() || *level_end != '\0' ||
        parsed_level > std::numeric_limits<unsigned>::max() ||
        parsed_size == 0 ||
        (parsed_level < result.level ||
         (parsed_level == result.level && parsed_size <= result.bytes))) {
      continue;
    }

    std::string shared_cpus;
    std::string cache_id;
    // shared_cpu_list identifies the physical cache instance across arbitrary
    // CPU numbering, chiplet and sub-NUMA topologies.  The kernel cache ID is
    // the next-best identity on systems that omit that list.
    if (ReadTextFile(prefix + "shared_cpu_list", &shared_cpus)) {
      result.identity = "L" + std::to_string(parsed_level) + ":" +
                        type + ":cpus=" + shared_cpus;
    } else if (ReadTextFile(prefix + "id", &cache_id)) {
      result.identity = "L" + std::to_string(parsed_level) + ":" +
                        type + ":id=" + cache_id;
    } else {
      result.identity.clear();
    }
    result.level = static_cast<unsigned>(parsed_level);
    result.bytes = parsed_size;
  }
#else
  (void)cpu;
#endif
  return result;
}

uint64_t UniqueLastLevelCacheBytes(const std::vector<int> &cpus) {
  std::map<std::string, uint64_t> caches;
  for (int cpu : cpus) {
    const CacheDescription cache = LastLevelCacheForCpu(cpu);
    if (cache.bytes == 0 || cache.identity.empty()) {
      return 0;  // use the portable per-domain fallback
    }
    caches[cache.identity] = cache.bytes;
  }
  uint64_t total = 0;
  for (const auto &cache : caches) {
    if (cache.second > std::numeric_limits<uint64_t>::max() - total) {
      return std::numeric_limits<uint64_t>::max();
    }
    total += cache.second;
  }
  return total;
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
  std::vector<int> all_usable_cpus;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, &allowed_cpus)) all_usable_cpus.push_back(cpu);
  }
  topology.unique_last_level_cache_bytes =
      UniqueLastLevelCacheBytes(all_usable_cpus);

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
          LastLevelCacheForCpu(usable_cpus.front()).bytes);
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
        first_cpu < 0 ? 0 : LastLevelCacheForCpu(first_cpu).bytes);
  }
  return topology;
}

unsigned ParsePositiveUnsignedEnvironment(const char *name,
                                          unsigned fallback) {
  const char *text = std::getenv(name);
  if (text == nullptr || *text == '\0') return fallback;
  char *end = nullptr;
  errno = 0;
  const unsigned long parsed = std::strtoul(text, &end, 10);
  return errno == 0 && end != text && *end == '\0' && parsed > 0 &&
                 parsed <= std::numeric_limits<unsigned>::max()
             ? static_cast<unsigned>(parsed)
             : fallback;
}

uint64_t ParsePositiveUint64Environment(const char *name) {
  const char *text = std::getenv(name);
  if (text == nullptr || *text == '\0') return 0;
  char *end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(text, &end, 10);
  return errno == 0 && end != text && *end == '\0' && parsed > 0
             ? static_cast<uint64_t>(parsed)
             : 0;
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
  if (unique_last_level_cache_bytes != 0) {
    return unique_last_level_cache_bytes;
  }
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

bool ConfigureProcessNumaMemoryPolicy(std::string *error_message) {
  const char *text = std::getenv("MEGAHIT_NUMA_NODE");
  if (text == nullptr || *text == '\0') {
    return true;
  }

#if defined(__linux__) && defined(SYS_set_mempolicy)
  char *end = nullptr;
  errno = 0;
  const long parsed = std::strtol(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || parsed < 0 ||
      parsed > std::numeric_limits<int>::max()) {
    if (error_message != nullptr) {
      *error_message = "invalid MEGAHIT_NUMA_NODE value '" +
                       std::string(text) + "'";
    }
    return false;
  }

  const int selected_node = static_cast<int>(parsed);
  const NumaTopology &topology = GetNumaTopology();
  if (std::find(topology.node_ids.begin(), topology.node_ids.end(),
                selected_node) == topology.node_ids.end()) {
    if (error_message != nullptr) {
      *error_message = "NUMA node " + std::to_string(selected_node) +
                       " is outside this process's CPU/memory allocation";
    }
    return false;
  }

  constexpr unsigned kBitsPerWord = sizeof(unsigned long) * 8u;
  const unsigned node = static_cast<unsigned>(selected_node);
  const unsigned num_words = node / kBitsPerWord + 1u;
  std::vector<unsigned long> mask(num_words, 0);
  mask[node / kBitsPerWord] |= 1ul << (node % kBitsPerWord);
  if (::syscall(SYS_set_mempolicy, MPOL_BIND, mask.data(),
                num_words * kBitsPerWord) == 0) {
    return true;
  }
  const int saved_errno = errno;
  if (error_message != nullptr) {
    *error_message = "set_mempolicy(MPOL_BIND, node " +
                     std::to_string(selected_node) + ") failed: " +
                     std::strerror(saved_errno);
  }
  return false;
#else
  if (error_message != nullptr) {
    *error_message =
        "explicit NUMA memory binding is unavailable on this platform";
  }
  return false;
#endif
}

const RuntimeResourcePolicy &GetRuntimeResourcePolicy() {
  static const RuntimeResourcePolicy policy = [] {
    RuntimeResourcePolicy value;
    value.jobs_per_node =
        ParsePositiveUnsignedEnvironment("MEGAHIT_JOBS_PER_NODE", 1u);
    value.memory_budget_per_job = ParsePositiveUint64Environment(
        "MEGAHIT_MEMORY_BUDGET_PER_JOB");
    return value;
  }();
  return policy;
}

uint64_t RuntimeResourcePolicy::last_level_cache_budget_bytes() const {
  const uint64_t total = GetNumaTopology().total_last_level_cache_bytes();
  return total == 0 ? 0 : std::max<uint64_t>(1u, total / jobs_per_node);
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

  constexpr unsigned kBitsPerWord = sizeof(unsigned long) * 8u;
  const NumaTopology &topology = GetNumaTopology();
  int max_node = -1;
  unsigned selected_nodes = 0;
  for (int node : topology.node_ids) {
    if (node >= 0) {
      max_node = std::max(max_node, node);
      ++selected_nodes;
    }
  }
  if (selected_nodes < 2u || max_node < 0) {
    return false;
  }
  const unsigned num_words =
      (static_cast<unsigned>(max_node) + kBitsPerWord) / kBitsPerWord;
  std::vector<unsigned long> selected(num_words, 0);
  for (int node : topology.node_ids) {
    if (node >= 0) {
      selected[static_cast<unsigned>(node) / kBitsPerWord] |=
          1ul << (static_cast<unsigned>(node) % kBitsPerWord);
    }
  }
  const unsigned maxnode_bits = num_words * kBitsPerWord;
  return ::syscall(SYS_mbind, aligned_address, aligned_bytes, MPOL_INTERLEAVE,
                   selected.data(), maxnode_bits, 0) == 0;
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
