//
// Created by vout on 6/23/19.
//

#include "base_engine.h"
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>

#include "utils/startup_affinity.h"
#endif

const unsigned BaseSequenceSortingEngine::kNumBuckets;

namespace {

void AdviseTransparentHugePages(void *address, size_t bytes) {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  const char *disable =
      std::getenv("MEGAHIT_DISABLE_DIRECT_HUGEPAGE_ADVICE");
  if (address == nullptr || bytes == 0 ||
      (disable != nullptr && std::strcmp(disable, "1") == 0)) {
    return;
  }

  const long page_size_long = ::sysconf(_SC_PAGESIZE);
  if (page_size_long <= 0) {
    return;
  }
  const std::uintptr_t page_size =
      static_cast<std::uintptr_t>(page_size_long);
  const std::uintptr_t allocation_begin =
      reinterpret_cast<std::uintptr_t>(address);
  if (bytes >
      std::numeric_limits<std::uintptr_t>::max() - allocation_begin) {
    return;
  }
  const std::uintptr_t allocation_end = allocation_begin + bytes;
  const std::uintptr_t begin_remainder = allocation_begin % page_size;
  const std::uintptr_t begin_advance =
      begin_remainder == 0 ? 0 : page_size - begin_remainder;
  // madvise operates on complete pages.  Exclude the allocator-owned partial
  // pages at both ends instead of rounding outwards into adjacent allocations.
  if (begin_advance > bytes) {
    return;
  }
  const std::uintptr_t aligned_begin =
      allocation_begin + begin_advance;
  const std::uintptr_t aligned_end =
      allocation_end - allocation_end % page_size;
  if (aligned_begin >= aligned_end) {
    return;
  }
  (void)::madvise(reinterpret_cast<void *>(aligned_begin),
                  aligned_end - aligned_begin, MADV_HUGEPAGE);
#else
  (void)address;
  (void)bytes;
#endif
}

/**
 * Try to fully use available memory to fit as much lv1 items as possible
 */
static std::pair<int64_t, int64_t> AdjustItemNumbers(int64_t mem_avail,
                                                     int64_t bytes_per_lv2_item,
                                                     int64_t bytes_per_lv1_item,
                                                     int64_t min_lv1_items,
                                                     int64_t min_lv2_items,
                                                     int64_t max_lv2_items) {
  int64_t num_lv1_items = 0;
  int64_t num_lv2_items = max_lv2_items;
  min_lv1_items = std::max(min_lv1_items, min_lv2_items);

  if (mem_avail < min_lv1_items * bytes_per_lv1_item +
                      min_lv2_items * bytes_per_lv2_item) {
    xfatal("Failed to adjust items number to fit in {} bytes\n", mem_avail);
  }

  // First test the requested Lv2 capacity as-is.  The old loop initialized
  // num_lv1_items to zero, so it unconditionally shrank Lv2 by 5% even when
  // all Lv1 offsets and the full scratch space already fit.  On skewed bucket
  // distributions that artificial shrink alone could force an extra complete
  // input scan.
  num_lv1_items =
      (mem_avail - bytes_per_lv2_item * num_lv2_items) /
      bytes_per_lv1_item;

  while (num_lv1_items < num_lv2_items || num_lv1_items < min_lv1_items ||
         num_lv2_items < min_lv2_items) {
    num_lv2_items = std::max(static_cast<int64_t>(lround(num_lv2_items * 0.95)),
                             min_lv2_items);
    num_lv1_items = (mem_avail - bytes_per_lv2_item * num_lv2_items) /
                    bytes_per_lv1_item;
    if (num_lv2_items == min_lv2_items && num_lv1_items < min_lv1_items) {
      xfatal("No enough memory during item adjustment. Impossible!\n");
    }
  }

  // --- adjust num_lv2_items to fit more lv1 item ---
  while (num_lv2_items * 4 > num_lv1_items) {
    if (num_lv2_items * 0.95 >= min_lv2_items) {
      num_lv2_items *= 0.95;
      num_lv1_items = (mem_avail - bytes_per_lv2_item * num_lv2_items) /
                      bytes_per_lv1_item;
    } else {
      break;
    }
  }

  return {num_lv1_items, num_lv2_items};
}
}  // namespace

void BaseSequenceSortingEngine::JoinPrefaultThreads() {
  for (auto &t : prefault_threads_) {
    if (t.joinable()) {
      t.join();
    }
  }
  prefault_threads_.clear();
}

void BaseSequenceSortingEngine::PreallocateDirectItems(
    int64_t estimated_items, int64_t estimated_memory_for_data) {
  if (prefault_direct_items_ || estimated_items <= 0 || mem_flag_ == 0 ||
      !Lv1SupportsDirectItems()) {
    return;
  }
  const int64_t words_per_item = Lv1DirectWordsPerItem();
  const int64_t aux_words_per_item = Lv1DirectAuxWordsPerItem();
  if (words_per_item <= 0 || aux_words_per_item < 0 ||
      aux_words_per_item >= words_per_item ||
      words_per_item > std::numeric_limits<int64_t>::max() /
                           static_cast<int64_t>(sizeof(uint32_t))) {
    return;
  }
  const int64_t bytes_per_item = words_per_item * sizeof(uint32_t);
  if (estimated_items > std::numeric_limits<int64_t>::max() / bytes_per_item) {
    return;
  }
  const int64_t direct_bytes = estimated_items * bytes_per_item;
  if (direct_bytes > Lv1DirectMemoryLimit()) {
    return;
  }
  // Mirror the AdjustMemory headroom policy against the estimated data size.
  const int64_t per_thread_bytes =
      static_cast<int64_t>(kNumBuckets) * sizeof(int64_t) * 3;
  const int64_t fixed_headroom = std::max<int64_t>(
      int64_t{256} << 20u, std::max<int64_t>(0, host_mem_ / 100));
  const int64_t headroom = fixed_headroom + per_thread_bytes * n_threads_;
  const int64_t mem_remained = host_mem_ - estimated_memory_for_data;
  if (headroom > mem_remained || direct_bytes > mem_remained - headroom) {
    return;
  }

  void *raw = nullptr;
  if (posix_memalign(&raw, size_t{2} << 20u,
                     static_cast<size_t>(direct_bytes)) != 0 ||
      raw == nullptr) {
    return;
  }
  prefault_direct_items_.reset(static_cast<uint32_t *>(raw));
  prefault_capacity_items_ = estimated_items;
  prefault_words_per_item_ = words_per_item;
  if (n_threads_ > 1) {
    InterleaveMemoryPages(raw, static_cast<size_t>(direct_bytes));
  }
  AdviseTransparentHugePages(raw, static_cast<size_t>(direct_bytes));

  // Fault every page with writer threads so kernel zeroing (and any hugepage
  // compaction) runs while the caller loads its input.  Joined before the
  // buffer is either adopted or freed, so nothing else touches it meanwhile.
  // Page faults include kernel allocation/zeroing and can sleep behind reclaim
  // or compaction after a write-heavy stage.  Use the requested worker team so
  // those latencies remain overlapped; bandwidth saturation naturally limits
  // the number that run simultaneously, while smaller jobs scale down.
  const unsigned n_threads =
      std::max<unsigned>(1, static_cast<unsigned>(n_threads_));
  char *base = static_cast<char *>(raw);
  const int64_t total = direct_bytes;
  for (unsigned t = 0; t < n_threads; ++t) {
    prefault_threads_.emplace_back([base, total, t, n_threads] {
      // Escape any narrow OpenMP-inherited affinity mask.
      ResetThreadAffinityToStartupMask();
      const int64_t begin = total * t / n_threads;
      const int64_t end = total * (t + 1) / n_threads;
      constexpr int64_t kPage = 4096;
      for (int64_t pos = begin - begin % kPage; pos < end; pos += kPage) {
        if (pos >= begin) {
          base[pos] = 0;
        }
      }
    });
  }
}

void BaseSequenceSortingEngine::AdjustMemory() {
  int64_t max_bucket_size =
      *std::max_element(bucket_sizes_.begin(), bucket_sizes_.end());
  int64_t total_bucket_size = 0;
  int num_non_empty = 0;
  for (unsigned i = 0; i < kNumBuckets; ++i) {
    total_bucket_size += bucket_sizes_[i];
    if (bucket_sizes_[i] > 0) {
      num_non_empty++;
    }
  }

  int64_t max_parallel_bucket_items = 0;
  for (unsigned i = 0; i < std::min(n_threads_, kNumBuckets); ++i) {
    max_parallel_bucket_items += bucket_sizes_[i];
  }
  int64_t est_lv2_items = std::max(
      3 * total_bucket_size / std::max(1, num_non_empty) * n_threads_,
      max_parallel_bucket_items);

  int64_t mem_remained = host_mem_ - meta_.memory_for_data;
  const int64_t min_lv1_items = std::max(
      static_cast<int64_t>(total_bucket_size / (kMaxLv1ScanTime - 0.5)),
      max_bucket_size);
  const int64_t lv2_bytes_per_item = meta_.words_per_lv2 * sizeof(uint32_t);
  std::pair<int64_t, int64_t> n_items;

  // Most producers deliberately restrict direct materialization to a complete
  // Lv1 pass.  A derived producer that has made its cross-bucket state
  // independent can explicitly reuse a bounded direct arena over several
  // scans.  Fundamental arrays are left uninitialized; worker stores provide
  // parallel first touch rather than serially clearing many GiB on the master
  // thread.
  const int64_t direct_words_per_item = Lv1DirectWordsPerItem();
  const int64_t direct_aux_words_per_item = Lv1DirectAuxWordsPerItem();
  if (Lv1SupportsDirectItems() && mem_flag_ != 0 && total_bucket_size > 0 &&
      direct_words_per_item > 0 && direct_aux_words_per_item >= 0 &&
      direct_aux_words_per_item < direct_words_per_item &&
      direct_words_per_item <=
          std::numeric_limits<int64_t>::max() /
              static_cast<int64_t>(sizeof(uint32_t))) {
    const int64_t direct_bytes_per_item =
        direct_words_per_item * sizeof(uint32_t);
    if (total_bucket_size <=
        std::numeric_limits<int64_t>::max() / direct_bytes_per_item) {
      const int64_t full_direct_bytes =
          total_bucket_size * direct_bytes_per_item;
      const int64_t direct_memory_limit = Lv1DirectMemoryLimit();
      // thread_meta_ retains two int64 tables per worker; direct fillers need
      // one cursor table. Preserve additional headroom for writer buffers,
      // stacks, allocator metadata and small input-dependent structures that
      // are intentionally outside MemoryStat.
      const int64_t per_thread_bytes =
          static_cast<int64_t>(kNumBuckets) * sizeof(int64_t) * 3;
      const int64_t fixed_headroom = std::max<int64_t>(
          int64_t{256} << 20u, std::max<int64_t>(0, host_mem_ / 100));
      const int64_t direct_headroom =
          fixed_headroom + per_thread_bytes * n_threads_;
      int64_t direct_capacity_items = total_bucket_size;
      if (direct_memory_limit < full_direct_bytes ||
          direct_headroom > mem_remained ||
          full_direct_bytes > mem_remained - direct_headroom) {
        if (Lv1AllowsPartialDirectItems() && direct_memory_limit > 0 &&
            direct_headroom <= mem_remained) {
          const int64_t usable_bytes = std::min(
              direct_memory_limit, mem_remained - direct_headroom);
          direct_capacity_items = usable_bytes / direct_bytes_per_item;
          if (direct_capacity_items < max_bucket_size) {
            direct_capacity_items = 0;
          }
        } else {
          direct_capacity_items = 0;
        }
      }
      if (direct_capacity_items == 0) {
        xinfo("Direct Lv1 items skipped: {} bytes exceeds usable policy "
              "capacity {}\n",
              full_direct_bytes,
              std::max<int64_t>(0, std::min(
                  direct_memory_limit, mem_remained - direct_headroom)));
      } else {
        const int64_t direct_bytes =
            direct_capacity_items * direct_bytes_per_item;
        do {
          JoinPrefaultThreads();
          if (prefault_direct_items_ &&
              prefault_words_per_item_ == direct_words_per_item &&
              prefault_capacity_items_ >= direct_capacity_items) {
            // Adopt the speculative buffer: its pages are already faulted.
            lv1_direct_items_ = std::move(prefault_direct_items_);
            prefault_capacity_items_ = 0;
            prefault_words_per_item_ = 0;
            lv1_direct_capacity_items_ = direct_capacity_items;
            lv1_direct_words_per_item_ = direct_words_per_item;
            lv1_direct_aux_words_per_item_ = direct_aux_words_per_item;
            lv1_direct_mode_ = true;
            substr_sort_ = SelectSortingFunc(
                direct_words_per_item - direct_aux_words_per_item,
                direct_aux_words_per_item, Lv2SortIgnoredLowBytes(),
                Lv2SortIgnoredHighBytes());
            xinfo("Direct Lv1 items: {} (prefaulted), words/item: {}\n",
                  direct_capacity_items, direct_words_per_item);
            return;
          }
          prefault_direct_items_.reset();
          // 2 MiB alignment satisfies the 64-byte requirement of the
          // non-temporal flushes and keeps hugepage-backed regions tidy.
          void *raw = nullptr;
          if (posix_memalign(&raw, size_t{2} << 20u,
                             static_cast<size_t>(direct_bytes)) != 0 ||
              raw == nullptr) {
            lv1_direct_items_.reset();
            lv1_direct_capacity_items_ = 0;
            lv1_direct_words_per_item_ = 0;
            lv1_direct_aux_words_per_item_ = 0;
            lv1_direct_mode_ = false;
            xwarn("Direct Lv1 allocation failed; using compact offsets\n");
            break;
          }
          lv1_direct_items_.reset(static_cast<uint32_t *>(raw));
          if (n_threads_ > 1 &&
              InterleaveMemoryPages(raw, static_cast<size_t>(direct_bytes))) {
            xinfo("Interleaved shared direct-record workspace across allowed NUMA nodes\n");
          }
          AdviseTransparentHugePages(lv1_direct_items_.get(),
                                     static_cast<size_t>(direct_bytes));
          lv1_direct_capacity_items_ = direct_capacity_items;
          lv1_direct_words_per_item_ = direct_words_per_item;
          lv1_direct_aux_words_per_item_ = direct_aux_words_per_item;
          lv1_direct_mode_ = true;
          substr_sort_ = SelectSortingFunc(
              direct_words_per_item - direct_aux_words_per_item,
              direct_aux_words_per_item, Lv2SortIgnoredLowBytes(),
              Lv2SortIgnoredHighBytes());
          xinfo("Direct Lv1 items: {}, words/item: {}, bytes: {}\n",
                direct_capacity_items, direct_words_per_item, direct_bytes);
          if (direct_capacity_items < total_bucket_size) {
            xinfo("Bounded direct Lv1 mode: {} total items in a {}-item "
                  "reusable arena\n",
                  total_bucket_size, direct_capacity_items);
          }
          return;
        } while (false);
      }
    }
  }

  // Direct mode was rejected (or the allocation failed): drop any
  // speculative prefaulted buffer before sizing the compact layout.
  JoinPrefaultThreads();
  prefault_direct_items_.reset();
  prefault_capacity_items_ = 0;
  prefault_words_per_item_ = 0;

  lv1_offset_bytes_ = Lv1BytesPerOffset();
  if (lv1_offset_bytes_ == 0 || lv1_offset_bytes_ > sizeof(uint32_t)) {
    xfatal("Invalid Lv1 locator width: {} bytes\n", lv1_offset_bytes_);
  }
  auto min_memory_required = meta_.memory_for_data +
                             min_lv1_items * lv1_offset_bytes_ +
                             max_bucket_size * +lv2_bytes_per_item;
  xinfo("Minimum memory required: {} bytes\n", min_memory_required);

  if (min_memory_required > host_mem_) {
    xwarn(
        "Memory available in less than memory required ({} < {}), "
        "still trying to perform sorting\n",
        host_mem_, min_memory_required);
    mem_remained = min_memory_required - meta_.memory_for_data;
  }

  if (mem_flag_ == 0) {
    // min memory
    int64_t est_lv1_items = total_bucket_size / (kMaxLv1ScanTime - 0.5);
    est_lv1_items = std::max(est_lv1_items, max_bucket_size);
    int64_t mem_needed =
        est_lv1_items * lv1_offset_bytes_ +
        est_lv2_items * lv2_bytes_per_item;

    if (mem_needed > mem_remained) {
      n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item,
                                  lv1_offset_bytes_, min_lv1_items,
                                  max_bucket_size, est_lv2_items);
    } else {
      n_items = AdjustItemNumbers(mem_needed, lv2_bytes_per_item,
                                  lv1_offset_bytes_, min_lv1_items,
                                  max_bucket_size, est_lv2_items);
    }
  } else if (mem_flag_ == 1) {
    // Automatic mode uses a bounded multi-pass working set for every thread
    // count.  A giant offset/sort array trades scans for page faults, TLB
    // pressure and memory-bandwidth contention; that trade depends on the
    // workload and memory hierarchy, not on a magic CPU-count boundary.
    // Explicit mem_flag > 1 remains the opt-in all-memory policy below.
    const int64_t auto_workspace_limit = Lv1AutoWorkspaceLimit();
    if (auto_workspace_limit > 0) {
      const int64_t minimum_workspace =
          min_memory_required - meta_.memory_for_data;
      const int64_t workspace = std::min(
          mem_remained, std::max(minimum_workspace, auto_workspace_limit));
      n_items = AdjustItemNumbers(workspace, lv2_bytes_per_item,
                                  lv1_offset_bytes_, min_lv1_items,
                                  max_bucket_size, est_lv2_items);
    } else {
      int64_t est_lv1_items =
          total_bucket_size / (kDefaultLv1ScanTime - 0.5);
      est_lv1_items = std::max(est_lv1_items, max_bucket_size);
      int64_t mem_needed = est_lv1_items * lv1_offset_bytes_ +
                           est_lv2_items * lv2_bytes_per_item;
      if (mem_needed > mem_remained) {
        n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item,
                                    lv1_offset_bytes_, min_lv1_items,
                                    max_bucket_size, est_lv2_items);
      } else {
        n_items = {est_lv1_items, est_lv2_items};
      }
    }
  } else {
    // At modest thread counts minimize complete Lv1 scans within the user's
    // hard memory budget.  Explicit all-memory modes take this path too.
    n_items = AdjustItemNumbers(mem_remained, lv2_bytes_per_item,
                                lv1_offset_bytes_, min_lv1_items,
                                max_bucket_size, est_lv2_items);
  }

  // A workspace between two scan-count boundaries only raises RSS: it cannot
  // remove another complete producer scan.  In automatic bounded mode, find
  // the smallest Lv1 capacity that preserves the exact number of greedy
  // bucket ranges selected by the policy above.  This uses the observed
  // bucket distribution (not an input-size special case), retains the same
  // full-scan count and Lv2 parallel-bucket capacity, and returns the otherwise
  // idle tail of the locator array to the system.
  if (mem_flag_ == 1 && Lv1AutoWorkspaceLimit() > 0 &&
      n_items.first < total_bucket_size) {
    const auto scan_count_for = [&](int64_t lv1_capacity) {
      unsigned start = 0;
      int scans = 0;
      while (start < kNumBuckets) {
        unsigned end = start;
        unsigned used_threads = 0;
        int64_t lv1_items = 0;
        int64_t lv2_items = 0;
        while (end < kNumBuckets) {
          int64_t next_lv2_items = lv2_items;
          if (used_threads < n_threads_) {
            next_lv2_items += bucket_sizes_[end];
          }
          if (next_lv2_items > n_items.second ||
              lv1_items + bucket_sizes_[end] > lv1_capacity) {
            break;
          }
          lv2_items = next_lv2_items;
          lv1_items += bucket_sizes_[end];
          if (used_threads < n_threads_) {
            ++used_threads;
          }
          ++end;
        }
        if (end == start) {
          return std::numeric_limits<int>::max();
        }
        ++scans;
        start = end;
      }
      return scans;
    };

    const int selected_scans = scan_count_for(n_items.first);
    if (selected_scans > 1 &&
        selected_scans < std::numeric_limits<int>::max()) {
      int64_t low = max_bucket_size;
      int64_t high = n_items.first;
      while (low < high) {
        const int64_t mid = low + (high - low) / 2;
        if (scan_count_for(mid) <= selected_scans) {
          high = mid;
        } else {
          low = mid + 1;
        }
      }
      if (low < n_items.first) {
        xinfo(
            "Trimmed automatic Lv1 locator capacity: {} -> {} items, "
            "preserving {} full scans\n",
            n_items.first, low, selected_scans);
        n_items.first = low;
      }
    }
  }

  if (n_items.first < min_lv1_items) {
    xfatal("No enough memory");
  }

  if (n_items.first > total_bucket_size) {
    n_items.first = total_bucket_size;
  }

  // Allocate uninitialized.  Zeroing is unnecessary (every consumed slot is
  // written before it is read), and a serial value-initialization would waste
  // seconds.  Unlike thread-private data, this workspace is written by scan
  // partitions and read by independently scheduled bucket workers, so a
  // stable interleaved NUMA policy avoids cross-round page migration.
  lv1_offsets_.reset();
  lv1_offsets_capacity_items_ = 0;
  lv2_items_.reset();
  lv2_capacity_items_ = 0;
  {
    const size_t bytes = std::max<size_t>(
        static_cast<size_t>(n_items.first) * lv1_offset_bytes_, 1);
    void *raw = nullptr;
    if (posix_memalign(&raw, size_t{2} << 20u, bytes) != 0 || raw == nullptr) {
      xfatal("Failed to allocate {} bytes for Lv1 offsets\n", bytes);
    }
    lv1_offsets_.reset(static_cast<uint8_t *>(raw));
    lv1_offsets_capacity_items_ = n_items.first;
    if (n_threads_ > 1 && InterleaveMemoryPages(raw, bytes)) {
      xinfo("Interleaved shared Lv1 locator workspace across allowed NUMA nodes\n");
    }
    // No explicit pre-touch and no hugepage advice here.  The kernel faults
    // pages lazily under the policy above.  An explicit parallel memset
    // measured strictly slower (fault-storm contention plus redundant
    // multi-GiB stores), and MADV_HUGEPAGE stalled in direct compaction under
    // defrag=madvise.
  }
  {
    const size_t words = std::max<size_t>(
        static_cast<size_t>(n_items.second) * meta_.words_per_lv2, 1);
    const size_t bytes = words * sizeof(uint32_t);
    void *raw = nullptr;
    if (posix_memalign(&raw, size_t{2} << 20u, bytes) != 0 || raw == nullptr) {
      xfatal("Failed to allocate {} bytes for Lv2 records\n", bytes);
    }
    lv2_items_.reset(static_cast<uint32_t *>(raw));
    lv2_capacity_items_ = n_items.second;
    if (n_threads_ > 1 && InterleaveMemoryPages(raw, bytes)) {
      xinfo("Interleaved shared Lv2 radix workspace across allowed NUMA nodes\n");
    }
  }
  substr_sort_ = SelectSortingFunc(
      meta_.words_per_lv2 - meta_.aux_words_per_lv2, meta_.aux_words_per_lv2,
      Lv2SortIgnoredLowBytes(), Lv2SortIgnoredHighBytes());

  xinfo("Lv1 items: {}, Lv2 items: {}\n", n_items.first, n_items.second);
  xinfo("Memory of derived class: {}, Memory for Lv1+Lv2: {}\n",
        meta_.memory_for_data,
        n_items.first * lv1_offset_bytes_ +
            n_items.second * lv2_bytes_per_item);
}

void BaseSequenceSortingEngine::Run() {
  // Set the requested team size before any derived initialization or OpenMP
  // region.  Previously the first bucket-size pass and the first Lv1 scan used
  // the process-wide OpenMP default; on large machines this could create far
  // more workers than requested by --num_cpu_threads.
  omp_set_num_threads(n_threads_);

  SimpleTimer lv0_timer;
  // read input & prepare
  lv0_timer.reset();
  lv0_timer.start();
  xinfo("Preparing data...\n");

  meta_ = Initialize();

  // A derived initializer may temporarily reserve one OpenMP worker for an
  // asynchronous producer (for example, SeqToSdbg::GenMercyEdges).  Restore
  // the requested team size before the sorting passes: the producer has
  // already finished by this point, and otherwise every Lv0/Lv1/Lv2 region
  // silently runs with one fewer worker.
  omp_set_num_threads(n_threads_);

  lv0_timer.stop();
  xinfo("Preparing data... Done. Time elapsed: {.4}\n", lv0_timer.elapsed());

  // Some producers can partition their input into exact, self-contained runs
  // without first materialising a global bucket histogram.  Give them the
  // opportunity here, after Initialize() has fixed the record layout but
  // before the generic Lv0/Lv1 machinery allocates its large workspaces.
  lv0_timer.reset();
  lv0_timer.start();
  if (RunSpecializedMainLoop()) {
    JoinPrefaultThreads();
    prefault_direct_items_.reset();
    prefault_capacity_items_ = 0;
    prefault_words_per_item_ = 0;
    lv0_timer.stop();
    xinfo("Specialized main loop done. Time elapsed: {.4}\n",
          lv0_timer.elapsed());
    lv0_timer.reset();
    lv0_timer.start();
    xinfo("Postprocessing...\n");
    Lv0Postprocess();
    lv0_timer.stop();
    xinfo("Postprocess done. Time elapsed: {.4}\n", lv0_timer.elapsed());
    return;
  }
  lv0_timer.reset();
  lv0_timer.start();
  xinfo("Preparing partitions and calculating bucket sizes...\n");

  // prepare rp bp and op
  Lv0PrepareThreadPartition();
  // calc bucket size
  Lv0CalcBucketSizeLaunchMt();
  Lv0ReorderBuckets();
  AdjustMemory();
  lv0_timer.stop();
  xinfo(
      "Preparing partitions and calculating bucket sizes... Done. Time "
      "elapsed: {.4}\n",
      lv0_timer.elapsed());

  lv0_timer.reset();
  lv0_timer.start();
  xinfo("Start main loop...\n");
  int lv1_iteration = 0;
  lv1_start_bucket_ = 0;

  while (lv1_start_bucket_ < kNumBuckets) {
    SimpleTimer lv1_timer;
    lv1_iteration++;
    // --- finds the bucket range for this iteration ---
    lv1_end_bucket_ = Lv1FindEndBuckets(lv1_start_bucket_);
    assert(lv1_start_bucket_ < lv1_end_bucket_);

    lv1_timer.reset();
    lv1_timer.start();
    xinfo("Lv1 scanning from bucket {} to {}\n", lv1_start_bucket_,
          lv1_end_bucket_);

    // --- scan to fill offset ---
    Lv1FillOffsetsLaunchMt();

    lv1_timer.stop();
    xinfo("Lv1 scanning done. Large diff: {}. Time elapsed: {.4}\n",
          lv1_special_offsets_.size(), lv1_timer.elapsed());
    lv1_timer.reset();
    lv1_timer.start();
    Lv1FetchAndSortLaunchMt();
    lv1_timer.stop();
    xinfo("Lv1 fetching & sorting done. Time elapsed: {.4}\n",
          lv1_timer.elapsed());
    lv1_start_bucket_ = lv1_end_bucket_;
  }

  if (Lv2DeferPostprocess()) {
    if (!lv1_direct_mode_ && lv1_num_items_ != 0) {
      xfatal("Deferred Lv2 postprocessing requires direct records\n");
    }
    Lv2PostprocessDeferred();
  }

  // Direct records are consumed completely by the per-bucket sorts above.
  // Release the potentially multi-GiB transient before derived
  // postprocessing writes candidates/statistics.
  lv1_direct_items_.reset();
  lv1_direct_capacity_items_ = 0;
  lv1_direct_words_per_item_ = 0;
  lv1_direct_aux_words_per_item_ = 0;
  lv1_direct_mode_ = false;

  lv0_timer.stop();
  xinfo("Main loop done. Time elapsed: {.4}\n", lv0_timer.elapsed());
  lv0_timer.reset();
  lv0_timer.start();
  xinfo("Postprocessing...\n");
  Lv0Postprocess();
  lv0_timer.stop();
  xinfo("Postprocess done. Time elapsed: {.4}\n", lv0_timer.elapsed());
}

void BaseSequenceSortingEngine::Lv0PrepareThreadPartition() {
  thread_meta_.resize(n_threads_);
  std::vector<std::pair<int64_t, int64_t>> balanced_ranges;
  const bool has_balanced_ranges =
      Lv0BuildBalancedRanges(&balanced_ranges) &&
      balanced_ranges.size() == n_threads_;
  for (unsigned t = 0; t < n_threads_; ++t) {
    ThreadMeta &meta = thread_meta_[t];
    if (has_balanced_ranges) {
      meta.seq_from = balanced_ranges[t].first;
      meta.seq_to = balanced_ranges[t].second;
    } else {
      // Distribute fixed-size reads by sequence count.
      const int64_t average = meta_.num_sequences / n_threads_;
      meta.seq_from = t * average;
      meta.seq_to =
          t < n_threads_ - 1 ? (t + 1) * average : meta_.num_sequences;
    }
    meta.offset_base = meta.seq_from < meta_.num_sequences ?
        Lv0EncodeDiffBase(meta.seq_from) : std::numeric_limits<int64_t>::max();
  }

  for (unsigned i = 0; i < kNumBuckets; ++i) {
    bucket_real_id[i] = i;
    bucket_rank_[i] = i;
  }
}

void BaseSequenceSortingEngine::Lv0ReorderBuckets() {
  std::vector<std::pair<int64_t, int>> size_and_id(kNumBuckets);

  for (unsigned i = 0; i < kNumBuckets; ++i) {
    size_and_id[i] = std::make_pair(bucket_sizes_[i], i);
  }

  std::sort(size_and_id.rbegin(), size_and_id.rend());

  for (unsigned i = 0; i < kNumBuckets; ++i) {
    bucket_sizes_[i] = size_and_id[i].first;
    bucket_real_id[i] = size_and_id[i].second;
    bucket_rank_[size_and_id[i].second] = i;
  }

  for (unsigned tid = 0; tid < n_threads_; ++tid) {
    auto old_bucket_sizes = thread_meta_[tid].bucket_sizes;
    for (unsigned i = 0; i < kNumBuckets; ++i) {
      thread_meta_[tid].bucket_sizes[i] = old_bucket_sizes[bucket_real_id[i]];
    }
  }
}

unsigned BaseSequenceSortingEngine::Lv1FindEndBuckets(unsigned start_bucket) {
  unsigned end_bucket = start_bucket;
  unsigned used_threads = 0;
  int64_t num_lv2 = 0;
  lv1_num_items_ = 0;

  std::fill(cur_lv1_buckets_.begin(), cur_lv1_buckets_.end(), 0);

  if (lv1_direct_mode_) {
    while (end_bucket < kNumBuckets) {
      if (lv1_num_items_ + bucket_sizes_[end_bucket] >
          lv1_direct_capacity_items_) {
        return end_bucket;
      }
      lv1_num_items_ += bucket_sizes_[end_bucket];
      const unsigned real_bucket = bucket_real_id[end_bucket];
      cur_lv1_buckets_[real_bucket >> 6u] |=
          uint64_t{1} << (real_bucket & 63u);
      ++end_bucket;
    }
    return end_bucket;
  }

  while (end_bucket < kNumBuckets) {
    if (used_threads < n_threads_) {
      num_lv2 += bucket_sizes_[end_bucket];
      ++used_threads;
    }

    if (num_lv2 > lv2_capacity_items_ ||
        lv1_num_items_ + bucket_sizes_[end_bucket] >
            lv1_offsets_capacity_items_) {
      return end_bucket;
    }

    lv1_num_items_ += bucket_sizes_[end_bucket];
    const unsigned real_bucket = bucket_real_id[end_bucket];
    cur_lv1_buckets_[real_bucket >> 6u] |= uint64_t{1} << (real_bucket & 63u);
    ++end_bucket;
  }

  return kNumBuckets;
}

void BaseSequenceSortingEngine::Lv1ComputeThreadBegin() {
  // set the bucket begin for the first thread
  thread_meta_[0].bucket_begin[lv1_start_bucket_] = 0;
  for (unsigned b = lv1_start_bucket_ + 1; b < lv1_end_bucket_; ++b) {
    thread_meta_[0].bucket_begin[b] = thread_meta_[0].bucket_begin[b - 1] +
                                      bucket_sizes_[b - 1];  // accumulate
  }

  // then for the remaining threads
  for (unsigned t = 1; t < n_threads_; ++t) {
    auto &current_thread_begin = thread_meta_[t].bucket_begin;
    auto &prev_thread_begin = thread_meta_[t - 1].bucket_begin;
    auto &prev_thread_sizes = thread_meta_[t - 1].bucket_sizes;
    for (unsigned b = lv1_start_bucket_; b < lv1_end_bucket_; ++b) {
      current_thread_begin[b] = prev_thread_begin[b] + prev_thread_sizes[b];
    }
  }
}

void BaseSequenceSortingEngine::Lv0CalcBucketSizeLaunchMt() {
#pragma omp parallel for
  for (unsigned t = 0; t < n_threads_; ++t) {
    auto &thread_meta = thread_meta_[t];
    Lv0CalcBucketSize(thread_meta.seq_from, thread_meta.seq_to,
                      &thread_meta.bucket_sizes);
  }
  std::fill(bucket_sizes_.begin(), bucket_sizes_.end(), 0);

  for (unsigned t = 0; t < n_threads_; ++t) {
    for (unsigned b = 0; b < kNumBuckets; ++b) {
      bucket_sizes_[b] += thread_meta_[t].bucket_sizes[b];
    }
  }
}

void BaseSequenceSortingEngine::Lv1FetchAndSortLaunchMt() {
  Lv2ThreadStatus thread_status{};
  thread_status.thread_offset.resize(n_threads_, -1);
  thread_status.rank.resize(n_threads_, 0);
  thread_status.profile =
      std::getenv("MEGAHIT_PROFILE_LV2_PHASES") != nullptr;
  if (thread_status.profile) {
    thread_status.extract_seconds.resize(n_threads_, 0.0);
    thread_status.sort_seconds.resize(n_threads_, 0.0);
    thread_status.postprocess_seconds.resize(n_threads_, 0.0);
  }
#pragma omp parallel for schedule(dynamic)
  for (auto i = lv1_start_bucket_; i < lv1_end_bucket_; ++i) {
    Lv2Sort(&thread_status, i, omp_get_thread_num());
  }
  if (thread_status.profile) {
    const auto sum = [](const std::vector<double> &values) {
      return std::accumulate(values.begin(), values.end(), 0.0);
    };
    xinfo("Lv2 aggregate CPU-seconds: extract={.4}, sort={.4}, "
          "postprocess={.4}\n",
          sum(thread_status.extract_seconds),
          sum(thread_status.sort_seconds),
          sum(thread_status.postprocess_seconds));
  }
}

void BaseSequenceSortingEngine::Lv2Sort(
    BaseSequenceSortingEngine::Lv2ThreadStatus *thread_status, unsigned b,
    int tid) {
  if (thread_status->thread_offset[tid] == -1) {
    std::lock_guard<std::mutex> lk(thread_status->mutex);
    thread_status->thread_offset[tid] = thread_status->acc;
    thread_status->acc += bucket_sizes_[b];
    thread_status->rank[tid] = thread_status->seen;
    thread_status->seen++;
  }

  if (bucket_sizes_[b] == 0) {
    const unsigned real_bucket = bucket_real_id[b];
    if (Lv2NeedsEmptyBucketPostprocess(real_bucket)) {
      Lv2Postprocess(0, 0, tid, nullptr, real_bucket);
    }
    return;
  }

  if (lv1_direct_mode_) {
    uint32_t *substr_ptr =
        lv1_direct_items_.get() +
        thread_meta_[0].bucket_begin[b] * lv1_direct_words_per_item_;
    const double sort_start =
        thread_status->profile ? omp_get_wtime() : 0.0;
    substr_sort_(substr_ptr, bucket_sizes_[b]);
    const double postprocess_start =
        thread_status->profile ? omp_get_wtime() : 0.0;
    if (thread_status->profile) {
      thread_status->sort_seconds[tid] += postprocess_start - sort_start;
    }
    if (!Lv2DeferPostprocess()) {
      Lv2Postprocess(0, bucket_sizes_[b], tid, substr_ptr,
                     bucket_real_id[b]);
    }
    if (thread_status->profile) {
      thread_status->postprocess_seconds[tid] +=
          omp_get_wtime() - postprocess_start;
    }
    return;
  }

  const size_t offset =
      thread_status->thread_offset[tid] * meta_.words_per_lv2;
  auto substr_ptr = lv2_items_.get() + offset;
  auto fetcher = OffsetFetcher(this, b, b + 1);
  const double extract_start =
      thread_status->profile ? omp_get_wtime() : 0.0;
  Lv2ExtractSubString(fetcher, substr_ptr);
  const double sort_start =
      thread_status->profile ? omp_get_wtime() : 0.0;
  substr_sort_(substr_ptr, bucket_sizes_[b]);
  const double postprocess_start =
      thread_status->profile ? omp_get_wtime() : 0.0;
  Lv2Postprocess(0, bucket_sizes_[b], tid, substr_ptr, bucket_real_id[b]);
  if (thread_status->profile) {
    thread_status->extract_seconds[tid] += sort_start - extract_start;
    thread_status->sort_seconds[tid] += postprocess_start - sort_start;
    thread_status->postprocess_seconds[tid] +=
        omp_get_wtime() - postprocess_start;
  }
}

void BaseSequenceSortingEngine::Lv1FillOffsetsLaunchMt() {
  lv1_special_offsets_.clear();
  Lv1ComputeThreadBegin();

#pragma omp parallel for
  for (unsigned t = 0; t < n_threads_; ++t) {
    OffsetFiller filler(this, lv1_start_bucket_, lv1_end_bucket_,
                        thread_meta_[t]);
    Lv1FillOffsets(filler, thread_meta_[t].seq_from, thread_meta_[t].seq_to);
  }
  std::sort(lv1_special_offsets_.begin(), lv1_special_offsets_.end(),
            [](const SpecialOffset &lhs, const SpecialOffset &rhs) {
              return lhs.item_index < rhs.item_index;
            });
}
