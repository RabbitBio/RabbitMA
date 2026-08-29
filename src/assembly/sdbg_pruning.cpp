/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include "sdbg_pruning.h"
#include <omp.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <queue>
#include <vector>

#include "kmlib/kmbitvector.h"
#include "utils/histgram.h"
#include "utils/utils.h"

namespace sdbg_pruning {

namespace {

using BitWord = AtomicBitVector::word_type;

inline BitWord ValidBitMask(uint64_t word_index, uint64_t bit_size) {
  const uint64_t bits_per_word = AtomicBitVector::bits_per_word();
  const uint64_t bit_begin = word_index * bits_per_word;
  const uint64_t valid_bits = std::min(bits_per_word, bit_size - bit_begin);
  if (valid_bits == bits_per_word) {
    return ~BitWord(0);
  }
  return (BitWord(1) << valid_bits) - 1;
}

inline unsigned TrailingZeroCount(BitWord word) {
  return __builtin_ctzl(static_cast<unsigned long>(word));
}

inline uint64_t BitWordIndex(uint64_t bit_index) {
  return bit_index / AtomicBitVector::bits_per_word();
}

inline void TrackRemovedWord(uint64_t edge_id,
                             AtomicBitVector *remove_word_seen) {
  if (remove_word_seen != nullptr) {
    remove_word_seen->set(BitWordIndex(edge_id));
  }
}

inline void TrackExposedWord(uint64_t edge_id,
                             AtomicBitVector *active_word_flags,
                             std::vector<uint64_t> *exposed_words,
                             std::vector<uint32_t> *exposed_edges) {
  if (exposed_edges != nullptr) {
    assert(edge_id <= std::numeric_limits<uint32_t>::max());
    exposed_edges->push_back(static_cast<uint32_t>(edge_id));
    return;
  }
  if (active_word_flags == nullptr || exposed_words == nullptr) {
    return;
  }
  const uint64_t word_index = BitWordIndex(edge_id);
  if (active_word_flags->try_lock(word_index)) {
    exposed_words->push_back(word_index);
  }
}

bool MarkBackwardTip(SDBG &dbg, int len, AtomicBitVector &ignored,
                     AtomicBitVector &to_remove, uint64_t id,
                     std::vector<uint64_t> &path,
                     bool combined_degree_query,
                     AtomicBitVector *remove_word_seen = nullptr,
                     AtomicBitVector *active_word_flags = nullptr,
                     std::vector<uint64_t> *exposed_words = nullptr,
                     std::vector<uint32_t> *exposed_edges = nullptr,
                     bool known_backward_endpoint = false,
                     bool track_every_exposure = false) {
  if (!known_backward_endpoint && !dbg.EdgeOutdegreeZero(id)) {
    return false;
  }

  uint64_t prev = SDBG::kNullID;
  uint64_t cur = id;
  bool is_tip = false;
  path.clear();
  path.push_back(id);

  for (int i = 1; i < len; ++i) {
    if (combined_degree_query) {
      const int indegree = dbg.UniqueIncomingEdge(cur, &prev);
      if (indegree != 1) {
        // UniquePrevEdge followed by EdgeIndegreeZero performs the same
        // BOSS backward/rank-select walk twice for zero- and branch-degree
        // vertices.  One complete degree query distinguishes both cases.
        is_tip = indegree == 0;
        prev = SDBG::kNullID;
        break;
      }
      if (dbg.UniqueNextEdge(prev) == SDBG::kNullID) {
        is_tip = true;
        break;
      }
    } else {
      prev = dbg.UniquePrevEdge(cur);
      if (prev == SDBG::kNullID) {
        is_tip = dbg.EdgeIndegreeZero(cur);
        break;
      } else if (dbg.UniqueNextEdge(prev) == SDBG::kNullID) {
        is_tip = true;
        break;
      }
    }
    path.push_back(prev);
    cur = prev;
  }
  if (!is_tip) {
    return false;
  }

  for (uint64_t edge_id : path) {
    to_remove.set(edge_id);
    TrackRemovedWord(edge_id, remove_word_seen);
  }
  ignored.set(id);
  ignored.set(path.back());
  if (prev != SDBG::kNullID) {
    const bool newly_exposed = ignored.try_unset(prev);
    if (newly_exposed || track_every_exposure) {
      TrackExposedWord(prev, active_word_flags, exposed_words, exposed_edges);
    }
  }

  return true;
}

bool MarkForwardTip(SDBG &dbg, int len, AtomicBitVector &ignored,
                    AtomicBitVector &to_remove, uint64_t id,
                    std::vector<uint64_t> &path,
                    bool combined_degree_query,
                    AtomicBitVector *remove_word_seen = nullptr,
                    AtomicBitVector *active_word_flags = nullptr,
                    std::vector<uint64_t> *exposed_words = nullptr,
                    std::vector<uint32_t> *exposed_edges = nullptr,
                    bool known_forward_endpoint = false,
                    bool track_every_exposure = false) {
  if (!known_forward_endpoint && !dbg.EdgeIndegreeZero(id)) {
    return false;
  }

  uint64_t next = SDBG::kNullID;
  uint64_t cur = id;
  bool is_tip = false;
  path.clear();
  path.push_back(id);

  for (int i = 1; i < len; ++i) {
    // The pre-pruning snapshot is only used for a still-valid non-zero simple
    // link.  Such a link cannot become ambiguous after deletions.  A missing,
    // overflow or stale entry falls through to the exact live degree query,
    // which also handles newly exposed simple paths.
    if (dbg.TryCachedNextSimplePathEdge(cur, &next)) {
      path.push_back(next);
      cur = next;
      continue;
    }
    if (combined_degree_query) {
      const int outdegree = dbg.UniqueOutgoingEdge(cur, &next);
      if (outdegree != 1) {
        is_tip = outdegree == 0;
        next = SDBG::kNullID;
        break;
      }
      if (dbg.UniquePrevEdge(next) == SDBG::kNullID) {
        is_tip = true;
        break;
      }
    } else {
      next = dbg.UniqueNextEdge(cur);
      if (next == SDBG::kNullID) {
        is_tip = dbg.EdgeOutdegreeZero(cur);
        break;
      } else if (dbg.UniquePrevEdge(next) == SDBG::kNullID) {
        is_tip = true;
        break;
      }
    }
    path.push_back(next);
    cur = next;
  }
  if (!is_tip) {
    return false;
  }

  for (uint64_t edge_id : path) {
    to_remove.set(edge_id);
    TrackRemovedWord(edge_id, remove_word_seen);
  }
  ignored.set(id);
  ignored.set(path.back());
  if (next != SDBG::kNullID) {
    const bool newly_exposed = ignored.try_unset(next);
    if (newly_exposed || track_every_exposure) {
      TrackExposedWord(next, active_word_flags, exposed_words, exposed_edges);
    }
  }
  return true;
}

// active_words is kept sorted.  active_word_flags makes additions idempotent
// without sorting the (usually much larger) existing frontier again.
void MergeExposedWords(std::vector<uint64_t> &active_words,
                       std::vector<std::vector<uint64_t>> &local_exposed) {
  size_t extra_size = 0;
  for (const auto &words : local_exposed) {
    extra_size += words.size();
  }
  if (extra_size == 0) {
    return;
  }

  std::vector<uint64_t> extra;
  extra.reserve(extra_size);
  for (auto &words : local_exposed) {
    extra.insert(extra.end(), words.begin(), words.end());
    words.clear();
  }
  std::sort(extra.begin(), extra.end());

  std::vector<uint64_t> merged;
  merged.reserve(active_words.size() + extra.size());
  std::merge(active_words.begin(), active_words.end(), extra.begin(),
             extra.end(), std::back_inserter(merged));
  active_words.swap(merged);
}

void FilterActiveWords(const AtomicBitVector &ignored, uint64_t edge_count,
                       AtomicBitVector &active_word_flags,
                       std::vector<uint64_t> &active_words) {
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint64_t>> local_active(max_threads);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const size_t begin = active_words.size() * tid / thread_count;
    const size_t end = active_words.size() * (tid + 1) / thread_count;
    auto &output = local_active[tid];
    output.reserve(end - begin);
    for (size_t i = begin; i < end; ++i) {
      const uint64_t word_index = active_words[i];
      const BitWord candidates =
          ~ignored.load_word(word_index) & ValidBitMask(word_index, edge_count);
      if (candidates != 0) {
        output.push_back(word_index);
      } else {
        active_word_flags.unset(word_index);
      }
    }
  }

  size_t output_size = 0;
  for (const auto &words : local_active) {
    output_size += words.size();
  }
  std::vector<uint64_t> filtered;
  filtered.reserve(output_size);
  for (auto &words : local_active) {
    filtered.insert(filtered.end(), words.begin(), words.end());
  }
  active_words.swap(filtered);
}

// Exact edge IDs are materially more compact than word IDs when endpoint
// candidates are sparse but distributed over most bit-vector words.  Keep the
// list sorted so the serial/tie semantics remain the same as an ascending
// ignored-bit scan.  A newly exposed endpoint can already be present in the
// retained frontier after a set/unset cycle, hence the final unique step.
void MergeExposedEdges(std::vector<uint32_t> &active_edges,
                       std::vector<std::vector<uint32_t>> &local_exposed) {
  size_t extra_size = 0;
  for (const auto &edges : local_exposed) extra_size += edges.size();
  if (extra_size == 0u) return;

  std::vector<uint32_t> extra;
  extra.reserve(extra_size);
  for (auto &edges : local_exposed) {
    extra.insert(extra.end(), edges.begin(), edges.end());
    edges.clear();
  }
  std::sort(extra.begin(), extra.end());
  extra.erase(std::unique(extra.begin(), extra.end()), extra.end());

  std::vector<uint32_t> merged;
  merged.reserve(active_edges.size() + extra.size());
  std::merge(active_edges.begin(), active_edges.end(), extra.begin(),
             extra.end(), std::back_inserter(merged));
  merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
  active_edges.swap(merged);
}

void FilterActiveEdges(const AtomicBitVector &ignored,
                       std::vector<uint32_t> &active_edges) {
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint32_t>> local_active(max_threads);
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const size_t begin = active_edges.size() * tid / thread_count;
    const size_t end = active_edges.size() * (tid + 1) / thread_count;
    auto &output = local_active[tid];
    output.reserve(end - begin);
    for (size_t i = begin; i < end; ++i) {
      if (!ignored.at(active_edges[i])) output.push_back(active_edges[i]);
    }
  }

  size_t output_size = 0;
  for (const auto &edges : local_active) output_size += edges.size();
  std::vector<uint32_t> filtered;
  filtered.reserve(output_size);
  for (auto &edges : local_active) {
    filtered.insert(filtered.end(), edges.begin(), edges.end());
  }
  active_edges.swap(filtered);
}

// Split the exact initial endpoint set into its two traversal directions.
// The input is already sorted and partitioned by edge ID, so a count/fill
// pair preserves that order without a global sort.  Degree is queried only
// here; later thresholds retain the direction until the endpoint is removed.
// An isolated edge is intentionally present in both lists, matching the old
// backward-then-forward passes.
void SplitDirectedEndpointFrontier(
    SDBG &dbg, const std::vector<uint32_t> &active_edges,
    std::vector<uint32_t> &backward_edges,
    std::vector<uint32_t> &forward_edges) {
  const int max_threads = omp_get_max_threads();
  std::vector<uint64_t> backward_counts(max_threads, 0u);
  std::vector<uint64_t> forward_counts(max_threads, 0u);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const size_t begin = active_edges.size() * tid / thread_count;
    const size_t end = active_edges.size() * (tid + 1) / thread_count;
    uint64_t backward_count = 0;
    uint64_t forward_count = 0;
    for (size_t i = begin; i < end; ++i) {
      const uint64_t id = active_edges[i];
      backward_count += dbg.EdgeOutdegreeZero(id);
      forward_count += dbg.EdgeIndegreeZero(id);
    }
    backward_counts[tid] = backward_count;
    forward_counts[tid] = forward_count;
  }

  std::vector<uint64_t> backward_offsets(max_threads + 1u, 0u);
  std::vector<uint64_t> forward_offsets(max_threads + 1u, 0u);
  for (int tid = 0; tid < max_threads; ++tid) {
    backward_offsets[tid + 1u] = backward_offsets[tid] + backward_counts[tid];
    forward_offsets[tid + 1u] = forward_offsets[tid] + forward_counts[tid];
  }
  backward_edges.resize(static_cast<size_t>(backward_offsets.back()));
  forward_edges.resize(static_cast<size_t>(forward_offsets.back()));

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thread_count = omp_get_num_threads();
    const size_t begin = active_edges.size() * tid / thread_count;
    const size_t end = active_edges.size() * (tid + 1) / thread_count;
    size_t backward_output = static_cast<size_t>(backward_offsets[tid]);
    size_t forward_output = static_cast<size_t>(forward_offsets[tid]);
    for (size_t i = begin; i < end; ++i) {
      const uint32_t id = active_edges[i];
      if (dbg.EdgeOutdegreeZero(id)) backward_edges[backward_output++] = id;
      if (dbg.EdgeIndegreeZero(id)) forward_edges[forward_output++] = id;
    }
    assert(backward_output == backward_offsets[tid + 1u]);
    assert(forward_output == forward_offsets[tid + 1u]);
  }
}

// Classify only endpoints exposed by the just-published deletion set.  The
// old unified frontier kept every exposed branch edge and retried both degree
// directions at every later threshold.  Dropping non-endpoints is exact:
// should a later deletion expose them again, track_every_exposure records the
// edge even when its ignored bit was already clear.
void ClassifyExposedEndpoints(
    SDBG &dbg, std::vector<std::vector<uint32_t>> &local_exposed,
    std::vector<std::vector<uint32_t>> &local_backward,
    std::vector<std::vector<uint32_t>> &local_forward) {
#pragma omp parallel for schedule(static)
  for (int tid = 0; tid < static_cast<int>(local_exposed.size()); ++tid) {
    auto &backward = local_backward[tid];
    auto &forward = local_forward[tid];
    for (uint32_t id : local_exposed[tid]) {
      if (!dbg.IsValidEdge(id)) continue;
      if (dbg.EdgeOutdegreeZero(id)) backward.push_back(id);
      if (dbg.EdgeIndegreeZero(id)) forward.push_back(id);
    }
    local_exposed[tid].clear();
  }
}

void ProfileCandidateDensity(const AtomicBitVector &ignored,
                             uint64_t edge_count, int len) {
  uint64_t candidate_words = 0;
  uint64_t candidate_edges = 0;
#pragma omp parallel for schedule(static) reduction(+ : candidate_words, candidate_edges)
  for (uint64_t word_index = 0; word_index < ignored.word_count();
       ++word_index) {
    const BitWord candidates =
        ~ignored.load_word(word_index) & ValidBitMask(word_index, edge_count);
    candidate_words += candidates != 0;
    candidate_edges += __builtin_popcountl(
        static_cast<unsigned long>(candidates));
  }
  xinfo("SDBG tip density after length {}: {} candidate words / {}, {} "
        "candidate edges / {}\n",
        len, candidate_words, ignored.word_count(), candidate_edges,
        edge_count);
}

void InvalidateTrackedEdges(SDBG &dbg, AtomicBitVector &to_remove,
                            AtomicBitVector &remove_word_seen) {
#pragma omp parallel for schedule(static)
  for (uint64_t summary_word = 0;
       summary_word < remove_word_seen.word_count(); ++summary_word) {
    BitWord touched = remove_word_seen.load_word(summary_word) &
                      ValidBitMask(summary_word, to_remove.word_count());
    while (touched != 0) {
      const uint64_t remove_word =
          summary_word * AtomicBitVector::bits_per_word() +
          TrailingZeroCount(touched);
      touched &= touched - 1;

      BitWord candidates = to_remove.load_word(remove_word) &
                           ValidBitMask(remove_word, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            remove_word * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        dbg.SetInvalidEdge(id);
      }
      // Both summaries are persistent across thresholds.  Only touched words
      // need clearing, avoiding a 1.4-billion-bit zero-fill on every pass.
      to_remove.store_word(remove_word, 0);
    }
    remove_word_seen.store_word(summary_word, 0);
  }
}

int64_t TrimDirectedEdgeFrontier(
    SDBG &dbg, int len, AtomicBitVector &ignored,
    AtomicBitVector &to_remove, AtomicBitVector &remove_word_seen,
    std::vector<uint32_t> &backward_edges,
    std::vector<uint32_t> &forward_edges, bool combined_degree_query) {
  int64_t number_tips = 0;
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint32_t>> backward_exposed(max_threads);
  std::vector<std::vector<uint32_t>> forward_exposed(max_threads);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = backward_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t i = 0; i < backward_edges.size(); ++i) {
      const uint64_t id = backward_edges[i];
      if (!ignored.at(id)) {
        number_tips += MarkBackwardTip(
            dbg, len, ignored, to_remove, id, path, combined_degree_query,
            &remove_word_seen, nullptr, nullptr, &exposed,
            /*known_backward_endpoint=*/true,
            /*track_every_exposure=*/true);
      }
    }
  }

  // Preserve the historical barrier semantics: a branch exposed by the
  // backward pass is allowed to participate in the current forward pass if
  // it was already a source before deletions are published.
  std::vector<std::vector<uint32_t>> backward_source(max_threads);
#pragma omp parallel for schedule(static)
  for (int tid = 0; tid < max_threads; ++tid) {
    auto &sources = backward_source[tid];
    for (uint32_t id : backward_exposed[tid]) {
      if (!ignored.at(id) && dbg.EdgeIndegreeZero(id)) sources.push_back(id);
    }
  }
  MergeExposedEdges(forward_edges, backward_source);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = forward_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t i = 0; i < forward_edges.size(); ++i) {
      const uint64_t id = forward_edges[i];
      if (!ignored.at(id)) {
        number_tips += MarkForwardTip(
            dbg, len, ignored, to_remove, id, path, combined_degree_query,
            &remove_word_seen, nullptr, nullptr, &exposed,
            /*known_forward_endpoint=*/true,
            /*track_every_exposure=*/true);
      }
    }
  }

  InvalidateTrackedEdges(dbg, to_remove, remove_word_seen);

  // Topology changes become visible only after the legacy two directional
  // passes.  Classify the small delta now and merge it into the retained,
  // sorted directional frontiers for the next threshold.
  std::vector<std::vector<uint32_t>> local_backward(max_threads);
  std::vector<std::vector<uint32_t>> local_forward(max_threads);
  ClassifyExposedEndpoints(dbg, backward_exposed, local_backward,
                           local_forward);
  ClassifyExposedEndpoints(dbg, forward_exposed, local_backward,
                           local_forward);
  MergeExposedEdges(backward_edges, local_backward);
  MergeExposedEdges(forward_edges, local_forward);
  FilterActiveEdges(ignored, backward_edges);
  FilterActiveEdges(ignored, forward_edges);
  return number_tips;
}

int64_t TrimFrontier(SDBG &dbg, int len, AtomicBitVector &ignored,
                     AtomicBitVector &to_remove,
                     AtomicBitVector &remove_word_seen,
                     AtomicBitVector &active_word_flags,
                     std::vector<uint64_t> &active_words,
                     bool combined_degree_query) {
  int64_t number_tips = 0;
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint64_t>> local_exposed(max_threads);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = local_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t active_index = 0; active_index < active_words.size();
         ++active_index) {
      const uint64_t word_index = active_words[active_index];
      BitWord candidates =
          ~ignored.load_word(word_index) & ValidBitMask(word_index, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            word_index * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        if (!ignored.at(id)) {
          number_tips += MarkBackwardTip(
              dbg, len, ignored, to_remove, id, path, combined_degree_query,
              &remove_word_seen, &active_word_flags, &exposed);
        }
      }
    }
  }

  // The legacy full scan sees predecessors exposed by the backward pass in
  // its following forward pass.  Merge them before starting that pass to
  // preserve the same directional semantics.
  MergeExposedWords(active_words, local_exposed);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = local_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t active_index = 0; active_index < active_words.size();
         ++active_index) {
      const uint64_t word_index = active_words[active_index];
      BitWord candidates =
          ~ignored.load_word(word_index) & ValidBitMask(word_index, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            word_index * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        if (!ignored.at(id)) {
          number_tips += MarkForwardTip(
              dbg, len, ignored, to_remove, id, path, combined_degree_query,
              &remove_word_seen, &active_word_flags, &exposed);
        }
      }
    }
  }

  MergeExposedWords(active_words, local_exposed);
  InvalidateTrackedEdges(dbg, to_remove, remove_word_seen);
  FilterActiveWords(ignored, dbg.size(), active_word_flags, active_words);
  return number_tips;
}

int64_t TrimEdgeFrontier(SDBG &dbg, int len, AtomicBitVector &ignored,
                         AtomicBitVector &to_remove,
                         AtomicBitVector &remove_word_seen,
                         std::vector<uint32_t> &active_edges,
                         bool combined_degree_query) {
  int64_t number_tips = 0;
  const int max_threads = omp_get_max_threads();
  std::vector<std::vector<uint32_t>> local_exposed(max_threads);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = local_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t active_index = 0; active_index < active_edges.size();
         ++active_index) {
      const uint64_t id = active_edges[active_index];
      if (!ignored.at(id)) {
        number_tips += MarkBackwardTip(
            dbg, len, ignored, to_remove, id, path, combined_degree_query,
            &remove_word_seen, nullptr, nullptr, &exposed);
      }
    }
  }

  // Match the historical two-pass barrier: endpoints exposed by the
  // backward walk participate in this threshold's forward walk, while those
  // exposed by the forward walk wait for the next threshold.
  MergeExposedEdges(active_edges, local_exposed);

#pragma omp parallel reduction(+ : number_tips)
  {
    std::vector<uint64_t> path;
    auto &exposed = local_exposed[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (size_t active_index = 0; active_index < active_edges.size();
         ++active_index) {
      const uint64_t id = active_edges[active_index];
      if (!ignored.at(id)) {
        number_tips += MarkForwardTip(
            dbg, len, ignored, to_remove, id, path, combined_degree_query,
            &remove_word_seen, nullptr, nullptr, &exposed);
      }
    }
  }

  MergeExposedEdges(active_edges, local_exposed);
  InvalidateTrackedEdges(dbg, to_remove, remove_word_seen);
  FilterActiveEdges(ignored, active_edges);
  return number_tips;
}

}  // namespace

double InferMinDepth(SDBG &dbg) {
  Histgram<mul_t> hist;
  const int num_threads = omp_get_max_threads();
  std::vector<std::vector<uint64_t>> local_counts(
      num_threads, std::vector<uint64_t>(kMaxMul + 1, 0));

#pragma omp parallel for
  for (uint64_t i = 0; i < dbg.size(); ++i) {
    if (dbg.IsValidEdge(i)) {
      ++local_counts[omp_get_thread_num()][dbg.EdgeMultiplicity(i)];
    }
  }

  for (int multiplicity = 0; multiplicity <= kMaxMul; ++multiplicity) {
    uint64_t count = 0;
    for (int tid = 0; tid < num_threads; ++tid) {
      count += local_counts[tid][multiplicity];
    }
    if (count != 0) {
      hist.insert(static_cast<mul_t>(multiplicity), count);
    }
  }

  double cov = hist.FirstLocalMinimum();
  for (int repeat = 1; repeat <= 100; ++repeat) {
    hist.TrimLow(static_cast<mul_t>(roundf(cov)));
    unsigned median = hist.median();
    double cov1 = sqrt(median);
    if (fabs(cov - cov1) < 1e-2) {
      return cov;
    }
    cov = cov1;
  }

  xwarn("Cannot detect min depth: unconverged");
  return 1;
}

int64_t Trim(SDBG &dbg, int len, AtomicBitVector &ignored,
             bool sparse_scan, bool combined_degree_query,
             AtomicBitVector *tracked_remove = nullptr,
             AtomicBitVector *remove_word_seen = nullptr) {
  int64_t number_tips = 0;
  AtomicBitVector local_remove;
  if (tracked_remove == nullptr) {
    local_remove.reset(dbg.size());
    tracked_remove = &local_remove;
  }
  AtomicBitVector &to_remove = *tracked_remove;
  std::vector<uint64_t> path;

  // A backward pass only clears a predecessor that still has the current
  // valid edge as an outgoing edge, so that predecessor cannot become an
  // outdegree-zero candidate during this pass.  The forward case is
  // symmetric.  It is therefore safe for a word snapshot to omit bits that
  // are cleared later in the same directional pass.  Bits set after the
  // snapshot are filtered by the per-edge ignored.at() recheck below.
  if (sparse_scan) {
#pragma omp parallel for reduction(+ : number_tips) private(path) schedule(static)
    for (uint64_t word_index = 0; word_index < ignored.word_count();
         ++word_index) {
      BitWord candidates =
          ~ignored.load_word(word_index) & ValidBitMask(word_index, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            word_index * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        if (!ignored.at(id)) {
          number_tips += MarkBackwardTip(dbg, len, ignored, to_remove, id, path,
                                         combined_degree_query,
                                         remove_word_seen);
        }
      }
    }
  } else {
#pragma omp parallel for reduction(+ : number_tips) private(path)
    for (uint64_t id = 0; id < dbg.size(); ++id) {
      if (!ignored.at(id)) {
        number_tips += MarkBackwardTip(dbg, len, ignored, to_remove, id, path,
                                       combined_degree_query,
                                       remove_word_seen);
      }
    }
  }

  if (sparse_scan) {
#pragma omp parallel for reduction(+ : number_tips) private(path) schedule(static)
    for (uint64_t word_index = 0; word_index < ignored.word_count();
         ++word_index) {
      BitWord candidates =
          ~ignored.load_word(word_index) & ValidBitMask(word_index, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            word_index * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        if (!ignored.at(id)) {
          number_tips += MarkForwardTip(dbg, len, ignored, to_remove, id, path,
                                        combined_degree_query,
                                        remove_word_seen);
        }
      }
    }
  } else {
#pragma omp parallel for reduction(+ : number_tips) private(path)
    for (uint64_t id = 0; id < dbg.size(); ++id) {
      if (!ignored.at(id)) {
        number_tips += MarkForwardTip(dbg, len, ignored, to_remove, id, path,
                                      combined_degree_query,
                                      remove_word_seen);
      }
    }
  }

  if (remove_word_seen != nullptr) {
    InvalidateTrackedEdges(dbg, to_remove, *remove_word_seen);
  } else if (sparse_scan) {
#pragma omp parallel for schedule(static)
    for (uint64_t word_index = 0; word_index < to_remove.word_count();
         ++word_index) {
      BitWord candidates =
          to_remove.load_word(word_index) & ValidBitMask(word_index, dbg.size());
      while (candidates != 0) {
        const uint64_t id =
            word_index * AtomicBitVector::bits_per_word() +
            TrailingZeroCount(candidates);
        candidates &= candidates - 1;
        dbg.SetInvalidEdge(id);
      }
    }
  } else {
#pragma omp parallel for
    for (uint64_t id = 0; id < dbg.size(); ++id) {
      if (to_remove.at(id)) {
        dbg.SetInvalidEdge(id);
      }
    }
  }
  return number_tips;
}

uint64_t RemoveTips(SDBG &dbg, int max_tip_len) {
  uint64_t number_tips = 0;
  SimpleTimer timer;
  timer.start();
  if (dbg.BuildSimplePathSnapshot(dbg.HostMemoryBudget())) {
    timer.stop();
    xinfo("Shared simple-path snapshot: {} bytes, {} overflow edges, "
          "built in {.4}s\n",
          dbg.SimplePathSnapshotBytes(),
          dbg.SimplePathSnapshotOverflowCount(), timer.elapsed());
  } else {
    timer.stop();
  }
  AtomicBitVector ignored(dbg.size());
  const bool sparse_scan =
      std::getenv("MEGAHIT_DISABLE_SPARSE_TIP_SCAN") == nullptr;
  const bool combined_degree_query =
      std::getenv("MEGAHIT_DISABLE_COMBINED_TIP_DEGREE") == nullptr;
  const bool frontier_candidate =
      std::getenv("MEGAHIT_DISABLE_TIP_FRONTIER") == nullptr;
  const bool force_frontier =
      std::getenv("MEGAHIT_FORCE_TIP_FRONTIER") != nullptr;
  bool frontier_scan = frontier_candidate;
  const bool profile_frontier =
      std::getenv("MEGAHIT_PROFILE_TIP_FRONTIER") != nullptr;
  const bool sentinel_endpoint_seed =
      std::getenv("MEGAHIT_DISABLE_SENTINEL_ENDPOINT_SEED") == nullptr;

  // The deletion summaries scale with edge_count / bits_per_word.  The
  // candidate frontier is selected below from exact 32-bit edge IDs or
  // 64-bit ignored-word IDs using their measured representation sizes.
  AtomicBitVector to_remove;
  AtomicBitVector remove_word_seen;
  AtomicBitVector active_word_flags;
  const int max_threads = omp_get_max_threads();
  std::vector<uint64_t> candidate_edges_per_partition(max_threads, 0u);
  std::vector<uint64_t> candidate_words_per_partition(max_threads, 0u);

  const double endpoint_seed_begin = omp_get_wtime();
  uint64_t seeded_endpoints = 0;
  if (sentinel_endpoint_seed) {
    seeded_endpoints = dbg.InitializeTipEndpointMask(&ignored);
    if (std::getenv("MEGAHIT_VALIDATE_SENTINEL_ENDPOINT_SEED") != nullptr) {
      uint64_t mismatches = 0;
#pragma omp parallel for schedule(static) reduction(+ : mismatches)
      for (uint64_t id = 0; id < dbg.size(); ++id) {
        const bool endpoint =
            dbg.EdgeIndegreeZero(id) || dbg.EdgeOutdegreeZero(id);
        mismatches += endpoint != !ignored.at(id);
      }
      if (mismatches != 0) {
        xfatal("Sentinel endpoint seed differs from exact degree scan at {} "
               "edges\n",
               mismatches);
      }
    }
  }

#pragma omp parallel for schedule(static)
  for (int partition = 0; partition < max_threads; ++partition) {
    const uint64_t word_begin =
        ignored.word_count() * static_cast<uint64_t>(partition) /
        max_threads;
    const uint64_t word_end =
        ignored.word_count() * static_cast<uint64_t>(partition + 1) /
        max_threads;
    uint64_t candidate_edges = 0;
    uint64_t candidate_words = 0;
    for (uint64_t word_idx = word_begin; word_idx < word_end; ++word_idx) {
      BitWord word;
      if (sentinel_endpoint_seed) {
        word = ignored.load_word(word_idx);
      } else {
        word = 0;
        const uint64_t id_begin = word_idx * ignored.bits_per_word();
        const uint64_t id_end = std::min<uint64_t>(
            id_begin + ignored.bits_per_word(), dbg.size());
        for (uint64_t id = id_begin; id < id_end; ++id) {
          if (!dbg.EdgeIndegreeZero(id) && !dbg.EdgeOutdegreeZero(id)) {
            word |= BitWord(1) << (id - id_begin);
          }
        }
        ignored.store_word(word_idx, word);
      }
      const BitWord candidates =
          (~word) & ValidBitMask(word_idx, dbg.size());
      if (candidates != 0) {
        ++candidate_words;
        candidate_edges += __builtin_popcountl(
            static_cast<unsigned long>(candidates));
      }
    }
    candidate_edges_per_partition[partition] = candidate_edges;
    candidate_words_per_partition[partition] = candidate_words;
  }

  uint64_t candidate_edges = 0;
  uint64_t candidate_words = 0;
  for (int partition = 0; partition < max_threads; ++partition) {
    candidate_edges += candidate_edges_per_partition[partition];
    candidate_words += candidate_words_per_partition[partition];
  }
  seeded_endpoints = candidate_edges;
  xinfo("SDBG tip endpoint seed: {} candidates in {.4}s ({s})\n",
        seeded_endpoints, omp_get_wtime() - endpoint_seed_begin,
        sentinel_endpoint_seed ? "sentinel" : "degree-scan");

  std::vector<uint64_t> active_words;
  std::vector<uint32_t> active_edges;
  std::vector<uint32_t> backward_edges;
  std::vector<uint32_t> forward_edges;
  bool edge_frontier = false;
  bool directional_edge_frontier = false;
  if (frontier_candidate) {
    const uint64_t word_flag_bytes =
        ((ignored.word_count() + AtomicBitVector::bits_per_word() - 1u) /
         AtomicBitVector::bits_per_word()) * sizeof(BitWord);
    const long double word_frontier_bytes =
        static_cast<long double>(candidate_words) * sizeof(uint64_t) +
        word_flag_bytes;
    const long double edge_frontier_bytes =
        static_cast<long double>(candidate_edges) * sizeof(uint32_t);
    edge_frontier =
        std::getenv("MEGAHIT_DISABLE_TIP_EDGE_FRONTIER") == nullptr &&
        dbg.size() <= std::numeric_limits<uint32_t>::max() &&
        (edge_frontier_bytes <= word_frontier_bytes ||
         std::getenv("MEGAHIT_FORCE_TIP_EDGE_FRONTIER") != nullptr);

    // A dense pass performs three complete contiguous bit streams.  The word
    // frontier is profitable only after candidate words contract, whereas an
    // exact edge frontier can be selected whenever it uses no more memory
    // than that word representation.  Both decisions depend on graph
    // structure, never on a fixed input size or socket count.
    frontier_scan = edge_frontier || force_frontier ||
                    candidate_words * uint64_t{2} < ignored.word_count();
    if (frontier_scan) {
      std::vector<uint64_t> partition_offsets(max_threads + 1u, 0u);
      const std::vector<uint64_t> &partition_counts =
          edge_frontier ? candidate_edges_per_partition
                        : candidate_words_per_partition;
      for (int partition = 0; partition < max_threads; ++partition) {
        partition_offsets[partition + 1u] =
            partition_offsets[partition] + partition_counts[partition];
      }
      if (partition_offsets.back() >
          static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        throw std::length_error("SDBG tip frontier exceeds address space");
      }
      if (edge_frontier) {
        active_edges.resize(static_cast<size_t>(partition_offsets.back()));
#pragma omp parallel for schedule(static)
        for (int partition = 0; partition < max_threads; ++partition) {
          size_t output = static_cast<size_t>(partition_offsets[partition]);
          const uint64_t word_begin =
              ignored.word_count() * static_cast<uint64_t>(partition) /
              max_threads;
          const uint64_t word_end =
              ignored.word_count() * static_cast<uint64_t>(partition + 1) /
              max_threads;
          for (uint64_t word_idx = word_begin; word_idx < word_end;
               ++word_idx) {
            BitWord candidates = ~ignored.load_word(word_idx) &
                                 ValidBitMask(word_idx, dbg.size());
            while (candidates != 0u) {
              const uint64_t id =
                  word_idx * AtomicBitVector::bits_per_word() +
                  TrailingZeroCount(candidates);
              candidates &= candidates - 1u;
              active_edges[output++] = static_cast<uint32_t>(id);
            }
          }
          assert(output == partition_offsets[partition + 1u]);
        }
      } else {
        active_words.resize(static_cast<size_t>(partition_offsets.back()));
#pragma omp parallel for schedule(static)
        for (int partition = 0; partition < max_threads; ++partition) {
          size_t output = static_cast<size_t>(partition_offsets[partition]);
          const uint64_t word_begin =
              ignored.word_count() * static_cast<uint64_t>(partition) /
              max_threads;
          const uint64_t word_end =
              ignored.word_count() * static_cast<uint64_t>(partition + 1) /
              max_threads;
          for (uint64_t word_idx = word_begin; word_idx < word_end;
               ++word_idx) {
            if (((~ignored.load_word(word_idx)) &
                 ValidBitMask(word_idx, dbg.size())) != 0u) {
              active_words[output++] = word_idx;
            }
          }
          assert(output == partition_offsets[partition + 1u]);
        }
      }
      to_remove.reset(dbg.size());
      remove_word_seen.reset(ignored.word_count());
      if (!edge_frontier) active_word_flags.reset(ignored.word_count());
#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < active_words.size(); ++i) {
        active_word_flags.set(active_words[i]);
      }

      if (edge_frontier &&
          std::getenv("MEGAHIT_DISABLE_DIRECTIONAL_TIP_FRONTIER") == nullptr) {
        const double split_begin = omp_get_wtime();
        SplitDirectedEndpointFrontier(dbg, active_edges, backward_edges,
                                      forward_edges);
        std::vector<uint32_t>().swap(active_edges);
        directional_edge_frontier = true;
        xinfo("Directional SDBG tip frontier: {} sink / {} source endpoints "
              "in {.4}s\n",
              backward_edges.size(), forward_edges.size(),
              omp_get_wtime() - split_begin);
      }
    }
  }
  if (!frontier_scan && sparse_scan) {
    // Reuse one deletion bitmap across all geometric thresholds.  A second,
    // 64x smaller summary marks only words that received removals, so neither
    // clearing nor publishing deletions scans the 1-bit-per-edge array again.
    to_remove.reset(dbg.size());
    remove_word_seen.reset(ignored.word_count());
  }

  auto trim_at_length = [&](int len) {
    if (profile_frontier && frontier_scan) {
      const uint64_t active_count =
          directional_edge_frontier
              ? backward_edges.size() + forward_edges.size()
              : (edge_frontier ? active_edges.size() : active_words.size());
      xinfo("SDBG tip frontier at length {}: {} active {s}\n", len,
            active_count,
            directional_edge_frontier
                ? "directed-edges"
                : (edge_frontier ? "edges" : "words"));
    }
    if (frontier_scan) {
      if (edge_frontier) {
        if (directional_edge_frontier) {
          return TrimDirectedEdgeFrontier(
              dbg, len, ignored, to_remove, remove_word_seen,
              backward_edges, forward_edges, combined_degree_query);
        }
        return TrimEdgeFrontier(dbg, len, ignored, to_remove,
                                remove_word_seen, active_edges,
                                combined_degree_query);
      }
      return TrimFrontier(dbg, len, ignored, to_remove, remove_word_seen,
                          active_word_flags, active_words,
                          combined_degree_query);
    }
    return Trim(dbg, len, ignored, sparse_scan, combined_degree_query,
                sparse_scan ? &to_remove : nullptr,
                sparse_scan ? &remove_word_seen : nullptr);
  };

  for (int len = 2; len < max_tip_len; len *= 2) {
    xinfo("Removing tips with length less than {}; ", len);
    timer.reset();
    timer.start();
    number_tips += trim_at_length(len);
    timer.stop();
    xinfoc("Accumulated tips removed: {}; time elapsed: {.4}\n", number_tips,
           timer.elapsed());
    if (profile_frontier && !frontier_scan) {
      ProfileCandidateDensity(ignored, dbg.size(), len);
    }
  }

  xinfo("Removing tips with length less than {}; ", max_tip_len);
  timer.reset();
  timer.start();
  number_tips += trim_at_length(max_tip_len);
  timer.stop();
  xinfoc("Accumulated tips removed: {}; time elapsed: {.4}\n", number_tips,
         timer.elapsed());
  if (profile_frontier && !frontier_scan) {
    ProfileCandidateDensity(ignored, dbg.size(), max_tip_len);
  }

  return number_tips;
}

}  // namespace sdbg_pruning
