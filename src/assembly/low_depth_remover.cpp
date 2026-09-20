//
// Created by vout on 11/21/18.
//

#include "low_depth_remover.h"

#include <cstdlib>
#include <limits>

#include "unitig_graph.h"
#include "utils/utils.h"

namespace {

double LocalDepth(UnitigGraph &graph, UnitigGraph::VertexAdapter &adapter,
                  uint32_t local_width) {
  double total_depth = 0;
  uint64_t num_added_edges = 0;

  for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
    UnitigGraph::VertexAdapter outs[4];
    int degree = graph.GetNextAdapters(adapter, outs);

    for (int i = 0; i < degree; ++i) {
      if (outs[i].GetLength() <= local_width) {
        num_added_edges += outs[i].GetLength();
        total_depth += outs[i].GetTotalDepth();
      } else {
        num_added_edges += local_width;
        total_depth += outs[i].GetAvgDepth() * local_width;
      }
    }
  }

  if (num_added_edges == 0) {
    return 0;
  } else {
    return total_depth / num_added_edges;
  }
}

// This is called only after a pruning pass removed nothing, so the graph is
// unchanged. Find the smallest unitig depth that can become removable as the
// global threshold increases. A unitig can cross a future threshold only
// when its depth is below its fixed local threshold.
bool FindNextLocalLowDepthEvent(UnitigGraph &graph, uint32_t max_len,
                                uint32_t local_width, double local_ratio,
                                double *event_depth) {
  double next_event = std::numeric_limits<double>::infinity();

#pragma omp parallel for reduction(min : next_event)
  for (UnitigGraph::size_type active_index = 0; active_index < graph.size();
       ++active_index) {
    auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
    if (adapter.IsStandalone() || adapter.GetLength() > max_len) continue;
    int indegree = graph.InDegree(adapter);
    int outdegree = graph.OutDegree(adapter);
    if (indegree + outdegree == 0) continue;
    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      const double depth = adapter.GetAvgDepth();
      const double local_threshold =
          LocalDepth(graph, adapter, local_width) * local_ratio;
      if (depth < local_threshold && depth < next_event) next_event = depth;
    }
  }

  *event_depth = next_event;
  return next_event < std::numeric_limits<double>::infinity();
}

// Keep the original sequence t, 1.1*t, 1.1^2*t, ... exactly.  Returning the
// first member strictly above event_depth preserves the pass on which the
// strict depth < threshold predicate first becomes true.
bool AdvanceToLocalLowDepthEvent(double current_threshold, double event_depth,
                                 double *next_threshold,
                                 uint32_t *num_noop_thresholds) {
  double next = current_threshold * 1.1;
  uint32_t skipped = 0;
  while (next < kMaxMul && next <= event_depth) {
    const double previous = next;
    next *= 1.1;
    ++skipped;
    if (next <= previous) {
      return false;
    }
  }
  *next_threshold = next;
  *num_noop_thresholds = skipped;
  return next < kMaxMul;
}

}  // namespace

bool RemoveLocalLowDepth(UnitigGraph &graph, double min_depth, uint32_t max_len,
                         uint32_t local_width, double local_ratio,
                         bool permanent_rm, uint32_t *num_removed) {
  bool need_refresh = false;
  uint32_t removed = 0;
  std::atomic_bool is_changed{false};

#pragma omp parallel for reduction(+ : removed) reduction(|| : need_refresh)
  for (UnitigGraph::size_type active_index = 0;
       active_index < graph.size(); ++active_index) {
    auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
    if (adapter.IsStandalone() || adapter.GetLength() > max_len) {
      continue;
    }
    int indegree = graph.InDegree(adapter);
    int outdegree = graph.OutDegree(adapter);
    if (indegree + outdegree == 0) {
      continue;
    }

    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      double depth = adapter.GetAvgDepth();
      if (is_changed.load(std::memory_order_relaxed) && depth > min_depth)
        continue;
      double mean = LocalDepth(graph, adapter, local_width);
      double threshold = min_depth;

      if (min_depth < mean * local_ratio)
        is_changed.store(true, std::memory_order_relaxed);
      else
        threshold = mean * local_ratio;

      if (depth < threshold) {
        is_changed.store(true, std::memory_order_relaxed);
        need_refresh = true;
        bool success = adapter.SetToDelete();
        assert(success);
        removed += success;
      }
    }
  }

  if (need_refresh) {
    bool set_changed = !permanent_rm;
    graph.Refresh(set_changed);
  }
  *num_removed = removed;
  return is_changed;
}

uint32_t IterateLocalLowDepth(UnitigGraph &graph, double min_depth,
                              uint32_t min_len, uint32_t local_width,
                              double local_ratio, bool permanent_rm) {
  uint32_t total_removed = 0;
  uint32_t iteration = 0;
  static const bool profile_iterations =
      std::getenv("MEGAHIT_PROFILE_LOW_DEPTH") != nullptr;
  static const bool event_skip_enabled =
      std::getenv("MEGAHIT_EXPERIMENTAL_LOW_DEPTH_EVENT_SKIP") != nullptr;
  while (min_depth < kMaxMul) {
    ++iteration;
    SimpleTimer timer;
    timer.start();
    uint32_t num_removed = 0;
    const bool changed =
        RemoveLocalLowDepth(graph, min_depth, min_len, local_width,
                            local_ratio, permanent_rm, &num_removed);
    timer.stop();
    if (profile_iterations) {
      xinfo("Low-depth iteration {}: threshold={.4}, removed={}, time={.4}\n",
            iteration, min_depth, num_removed, timer.elapsed());
    }
    if (!changed) {
      break;
    }
    total_removed += num_removed;
    if (event_skip_enabled && num_removed == 0) {
      double event_depth = 0;
      if (!FindNextLocalLowDepthEvent(graph, min_len, local_width, local_ratio,
                                      &event_depth)) {
        if (profile_iterations) {
          xinfo(
              "Low-depth event skip: no future deletion event; stop after "
              "no-op iteration {}\n",
              iteration);
        }
        break;
      }
      double next_threshold = min_depth;
      uint32_t skipped_thresholds = 0;
      if (!AdvanceToLocalLowDepthEvent(min_depth, event_depth, &next_threshold,
                                       &skipped_thresholds)) {
        if (profile_iterations) {
          xinfo(
              "Low-depth event skip: next deletion event is outside the "
              "threshold range\n");
        }
        break;
      }
      if (profile_iterations && skipped_thresholds > 0) {
        xinfo(
            "Low-depth event skip: skipped={} no-op thresholds, "
            "next={.4}, trigger_depth={.4}\n",
            skipped_thresholds, next_threshold, event_depth);
      }
      min_depth = next_threshold;
    } else {
      min_depth *= 1.1;
    }
  }
  return total_removed;
}

uint32_t RemoveLowDepth(UnitigGraph &graph, double min_depth) {
  uint32_t num_removed = 0;
#pragma omp parallel for reduction(+ : num_removed)
  for (UnitigGraph::size_type active_index = 0;
       active_index < graph.size(); ++active_index) {
    auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
    if (adapter.GetAvgDepth() < min_depth) {
      bool success = adapter.SetToDelete();
      assert(success);
      num_removed += success;
    }
  }
  if (num_removed != 0) {
    graph.Refresh(false);
  }
  return num_removed;
}
