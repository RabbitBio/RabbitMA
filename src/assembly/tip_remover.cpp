//
// Created by vout on 11/21/18.
//

#include "tip_remover.h"
#include "unitig_graph.h"

#include <algorithm>

uint64_t RemoveInitialTips(UnitigGraph &graph, uint32_t max_tip_len) {
  uint64_t removed_edge_length = 0;
  if (max_tip_len == 0) return 0;

  for (uint32_t threshold = 2;;
       threshold = std::min<uint32_t>(threshold * 2u, max_tip_len)) {
    uint64_t removed_this_pass = 0;
#pragma omp parallel for reduction(+ : removed_this_pass)
    for (UnitigGraph::size_type active_index = 0;
         active_index < graph.size(); ++active_index) {
      auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
      if (adapter.GetLength() >= threshold || adapter.IsLoop()) continue;

      UnitigGraph::VertexAdapter nexts[4], prevs[4];
      const int out_degree = graph.GetNextAdapters(adapter, nexts);
      const int in_degree = graph.GetPrevAdapters(adapter, prevs);
      const bool is_tip =
          in_degree + out_degree == 0 ||
          (in_degree == 0 && out_degree == 1) ||
          (out_degree == 0 && in_degree == 1);
      if (is_tip) {
        const uint32_t edge_length = adapter.GetLength();
        if (adapter.SetToDelete()) removed_this_pass += edge_length;
      }
    }

    removed_edge_length += removed_this_pass;
    if (removed_this_pass != 0) graph.Refresh(false);
    if (threshold >= max_tip_len) break;
  }
  return removed_edge_length;
}

uint32_t RemoveTips(UnitigGraph &graph, uint32_t max_tip_len) {
  uint32_t num_removed = 0;
  // If the previous pass removed nothing, the graph is unchanged, so any
  // vertex already evaluated at a smaller threshold would produce the same
  // (negative) result: the tip test only depends on the length threshold and
  // the neighbor degrees.  Re-examine only the newly eligible length range,
  // which skips the expensive SDBG degree queries for everything else.
  uint32_t evaluated_below = 0;
  for (uint32_t thre = 2; thre < max_tip_len;
       thre = std::min(thre * 2, max_tip_len)) {
    const uint32_t skip_below = evaluated_below;
    uint32_t removed_this_pass = 0;
#pragma omp parallel for reduction(+ : removed_this_pass)
    for (UnitigGraph::size_type active_index = 0;
         active_index < graph.size(); ++active_index) {
      auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
      if (adapter.GetLength() >= thre || adapter.GetLength() < skip_below) {
        continue;
      }
      if (adapter.IsStandalone()) {
        bool success = adapter.SetToDelete();
        assert(success);
        removed_this_pass += success;
      } else {
        UnitigGraph::VertexAdapter nexts[4], prevs[4];
        int outd = graph.GetNextAdapters(adapter, nexts);
        int ind = graph.GetPrevAdapters(adapter, prevs);

        if (ind + outd == 0) {
          bool success = adapter.SetToDelete();
          assert(success);
          removed_this_pass += success;
        } else if (outd == 1 && ind == 0) {
          if (nexts[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
            bool success = adapter.SetToDelete();
            assert(success);
            removed_this_pass += success;
          }
        } else if (outd == 0 && ind == 1) {
          if (prevs[0].GetAvgDepth() > 8 * adapter.GetAvgDepth()) {
            bool success = adapter.SetToDelete();
            assert(success);
            removed_this_pass += success;
          }
        }
      }
    }

    num_removed += removed_this_pass;
    if (removed_this_pass != 0) {
      graph.Refresh(false);
      evaluated_below = 0;  // merged unitigs invalidate previous evaluations
    } else {
      evaluated_below = thre;
    }
    if (thre >= max_tip_len) {
      break;
    }
  }
  return num_removed;
}
