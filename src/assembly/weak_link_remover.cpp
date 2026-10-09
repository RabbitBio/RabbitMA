//
// Created by vout on 11/21/18.
//

#include "weak_link_remover.h"
#include "unitig_graph.h"

uint32_t DisconnectWeakLinks(UnitigGraph &graph, double local_ratio = 0.1) {
  const bool profile = std::getenv("MEGAHIT_PROFILE_GRAPH_REFRESH") != nullptr;
  const double detect_begin = profile ? omp_get_wtime() : 0;
  const size_t active_before = graph.size();
  uint32_t num_disconnected = 0;
#pragma omp parallel for reduction(+ : num_disconnected)
  for (UnitigGraph::size_type active_index = 0;
       active_index < graph.size(); ++active_index) {
    auto adapter = graph.MakeVertexAdapter(graph.active_id(active_index));
    if (adapter.IsStandalone() || adapter.IsPalindrome()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      UnitigGraph::VertexAdapter nexts[4];
      double depths[4];
      double total_depth = 0;
      int degree = graph.GetNextAdapters(adapter, nexts);
      if (degree <= 1) {
        continue;
      }
      for (int j = 0; j < degree; ++j) {
        depths[j] = nexts[j].GetAvgDepth();
        total_depth += depths[j];
      }
      for (int j = 0; j < degree; ++j) {
        if (depths[j] <= local_ratio * total_depth) {
          num_disconnected += nexts[j].SetToDisconnect();
        }
      }
    }
  }
  const double detect_time = profile ? omp_get_wtime() - detect_begin : 0;
  const double refresh_begin = profile ? omp_get_wtime() : 0;
  if (num_disconnected != 0) {
    graph.Refresh(false);
  }
  if (profile) {
    xinfo("Weak-link profile: active-before={}, active-after={}, "
          "disconnected={}, refresh-calls={}, detect={.4}, refresh={.4}\n",
          active_before, graph.size(), num_disconnected,
          num_disconnected != 0 ? 1 : 0, detect_time,
          omp_get_wtime() - refresh_begin);
  }
  return num_disconnected;
}
