/**
 * @file contig_graph_branch_group.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2011-12-27
 */

#include "contig_graph_branch_group.h"

#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;

bool ContigGraphBranchGroup::Search() {
  int kmer_size = contig_graph_->kmer_size();
  branches_.reserve(max_branches_);

  if (branches_.empty()) branches_.emplace_back();
  branches_[0].clear();
  branches_[0].Append(begin_, 0);
  active_branches_ = 1;

  if ((int)begin_.out_edges().size() <= 1 ||
      (int)begin_.out_edges().size() > max_branches_ ||
      (int)begin_.contig_size() == kmer_size)
    return false;

  bool is_converge = false;
  for (int k = 1; k < max_length_; ++k) {
    int num_branches = static_cast<int>(active_branches_);
    bool is_extend = false;
    for (int i = 0; i < num_branches; ++i) {
      if ((int)branches_[i].internal_size(kmer_size) >= max_length_) continue;

      ContigGraphVertexAdaptor current = branches_[i].back();

      if (current.out_edges().size() == 0) return false;

      bool is_first = true;
      // The old code copied the complete prefix before every one-edge
      // extension even though that copy is used only when the vertex really
      // branches.  Snapshot lazily at an actual branch and append directly to
      // the newly created path.  Branch order remains the original base order.
      if (current.out_edges().size() > 1) prefix_workspace_ = branches_[i];
      uint8_t edges = static_cast<uint8_t>(current.out_edges());
      while (edges != 0u) {
        const int x = __builtin_ctz(edges);
        edges &= static_cast<uint8_t>(edges - 1u);
          ContigGraphVertexAdaptor next =
              contig_graph_->GetNeighbor(current, x);

          if (next.status().IsDead()) return false;

          if (is_first) {
            branches_[i].Append(next, -kmer_size + 1);
            is_first = false;
          } else {
            if (static_cast<int>(active_branches_) == max_branches_)
              return false;
            if (active_branches_ == branches_.size()) {
              branches_.emplace_back();
            }
            branches_[active_branches_] = prefix_workspace_;
            branches_[active_branches_].Append(next, -kmer_size + 1);
            ++active_branches_;
          }
        is_extend = true;
      }
    }

    end_ = branches_[0].back();

    if ((int)end_.contig_size() > kmer_size) {
      is_converge = true;
      for (size_t i = 0; i < active_branches_; ++i) {
        if (branches_[i].back() != end_ ||
            (int)branches_[i].internal_size(kmer_size) != max_length_) {
          is_converge = false;
          break;
        }
      }

      if (is_converge) break;
    }

    if (!is_extend) break;
  }

  return is_converge && begin_ != end_;
}

void ContigGraphBranchGroup::Merge() {
  unsigned best = 0;
  for (size_t i = 1; i < active_branches_; ++i) {
    if (branches_[i].kmer_count() > branches_[best].kmer_count()) best = i;
  }

  for (size_t i = 0; i < active_branches_; ++i) {
    ContigGraphPath &path = branches_[i];
    path.front().out_edges() = 0;
    path.back().in_edges() = 0;
    for (unsigned j = 1; j + 1 < path.num_nodes(); ++j) {
      path[j].in_edges() = 0;
      path[j].out_edges() = 0;
      path[j].status().SetDeadFlag();
    }
  }

  ContigGraphPath &path = branches_[best];
  for (unsigned j = 1; j + 1 < path.num_nodes(); ++j)
    path[j].status().ResetDeadFlag();

  for (unsigned j = 0; j + 1 < path.num_nodes(); ++j)
    contig_graph_->AddEdge(path[j], path[j + 1]);

}
