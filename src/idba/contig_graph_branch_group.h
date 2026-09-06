/**
 * @file contig_graph_branch_group.h
 * @brief ContigGraphBranchGroup Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2011-12-27
 */

#ifndef __GRAPH_CONTIG_GRAPH_BRANCH_GROUP_H_

#define __GRAPH_CONTIG_GRAPH_BRANCH_GROUP_H_

#include <vector>

#include "idba/contig_graph.h"
#include "idba/contig_graph_path.h"

/**
 * @brief It is used to contain a branch group in ContigGraph.
 */
class ContigGraphBranchGroup {
 public:
  ContigGraphBranchGroup()
      : contig_graph_(nullptr),
        max_branches_(0),
        max_length_(0),
        active_branches_(0) {}

  ContigGraphBranchGroup(ContigGraph *graph, ContigGraphVertexAdaptor begin,
                         int max_branches = 2, int max_length = 0)
      : active_branches_(0) {
    Reset(graph, begin, max_branches, max_length);
  }

  void Reset(ContigGraph *graph, ContigGraphVertexAdaptor begin,
             int max_branches = 2, int max_length = 0) {
    contig_graph_ = graph;
    begin_ = begin;
    max_branches_ = max_branches;
    max_length_ = max_length;
    if (max_length_ == 0) max_length_ = 2 * contig_graph_->kmer_size() + 2;
    for (ContigGraphPath &path : branches_) path.clear();
    prefix_workspace_.clear();
    active_branches_ = 0;
  }

  bool Search();
  void Merge();

  ContigGraphVertexAdaptor begin() { return begin_; }
  ContigGraphVertexAdaptor end() { return end_; }

 private:
  ContigGraph *contig_graph_;
  ContigGraphVertexAdaptor begin_;
  ContigGraphVertexAdaptor end_;
  std::vector<ContigGraphPath> branches_;
  ContigGraphPath prefix_workspace_;
  int max_branches_;
  int max_length_;
  size_t active_branches_;
};

#endif
