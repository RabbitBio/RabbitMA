/**
 * @file contig_graph.h
 * @brief ContigGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-16
 */

#ifndef __GRAPH_CONTIG_GRAPH_H_

#define __GRAPH_CONTIG_GRAPH_H_

#include <algorithm>
#include <deque>
#include <map>
#include <vector>
#include "bit_operation.h"
#include "idba/compact_endpoint_map.h"
#include "idba/contig_graph_path.h"
#include "idba/contig_graph_vertex.h"
#include "idba/contig_info.h"
#include "idba/hash.h"
#include "idba/hash_graph.h"
#include "idba/kmer.h"
#include "idba/sequence.h"

/**
 * @brief It is compact version de Bruijn graph in which each vertex is a contig
 * and each edge between contigs means they are connected in de Bruijn graph.
 */
class ContigGraph {
 public:
  explicit ContigGraph(uint32_t kmer_size = 0)
      : num_edges_(0),
        kmer_size_(kmer_size),
        endpoint_keys_unique_(false),
        explicit_adjacency_(false) {}

  ~ContigGraph() { clear(); }

  void Initialize(std::deque<Sequence> &contigs,
                  std::deque<ContigInfo> &contig_infos);
  void Initialize(std::vector<ContigGraphVertex> &vertices);
  void InitializeWithAdjacency(std::vector<ContigGraphVertex> &vertices,
                               std::vector<uint32_t> &neighbor_codes,
                               uint64_t num_edges);

  void Refresh();
  void RefreshVertices();
  void RefreshEdges();

  void AddEdge(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to) {
    tip_candidates_valid_ = false;
    bubble_candidates_valid_ = false;
    coverage_candidates_valid_ = false;
    const uint8_t forward_base = to.get_base(kmer_size_ - 1u);
    from.out_edges().Add(forward_base);
    if (explicit_adjacency_) {
      neighbor_codes_[NeighborSlot(from.id(), from.is_reverse(), forward_base)] =
          EncodeEndpoint(to.id(), to.is_reverse());
    }
    from.ReverseComplement();
    to.ReverseComplement();
    std::swap(from, to);
    const uint8_t reverse_base = to.get_base(kmer_size_ - 1u);
    from.out_edges().Add(reverse_base);
    if (explicit_adjacency_) {
      neighbor_codes_[NeighborSlot(from.id(), from.is_reverse(), reverse_base)] =
          EncodeEndpoint(to.id(), to.is_reverse());
    }
  }

  void RemoveEdge(ContigGraphVertexAdaptor current, int x) {
    tip_candidates_valid_ = false;
    bubble_candidates_valid_ = false;
    coverage_candidates_valid_ = false;
    current.out_edges().Remove(x);
    ContigGraphVertexAdaptor next = GetNeighbor(current, x);
    if (explicit_adjacency_) {
      neighbor_codes_[NeighborSlot(current.id(), current.is_reverse(), x)] =
          ContigGraphVertex::kNoCachedNeighbor;
    }
    next.ReverseComplement();
    const uint8_t reverse_base = 3u - current.get_base(0);
    next.out_edges().Remove(reverse_base);
    if (explicit_adjacency_) {
      neighbor_codes_[NeighborSlot(next.id(), next.is_reverse(), reverse_base)] =
          ContigGraphVertex::kNoCachedNeighbor;
    }
  }

  void ClearStatus();

  void MergeSimplePaths();

  int64_t Trim(int min_length, int *next_tip_length = nullptr);

  int64_t RemoveDeadEnd(int min_length);
  int64_t RemoveBubble();

  double IterateCoverage(int min_length, double min_cover, double max_cover,
                         double factor = 1.1,
                         bool *changed_out = nullptr);

  bool RemoveLowCoverage(double min_cover, int min_length,
                         double *next_coverage = nullptr);

  int64_t Assemble(std::vector<Sequence> &contigs,
                   std::vector<ContigInfo> &contig_infos);

  ContigGraphVertexAdaptor GetNeighbor(const ContigGraphVertexAdaptor &current,
                                       int x) {
    if (explicit_adjacency_) {
      const uint32_t code =
          neighbor_codes_[NeighborSlot(current.id(), current.is_reverse(), x)];
      if (code < ContigGraphVertex::kUncacheableNeighbor) {
        return ContigGraphVertexAdaptor(&vertices_[code >> 1u], code & 1u);
      }
      return ContigGraphVertexAdaptor();
    }
    if (current.out_edges().size() == 1 && current.out_edges()[x]) {
      const uint32_t cached = current.next_neighbor();
      if (cached < ContigGraphVertex::kUncacheableNeighbor) {
        return ContigGraphVertexAdaptor(&vertices_[cached >> 1u], cached & 1u);
      }
    }
    IdbaKmer kmer = current.end_kmer(kmer_size_);
    kmer.ShiftAppend(x);
    return FindVertexAdaptorByBeginIdbaKmer(kmer);
  }

  void GetNeighbors(const ContigGraphVertexAdaptor &current,
                    std::deque<ContigGraphVertexAdaptor> &neighbors) {
    neighbors.clear();
    for (int x = 0; x < 4; ++x) {
      if (current.out_edges()[x]) neighbors.push_back(GetNeighbor(current, x));
    }
  }

  void swap(ContigGraph &contig_graph) {
    begin_kmer_map_.swap(contig_graph.begin_kmer_map_);
    vertices_.swap(contig_graph.vertices_);
    neighbor_codes_.swap(contig_graph.neighbor_codes_);
    tip_candidates_.swap(contig_graph.tip_candidates_);
    bubble_candidates_.swap(contig_graph.bubble_candidates_);
    coverage_candidates_.swap(contig_graph.coverage_candidates_);
    bubble_verified_candidates_.swap(
        contig_graph.bubble_verified_candidates_);
    merge_vertices_workspace_.swap(contig_graph.merge_vertices_workspace_);
    merge_neighbors_workspace_.swap(contig_graph.merge_neighbors_workspace_);
    merge_terminals_workspace_.swap(contig_graph.merge_terminals_workspace_);
    merge_payload_sources_workspace_.swap(
        contig_graph.merge_payload_sources_workspace_);
    merge_path_workspace_.swap(contig_graph.merge_path_workspace_);
    assemble_path_workspace_.swap(contig_graph.assemble_path_workspace_);
    std::swap(num_edges_, contig_graph.num_edges_);
    std::swap(kmer_size_, contig_graph.kmer_size_);
    std::swap(endpoint_keys_unique_, contig_graph.endpoint_keys_unique_);
    std::swap(explicit_adjacency_, contig_graph.explicit_adjacency_);
    std::swap(tip_candidates_valid_, contig_graph.tip_candidates_valid_);
    std::swap(bubble_candidates_valid_,
              contig_graph.bubble_candidates_valid_);
    std::swap(coverage_candidates_valid_,
              contig_graph.coverage_candidates_valid_);
    std::swap(coverage_candidate_min_length_,
              contig_graph.coverage_candidate_min_length_);
    std::swap(pending_dead_vertices_, contig_graph.pending_dead_vertices_);
    trim_dead_workspace_.swap(contig_graph.trim_dead_workspace_);
    touched_vertex_workspace_.swap(contig_graph.touched_vertex_workspace_);
  }

  uint32_t kmer_size() const { return kmer_size_; }
  void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

  void clear() {
    num_edges_ = 0;
    vertices_.clear();
    // Retain both capacity and logical size as a persistent per-thread
    // workspace. RefreshEdges overwrites every slot whose edge bit is live;
    // stale slots for absent edges are never queried. This avoids zeroing
    // eight uint32 entries per vertex for every endpoint/k round.
    begin_kmer_map_.clear();
    endpoint_keys_unique_ = false;
    explicit_adjacency_ = false;
    tip_candidates_.clear();
    bubble_candidates_.clear();
    coverage_candidates_.clear();
    bubble_verified_candidates_.clear();
    merge_vertices_workspace_.clear();
    // Keep the logical high-water size as well as capacity. During a remap,
    // only slots selected by the new vertex edge masks are observable and
    // every such slot is overwritten. Avoid zero-filling 8*V uint32 values
    // again for every cleaning round and endpoint.
    merge_terminals_workspace_.clear();
    merge_payload_sources_workspace_.clear();
    merge_path_workspace_.clear();
    assemble_path_workspace_.clear();
    tip_candidates_valid_ = false;
    bubble_candidates_valid_ = false;
    coverage_candidates_valid_ = false;
    pending_dead_vertices_ = false;
  }

 private:
  ContigGraph(const ContigGraph &);
  const ContigGraph &operator=(const ContigGraph &);

  void BuildBeginIdbaKmerMap();
  void DisconnectDeadEdges();
  bool DisconnectDeadEdgesIncremental(
      const std::vector<uint32_t> &newly_dead);
  void EnsureTipCandidates();
  void EnsureBubbleCandidates();
  void EnsureCoverageCandidates(int min_length);

  static uint32_t EncodeEndpoint(uint32_t vertex_id, bool strand) {
    if (vertex_id >= (UINT32_MAX - 1u) / 2u) {
      return ContigGraphVertex::kUncacheableNeighbor;
    }
    return (vertex_id << 1u) | uint32_t(strand);
  }

  static size_t NeighborSlot(uint32_t vertex_id, bool strand, int base) {
    return size_t(vertex_id) * 8u + size_t(strand) * 4u +
           static_cast<uint8_t>(base);
  }

  bool GetNextVertexAdaptor(ContigGraphVertexAdaptor &current,
                            ContigGraphVertexAdaptor &next) {
    if (current.out_edges().size() != 1) return false;

    const uint32_t cached = current.next_neighbor();
    if (cached < ContigGraphVertex::kUncacheableNeighbor) {
      next.set_vertex(&vertices_[cached >> 1u], cached & 1u);
    } else {
      next =
          GetNeighbor(current, bit_operation::BitToIndex(current.out_edges()));
    }
    return next.in_edges().size() == 1 &&
           ((kmer_size_ & 1u) != 0u || next.contig_size() != kmer_size_ ||
            !next.contig().IsPalindrome());
  }

  bool IsLoop(const ContigGraphPath &path,
              const ContigGraphVertexAdaptor &next) {
    return path.front().id() == next.id();
  }

  bool IsPalindromeLoop(const ContigGraphPath &path,
                        const ContigGraphVertexAdaptor &next) {
    return path.back().id() == next.id();
  }

  ContigGraphVertexAdaptor FindVertexAdaptorByBeginIdbaKmer(
      const IdbaKmer &begin_kmer) {
    IdbaKmer key = begin_kmer.unique_format();
    const uint64_t *found = begin_kmer_map_.find(key);
    if (found != nullptr) {
      const uint64_t code = *found;
      const uint32_t vertex_id = static_cast<uint32_t>(code >> 4u);
      const uint8_t endpoint_mask = code & 3u;
      const uint8_t canonical_mask = (code >> 2u) & 3u;
      const bool query_is_canonical = begin_kmer == key;
      for (uint8_t strand = 0; strand < 2; ++strand) {
        if ((endpoint_mask & (1u << strand)) != 0 &&
            (((canonical_mask >> strand) & 1u) != 0) ==
                query_is_canonical) {
          return ContigGraphVertexAdaptor(&vertices_[vertex_id], strand);
        }
      }
    }

    return ContigGraphVertexAdaptor();
  }

  bool CycleDetect(ContigGraphVertexAdaptor current,
                   std::map<int, int> &status);
  void TopSortDFS(std::deque<ContigGraphVertexAdaptor> &order,
                  ContigGraphVertexAdaptor current, std::map<int, int> &status);
  int GetDepth(ContigGraphVertexAdaptor current, int length, int &maximum,
               int min_length);

  CompactEndpointMap begin_kmer_map_;
  // A cleaning round only compacts this container; it never grows while
  // vertex adaptors are live.  Contiguous storage therefore preserves IDs
  // and traversal order while making the repeated tip/bubble/coverage scans
  // cache- and prefetch-friendly.
  std::vector<ContigGraphVertex> vertices_;
  // The dense neighbor array is enabled only when BuildBeginIdbaKmerMap has
  // proved that no two different vertices share a canonical endpoint key.
  // In that case endpoint lookup is a true one-to-one relation and component
  // compression can remap topology exactly without rebuilding the hash map.
  // Graphs with duplicate endpoint keys retain the historical map semantics.
  std::vector<uint32_t> neighbor_codes_;
  // Topology-only candidate set reused by the geometrically increasing tip
  // thresholds.  A threshold change cannot create a new candidate; graph
  // mutation invalidates the set and the next Trim rebuilds it exactly.
  std::vector<uint32_t> tip_candidates_;
  std::vector<uint32_t> bubble_candidates_;
  std::vector<uint32_t> coverage_candidates_;
  std::vector<uint32_t> bubble_verified_candidates_;
  // Geometric tip thresholds frequently delete leaves without making their
  // high-degree attachment points mergeable.  Keep those deletions logical
  // until a simple path is actually exposed, then perform the same stable
  // whole-graph merge once.  These tapes are persistent per-worker storage.
  std::vector<uint32_t> trim_dead_workspace_;
  std::vector<uint32_t> touched_vertex_workspace_;
  // MergeSimplePaths is called repeatedly as tip/bubble/coverage cleaning
  // progresses. Keep its destination arrays and path tape as graph-owned
  // workspace, so each refresh destroys only logical payloads and reuses the
  // worker's high-water allocations on the next round/endpoint.
  std::vector<ContigGraphVertex> merge_vertices_workspace_;
  std::vector<uint32_t> merge_neighbors_workspace_;
  std::vector<uint32_t> merge_terminals_workspace_;
  std::vector<uint32_t> merge_payload_sources_workspace_;
  ContigGraphPath merge_path_workspace_;
  ContigGraphPath assemble_path_workspace_;
  uint64_t num_edges_;
  uint32_t kmer_size_;
  bool endpoint_keys_unique_;
  bool explicit_adjacency_;
  bool tip_candidates_valid_{false};
  bool bubble_candidates_valid_{false};
  bool coverage_candidates_valid_{false};
  int coverage_candidate_min_length_{0};
  bool pending_dead_vertices_{false};

};

#endif
