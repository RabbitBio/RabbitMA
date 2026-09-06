/**
 * @file contig_graph.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "contig_graph.h"

#include <cmath>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "contig_graph_branch_group.h"
#include "sequence.h"

using namespace std;

void ContigGraph::Initialize(deque<Sequence> &contigs,
                             deque<ContigInfo> &contig_infos) {
  pending_dead_vertices_ = false;
  vertices_.clear();
  vertices_.resize(contigs.size());

  for (int64_t i = 0; i < (int64_t)contigs.size(); ++i) {
    vertices_[i].clear();
    vertices_[i].take_contig(contigs[i]);
    vertices_[i].take_contig_info(contig_infos[i]);
    vertices_[i].set_id(i);
  }
  RefreshEdges();
}

void ContigGraph::Initialize(std::vector<ContigGraphVertex> &vertices) {
  pending_dead_vertices_ = false;
  vertices_.swap(vertices);
  vertices.clear();
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    vertices_[i].set_id(i);
  }
  RefreshEdges();
}

void ContigGraph::InitializeWithAdjacency(
    std::vector<ContigGraphVertex> &vertices,
    std::vector<uint32_t> &neighbor_codes, uint64_t num_edges) {
  vertices_.swap(vertices);
  vertices.clear();
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    vertices_[i].set_id(i);
  }
  neighbor_codes_.swap(neighbor_codes);
  begin_kmer_map_.clear();
  num_edges_ = num_edges;
  endpoint_keys_unique_ = true;
  explicit_adjacency_ = true;
  pending_dead_vertices_ = false;
  tip_candidates_valid_ = false;
  bubble_candidates_valid_ = false;
  coverage_candidates_valid_ = false;
}

void ContigGraph::Refresh() {
  RefreshVertices();
  RefreshEdges();
}

void ContigGraph::RefreshVertices() {
  tip_candidates_valid_ = false;
  bubble_candidates_valid_ = false;
  coverage_candidates_valid_ = false;
  uint64_t index = 0;
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    if (!vertices_[i].status().IsDead()) {
      vertices_[index].swap(vertices_[i]);
      vertices_[index].set_id(index);
      ++index;
    }
  }
  vertices_.resize(index);
  pending_dead_vertices_ = false;
}

void ContigGraph::RefreshEdges() {
  tip_candidates_valid_ = false;
  bubble_candidates_valid_ = false;
  coverage_candidates_valid_ = false;
  explicit_adjacency_ = false;
  BuildBeginIdbaKmerMap();
  if (endpoint_keys_unique_) {
    const size_t required = vertices_.size() * 8u;
    if (neighbor_codes_.size() < required) neighbor_codes_.resize(required);
  }

  uint64_t total_degree = 0;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);
      const IdbaKmer end_kmer = current.end_kmer(kmer_size_);
      uint32_t neighbor_code = ContigGraphVertex::kNoCachedNeighbor;

      for (int x = 0; x < 4; ++x) {
        if (current.out_edges()[x]) {
          IdbaKmer kmer = end_kmer;
          kmer.ShiftAppend(x);
          ContigGraphVertexAdaptor next =
              FindVertexAdaptorByBeginIdbaKmer(kmer);
          if (next.is_null()) {
            current.out_edges().Remove(x);
          } else if (next.id() < (UINT32_MAX - 1u) / 2u) {
            neighbor_code =
                (next.id() << 1u) | uint32_t(next.is_reverse());
            if (endpoint_keys_unique_) {
              neighbor_codes_[NeighborSlot(
                  static_cast<uint32_t>(i), strand != 0, x)] = neighbor_code;
            }
          } else {
            neighbor_code = ContigGraphVertex::kUncacheableNeighbor;
          }
        }
      }

      current.next_neighbor() =
          current.out_edges().size() == 1
              ? neighbor_code
              : ContigGraphVertex::kUncacheableNeighbor;

      total_degree += current.out_edges().size();
    }

    if ((kmer_size_ & 1u) == 0u &&
        vertices_[i].contig_size() == kmer_size_ &&
        vertices_[i].contig().IsPalindrome()) {
      vertices_[i].in_edges() =
          vertices_[i].out_edges() | vertices_[i].out_edges();
      vertices_[i].out_edges() = vertices_[i].in_edges();
      vertices_[i].next_neighbor(false) =
          ContigGraphVertex::kUncacheableNeighbor;
      vertices_[i].next_neighbor(true) =
          ContigGraphVertex::kUncacheableNeighbor;
    }

  }

  num_edges_ = total_degree / 2;
  explicit_adjacency_ = endpoint_keys_unique_;
}

void ContigGraph::ClearStatus() {
  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    vertices_[i].status().clear();
}

void ContigGraph::EnsureTipCandidates() {
  if (tip_candidates_valid_) return;
  tip_candidates_.clear();
  tip_candidates_.reserve(vertices_.size() / 8u + 1u);
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;
    const uint8_t in_edges = vertices_[i].in_edges();
    const uint8_t out_edges = vertices_[i].out_edges();
    if ((in_edges == 0u || out_edges == 0u) &&
        kmlib::bit::Popcount(in_edges) + kmlib::bit::Popcount(out_edges) <= 1) {
      if ((kmer_size_ & 1u) == 0u &&
          vertices_[i].contig_size() == kmer_size_ &&
          vertices_[i].contig().IsPalindrome()) {
        continue;
      }
      tip_candidates_.push_back(i);
    }
  }
  tip_candidates_valid_ = true;
}

void ContigGraph::EnsureBubbleCandidates() {
  if (bubble_candidates_valid_) return;
  bubble_candidates_.clear();
  bubble_candidates_.reserve(vertices_.size() / 8u + 1u);
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;
    if (vertices_[i].contig_size() <= kmer_size_) continue;
    if (kmlib::bit::Popcount(static_cast<uint8_t>(vertices_[i].out_edges())) >
        1) {
      bubble_candidates_.push_back(i << 1u);
    }
    if (kmlib::bit::Popcount(static_cast<uint8_t>(vertices_[i].in_edges())) >
        1) {
      bubble_candidates_.push_back((i << 1u) | 1u);
    }
  }
  bubble_candidates_valid_ = true;
}

void ContigGraph::EnsureCoverageCandidates(int min_length) {
  if (coverage_candidates_valid_ &&
      coverage_candidate_min_length_ == min_length) {
    return;
  }
  coverage_candidates_.clear();
  coverage_candidates_.reserve(vertices_.size() / 4u + 1u);
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;
    const uint8_t in_edges = vertices_[i].in_edges();
    const uint8_t out_edges = vertices_[i].out_edges();
    if (vertices_[i].contig_size() < min_length + kmer_size_ - 1 &&
        ((kmlib::bit::Popcount(in_edges) <= 1 &&
          kmlib::bit::Popcount(out_edges) <= 1) ||
         in_edges == 0u || out_edges == 0u)) {
      coverage_candidates_.push_back(i);
    }
  }
  coverage_candidate_min_length_ = min_length;
  coverage_candidates_valid_ = true;
}

int64_t ContigGraph::Trim(int min_length, int *next_tip_length) {
  std::vector<uint32_t> &newly_dead = trim_dead_workspace_;
  newly_dead.clear();
  if (next_tip_length != nullptr) *next_tip_length = INT_MAX;
  EnsureTipCandidates();

  for (uint32_t i : tip_candidates_) {
    if (vertices_[i].status().IsDead()) continue;
    const int num_kmers = static_cast<int>(vertices_[i].contig_size()) -
                          static_cast<int>(kmer_size_) + 1;
    if (num_kmers < min_length) {
      vertices_[i].status().SetDeadFlag();
      newly_dead.push_back(i);
    } else if (next_tip_length != nullptr) {
      *next_tip_length = std::min(*next_tip_length, num_kmers + 1);
    }
  }
  // Refresh performs two endpoint-map rebuilds (directly and again after
  // MergeSimplePaths).  With hundreds of thousands of tiny local graphs most
  // threshold probes remove nothing, in which case every observable graph
  // field is already unchanged.
  if (newly_dead.empty()) return 0;
  // The endpoint map and neighbor cache still describe the pre-deletion
  // graph.  Use them to cut only edges incident to marked vertices, then
  // compress the surviving graph directly.  Refresh()+MergeSimplePaths()
  // rebuilt the complete endpoint table both before and after compression.
  const bool use_incremental_disconnect =
      explicit_adjacency_ &&
      std::getenv("MEGAHIT_DISABLE_LOCAL_DELTA_TIP") == nullptr;
  if (use_incremental_disconnect) {
    const bool merge_exposed = DisconnectDeadEdgesIncremental(newly_dead);
    pending_dead_vertices_ = true;
    if (merge_exposed) MergeSimplePaths();
  } else {
    DisconnectDeadEdges();
    MergeSimplePaths();
  }

  // RemoveDeadEnd only needs a non-zero value to preserve its historical
  // second pass.  The public caller ignores the exact compaction count, while
  // biological deletion and traversal order are unchanged.
  return static_cast<int64_t>(newly_dead.size());
}

int64_t ContigGraph::RemoveDeadEnd(int min_length) {
  uint64_t num_deadend = 0;
  int l = std::min(2, min_length);
  int64_t last_removed = 0;
  while (true) {
    int next_tip_length = INT_MAX;
    last_removed = Trim(l, &next_tip_length);
    num_deadend += last_removed;

    if (l == min_length) break;

    int next_l = min(2 * l, min_length);
    if (last_removed == 0) {
      // With unchanged topology, a tip whose k-mer length is not yet below
      // the threshold cannot be removed by any intervening probe.  Jump to
      // the first historical power-of-two threshold that can affect one.
      if (next_tip_length > min_length) break;
      while (next_l < next_tip_length && next_l < min_length) {
        next_l = min(2 * next_l, min_length);
      }
    }
    l = next_l;
  }
  // The historical second pass exposes tips created by the first pass.  If
  // the first pass changed nothing, the second pass is provably identical.
  if (last_removed != 0) num_deadend += Trim(min_length);
  if (pending_dead_vertices_) MergeSimplePaths();
  return num_deadend;
}

int64_t ContigGraph::RemoveBubble() {
  std::vector<uint32_t> &candidates = bubble_verified_candidates_;
  candidates.clear();
  EnsureBubbleCandidates();
  // RemoveBubble is endpoint-local and never recursive. Persist the branch
  // path tapes per worker so hundreds of thousands of k rounds reset vector
  // sizes instead of allocating the same small paths again.
  static thread_local ContigGraphBranchGroup branch_group;
  static thread_local ContigGraphBranchGroup rev_branch_group;

  for (uint32_t code : bubble_candidates_) {
    ContigGraphVertexAdaptor current(&vertices_[code >> 1u], code & 1u);
    branch_group.Reset(this, current, 4, kmer_size_ + 2);

    if (branch_group.Search()) {
      ContigGraphVertexAdaptor begin = branch_group.begin();
      ContigGraphVertexAdaptor end = branch_group.end();

      begin.ReverseComplement();
      end.ReverseComplement();
      std::swap(begin, end);
      rev_branch_group.Reset(this, begin, 4, kmer_size_ + 2);

      if (rev_branch_group.Search() && rev_branch_group.end() == end) {
        candidates.push_back(code);
      }
    }
  }

  int64_t bubble = 0;
  for (uint32_t code : candidates) {
    ContigGraphVertexAdaptor current(&vertices_[code >> 1u], code & 1u);

    if (current.out_edges().size() > 1) {
      branch_group.Reset(this, current, 4, kmer_size_ + 2);

      if (branch_group.Search()) {
        ContigGraphVertexAdaptor begin = branch_group.begin();
        ContigGraphVertexAdaptor end = branch_group.end();

        begin.ReverseComplement();
        end.ReverseComplement();
        std::swap(begin, end);
        rev_branch_group.Reset(this, begin, 4, kmer_size_ + 2);

        if (rev_branch_group.Search() && rev_branch_group.end() == end) {
          branch_group.Merge();
          ++bubble;
        }
      }
    }
  }

  if (bubble != 0) {
    // BranchGroup::Merge has completed the historical whole-round mutation
    // before this barrier. Rebuilding the endpoint map in Refresh() is
    // redundant: MergeSimplePaths immediately compresses the same surviving
    // graph and rebuilds/remaps topology again. Cut edges incident to marked
    // vertices in the still-valid pre-compaction adjacency and compress once.
    MergeSimplePaths();
  }

  return bubble;
}

double ContigGraph::IterateCoverage(int min_length, double min_cover,
                                    double max_cover, double factor,
                                    bool *changed_out) {
  bool any_changed = false;
  min_cover = min(min_cover, max_cover);
  while (true) {
    double next_candidate = std::numeric_limits<double>::infinity();
    const bool changed =
        RemoveLowCoverage(min_cover, min_length, &next_candidate);
    any_changed = any_changed || changed;
    min_cover *= factor;
    if (!changed) {
      // With unchanged topology every candidate retains its coverage.  Replay
      // the exact geometric threshold sequence, but omit scans until the
      // strict `coverage < threshold` predicate can first become true.
      if (!std::isfinite(next_candidate)) {
        while (min_cover < max_cover) min_cover *= factor;
        break;
      }
      while (min_cover <= next_candidate && min_cover < max_cover) {
        min_cover *= factor;
      }
    }
    if (min_cover >= max_cover) break;
  }
  if (changed_out != nullptr) *changed_out = any_changed;
  return min_cover;
}

bool ContigGraph::RemoveLowCoverage(double min_cover, int min_length,
                                    double *next_coverage) {
  bool is_changed = false;
  if (next_coverage != nullptr) {
    *next_coverage = std::numeric_limits<double>::infinity();
  }
  EnsureCoverageCandidates(min_length);

  for (uint32_t i : coverage_candidates_) {
    ContigGraphVertexAdaptor current(&vertices_[i]);
    const double coverage = current.coverage();
    if (coverage < min_cover) {
      is_changed = true;
      current.status().SetDeadFlag();
    } else if (next_coverage != nullptr) {
      *next_coverage = std::min(*next_coverage, coverage);
    }
  }

  if (is_changed) {
    // Keep the original batch-delete barrier, but avoid RefreshEdges followed
    // immediately by another complete rebuild in MergeSimplePaths. This is
    // the same exact topology-preserving path used by Trim().
    DisconnectDeadEdges();
    MergeSimplePaths();
  }

  return is_changed;
}

void ContigGraph::MergeSimplePaths() {
  tip_candidates_valid_ = false;
  bubble_candidates_valid_ = false;
  coverage_candidates_valid_ = false;
  const bool can_remap_topology = explicit_adjacency_;
  bool remap_topology = can_remap_topology;
  constexpr uint32_t kComponentMarker = uint32_t{1} << 31u;
  std::vector<uint32_t> &old_terminal = merge_terminals_workspace_;
  old_terminal.clear();
  if (can_remap_topology) {
    old_terminal.reserve(vertices_.size() * 2u);
  }

  std::vector<ContigGraphVertex> &merged_vertices = merge_vertices_workspace_;
  merged_vertices.clear();
  merged_vertices.reserve(vertices_.size());
  std::vector<uint32_t> &payload_sources =
      merge_payload_sources_workspace_;
  payload_sources.clear();
  payload_sources.reserve(vertices_.size());
  // Preserve Assemble()'s historical palindrome-first output order.  A DNA
  // reverse-complement palindrome has even length, so the normal odd-k local
  // path skips this entire full-vertex sequence scan.
  if ((kmer_size_ & 1u) == 0u) {
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
      if (vertices_[i].status().IsDead()) continue;
      if (vertices_[i].contig_size() == kmer_size_ &&
          vertices_[i].contig().IsPalindrome()) {
        vertices_[i].status().Lock(1);

        merged_vertices.emplace_back(vertices_[i].contig(),
                                     vertices_[i].contig_info());
        payload_sources.push_back(UINT32_MAX);
        merged_vertices.back().set_id(
            static_cast<uint32_t>(merged_vertices.size() - 1u));
        if (can_remap_topology) remap_topology = false;
      }
    }
  }

  ContigGraphPath &path = merge_path_workspace_;
  path.clear();
  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;
    if (!vertices_[i].status().Lock(0)) continue;

    path.clear();
    path.Append(ContigGraphVertexAdaptor(&vertices_[i]), 0);
    bool failed = false;

    for (int strand = 0; strand < 2 && !failed; ++strand) {
      while (true) {
        ContigGraphVertexAdaptor current = path.back();
        ContigGraphVertexAdaptor next;

        if (!GetNextVertexAdaptor(current, next)) break;
        if (IsPalindromeLoop(path, next)) break;
        if (IsLoop(path, next) || !next.status().LockPreempt(0)) {
          failed = true;
          break;
        }

        path.Append(next, -kmer_size_ + 1);
      }

      if (!failed) path.ReverseComplement();
    }

    if (!failed) {
      if (path.num_nodes() == 1u && !path[0].is_reverse()) {
        // The overwhelming common case after a sparse cleaning mutation is
        // an unchanged, already-compressed unitig.  Defer moving its payload
        // until component discovery has finished, because later walks still
        // consult the old vertex's degree/status fields.
        merged_vertices.emplace_back();
        payload_sources.push_back(path[0].id());
      } else {
        Sequence contig;
        ContigInfo contig_info;
        path.Assemble(contig, contig_info);
        merged_vertices.emplace_back();
        merged_vertices.back().take_contig(contig);
        merged_vertices.back().take_contig_info(contig_info);
        payload_sources.push_back(UINT32_MAX);
      }
      merged_vertices.back().set_id(
          static_cast<uint32_t>(merged_vertices.size() - 1u));

      if (can_remap_topology) {
        const uint32_t new_id =
            static_cast<uint32_t>(merged_vertices.size() - 1u);
        for (uint32_t j = 0; j < path.num_nodes(); ++j) {
          const uint32_t old_code =
              (path[j].id() << 1u) | uint32_t(path[j].is_reverse());
          const uint32_t new_code = new_id << 1u;
          uint32_t &forward_map =
              vertices_[old_code >> 1u].next_neighbor(old_code & 1u);
          uint32_t &reverse_map =
              vertices_[old_code >> 1u].next_neighbor((old_code & 1u) ^ 1u);
          if ((forward_map < ContigGraphVertex::kUncacheableNeighbor &&
               (forward_map & kComponentMarker) != 0u &&
               (forward_map & ~kComponentMarker) != new_code) ||
              (reverse_map < ContigGraphVertex::kUncacheableNeighbor &&
               (reverse_map & kComponentMarker) != 0u &&
               (reverse_map & ~kComponentMarker) != (new_code | 1u))) {
            remap_topology = false;
          }
          forward_map = kComponentMarker | new_code;
          reverse_map = kComponentMarker | new_code | 1u;
        }
        const uint32_t forward_terminal =
            (path.back().id() << 1u) | uint32_t(path.back().is_reverse());
        const uint32_t reverse_terminal =
            ((path.front().id() << 1u) |
             uint32_t(path.front().is_reverse())) ^
            1u;
        old_terminal.push_back(forward_terminal);
        old_terminal.push_back(reverse_terminal);
      }
    } else if (can_remap_topology) {
      remap_topology = false;
    }
  }

  assert(payload_sources.size() == merged_vertices.size());
  for (uint32_t new_id = 0; new_id < merged_vertices.size(); ++new_id) {
    const uint32_t source = payload_sources[new_id];
    if (source != UINT32_MAX) {
      merged_vertices[new_id].take_payload(vertices_[source]);
      merged_vertices[new_id].set_id(new_id);
    }
  }

  std::vector<uint32_t> &merged_neighbors = merge_neighbors_workspace_;
  uint64_t merged_degree = 0;
  if (remap_topology && old_terminal.size() == merged_vertices.size() * 2u) {
    const size_t required_neighbors = merged_vertices.size() * 8u;
    if (merged_neighbors.size() < required_neighbors) {
      merged_neighbors.resize(required_neighbors);
    }
    for (uint32_t i = 0; i < merged_vertices.size() && remap_topology; ++i) {
      for (uint32_t strand = 0; strand < 2u; ++strand) {
        ContigGraphVertexAdaptor current(&merged_vertices[i], strand != 0);
        const uint32_t terminal = old_terminal[i * 2u + strand];
        uint32_t sole_neighbor = ContigGraphVertex::kNoCachedNeighbor;
        uint8_t edges = static_cast<uint8_t>(current.out_edges());
        while (edges != 0u) {
          const uint32_t x = static_cast<uint32_t>(__builtin_ctz(edges));
          edges &= static_cast<uint8_t>(edges - 1u);
          const uint32_t old_neighbor =
              neighbor_codes_[NeighborSlot(terminal >> 1u, terminal & 1u, x)];
          if (old_neighbor >= ContigGraphVertex::kUncacheableNeighbor ||
              (old_neighbor >> 1u) >= vertices_.size()) {
            remap_topology = false;
            break;
          }
          const uint32_t marked_neighbor =
              vertices_[old_neighbor >> 1u].next_neighbor(old_neighbor & 1u);
          if (marked_neighbor >= ContigGraphVertex::kUncacheableNeighbor ||
              (marked_neighbor & kComponentMarker) == 0u) {
            remap_topology = false;
            break;
          }
          const uint32_t new_neighbor = marked_neighbor & ~kComponentMarker;
          merged_neighbors[NeighborSlot(i, strand != 0, x)] = new_neighbor;
          sole_neighbor = new_neighbor;
        }
        current.next_neighbor() =
            current.out_edges().size() == 1
                ? sole_neighbor
                : ContigGraphVertex::kUncacheableNeighbor;
        merged_degree += current.out_edges().size();
      }
    }
  } else {
    remap_topology = false;
  }

  vertices_.swap(merged_vertices);
  if (remap_topology) {
    neighbor_codes_.swap(merged_neighbors);
    num_edges_ = merged_degree / 2u;
    endpoint_keys_unique_ = true;
    explicit_adjacency_ = true;
    // The old endpoint table contains stale vertex IDs.  No query uses it
    // while explicit adjacency is active; clearing keys also releases their
    // owned k-mer state before the next graph round reuses the capacity.
    begin_kmer_map_.clear();
  } else {
    explicit_adjacency_ = false;
    RefreshEdges();
  }
  pending_dead_vertices_ = false;
}

bool ContigGraph::DisconnectDeadEdgesIncremental(
    const std::vector<uint32_t> &newly_dead) {
  assert(explicit_adjacency_);
  std::vector<uint32_t> &touched = touched_vertex_workspace_;
  touched.clear();
  touched.reserve(newly_dead.size() * 2u);

  // Every physical edge is represented in both orientations.  Enumerating
  // the two oriented edge masks of newly dead vertices therefore identifies
  // every live vertex whose reciprocal mask can have changed.
  for (uint32_t dead_id : newly_dead) {
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      const ContigGraphVertexAdaptor dead(
          &vertices_[dead_id], strand != 0u);
      uint8_t edges = static_cast<uint8_t>(dead.out_edges());
      while (edges != 0u) {
        const uint32_t base = static_cast<uint32_t>(__builtin_ctz(edges));
        edges &= static_cast<uint8_t>(edges - 1u);
        const uint32_t code =
            neighbor_codes_[NeighborSlot(dead_id, strand != 0u, base)];
        if (code < ContigGraphVertex::kUncacheableNeighbor &&
            (code >> 1u) < vertices_.size() &&
            !vertices_[code >> 1u].status().IsDead()) {
          touched.push_back(code >> 1u);
        }
      }
    }
  }
  std::sort(touched.begin(), touched.end());
  touched.erase(std::unique(touched.begin(), touched.end()), touched.end());

  for (uint32_t vertex_id : touched) {
    ContigGraphVertex &vertex = vertices_[vertex_id];
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      BitEdges &edge_set = strand == 0u ? vertex.out_edges()
                                        : vertex.in_edges();
      uint8_t live_edges = static_cast<uint8_t>(edge_set);
      uint8_t scan = live_edges;
      uint32_t sole_neighbor = ContigGraphVertex::kNoCachedNeighbor;
      while (scan != 0u) {
        const uint32_t base = static_cast<uint32_t>(__builtin_ctz(scan));
        scan &= static_cast<uint8_t>(scan - 1u);
        const size_t slot = NeighborSlot(vertex_id, strand != 0u, base);
        const uint32_t code = neighbor_codes_[slot];
        if (code >= ContigGraphVertex::kUncacheableNeighbor ||
            (code >> 1u) >= vertices_.size() ||
            vertices_[code >> 1u].status().IsDead()) {
          live_edges &= static_cast<uint8_t>(~(uint8_t{1} << base));
          neighbor_codes_[slot] = ContigGraphVertex::kNoCachedNeighbor;
        } else {
          sole_neighbor = code;
        }
      }
      edge_set = live_edges;
      vertex.next_neighbor(strand != 0u) =
          live_edges != 0u && (live_edges & (live_edges - 1u)) == 0u
              ? sole_neighbor
              : ContigGraphVertex::kUncacheableNeighbor;
    }
  }

  tip_candidates_valid_ = false;
  bubble_candidates_valid_ = false;
  coverage_candidates_valid_ = false;

  // The preceding graph was already maximally unitig-compressed.  Therefore
  // a new simple link can only touch a vertex adjacent to this deletion.
  // When none exists, stable compaction and sequence materialization can be
  // deferred across geometric thresholds with no semantic effect.
  for (uint32_t vertex_id : touched) {
    if (vertices_[vertex_id].status().IsDead()) continue;
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[vertex_id], strand != 0u);
      ContigGraphVertexAdaptor next;
      if (GetNextVertexAdaptor(current, next) &&
          !next.status().IsDead()) {
        return true;
      }
    }
  }
  return false;
}

void ContigGraph::DisconnectDeadEdges() {
  if (explicit_adjacency_) {
    uint64_t total_degree = 0;
    for (uint32_t i = 0; i < vertices_.size(); ++i) {
      ContigGraphVertex &vertex = vertices_[i];
      if (vertex.status().IsDead()) continue;

      for (uint32_t strand = 0; strand < 2u; ++strand) {
        BitEdges &edge_set =
            strand == 0u ? vertex.out_edges() : vertex.in_edges();
        uint8_t live_edges = static_cast<uint8_t>(edge_set);
        uint8_t scan = live_edges;
        uint32_t sole_neighbor = ContigGraphVertex::kNoCachedNeighbor;
        while (scan != 0u) {
          const uint32_t base = static_cast<uint32_t>(__builtin_ctz(scan));
          scan &= static_cast<uint8_t>(scan - 1u);
          const uint32_t code =
              neighbor_codes_[NeighborSlot(i, strand != 0u, base)];
          if (code >= ContigGraphVertex::kUncacheableNeighbor ||
              (code >> 1u) >= vertices_.size() ||
              vertices_[code >> 1u].status().IsDead()) {
            live_edges &= static_cast<uint8_t>(~(uint8_t{1} << base));
          } else {
            sole_neighbor = code;
          }
        }
        edge_set = live_edges;
        vertex.next_neighbor(strand != 0u) =
            live_edges != 0u && (live_edges & (live_edges - 1u)) == 0u
                ? sole_neighbor
                : ContigGraphVertex::kUncacheableNeighbor;
        total_degree += kmlib::bit::Popcount(live_edges);
      }

      if ((kmer_size_ & 1u) == 0u &&
          vertex.contig_size() == kmer_size_ &&
          vertex.contig().IsPalindrome()) {
        vertex.in_edges() = vertex.out_edges() | vertex.out_edges();
        vertex.out_edges() = vertex.in_edges();
        vertex.next_neighbor(false) =
            ContigGraphVertex::kUncacheableNeighbor;
        vertex.next_neighbor(true) =
            ContigGraphVertex::kUncacheableNeighbor;
      }
    }
    num_edges_ = total_degree / 2u;
    return;
  }

  uint64_t total_degree = 0;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;

    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);
      uint32_t neighbor_code = ContigGraphVertex::kNoCachedNeighbor;

      for (int x = 0; x < 4; ++x) {
        if (!current.out_edges()[x]) continue;

        ContigGraphVertexAdaptor next = GetNeighbor(current, x);
        if (next.is_null() || next.status().IsDead()) {
          current.out_edges().Remove(x);
        } else if (next.id() < (UINT32_MAX - 1u) / 2u) {
          neighbor_code =
              (next.id() << 1u) | uint32_t(next.is_reverse());
        } else {
          neighbor_code = ContigGraphVertex::kUncacheableNeighbor;
        }
      }

      current.next_neighbor() =
          current.out_edges().size() == 1
              ? neighbor_code
              : ContigGraphVertex::kUncacheableNeighbor;
      total_degree += current.out_edges().size();
    }

    if ((kmer_size_ & 1u) == 0u &&
        vertices_[i].contig_size() == kmer_size_ &&
        vertices_[i].contig().IsPalindrome()) {
      vertices_[i].in_edges() =
          vertices_[i].out_edges() | vertices_[i].out_edges();
      vertices_[i].out_edges() = vertices_[i].in_edges();
      vertices_[i].next_neighbor(false) =
          ContigGraphVertex::kUncacheableNeighbor;
      vertices_[i].next_neighbor(true) =
          ContigGraphVertex::kUncacheableNeighbor;
    }
  }

  num_edges_ = total_degree / 2;
}

int64_t ContigGraph::Assemble(vector<Sequence> &contigs,
                              vector<ContigInfo> &contig_infos) {
  contigs.clear();
  contig_infos.clear();
  contigs.reserve(vertices_.size());
  contig_infos.reserve(vertices_.size());

  if ((kmer_size_ & 1u) == 0u) {
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
      if (vertices_[i].contig_size() == kmer_size_ &&
          vertices_[i].contig().IsPalindrome()) {
        vertices_[i].status().Lock(1);

        Sequence contig = vertices_[i].contig();
        ContigInfo contig_info;
        contig_info.set_kmer_count(vertices_[i].kmer_count());
        contig_info.in_edges() = vertices_[i].in_edges();
        contig_info.out_edges() = vertices_[i].out_edges();

        contigs.push_back(contig);
        contig_infos.push_back(contig_info);
      }
    }
  }

  ContigGraphPath &path = assemble_path_workspace_;
  path.clear();
  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    if (vertices_[i].status().IsDead()) continue;
    if (!vertices_[i].status().Lock(0)) continue;

    path.clear();
    path.Append(ContigGraphVertexAdaptor(&vertices_[i]), 0);

    Sequence contig;
    ContigInfo contig_info;
    for (int strand = 0; strand < 2; ++strand) {
      while (true) {
        ContigGraphVertexAdaptor current = path.back();
        ContigGraphVertexAdaptor next;

        if (!GetNextVertexAdaptor(current, next)) break;

        if (IsPalindromeLoop(path, next)) break;

        if (IsLoop(path, next)) goto FAIL;

        if (!next.status().LockPreempt(0)) goto FAIL;

        path.Append(next, -kmer_size_ + 1);
      }

      path.ReverseComplement();
    }

    path.Assemble(contig, contig_info);
    contigs.emplace_back();
    contig_infos.emplace_back();
    contigs.back().swap(contig);
    contig_infos.back().swap(contig_info);
  FAIL:;
  }

  return contigs.size();
}

struct SearchNode {
  ContigGraphVertexAdaptor node;
  int distance;
  int label;
};

void ContigGraph::BuildBeginIdbaKmerMap() {
  begin_kmer_map_.clear();
  begin_kmer_map_.reserve(vertices_.size() * 2u);
  endpoint_keys_unique_ = true;

  for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i) {
    for (int strand = 0; strand < 2; ++strand) {
      ContigGraphVertexAdaptor current(&vertices_[i], strand);
      IdbaKmer kmer = current.begin_kmer(kmer_size_);

      IdbaKmer key = kmer.unique_format();
      const uint64_t endpoint_bit = uint64_t{1} << strand;
      const uint64_t canonical_bit =
          (kmer == key ? uint64_t{1} : uint64_t{0}) << (2u + strand);
      const uint64_t code = (uint64_t(i) << 4u) | endpoint_bit | canonical_bit;
      uint64_t *existing = begin_kmer_map_.find(key);
      if (existing == nullptr) {
        begin_kmer_map_.insert(key, code);
      } else if ((*existing >> 4u) == static_cast<uint64_t>(i)) {
        // Both oriented endpoints of a k-long/palindromic vertex may share a
        // canonical key.  Preserve both orientations while retaining the
        // historical strand-0-first lookup rule.
        *existing |= endpoint_bit | canonical_bit;
      } else {
        // Preserve operator[] assignment semantics: a later vertex with the
        // same endpoint key replaces the earlier one.
        endpoint_keys_unique_ = false;
        *existing = code;
      }
    }
  }
}

bool ContigGraph::CycleDetect(ContigGraphVertexAdaptor current,
                              map<int, int> &status) {
  if (status[current.id()] == 0) {
    bool flag = false;
    status[current.id()] = 1;
    for (int x = 0; x < 4; ++x) {
      if (current.out_edges()[x]) {
        if (CycleDetect(GetNeighbor(current, x), status)) flag = true;
      }
    }
    status[current.id()] = 2;
    return flag;
  } else if (status[current.id()] == 1)
    return true;
  else
    return false;
}

void ContigGraph::TopSortDFS(deque<ContigGraphVertexAdaptor> &order,
                             ContigGraphVertexAdaptor current,
                             map<int, int> &status) {
  if (status[current.id()] == 0) {
    status[current.id()] = 1;
    for (int x = 0; x < 4; ++x) {
      if (current.out_edges()[x])
        TopSortDFS(order, GetNeighbor(current, x), status);
    }
    order.push_back(current);
  }
}

int ContigGraph::GetDepth(ContigGraphVertexAdaptor current, int depth,
                          int &maximum, int min_length) {
  if (depth > maximum) maximum = depth;

  if (maximum >= min_length) return min_length;

  deque<ContigGraphVertexAdaptor> neighbors;
  GetNeighbors(current, neighbors);
  for (unsigned i = 0; i < neighbors.size(); ++i) {
    if (neighbors[i].status().IsDead()) continue;

    GetDepth(neighbors[i], depth - kmer_size_ + 1 + neighbors[i].contig_size(),
             maximum, min_length);
  }

  return min(maximum, min_length);
}
