//
// Created by vout on 11/10/18.
//

#ifndef MEGAHIT_UNITIG_GRAPH_H
#define MEGAHIT_UNITIG_GRAPH_H

#include <array>
#include <limits>
#include <memory>
#include <vector>
#include "parallel_hashmap/phmap.h"
#include "sdbg/sdbg.h"
#include "unitig_graph_vertex.h"
#include "utils/mutex.h"
#include "utils/utils.h"

class UnitigGraph {
 public:
  using Vertex = UnitigGraphVertex;
  using VertexAdapter = UnitigGraphVertex::Adapter;
  using size_type = VertexAdapter::size_type;
  static const size_type kMaxNumVertices =
      std::numeric_limits<size_type>::max() - 1;
  static const size_type kNullVertexID = kMaxNumVertices + 1;

 public:
  explicit UnitigGraph(SDBG *sdbg);
  UnitigGraph(const UnitigGraph &) = delete;
  UnitigGraph(const UnitigGraph &&) = delete;
  ~UnitigGraph() = default;
  // `vertices_` is a stable slot arena.  Cleaning removes IDs from the dense
  // active list instead of moving 40-byte vertices and invalidating every
  // graph handle.  Public iteration remains dense and in the same stable
  // order as the historical stable compaction.
  size_type size() const { return active_ids_.size(); }
  size_type active_id(size_type active_index) const {
    return active_ids_[active_index];
  }
  size_t k() const { return sdbg_->k(); }

 public:
  void Refresh(bool mark_changed = false);
  std::string VertexToDNAString(VertexAdapter adapter);

 public:
  /*
   * Function for VertexAdapter obtaining & traversal
   */
  VertexAdapter MakeVertexAdapter(size_type id, int strand = 0) {
    return adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(VertexAdapter &adapter, VertexAdapter *out) {
    return adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(VertexAdapter &adapter) {
    return adapter_impl_.OutDegree(adapter);
  }
  int InDegree(VertexAdapter &adapter) {
    return adapter_impl_.InDegree(adapter);
  }

 private:
  /*
   * Function for SudoVertexAdapter obtaining & traversal
   */
  using SudoVertexAdapter = UnitigGraphVertex::SudoAdapter;
  SudoVertexAdapter MakeSudoAdapter(size_type id, int strand = 0) {
    return sudo_adapter_impl_.MakeVertexAdapter(id, strand);
  }
  int GetNextAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetNextAdapters(adapter, out);
  }
  int GetPrevAdapters(SudoVertexAdapter &adapter, SudoVertexAdapter *out) {
    return sudo_adapter_impl_.GetPrevAdapters(adapter, out);
  }
  int OutDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.OutDegree(adapter);
  }
  int InDegree(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.InDegree(adapter);
  }
  SudoVertexAdapter NextSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.NextSimplePathAdapter(adapter);
  }
  SudoVertexAdapter PrevSimplePathAdapter(SudoVertexAdapter &adapter) {
    return sudo_adapter_impl_.PrevSimplePathAdapter(adapter);
  }

 private:
  /**
   * A wrapper for operating different types of adapters
   * @tparam AdapterType type of the vertex adapter
   */
  template <class AdapterType>
  class AdapterImpl {
   public:
    AdapterImpl(UnitigGraph *graph) : graph_(graph) {}

   public:
    AdapterType MakeVertexAdapter(size_type id, int strand = 0) {
      return {graph_->vertices_[id], strand, id};
    }
    int GetNextAdapters(AdapterType &adapter, AdapterType *out) {
      if (graph_->use_direct_adjacency_ &&
          !graph_->direct_adjacency_suspended_) {
        if (adapter.UnitigId() >= graph_->vertices_.size()) {
          xfatal("Invalid unitig ID {} in direct adjacency traversal\n",
                 adapter.UnitigId());
        }
        graph_->EnsureDirectAdjacency(adapter.UnitigId());
        const auto &adjacency =
            graph_->direct_adjacency_[adapter.UnitigId()];
        const unsigned strand = adapter.strand();
        const uint32_t *targets = adjacency.targets + strand * 4u;
        int degree = 0;
        while (degree < 4 &&
               targets[degree] != std::numeric_limits<uint32_t>::max()) {
          ++degree;
        }
        if (out) {
          for (int i = 0; i < degree; ++i) {
            const uint32_t handle = targets[i];
            if ((handle >> 1u) >= graph_->vertices_.size()) {
              xfatal("Corrupt direct adjacency handle {} for unitig {}:{}\n",
                     handle, adapter.UnitigId(), strand);
            }
            out[i] = MakeVertexAdapter(handle >> 1u, handle & 1u);
          }
        }
        return degree;
      }
      uint64_t next_starts[4];
      int degree = graph_->sdbg_->OutgoingEdges(adapter.e(), next_starts);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i] = MakeVertexAdapterWithSdbgId(next_starts[i]);
        }
      }
      return degree;
    }
    int GetPrevAdapters(AdapterType &adapter, AdapterType *out) {
      adapter.ReverseComplement();
      int degree = GetNextAdapters(adapter, out);
      if (out) {
        for (int i = 0; i < degree; ++i) {
          out[i].ReverseComplement();
        }
      }
      adapter.ReverseComplement();
      return degree;
    }
    int OutDegree(AdapterType &adapter) {
      return GetNextAdapters(adapter, nullptr);
    }
    int InDegree(AdapterType &adapter) {
      adapter.ReverseComplement();
      int degree = OutDegree(adapter);
      adapter.ReverseComplement();
      return degree;
    }
    AdapterType NextSimplePathAdapter(AdapterType &adapter) {
      if (graph_->use_direct_adjacency_ &&
          !graph_->direct_adjacency_suspended_) {
        AdapterType next[4];
        if (GetNextAdapters(adapter, next) != 1) {
          return AdapterType{};
        }
        AdapterType reverse_next = next[0];
        reverse_next.ReverseComplement();
        if (GetNextAdapters(reverse_next, nullptr) != 1) {
          return AdapterType{};
        }
        return next[0];
      }
      uint64_t next_sdbg_id = graph_->sdbg_->NextSimplePathEdge(adapter.e());
      if (next_sdbg_id != SDBG::kNullID) {
        return MakeVertexAdapterWithSdbgId(next_sdbg_id);
      } else {
        return AdapterType{};
      }
    }
    AdapterType PrevSimplePathAdapter(AdapterType &adapter) {
      adapter.ReverseComplement();
      AdapterType ret = NextSimplePathAdapter(adapter);
      ret.ReverseComplement();
      adapter.ReverseComplement();
      return ret;
    }

   private:
    AdapterType MakeVertexAdapterWithSdbgId(uint64_t sdbg_id) {
      uint32_t id = graph_->VertexIdForSdbgId(sdbg_id);
      AdapterType adapter(graph_->vertices_[id], 0, id);
      if (adapter.b() != sdbg_id) {
        adapter.ReverseComplement();
      }
      return adapter;
    }

   private:
    UnitigGraph *graph_;
  };

  std::vector<size_type> CollectRefreshAffected();
  std::vector<size_type> CollectRefreshSeeds(
      const std::vector<size_type> &affected);
  void RefreshDisconnected(const std::vector<size_type> &affected);
  void RefreshDelta(bool set_changed,
                    const std::vector<size_type> &affected,
                    const std::vector<size_type> &seeds);
  void EnsureDirectAdjacency(size_type id);
  void BuildDirectAdjacency(size_type id);
  void InvalidateDirectAdjacency();
  uint64_t SimpleNextForMaterialization(uint64_t edge) const;
  void ReplaceVertexIdForSdbgId(uint64_t old_sdbg_id,
                                uint64_t new_sdbg_id,
                                size_type unitig_id) {
    if (use_dense_id_map_) {
      assert(old_sdbg_id < dense_id_map_.size());
      assert(new_sdbg_id < dense_id_map_.size());
      dense_id_map_[old_sdbg_id] = kNullVertexID;
      dense_id_map_[new_sdbg_id] = unitig_id;
    } else {
      id_map_.ReplaceConcurrent(old_sdbg_id, new_sdbg_id, unitig_id);
    }
  }
  size_type VertexIdForSdbgId(uint64_t sdbg_id) const {
    if (use_dense_id_map_) {
      assert(sdbg_id < dense_id_map_.size());
      size_type id = dense_id_map_[sdbg_id];
      assert(id != kNullVertexID);
      return id;
    }
    return id_map_.at(sdbg_id);
  }

  // The legacy one-thread constructor visits non-loop unitigs by ascending
  // terminal SDBG edge and appends the initial loop components afterwards.
  // Parallel multistream construction deliberately assigns physical slots in
  // completion order, so a slot ID must never be used as an algorithmic
  // tie-break.  This immutable key retains the historical order even after a
  // vertex's begin/end edges are changed by Refresh.
  uint64_t LegacyOrderKey(size_type id) const {
    assert(id < legacy_order_keys_.size());
    return legacy_order_keys_[id];
  }
  bool LegacyOrderLess(size_type lhs, size_type rhs) const {
    const uint64_t lhs_key = LegacyOrderKey(lhs);
    const uint64_t rhs_key = LegacyOrderKey(rhs);
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs < rhs;
  }
  void SetVertexIdForSdbgId(uint64_t sdbg_id, size_type unitig_id) {
    if (use_dense_id_map_) {
      assert(sdbg_id < dense_id_map_.size());
      dense_id_map_[sdbg_id] = unitig_id;
    } else {
      id_map_[sdbg_id] = unitig_id;
    }
  }
  void UpdateVertexIdForSdbgIdConcurrent(uint64_t sdbg_id,
                                         size_type unitig_id) {
    if (use_dense_id_map_) {
      assert(sdbg_id < dense_id_map_.size());
      dense_id_map_[sdbg_id] = unitig_id;
    } else {
      id_map_.UpdateExistingConcurrent(sdbg_id, unitig_id);
    }
  }
  void EraseVertexIdForSdbgId(uint64_t sdbg_id) {
    if (use_dense_id_map_) {
      assert(sdbg_id < dense_id_map_.size());
      dense_id_map_[sdbg_id] = kNullVertexID;
    } else {
      id_map_.erase(sdbg_id);
    }
  }

  // Endpoint IDs are immutable during graph traversals and are updated only
  // in RefreshDisconnected's dedicated mutation phase.  Sharding retains
  // lock-free phmap lookups while allowing independent refresh workers to
  // replace endpoints without a graph-wide critical section.
  class EndpointMap {
   private:
    struct Shard {
      phmap::flat_hash_map<uint32_t, size_type> compact_entries;
      phmap::flat_hash_map<uint64_t, size_type> wide_entries;
      SpinLock lock;
    };

    size_t ShardFor(uint64_t key) const {
      key ^= key >> 30u;
      key *= UINT64_C(0xbf58476d1ce4e5b9);
      key ^= key >> 27u;
      key *= UINT64_C(0x94d049bb133111eb);
      key ^= key >> 31u;
      return static_cast<size_t>(key) & shard_mask_;
    }

   public:
    struct BufferedEntry {
      uint64_t key;
      size_type value;
    };

    // Construction inserts two immutable endpoints per unitig.  Buffer a
    // short run for each shard so one lock handoff publishes many phmap
    // entries; later sparse Refresh mutations keep using the fine-grained
    // single-key methods below.  Buffers are thread-local and allocated lazily.
    class BufferedInserter {
     public:
      explicit BufferedInserter(EndpointMap *map)
          : map_(map),
            batch_size_(map->batch_size_),
            buffers_(map->num_shards_) {}
      BufferedInserter(const BufferedInserter &) = delete;
      BufferedInserter &operator=(const BufferedInserter &) = delete;

      void Insert(uint64_t key, size_type value) {
        const size_t shard = map_->ShardFor(key);
        auto &buffer = buffers_[shard];
        if (buffer.empty()) {
          buffer.reserve(batch_size_);
        }
        buffer.push_back({key, value});
        if (buffer.size() == batch_size_) {
          Flush(shard);
        }
      }

      void FlushAll() {
        for (size_t shard = 0; shard < map_->num_shards_; ++shard) {
          if (!buffers_[shard].empty()) {
            Flush(shard);
          }
        }
      }

     private:
      void Flush(size_t shard) {
        auto &buffer = buffers_[shard];
        Shard &target = map_->shards_[shard];
        std::lock_guard<SpinLock> guard(target.lock);
        for (const BufferedEntry &entry : buffer) {
          map_->SetUnlocked(target, entry.key, entry.value);
        }
        buffer.clear();
      }

      EndpointMap *map_;
      size_t batch_size_;
      std::vector<std::vector<BufferedEntry>> buffers_;
    };

    EndpointMap() : shards_(new Shard[1]) {}

    void clear() {
      for (size_t i = 0; i < num_shards_; ++i) {
        shards_[i].compact_entries.clear();
        shards_[i].wide_entries.clear();
      }
    }
    void reserve(size_t size, size_t concurrency,
                 uint64_t key_space_size) {
      concurrency = std::max<size_t>(1, concurrency);
      compact_keys_ =
          key_space_size <= uint64_t{std::numeric_limits<uint32_t>::max()} + 1;

      // Keep several independently locked shards per active worker, but do
      // not create more shards than the input can populate efficiently.  The
      // latter bound also prevents the thread-local batching matrix from
      // growing quadratically on future many-core machines with small graphs.
      const size_t concurrency_target =
          concurrency > std::numeric_limits<size_t>::max() / 4
              ? std::numeric_limits<size_t>::max()
              : concurrency * 4;
      const size_t workload_target =
          std::max<size_t>(1, size / concurrency / 32);
      const size_t shard_target =
          std::min(concurrency_target, workload_target);
      size_t selected_shards = 1;
      while (selected_shards <= shard_target / 2) {
        selected_shards *= 2;
      }
      if (selected_shards != num_shards_) {
        shards_.reset(new Shard[selected_shards]);
        num_shards_ = selected_shards;
        shard_mask_ = selected_shards - 1;
      }

      const size_t mean_per_worker_shard =
          size / concurrency / num_shards_;
      batch_size_ = std::max<size_t>(
          4, std::min<size_t>(64, mean_per_worker_shard / 4 + 1));
      const size_t per_shard =
          size / num_shards_ + (size % num_shards_ != 0);
      for (size_t i = 0; i < num_shards_; ++i) {
        if (compact_keys_) {
          shards_[i].compact_entries.reserve(per_shard);
        } else {
          shards_[i].wide_entries.reserve(per_shard);
        }
      }
    }
    size_t num_shards() const { return num_shards_; }
    size_t batch_size() const { return batch_size_; }
    unsigned key_bits() const { return compact_keys_ ? 32u : 64u; }
    size_t size() const {
      size_t total = 0;
      for (size_t i = 0; i < num_shards_; ++i) {
        total += compact_keys_ ? shards_[i].compact_entries.size()
                               : shards_[i].wide_entries.size();
      }
      return total;
    }
    size_type at(uint64_t key) const {
      const Shard &shard = shards_[ShardFor(key)];
      return compact_keys_
                 ? shard.compact_entries.at(static_cast<uint32_t>(key))
                 : shard.wide_entries.at(key);
    }
    size_type &operator[](uint64_t key) {
      Shard &shard = shards_[ShardFor(key)];
      return compact_keys_
                 ? shard.compact_entries[static_cast<uint32_t>(key)]
                 : shard.wide_entries[key];
    }
    void InsertConcurrent(uint64_t key, size_type value) {
      Shard &shard = shards_[ShardFor(key)];
      std::lock_guard<SpinLock> guard(shard.lock);
      SetUnlocked(shard, key, value);
    }
    void UpdateExistingConcurrent(uint64_t key, size_type value) {
      // Refresh changes only mapped IDs.  The table was fully built and
      // reserved before graph cleaning, so distinct existing elements can be
      // updated concurrently without mutating hash-table structure.
      Shard &shard = shards_[ShardFor(key)];
      if (compact_keys_) {
        auto it = shard.compact_entries.find(static_cast<uint32_t>(key));
        assert(it != shard.compact_entries.end());
        it->second = value;
      } else {
        auto it = shard.wide_entries.find(key);
        assert(it != shard.wide_entries.end());
        it->second = value;
      }
    }
    size_t erase(uint64_t key) {
      Shard &shard = shards_[ShardFor(key)];
      return compact_keys_
                 ? shard.compact_entries.erase(static_cast<uint32_t>(key))
                 : shard.wide_entries.erase(key);
    }
    void ReplaceConcurrent(uint64_t old_key, uint64_t new_key,
                           size_type value) {
      const size_t old_shard = ShardFor(old_key);
      const size_t new_shard = ShardFor(new_key);
      if (old_shard == new_shard) {
        std::lock_guard<SpinLock> guard(shards_[old_shard].lock);
        EraseUnlocked(shards_[old_shard], old_key);
        SetUnlocked(shards_[old_shard], new_key, value);
        return;
      }
      {
        std::lock_guard<SpinLock> guard(shards_[old_shard].lock);
        EraseUnlocked(shards_[old_shard], old_key);
      }
      {
        std::lock_guard<SpinLock> guard(shards_[new_shard].lock);
        SetUnlocked(shards_[new_shard], new_key, value);
      }
    }

   private:
    void SetUnlocked(Shard &shard, uint64_t key, size_type value) {
      if (compact_keys_) {
        shard.compact_entries[static_cast<uint32_t>(key)] = value;
      } else {
        shard.wide_entries[key] = value;
      }
    }
    void EraseUnlocked(Shard &shard, uint64_t key) {
      if (compact_keys_) {
        shard.compact_entries.erase(static_cast<uint32_t>(key));
      } else {
        shard.wide_entries.erase(key);
      }
    }

    std::unique_ptr<Shard[]> shards_;
    size_t num_shards_{1};
    size_t shard_mask_{0};
    size_t batch_size_{4};
    bool compact_keys_{false};
  };

 private:
  // Two directed orientations, four outgoing branches each.  Targets are
  // stable `(slot_id << 1) | strand` handles, so graph traversal is a pair of
  // contiguous array reads instead of an SDBG rank/select query followed by a
  // sparse endpoint-map probe.
  struct DirectAdjacency {
    uint32_t targets[8];
  };

  SDBG *sdbg_{};
  std::vector<UnitigGraphVertex> vertices_;
  std::vector<size_type> active_ids_;
  std::vector<uint64_t> legacy_order_keys_;
  // The block-compressed simple-successor snapshot is built anyway while
  // finding unitigs.  Retaining it lets sequence materialization reuse that
  // work instead of issuing another dependent rank/select traversal for
  // every base.  Missing/stale entries always fall back to the live SDBG, so
  // graph cleaning may still expose new simple links without changing
  // semantics.
  std::unique_ptr<uint8_t[]> materialization_simple_codes_;
  std::unique_ptr<uint32_t[]> materialization_simple_bases_;
  uint64_t materialization_simple_edge_count_{0};
  std::unique_ptr<DirectAdjacency[]> direct_adjacency_;
  std::vector<AtomicWrapper<uint32_t>> direct_adjacency_epoch_;
  uint32_t direct_adjacency_generation_{1};
  bool use_direct_adjacency_{false};
  bool direct_adjacency_suspended_{false};
  std::vector<size_type> dense_id_map_;
  EndpointMap id_map_;
  bool use_dense_id_map_{false};
  AdapterImpl<VertexAdapter> adapter_impl_;
  AdapterImpl<SudoVertexAdapter> sudo_adapter_impl_;
};

#endif  // MEGAHIT_UNITIG_GRAPH_H
