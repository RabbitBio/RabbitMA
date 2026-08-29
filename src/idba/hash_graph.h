/**
 * @file hash_graph.h
 * @brief HashGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_HASH_GRAPH_H_

#define __GRAPH_HASH_GRAPH_H_

#include <deque>
#include <istream>
#include <ostream>
#include <vector>

#include "bit_operation.h"
#include "idba/contig_info.h"
#include "idba/hash_graph_vertex.h"
#include "idba/hash_table.h"
#include "idba/kmer.h"
#include "idba/sequence.h"
#include "utils/histgram.h"

class IdbaKmer;
class Sequence;

/**
 * @brief It is a hash table based de Bruijn graph implementation.
 */
class HashGraph {
 public:
  typedef HashTableST<HashGraphVertex, IdbaKmer> vertex_table_type;

  explicit HashGraph(uint32_t kmer_size = 0) { set_kmer_size(kmer_size); }
  ~HashGraph() {}

  HashGraphVertexAdaptor FindVertexAdaptor(const IdbaKmer &kmer) {
    IdbaKmer key = kmer.unique_format();
    auto p = vertex_table_.find(key);
    return ((p != vertex_table_.end())
                ? HashGraphVertexAdaptor(&*p, kmer != key)
                : HashGraphVertexAdaptor(NULL));
  }

  int64_t InsertKmers(const Sequence &seq);
  int64_t InsertUncountKmers(const Sequence &seq);

  void ClearStatus() {
    ClearStatusFunc func;
    vertex_table_.for_each(func);
  }

  int64_t Assemble(std::deque<Sequence> &contigs,
                   std::deque<ContigInfo> &contig_infos);

  void reserve(uint64_t capacity) { vertex_table_.reserve(capacity); }

  uint32_t kmer_size() const { return kmer_size_; }
  void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

  int coverage_percentile(double percentile) {
    CoveragePercentileFunc func;
    vertex_table_.for_each(func);
    return func.percentile(percentile);
  }

  void swap(HashGraph &hash_graph) {
    if (this != &hash_graph) {
      vertex_table_.swap(hash_graph.vertex_table_);
      std::swap(kmer_size_, hash_graph.kmer_size_);
    }
  }

  uint64_t num_vertices() const { return vertex_table_.size(); }
  void clear() { vertex_table_.clear(); }

 private:
#if __cplusplus >= 201103L
  HashGraph(const HashGraph &) = delete;
  const HashGraph &operator=(const HashGraph &) = delete;
#else
  HashGraph(const HashGraph &);
  const HashGraph &operator=(const HashGraph &);
#endif

  bool GetNextVertexAdaptor(const HashGraphVertexAdaptor &current,
                            HashGraphVertexAdaptor &next) {
    if (current.out_edges().size() != 1) return false;

    IdbaKmer kmer = current.kmer();
    kmer.ShiftAppend(bit_operation::BitToIndex(current.out_edges()));
    next = FindVertexAdaptor(kmer);

    return !kmer.IsPalindrome() && next.in_edges().size() == 1;
  }

  bool IsLoop(const Sequence &contig, HashGraphVertexAdaptor &next) {
    IdbaKmer kmer = next.kmer();
    IdbaKmer rev_comp = kmer;
    rev_comp.ReverseComplement();

    return contig.GetIdbaKmer(0, kmer_size_) == kmer;
  }

  bool IsPalindromeLoop(const Sequence &contig, HashGraphVertexAdaptor &next) {
    IdbaKmer kmer = next.kmer();
    IdbaKmer rev_comp = kmer;
    rev_comp.ReverseComplement();

    return contig.GetIdbaKmer(contig.size() - kmer_size_, kmer_size_) ==
           rev_comp;
  }

  class ClearStatusFunc {
   public:
    ClearStatusFunc() {}

    void operator()(HashGraphVertex &vertex) { vertex.status().clear(); }
  };

  class AssembleFunc {
   public:
    AssembleFunc(HashGraph *hash_graph) : hash_graph_(hash_graph) {}
    ~AssembleFunc() {}

    void operator()(HashGraphVertex &vertex);

    std::deque<Sequence> &contigs() { return contigs_; }
    std::deque<ContigInfo> &contig_infos() { return contig_infos_; }

   private:
    HashGraph *hash_graph_;
    std::deque<Sequence> contigs_;
    std::deque<ContigInfo> contig_infos_;
  };

  class CoveragePercentileFunc {
   public:
    void operator()(HashGraphVertex &vertex) {
      const size_t coverage = static_cast<size_t>(vertex.count());
      if (coverage >= bins_.size()) {
        bins_.resize(coverage + 1u, 0);
      }
      ++bins_[coverage];
      ++size_;
    }

    int percentile(double p) const {
      if (size_ == 0) {
        return 0;
      }
      const size_t rank = static_cast<size_t>(size_ * p);
      size_t cumulative = 0;
      for (size_t coverage = 0; coverage < bins_.size(); ++coverage) {
        cumulative += bins_[coverage];
        if (cumulative > rank) {
          return static_cast<int>(coverage);
        }
      }
      return 0;
    }

   private:
    std::vector<size_t> bins_;
    size_t size_{0};
  };

  vertex_table_type vertex_table_;
  uint32_t kmer_size_;
};

namespace std {
inline void swap(HashGraph &x, HashGraph &y) { x.swap(y); }
}  // namespace std

#endif
