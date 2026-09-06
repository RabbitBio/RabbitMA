/**
 * @file contig_builder.h
 * @brief Contig Build Class which builds contig and related contig info.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2011-12-27
 */

#ifndef __GRAPH_CONTIG_BUILDER_H_

#define __GRAPH_CONTIG_BUILDER_H_

#include "idba/contig_graph_vertex.h"
#include "idba/contig_info.h"
#include "idba/hash_graph_vertex.h"
#include "idba/sequence.h"

/**
 * @brief It is a builder class for building contigs.
 */
class ContigBuilder {
 public:
  explicit ContigBuilder(bool materialize_position_counts = true)
      : materialize_position_counts_(materialize_position_counts) {}

  explicit ContigBuilder(HashGraphVertexAdaptor x)
      : materialize_position_counts_(true) { Append(x); }

  explicit ContigBuilder(ContigGraphVertexAdaptor x)
      : materialize_position_counts_(true) { Append(x, 0); }

  // Hash-graph unitig traversal knows its final number of vertices before it
  // materializes bases. Reserve both payloads once so the final, semantic
  // order can be streamed without geometric reallocations.
  void Reserve(size_t sequence_size, size_t count_size) {
    contig_.reserve(sequence_size);
    if (materialize_position_counts_) {
      contig_info_.mutable_counts().reserve(count_size);
    }
  }

  void ReserveHashUnitig(uint32_t kmer_size, size_t num_vertices) {
    if (num_vertices == 0) return;
    Reserve(size_t(kmer_size) + num_vertices - 1u, num_vertices);
  }

  void Append(HashGraphVertexAdaptor x) {
    if (contig_.size() == 0) {
      const IdbaKmer kmer = x.kmer();
      contig_.Assign(kmer);
      contig_info_.in_edges_ = x.in_edges();
      contig_info_.out_edges_ = x.out_edges();
      contig_info_.kmer_size_ = kmer.size();
      contig_info_.kmer_count_ = x.count();
      if (materialize_position_counts_) {
        SequenceCount &counts = contig_info_.mutable_counts();
        counts.resize(1);
        counts[0] = x.count();
      }
    } else {
      // A unitig extension contributes only its terminal oriented base.  The
      // historical path copied (and on the reverse strand transformed) the
      // complete fixed-capacity IdbaKmer merely to read this one symbol.
      contig_ += x.last_base();
      contig_info_.out_edges_ = x.out_edges();
      contig_info_.kmer_count_ += x.count();
      if (materialize_position_counts_) {
        contig_info_.mutable_counts().push_back(x.count());
      }
    }
  }

  void Append(ContigGraphVertexAdaptor x, int d) {
    const Sequence &source_contig = x.stored_contig();
    const SequenceCount &source_counts = x.stored_counts();
    if (contig_.size() == 0) {
      if (!x.is_reverse())
        contig_ = source_contig;
      else
        contig_.AssignReverseComplement(source_contig);
      contig_info_.in_edges_ = x.in_edges();
      contig_info_.out_edges_ = x.out_edges();
      contig_info_.kmer_size_ = x.kmer_size();
      contig_info_.kmer_count_ = x.kmer_count();
      if (materialize_position_counts_) {
        SequenceCount &counts = contig_info_.mutable_counts();
        if (!x.is_reverse())
          counts = source_counts;
        else
          counts.assign(source_counts.rbegin(), source_counts.rend());
      }
    } else {
      if (d <= 0) {
        const int offset = std::min(-d, (int)x.contig_size());
        if (!x.is_reverse())
          contig_.Append(source_contig, offset);
        else
          contig_.AppendReverseComplement(source_contig, offset);
        contig_info_.out_edges_ = x.out_edges();
        contig_info_.kmer_count_ += x.kmer_count();
        if (materialize_position_counts_) {
          SequenceCount &counts = contig_info_.mutable_counts();
          int start = std::min(-d - (int)contig_info_.kmer_size_ + 1,
                               (int)source_counts.size());
          if (!x.is_reverse()) {
            counts.insert(counts.end(), source_counts.begin() + start,
                          source_counts.end());
          } else {
            counts.insert(counts.end(), source_counts.rbegin() + start,
                          source_counts.rend());
          }
        }
      } else {
        contig_.Append(d, 4);
        if (!x.is_reverse())
          contig_.Append(source_contig);
        else
          contig_.AppendReverseComplement(source_contig);
        contig_info_.out_edges_ = x.out_edges();
        contig_info_.kmer_count_ += x.kmer_count();
        if (materialize_position_counts_) {
          SequenceCount &counts = contig_info_.mutable_counts();
          counts.insert(counts.end(), d, 0);

          if (!x.is_reverse()) {
            counts.insert(counts.end(), source_counts.begin(),
                          source_counts.end());
          } else {
            counts.insert(counts.end(), source_counts.rbegin(),
                          source_counts.rend());
          }
        }
      }
    }
  }

  const ContigBuilder &ReverseComplement() {
    contig_.ReverseComplement();
    contig_info_.ReverseComplement();
    return *this;
  }

  const Sequence &contig() const { return contig_; }

  const ContigInfo &contig_info() const { return contig_info_; }

  // Transfer the completed buffers to their consumer.  Local assembly builds
  // millions of short-lived unitigs; copying both the sequence string and the
  // per-k-mer count vector here used to materialize every unitig one extra
  // time before ContigGraph immediately consumed it.
  void Release(Sequence &contig, ContigInfo &contig_info) {
    contig.swap(contig_);
    contig_info.swap(contig_info_);
  }

  void clear() {
    contig_.clear();
    contig_info_.clear();
  }

 private:
  Sequence contig_;
  ContigInfo contig_info_;
  bool materialize_position_counts_;
};

#endif
