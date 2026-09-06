/**
 * @file contig_info.h
 * @brief ContigInfo Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#ifndef __GRAPH_CONTIG_INFO_H_

#define __GRAPH_CONTIG_INFO_H_

#include "bit_edges.h"

#include <algorithm>
#include <deque>
#include <memory>
#include <vector>
#include <istream>
#include <ostream>
#include <string>

typedef uint32_t SequenceCountUnitType;
typedef std::vector<SequenceCountUnitType> SequenceCount;

class ContigBuilder;

/**
 * @brief It is used to store information of contigs, like k-mer counts,
 * in-edges,
 * out-edges, etc.
 */
class ContigInfo {
  friend class ContigBuilder;

 public:
  ContigInfo() {
    kmer_count_ = 0;
    kmer_size_ = 0;
  }

  ContigInfo(const ContigInfo &contig_info) {
    in_edges_ = contig_info.in_edges_;
    out_edges_ = contig_info.out_edges_;
    kmer_count_ = contig_info.kmer_count_;
    kmer_size_ = contig_info.kmer_size_;
    if (contig_info.counts_) {
      counts_.reset(new SequenceCount(*contig_info.counts_));
    }
  }

  const ContigInfo &operator=(const ContigInfo &contig_info) {
    in_edges_ = contig_info.in_edges_;
    out_edges_ = contig_info.out_edges_;
    kmer_count_ = contig_info.kmer_count_;
    kmer_size_ = contig_info.kmer_size_;
    if (contig_info.counts_) {
      if (!counts_) counts_.reset(new SequenceCount());
      *counts_ = *contig_info.counts_;
    } else {
      counts_.reset();
    }
    return *this;
  }

  const ContigInfo &ReverseComplement() {
    std::swap(in_edges_, out_edges_);
    if (counts_) std::reverse(counts_->begin(), counts_->end());
    return *this;
  }

  BitEdges &in_edges() { return in_edges_; }
  const BitEdges &in_edges() const { return in_edges_; }

  BitEdges &out_edges() { return out_edges_; }
  const BitEdges &out_edges() const { return out_edges_; }

  uint32_t kmer_size() const { return kmer_size_; }
  void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

  uint32_t kmer_count() const { return kmer_count_; }
  void set_kmer_count(uint32_t kmer_count) { kmer_count_ = kmer_count; }

  const SequenceCount &counts() const {
    static const SequenceCount empty;
    return counts_ ? *counts_ : empty;
  }
  void set_counts(const SequenceCount &counts) {
    if (counts.empty()) {
      counts_.reset();
    } else {
      mutable_counts() = counts;
    }
  }

  void swap(ContigInfo &contig_info) {
    if (this != &contig_info) {
      std::swap(in_edges_, contig_info.in_edges_);
      std::swap(out_edges_, contig_info.out_edges_);
      std::swap(kmer_size_, contig_info.kmer_size_);
      std::swap(kmer_count_, contig_info.kmer_count_);
      counts_.swap(contig_info.counts_);
    }
  }

  void clear() {
    in_edges_ = 0;
    out_edges_ = 0;
    kmer_size_ = 0;
    kmer_count_ = 0;
    if (counts_) counts_->clear();
  }

 private:
  BitEdges in_edges_;
  BitEdges out_edges_;
  uint16_t kmer_size_;
  uint32_t kmer_count_;
  SequenceCount &mutable_counts() {
    if (!counts_) counts_.reset(new SequenceCount());
    return *counts_;
  }

  // Local assembly deliberately carries only aggregate coverage.  Keeping
  // the unused vector object inline cost 24 bytes in every transient unitig;
  // allocate it only for legacy callers that actually materialize per-k-mer
  // counts.
  std::unique_ptr<SequenceCount> counts_;
};

namespace std {
template <>
inline void swap(ContigInfo &x, ContigInfo &y) {
  x.swap(y);
}
}  // namespace std

#endif
