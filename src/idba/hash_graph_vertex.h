/**
 * @file hash_graph_vertex.h
 * @brief HashGraphVertex Class and HashGraphVertexAdaptor Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __GRAPH_HASH_GRAPH_VERTEX_H_

#define __GRAPH_HASH_GRAPH_VERTEX_H_

#include <algorithm>

#include "idba/bit_edges.h"
#include "idba/bit_operation.h"
#include "idba/kmer.h"
#include "idba/vertex_status.h"

/**
 * @brief It is the vertex class used in HashGraph.
 */
class HashGraphVertex {
 public:
  enum : uint32_t {
    kNoCachedNeighbor = UINT32_MAX,
    kUncacheableNeighbor = UINT32_MAX - 1u
  };
  enum : uint16_t { kNoCompactBucketNext = UINT16_MAX };

  explicit HashGraphVertex(const uint64_t *key_words = NULL,
                           uint32_t kmer_size = 0,
                           uint32_t bucket_hash = 0)
      : key_storage_(kmer_size <= 32u && key_words != NULL
                         ? key_words[0]
                         : static_cast<uint64_t>(reinterpret_cast<uintptr_t>(
                               key_words))),
        count_(0),
        kmer_size_(static_cast<uint16_t>(kmer_size)),
        bucket_next16_(kNoCompactBucketNext),
        next_neighbor_{kNoCachedNeighbor, kNoCachedNeighbor},
        bucket_hash_(bucket_hash) {}
  HashGraphVertex(const HashGraphVertex &x)
      : key_storage_(x.key_storage_),
        count_(x.count_),
        kmer_size_(x.kmer_size_),
        status_(x.status_),
        in_edges_(x.in_edges_),
        out_edges_(x.out_edges_),
        bucket_next16_(x.bucket_next16_),
        next_neighbor_{x.next_neighbor_[0], x.next_neighbor_[1]},
        bucket_hash_(x.bucket_hash_) {}

  const HashGraphVertex &operator=(const HashGraphVertex &x) {
    key_storage_ = x.key_storage_;
    count_ = x.count_;
    kmer_size_ = x.kmer_size_;
    status_ = x.status_;
    in_edges_ = x.in_edges_;
    out_edges_ = x.out_edges_;
    bucket_next16_ = x.bucket_next16_;
    next_neighbor_[0] = x.next_neighbor_[0];
    next_neighbor_[1] = x.next_neighbor_[1];
    bucket_hash_ = x.bucket_hash_;
    return *this;
  }

  void FixPalindromeEdges() {
    if (kmer().IsPalindrome()) out_edges_ = in_edges_ = (in_edges_ | out_edges_);
  }

  IdbaKmer key() const { return kmer(); }

  bool EqualsKey(const IdbaKmer &key) const {
    if (key.size() != kmer_size_) return false;
    const uint32_t words = (kmer_size_ + 31u) >> 5u;
    for (uint32_t word = 0; word < words; ++word) {
      if (key_words()[word] != key.word(word)) return false;
    }
    return true;
  }

  template <unsigned Words>
  bool EqualsPackedWords(const uint64_t *words,
                         uint32_t kmer_size) const {
    if (kmer_size_ != kmer_size) return false;
    for (unsigned word = 0; word < Words; ++word) {
      if (key_words()[word] != words[word]) return false;
    }
    return true;
  }

  IdbaKmer kmer() const {
    IdbaKmer result;
    result.AssignWords(key_words(), kmer_size_);
    return result;
  }

  uint64_t hash() const { return kmer().hash(); }
  uint32_t bucket_hash() const { return bucket_hash_; }

  uint8_t base(uint32_t index) const {
    return static_cast<uint8_t>(
        (key_words()[index >> 5u] >> ((index & 31u) << 1u)) & 3u);
  }

  uint8_t last_base(bool is_reverse) const {
    return !is_reverse ? base(kmer_size_ - 1u) : 3u - base(0);
  }

  int32_t &count() { return count_; }
  const int32_t &count() const { return count_; }

  VertexStatus &status() { return status_; }
  const VertexStatus &status() const { return status_; }

  BitEdges &in_edges() { return in_edges_; }
  const BitEdges &in_edges() const { return in_edges_; }

  BitEdges &out_edges() { return out_edges_; }
  const BitEdges &out_edges() const { return out_edges_; }

  uint16_t bucket_next16() const { return bucket_next16_; }
  void set_bucket_next16(uint16_t next) { bucket_next16_ = next; }

  uint32_t &next_neighbor(bool is_reverse) {
    return next_neighbor_[is_reverse ? 1 : 0];
  }
  uint32_t next_neighbor(bool is_reverse) const {
    return next_neighbor_[is_reverse ? 1 : 0];
  }

  void swap(HashGraphVertex &x) {
    if (this != &x) {
      std::swap(key_storage_, x.key_storage_);
      std::swap(count_, x.count_);
      std::swap(kmer_size_, x.kmer_size_);
      status_.swap(x.status_);
      in_edges_.swap(x.in_edges_);
      out_edges_.swap(x.out_edges_);
      std::swap(bucket_next16_, x.bucket_next16_);
      std::swap(next_neighbor_[0], x.next_neighbor_[0]);
      std::swap(next_neighbor_[1], x.next_neighbor_[1]);
      std::swap(bucket_hash_, x.bucket_hash_);
    }
  }

  uint32_t kmer_size() const { return kmer_size_; }
  const uint64_t *key_words() const {
    if (kmer_size_ <= 32u) return &key_storage_;
    return reinterpret_cast<const uint64_t *>(
        static_cast<uintptr_t>(key_storage_));
  }

  void clear() {
    in_edges_.clear();
    out_edges_.clear();
    status_.clear();
    count_ = 0;
    bucket_next16_ = kNoCompactBucketNext;
    next_neighbor_[0] = kNoCachedNeighbor;
    next_neighbor_[1] = kNoCachedNeighbor;
  }

 private:
  // Canonical bases live in HashGraphVertexTable's stable, active-width key
  // arena.  Keeping only a reference here avoids 72 fixed bytes per vertex
  // when a small k uses one or two machine words.
  // A one-word key occupies the same eight bytes as the pointer used by wider
  // keys. kmer_size_ is the tag, so k<=32 avoids an arena allocation and a
  // dependent pointer load without increasing the vertex footprint.
  uint64_t key_storage_;
  int32_t count_;
  uint16_t kmer_size_;
  VertexStatus status_;
  BitEdges in_edges_;
  BitEdges out_edges_;
  // Most endpoint graphs have far fewer than 65,535 vertices.  Their bucket
  // chain link fits in the two bytes that were padding here, keeping the hash
  // key and its next link on the same cache line.  HashGraphVertexTable
  // promotes losslessly to a side uint32_t array before this range is
  // exceeded, so large samples retain the full historical index domain.
  uint16_t bucket_next16_;
  uint32_t next_neighbor_[2];
  uint32_t bucket_hash_;
};

/**
 * @brief It is adaptor class used for accessing HashGraphVertex. Because
 * a k-mer and its reverse complemtn share the same vertex, using adaptor
 * makes sure the access to vertex consistant.
 */
class HashGraphVertexAdaptor {
 public:
  explicit HashGraphVertexAdaptor(HashGraphVertex *vertex = NULL,
                                  bool is_reverse = false) {
    vertex_ = vertex;
    is_reverse_ = is_reverse;
  }
  HashGraphVertexAdaptor(const HashGraphVertexAdaptor &x) {
    vertex_ = x.vertex_;
    is_reverse_ = x.is_reverse_;
  }

  const HashGraphVertexAdaptor &operator=(const HashGraphVertexAdaptor &x) {
    vertex_ = x.vertex_;
    is_reverse_ = x.is_reverse_;
    return *this;
  }

  bool operator<(const HashGraphVertexAdaptor &x) const {
    return (vertex_ != x.vertex_) ? (vertex_ < x.vertex_)
                                  : (is_reverse_ < x.is_reverse_);
  }
  bool operator>(const HashGraphVertexAdaptor &x) const {
    return (vertex_ != x.vertex_) ? (vertex_ > x.vertex_)
                                  : (is_reverse_ > x.is_reverse_);
  }

  bool operator==(const HashGraphVertexAdaptor &x) const {
    return vertex_ == x.vertex_ && is_reverse_ == x.is_reverse_;
  }
  bool operator!=(const HashGraphVertexAdaptor &x) const {
    return vertex_ != x.vertex_ || is_reverse_ != x.is_reverse_;
  }

  const HashGraphVertexAdaptor &ReverseComplement() {
    is_reverse_ = !is_reverse_;
    return *this;
  }

  IdbaKmer kmer() const {
    IdbaKmer kmer = vertex_->kmer();
    return !is_reverse_ ? kmer : kmer.ReverseComplement();
  }

  HashGraphVertex &vertex() { return *vertex_; }
  const HashGraphVertex &vertex() const { return *vertex_; }
  void set_vertex(HashGraphVertex *vertex, bool is_reverse = false) {
    vertex_ = vertex;
    is_reverse_ = is_reverse;
  }

  int32_t &count() { return vertex_->count(); }
  const int32_t &count() const { return vertex_->count(); }

  VertexStatus &status() { return vertex_->status(); }
  const VertexStatus &status() const { return vertex_->status(); }

  BitEdges &in_edges() {
    return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges();
  }
  const BitEdges &in_edges() const {
    return !is_reverse_ ? vertex_->in_edges() : vertex_->out_edges();
  }

  BitEdges &out_edges() {
    return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges();
  }
  const BitEdges &out_edges() const {
    return !is_reverse_ ? vertex_->out_edges() : vertex_->in_edges();
  }

  uint32_t &next_neighbor() {
    return vertex_->next_neighbor(is_reverse_);
  }
  uint32_t next_neighbor() const {
    return vertex_->next_neighbor(is_reverse_);
  }


  void swap(HashGraphVertexAdaptor &x) {
    if (this != &x) {
      std::swap(vertex_, x.vertex_);
      std::swap(is_reverse_, x.is_reverse_);
    }
  }

  bool is_null() const { return vertex_ == NULL; }
  bool is_reverse() const { return is_reverse_; }

  uint32_t kmer_size() const { return vertex_->kmer_size(); }

  uint8_t last_base() const {
    return vertex_->last_base(is_reverse_);
  }

  void clear() { vertex_->clear(); }

 private:
  HashGraphVertex *vertex_;
  bool is_reverse_;
};

namespace std {
template <>
inline void swap(HashGraphVertex &x, HashGraphVertex &y) {
  x.swap(y);
}
template <>
inline void swap(HashGraphVertexAdaptor &x, HashGraphVertexAdaptor &y) {
  x.swap(y);
}
}  // namespace std

#endif
