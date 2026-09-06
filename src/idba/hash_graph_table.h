#ifndef MEGAHIT_IDBA_HASH_GRAPH_TABLE_H
#define MEGAHIT_IDBA_HASH_GRAPH_TABLE_H

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "idba/hash.h"
#include "idba/hash_graph_vertex.h"
#include "idba/kmer.h"

// A compact, single-threaded table for the per-endpoint IDBA graph.
//
// The historical HashTableST stores every vertex in a separately addressed
// chained node.  Local assembly creates hundreds of thousands of small graphs,
// so those pointers and allocator chunks dominate cache misses.  This table
// keeps vertices and next links in flat arrays, but deliberately reproduces
// HashTableST's bucket count, head insertion and rehash traversal order.  The
// graph algorithms therefore observe the same vertex order and tie behavior.
class HashGraphVertexTable {
 public:
  using size_type = size_t;
  static constexpr uint32_t kNull = std::numeric_limits<uint32_t>::max();
  static constexpr uint32_t kDefaultNumBuckets = 1u << 12u;

  class iterator {
   public:
    iterator(HashGraphVertexTable *owner = nullptr, uint32_t index = kNull)
        : owner_(owner), index_(index) {}

    bool operator==(const iterator &other) const {
      return index_ == other.index_;
    }
    bool operator!=(const iterator &other) const {
      return index_ != other.index_;
    }
    HashGraphVertex &operator*() const { return owner_->values_[index_]; }
    HashGraphVertex *operator->() const { return &owner_->values_[index_]; }

   private:
    HashGraphVertexTable *owner_;
    uint32_t index_;
  };

  HashGraphVertexTable()
      : bucket_heads_(kDefaultNumBuckets, kNull),
        occupied_words_((kDefaultNumBuckets + 63u) / 64u, 0u) {}

  iterator find(const IdbaKmer &key) {
    const uint64_t hash = hasher_(key);
    uint32_t index = bucket_heads_[BucketIndex(hash)];
    while (index != kNull) {
      if (values_[index].EqualsKey(key)) return iterator(this, index);
      index = BucketNext(index);
    }
    return end();
  }

  const HashGraphVertex *find_value(const IdbaKmer &key) const {
    const uint64_t hash = hasher_(key);
    uint32_t index = bucket_heads_[BucketIndex(hash)];
    while (index != kNull) {
      if (values_[index].EqualsKey(key)) return &values_[index];
      index = BucketNext(index);
    }
    return nullptr;
  }

  HashGraphVertex &find_or_insert_key(const IdbaKmer &key,
                                      uint32_t *index_out = nullptr) {
    ApplyPendingRehash();
    const uint64_t hash = hasher_(key);
    const size_type bucket = BucketIndex(hash);
    uint32_t index = bucket_heads_[bucket];
    while (index != kNull) {
      if (values_[index].EqualsKey(key)) {
        if (index_out != nullptr) *index_out = index;
        return values_[index];
      }
      index = BucketNext(index);
    }

    if (values_.size() >= kNull) {
      throw std::length_error("local hash graph exceeds compact index range");
    }
    const uint32_t new_index = static_cast<uint32_t>(values_.size());
    EnsureBucketLinkWidth(new_index);
    const uint32_t words = (key.size() + 31u) >> 5u;
    const uint64_t *stored_key = key.words();
    if (words > 1u) {
      uint64_t *arena_key = key_arena_.Allocate(words);
      for (uint32_t word = 0; word < words; ++word) {
        arena_key[word] = key.word(word);
      }
      stored_key = arena_key;
    }
    values_.emplace_back(stored_key, key.size(), static_cast<uint32_t>(hash));
    AppendBucketNext(bucket_heads_[bucket]);
    if (bucket_heads_[bucket] == kNull) RecordOccupied(bucket);
    bucket_heads_[bucket] = new_index;
    MarkPendingRehash();
    if (index_out != nullptr) *index_out = new_index;
    return values_.back();
  }

  // Lookup directly from a compact rolling key.  Duplicate occurrences avoid
  // constructing/copying the fixed-capacity IdbaKmer; only a genuinely new
  // vertex is materialized.  HashPackedWords returns the exact historical
  // full-capacity XXH3 value, so bucket chains and traversal order are
  // unchanged.
  template <unsigned Words>
  HashGraphVertex &find_or_insert_packed(const uint64_t *key_words,
                                         uint32_t kmer_size,
                                         uint64_t fixed_hash_accumulator,
                                         uint32_t *index_out = nullptr) {
    ApplyPendingRehash();
    const uint64_t hash = IdbaKmer::HashPackedWords<Words>(
        key_words, kmer_size, fixed_hash_accumulator);
    const size_type bucket = BucketIndex(hash);
    uint32_t index = bucket_heads_[bucket];
    while (index != kNull) {
      if (values_[index].EqualsPackedWords<Words>(key_words, kmer_size)) {
        if (index_out != nullptr) *index_out = index;
        return values_[index];
      }
      index = BucketNext(index);
    }

    if (values_.size() >= kNull) {
      throw std::length_error("local hash graph exceeds compact index range");
    }
    const uint32_t new_index = static_cast<uint32_t>(values_.size());
    EnsureBucketLinkWidth(new_index);
    const uint64_t *stored_key = key_words;
    if (Words > 1u) {
      uint64_t *arena_key = key_arena_.Allocate(Words);
      std::memcpy(arena_key, key_words, Words * sizeof(uint64_t));
      stored_key = arena_key;
    }
    values_.emplace_back(stored_key, kmer_size, static_cast<uint32_t>(hash));
    AppendBucketNext(bucket_heads_[bucket]);
    if (bucket_heads_[bucket] == kNull) RecordOccupied(bucket);
    bucket_heads_[bucket] = new_index;
    MarkPendingRehash();
    if (index_out != nullptr) *index_out = new_index;
    return values_.back();
  }

  HashGraphVertex &value_at(uint32_t index) { return values_[index]; }
  const HashGraphVertex &value_at(uint32_t index) const {
    return values_[index];
  }

  // A duplicate occurrence reached through the graph's exact transition
  // cache still has to observe a pending rehash at the same logical point as
  // find_or_insert_packed().  Rehashing does not change compact vertex IDs,
  // so callers may then dereference the cached ID without another hash/probe.
  void prepare_cached_occurrence() { ApplyPendingRehash(); }

  uint32_t index_of(const HashGraphVertex &vertex) const {
    return static_cast<uint32_t>(&vertex - values_.data());
  }

  iterator end() { return iterator(); }

  template <typename UnaryProc>
  UnaryProc &for_each(UnaryProc &op) {
    // Historical traversal is ascending bucket order.  A 1-bit occupancy
    // directory produces that order directly and avoids sorting a touched
    // bucket vector millions of times across endpoint/k graphs.
    for (size_type word = 0; word < occupied_words_.size(); ++word) {
      uint64_t bits = occupied_words_[word];
      while (bits != 0u) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
        const size_type bucket = (word << 6u) + bit;
        uint32_t index = bucket_heads_[bucket];
        while (index != kNull) {
          op(values_[index]);
          index = BucketNext(index);
        }
        bits &= bits - 1u;
      }
    }
    return op;
  }

  // Reductions such as coverage histograms have no traversal-order
  // semantics.  Stream the compact vertex array instead of following legacy
  // bucket chains; Assemble() continues to use for_each() above for exact tie
  // behavior.
  template <typename UnaryProc>
  UnaryProc &for_each_value(UnaryProc &op) {
    for (HashGraphVertex &value : values_) op(value);
    return op;
  }

  void reserve(size_type capacity) {
    values_.reserve(capacity);
    if (wide_bucket_links_) next_.reserve(capacity);
    RehashIfNeeded(capacity);
  }

  void clear() {
    values_.clear();
    next_.clear();
    wide_bucket_links_ = false;
    key_arena_.clear();
    for (size_type word = 0; word < occupied_words_.size(); ++word) {
      uint64_t bits = occupied_words_[word];
      while (bits != 0u) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
        bucket_heads_[(word << 6u) + bit] = kNull;
        bits &= bits - 1u;
      }
      occupied_words_[word] = 0u;
    }
    rehash_pending_ = false;
  }

  size_type size() const { return values_.size(); }

  size_type bucket_count() const { return bucket_heads_.size(); }
  bool rehash_pending() const { return rehash_pending_; }

  void reset_empty_bucket_count(size_type count) {
    if (!values_.empty() || count == 0u || (count & (count - 1u)) != 0u) {
      throw std::logic_error("invalid empty local graph bucket reset");
    }
    bucket_heads_.assign(count, std::numeric_limits<uint32_t>::max());
    occupied_words_.assign((count + 63u) / 64u, 0u);
    rehash_pending_ = false;
  }

  void traversal_indices(std::vector<uint32_t> *indices) const {
    indices->clear();
    indices->reserve(values_.size());
    for (size_type word = 0; word < occupied_words_.size(); ++word) {
      uint64_t bits = occupied_words_[word];
      while (bits != 0u) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
        const size_type bucket = (word << 6u) + bit;
        uint32_t index = bucket_heads_[bucket];
        while (index != kNull) {
          indices->push_back(index);
          index = BucketNext(index);
        }
        bits &= bits - 1u;
      }
    }
  }

  // Retain allocated vectors between endpoint tasks while restoring the
  // exact logical state of a freshly constructed historical table.  In
  // particular, carrying an expanded bucket count into another endpoint
  // changes bucket traversal/tie order even when the table is empty.
  void reset_for_endpoint() {
    clear();
    bucket_heads_.resize(kDefaultNumBuckets);
    occupied_words_.resize((kDefaultNumBuckets + 63u) / 64u);
  }

  // A bulk reducer replays only first occurrences.  If its final unique key
  // was followed by duplicate occurrences, the historical table would check
  // and perform a pending rehash on the very next duplicate.  Reproduce that
  // occurrence boundary without probing the duplicate key again.
  void finish_bulk_occurrences(bool had_occurrence_after_last_unique) {
    if (had_occurrence_after_last_unique) ApplyPendingRehash();
  }

  void swap(HashGraphVertexTable &other) {
    values_.swap(other.values_);
    next_.swap(other.next_);
    bucket_heads_.swap(other.bucket_heads_);
    occupied_words_.swap(other.occupied_words_);
    key_arena_.swap(other.key_arena_);
    std::swap(wide_bucket_links_, other.wide_bucket_links_);
    std::swap(rehash_pending_, other.rehash_pending_);
  }

 private:
  // Geometrically growing blocks keep key addresses stable while storing only
  // the words active at the graph's current k.  clear() rewinds the arena and
  // retains its high-water allocation for the next endpoint task.
  class CompactKeyArena {
   public:
    uint64_t *Allocate(size_type words) {
      while (block_index_ < blocks_.size()) {
        Block &block = blocks_[block_index_];
        if (block.capacity - offset_ >= words) {
          uint64_t *result = block.data.get() + offset_;
          offset_ += words;
          return result;
        }
        ++block_index_;
        offset_ = 0;
      }

      const size_type previous =
          blocks_.empty() ? 0u : blocks_.back().capacity;
      size_type capacity = previous == 0u ? 64u : previous * 2u;
      if (capacity < words) capacity = words;
      if (previous != 0u && capacity < previous) {
        throw std::length_error("local compact key arena overflow");
      }
      blocks_.emplace_back(capacity);
      uint64_t *result = blocks_.back().data.get();
      offset_ = words;
      return result;
    }

    void clear() {
      block_index_ = 0;
      offset_ = 0;
    }

    void swap(CompactKeyArena &other) {
      blocks_.swap(other.blocks_);
      std::swap(block_index_, other.block_index_);
      std::swap(offset_, other.offset_);
    }

   private:
    struct Block {
      explicit Block(size_type word_capacity)
          : data(new uint64_t[word_capacity]), capacity(word_capacity) {}
      Block(Block &&other)
          : data(std::move(other.data)), capacity(other.capacity) {}
      Block &operator=(Block &&other) {
        data = std::move(other.data);
        capacity = other.capacity;
        return *this;
      }
      Block(const Block &) = delete;
      Block &operator=(const Block &) = delete;

      std::unique_ptr<uint64_t[]> data;
      size_type capacity;
    };

    std::vector<Block> blocks_;
    size_type block_index_{0};
    size_type offset_{0};
  };

  size_type BucketIndex(uint64_t hash) const {
    return static_cast<size_type>(hash) & (bucket_heads_.size() - 1u);
  }

  void RecordOccupied(size_type bucket) {
    occupied_words_[bucket >> 6u] |= uint64_t{1} << (bucket & 63u);
  }

  uint32_t BucketNext(uint32_t index) const {
    return wide_bucket_links_
               ? next_[index]
               : (values_[index].bucket_next16() ==
                          HashGraphVertex::kNoCompactBucketNext
                      ? kNull
                      : static_cast<uint32_t>(
                            values_[index].bucket_next16()));
  }

  void SetBucketNext(uint32_t index, uint32_t next) {
    if (wide_bucket_links_) {
      next_[index] = next;
    } else {
      assert(next == kNull ||
             next < HashGraphVertex::kNoCompactBucketNext);
      values_[index].set_bucket_next16(
          next == kNull
              ? HashGraphVertex::kNoCompactBucketNext
              : static_cast<uint16_t>(next));
    }
  }

  void AppendBucketNext(uint32_t next) {
    if (wide_bucket_links_) {
      next_.push_back(next);
    } else {
      assert(next == kNull ||
             next < HashGraphVertex::kNoCompactBucketNext);
      values_.back().set_bucket_next16(
          next == kNull
              ? HashGraphVertex::kNoCompactBucketNext
              : static_cast<uint16_t>(next));
    }
  }

  void EnsureBucketLinkWidth(uint32_t new_index) {
    if (wide_bucket_links_ ||
        new_index < HashGraphVertex::kNoCompactBucketNext) {
      return;
    }
    next_.clear();
    next_.reserve(values_.capacity());
    for (const HashGraphVertex &value : values_) {
      next_.push_back(
          value.bucket_next16() == HashGraphVertex::kNoCompactBucketNext
              ? kNull
              : static_cast<uint32_t>(value.bucket_next16()));
    }
    wide_bucket_links_ = true;
  }

  void RehashIfNeeded(size_type capacity) {
    if (capacity <= bucket_heads_.size() * 2u) return;
    size_type new_bucket_count = bucket_heads_.size();
    while (capacity > new_bucket_count * 2u) {
      if (new_bucket_count >
          std::numeric_limits<size_type>::max() / 2u) {
        throw std::length_error("local hash graph bucket count overflow");
      }
      new_bucket_count *= 2u;
    }
    Rehash(new_bucket_count);
    rehash_pending_ = false;
  }

  void MarkPendingRehash() {
    rehash_pending_ = values_.size() > bucket_heads_.size() * 2u;
  }

  void ApplyPendingRehash() {
    if (!rehash_pending_) return;
    RehashIfNeeded(values_.size());
  }

  void Rehash(size_type new_bucket_count) {
    if (new_bucket_count == bucket_heads_.size()) return;
    if ((new_bucket_count & (new_bucket_count - 1u)) != 0) {
      throw std::logic_error("invalid local hash graph bucket count");
    }

    std::vector<uint32_t> new_heads(new_bucket_count, kNull);
    std::vector<uint64_t> new_occupied_words(
        (new_bucket_count + 63u) / 64u, 0u);
    // This loop and head insertion intentionally match HashTableST::rehash.
    // Walking the occupancy directory is exactly the same ascending bucket
    // order as the historical full scan, but skips empty buckets.  Build the
    // new directory while scattering so no second full-bucket scan is needed.
    for (size_type word = 0; word < occupied_words_.size(); ++word) {
      uint64_t bits = occupied_words_[word];
      while (bits != 0u) {
        const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
        const size_type old_bucket = (word << 6u) + bit;
        uint32_t index = bucket_heads_[old_bucket];
        while (index != kNull) {
          const uint32_t old_next = BucketNext(index);
          const size_type new_bucket =
              static_cast<size_type>(values_[index].bucket_hash()) &
              (new_bucket_count - 1u);
          SetBucketNext(index, new_heads[new_bucket]);
          if (new_heads[new_bucket] == kNull) {
            new_occupied_words[new_bucket >> 6u] |=
                uint64_t{1} << (new_bucket & 63u);
          }
          new_heads[new_bucket] = index;
          index = old_next;
        }
        bits &= bits - 1u;
      }
    }
    bucket_heads_.swap(new_heads);
    occupied_words_.swap(new_occupied_words);
  }

  Hash<IdbaKmer> hasher_;
  std::vector<HashGraphVertex> values_;
  std::vector<uint32_t> next_;
  std::vector<uint32_t> bucket_heads_;
  std::vector<uint64_t> occupied_words_;
  CompactKeyArena key_arena_;
  bool wide_bucket_links_{false};
  bool rehash_pending_{false};
};

#endif  // MEGAHIT_IDBA_HASH_GRAPH_TABLE_H
