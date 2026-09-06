#ifndef MEGAHIT_IDBA_COMPACT_ENDPOINT_MAP_H
#define MEGAHIT_IDBA_COMPACT_ENDPOINT_MAP_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <utility>
#include <vector>

#include "idba/kmer.h"
#include "parallel_hashmap/phmap.h"
#include "xxhash/xxh3.h"

// Exact endpoint lookup with active-width keys.  IdbaKmer reserves 72 bytes
// for the global maximum k, while local assembly normally uses only 1--5
// words.  The flat table keeps a compact stable reference and the arena stores
// only those active words.  Map iteration is never semantic; insertion order
// and duplicate replacement remain controlled by ContigGraph.
class CompactEndpointMap {
 private:
  struct Key {
    Key(const uint64_t *key_words = nullptr, uint16_t active_words = 0)
        : words(key_words), word_count(active_words) {}
    const uint64_t *words;
    uint16_t word_count;
  };

  struct KeyHash {
    size_t operator()(const Key &key) const {
      return static_cast<size_t>(
          XXH3_64bits(key.words, size_t(key.word_count) * sizeof(uint64_t)));
    }
  };

  struct KeyEqual {
    bool operator()(const Key &lhs, const Key &rhs) const {
      return lhs.word_count == rhs.word_count &&
             std::memcmp(lhs.words, rhs.words,
                         size_t(lhs.word_count) * sizeof(uint64_t)) == 0;
    }
  };

  using Map = phmap::flat_hash_map<Key, uint64_t, KeyHash, KeyEqual>;

  class KeyArena {
   public:
    uint64_t *Allocate(size_t words) {
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

      const size_t previous = blocks_.empty() ? 0u : blocks_.back().capacity;
      size_t capacity = previous == 0u ? 64u : previous * 2u;
      if (capacity < words) capacity = words;
      blocks_.emplace_back(capacity);
      uint64_t *result = blocks_.back().data.get();
      offset_ = words;
      return result;
    }

    void clear() {
      block_index_ = 0;
      offset_ = 0;
    }

    void swap(KeyArena &other) {
      blocks_.swap(other.blocks_);
      std::swap(block_index_, other.block_index_);
      std::swap(offset_, other.offset_);
    }

   private:
    struct Block {
      explicit Block(size_t words)
          : data(new uint64_t[words]), capacity(words) {}
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
      size_t capacity;
    };

    std::vector<Block> blocks_;
    size_t block_index_{0};
    size_t offset_{0};
  };

 public:
  uint64_t *find(const IdbaKmer &key) {
    const uint16_t words = static_cast<uint16_t>((key.size() + 31u) >> 5u);
    auto iter = map_.find(Key(key.words(), words));
    return iter == map_.end() ? nullptr : &iter->second;
  }

  void insert(const IdbaKmer &key, uint64_t value) {
    const uint16_t words = static_cast<uint16_t>((key.size() + 31u) >> 5u);
    uint64_t *stored = arena_.Allocate(words);
    std::memcpy(stored, key.words(), size_t(words) * sizeof(uint64_t));
    map_.emplace(Key(stored, words), value);
  }

  void reserve(size_t size) { map_.reserve(size); }

  void clear() {
    map_.clear();
    arena_.clear();
  }

  void swap(CompactEndpointMap &other) {
    map_.swap(other.map_);
    arena_.swap(other.arena_);
  }

 private:
  Map map_;
  KeyArena arena_;
};

#endif  // MEGAHIT_IDBA_COMPACT_ENDPOINT_MAP_H
