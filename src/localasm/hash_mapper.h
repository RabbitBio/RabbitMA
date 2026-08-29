//
// Created by vout on 7/10/19.
//

#ifndef MEGAHIT_LOCALASM_HASH_MAPPER_H
#define MEGAHIT_LOCALASM_HASH_MAPPER_H

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "parallel_hashmap/phmap.h"
#include "sequence/kmer_plus.h"
#include "sequence/sequence_package.h"

struct MappingRecord {
  uint32_t contig_id;
  int32_t contig_from;
  int32_t contig_to;
  uint64_t query_id;
  int32_t query_from;
  int32_t query_to;
  uint32_t mismatch;
  uint8_t strand;
  bool valid;

  bool operator<(const MappingRecord &rhs) const {
    if (contig_id != rhs.contig_id) return contig_id < rhs.contig_id;
    if (contig_from != rhs.contig_from) return contig_from < rhs.contig_from;
    if (contig_to != rhs.contig_to) return contig_to < rhs.contig_to;
    if (query_id != rhs.query_id) return query_id < rhs.query_id;
    if (query_from != rhs.query_from) return query_from < rhs.query_from;
    if (query_to != rhs.query_to) return query_to < rhs.query_to;
    return strand < rhs.strand;
  }

  bool operator==(const MappingRecord &rhs) const {
    return contig_id == rhs.contig_id && contig_from == rhs.contig_from &&
           contig_to == rhs.contig_to && query_id == rhs.query_id &&
           query_from == rhs.query_from && query_to == rhs.query_to &&
           strand == rhs.strand;
  };
};

class HashMapper {
 public:
  using TKmer = Kmer<2, uint32_t>;
  using TSeedKey = uint64_t;
  struct IndexHasher {
    size_t operator()(TSeedKey key) const {
      // The mapper's seed is at most 32 bases, so its complete key is already
      // two uint32 words.  phmap applies its own strong multiplicative mix;
      // running XXH3 over these same eight bytes first only duplicates work in
      // the several-billion-lookup mapping loop.
      return static_cast<size_t>(key ^ (key >> (sizeof(size_t) * 4u)));
    }
  };

  using TIndexShard = phmap::flat_hash_map<TSeedKey, uint64_t, IndexHasher>;

  struct EndpointSeedWitness {
    uint64_t index_value{0};
    int32_t end_position{0};
    uint8_t query_strand{0};
    bool valid{false};
  };

  void LoadAndBuild(const std::string &contig_file, int32_t min_len,
                    int32_t seed_kmer_size, int32_t sparsity,
                    bool build_seed_filter = true);

  // Return the exact set of non-repetitive sampled contig seeds whose sole
  // occurrence lies in an endpoint range.  The persistent read index uses
  // this to reverse the local-mapping join (endpoint -> candidate reads),
  // while the original mapper remains the final semantic verifier.
  std::vector<TSeedKey> CollectUniqueEndpointSeeds(
      int32_t endpoint_range) const;

  void SetMappingThreshold(int32_t mapping_len, double similarity) {
    min_mapped_len_ = mapping_len;
    similarity_ = similarity;
  }

  // Build an exact-candidate prefilter for the only contig regions consumed
  // by local assembly. False Bloom positives are verified in the full index;
  // false negatives for a retained mapping are impossible.
  void BuildEndpointSeedFilter(int32_t endpoint_range);
  bool MayMapToEndpoint(
      const SeqPackage::SeqView &seq_view,
      EndpointSeedWitness *witness = nullptr) const;
  bool MayMapToEndpoint(const uint32_t *packed_words, unsigned length,
                        EndpointSeedWitness *witness = nullptr) const;

  MappingRecord TryMap(
      const SeqPackage ::SeqView &seq_view,
      const EndpointSeedWitness *witness = nullptr) const;
  MappingRecord TryMap(const uint32_t *packed_words, unsigned length,
                       uint64_t query_id,
                       const EndpointSeedWitness *witness = nullptr) const;
  const SeqPackage &refseq() const { return refseq_; }


  // Mapping is the sole consumer of the multi-GiB seed index and Bloom
  // filter.  Keep refseq for endpoint assembly, but return the transient
  // lookup storage before local graphs begin allocating their workspaces.
  void ReleaseIndex() {
    std::vector<TIndexShard>().swap(index_);
    std::vector<uint64_t>().swap(seed_filter_);
    std::vector<std::vector<uint64_t>>().swap(seed_filter_replicas_);
    std::vector<uint64_t>().swap(endpoint_seed_filter_);
    std::vector<std::vector<uint64_t>>().swap(
        endpoint_seed_filter_replicas_);
    seed_filter_mask_ = 0;
    endpoint_seed_filter_mask_ = 0;
    endpoint_range_ = 0;
    index_shard_shift_ = sizeof(size_t) * 8u;
    index_shards_power_of_two_ = true;
  }

 private:
  static TSeedKey SeedKey(const TKmer &kmer) {
    return (static_cast<uint64_t>(kmer.data()[0]) << 32u) | kmer.data()[1];
  }
  static size_t IndexHash(TSeedKey key) {
    // phmap's explicit-hash lookup expects its internally mixed hash, not the
    // raw user hasher result.
    return phmap::phmap_mix<sizeof(size_t)>()(IndexHasher{}(key));
  }
  unsigned IndexShard(size_t hash) const {
    // Use high bits for the outer shard.  phmap's probe sequence consumes the
    // low bits; conditioning those bits on a shard would make every key in a
    // shard start in the same small set of groups.
    if (index_shards_power_of_two_) {
      return index_shard_shift_ == sizeof(size_t) * 8u
                 ? 0u
                 : static_cast<unsigned>(hash >> index_shard_shift_);
    }
    if (sizeof(size_t) == sizeof(uint64_t)) {
      const uint64_t high = static_cast<uint32_t>(hash >> 32u);
      return static_cast<unsigned>((high * index_.size()) >> 32u);
    }
    return static_cast<unsigned>(
        (uint64_t(static_cast<uint32_t>(hash)) * index_.size()) >> 32u);
  }
  void BuildSeedFilter(size_t index_size);
  const uint64_t *LocalSeedFilter() const;
  const uint64_t *LocalEndpointSeedFilter() const;
  bool SeedMayContain(size_t hash, const uint64_t *filter) const {
    if (filter == nullptr) {
      return false;
    }
    const uint64_t bits = SeedFilterBits(hash);
    return (filter[hash & seed_filter_mask_] & bits) == bits;
  }
  static uint64_t SeedFilterBits(uint64_t hash) {
    return (uint64_t{1} << ((hash >> 32u) & 63u)) |
           (uint64_t{1} << ((hash >> 40u) & 63u)) |
           (uint64_t{1} << ((hash >> 48u) & 63u)) |
           (uint64_t{1} << ((hash >> 56u) & 63u));
  }

  MappingRecord TryMapRaw(const uint32_t *query_words,
                          unsigned query_shift, unsigned length,
                          uint64_t query_id,
                          const EndpointSeedWitness *witness) const;
  bool MayMapToEndpointRaw(const uint32_t *query_words,
                           unsigned query_shift, unsigned length,
                           EndpointSeedWitness *witness) const;
  int32_t Match(const uint32_t *query_words, unsigned query_shift,
                int query_from, int query_to, size_t contig_id, int ref_from,
                int ref_to, bool strand) const;

 private:
  std::vector<TIndexShard> index_;
  unsigned index_shard_shift_{sizeof(size_t) * 8u};
  bool index_shards_power_of_two_{true};
  // One cache-line-local probe rejects absent read seeds before they touch the
  // much larger contig index.  Eight bits/key and four bits/probe retain a low
  // false-positive rate; false negatives are impossible.
  std::vector<uint64_t> seed_filter_;
  // Optional topology-sized replicas.  They contain identical bits and only
  // change physical placement, never membership semantics.
  std::vector<std::vector<uint64_t>> seed_filter_replicas_;
  size_t seed_filter_mask_{0};
  std::vector<uint64_t> endpoint_seed_filter_;
  std::vector<std::vector<uint64_t>> endpoint_seed_filter_replicas_;
  size_t endpoint_seed_filter_mask_{0};
  int32_t endpoint_range_{0};
  SeqPackage refseq_;

  int32_t seed_kmer_size_{31};
  int32_t index_sparsity_{1};
  int32_t min_mapped_len_{50};
  double similarity_{0.95};
};

#endif  // MEGAHIT_LOCALASM_HASH_MAPPER_H
