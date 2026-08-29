//
// Created by vout on 7/10/19.
//

#include "hash_mapper.h"
#include <array>
#include <cstdlib>
#include <mutex>
#include <omp.h>
#include "sequence/io/contig/contig_reader.h"
#include "utils/startup_affinity.h"
#include "utils/utils.h"

namespace {

inline uint64_t EncodeContigOffset(uint32_t contig_id, uint32_t contig_offset,
                                   uint8_t strand) {
  return (uint64_t(contig_id) << 32) | (contig_offset << 1) | strand;
}

inline void DecodeContigOffset(uint64_t code, uint32_t &contig_id,
                               uint32_t &contig_offset, uint8_t &strand) {
  contig_id = code >> 32;
  contig_offset = (code & 0xFFFFFFFFULL) >> 1;
  strand = code & 1ULL;
}

inline uint32_t GetWord(const uint32_t *first_word, uint32_t first_shift,
                        int from, int len, bool strand) {
  int from_word_idx = (first_shift + from) / 16;
  int from_word_shift = (first_shift + from) % 16;
  uint32_t ret = *(first_word + from_word_idx) << from_word_shift * 2;
  assert(len <= 16);

  if (16 - from_word_shift < len) {
    ret |= *(first_word + from_word_idx + 1) >> (16 - from_word_shift) * 2;
  }

  if (len < 16) {
    ret >>= (16 - len) * 2;
    ret <<= (16 - len) * 2;
  }

  if (strand == 1) {
    ret = kmlib::bit::ReverseComplement<2>(ret);
    ret <<= (16 - len) * 2;
  }

  return ret;
}

inline int Mismatch(uint32_t x, uint32_t y) {
  x ^= y;
  x |= x >> 1;
  x &= 0x55555555U;
  return kmlib::bit::Popcount(x);
}

inline uint64_t GetWord64(const uint32_t *first_word, uint32_t first_shift,
                          int from, int len, bool strand) {
  const int absolute = static_cast<int>(first_shift) + from;
  const int word_index = absolute / 16;
  const int base_shift = absolute % 16;
  const unsigned bit_shift = static_cast<unsigned>(base_shift * 2);

  uint64_t joined = uint64_t(first_word[word_index]) << 32u;
  if (base_shift + len > 16) joined |= first_word[word_index + 1];
  uint64_t result = joined << bit_shift;
  if (base_shift + len > 32) {
    result |= uint64_t(first_word[word_index + 2]) >> (32u - bit_shift);
  }
  if (len < 32) {
    const unsigned padding = static_cast<unsigned>((32 - len) * 2);
    result = (result >> padding) << padding;
  }
  if (strand) {
    result = kmlib::bit::ReverseComplement<2>(result);
    result <<= static_cast<unsigned>((32 - len) * 2);
  }
  return result;
}

inline int Mismatch64(uint64_t x, uint64_t y) {
  x ^= y;
  x |= x >> 1u;
  x &= UINT64_C(0x5555555555555555);
  return kmlib::bit::Popcount(x);
}

}  // namespace

void HashMapper::LoadAndBuild(const std::string &contig_file, int32_t min_len,
                              int32_t seed_kmer_size, int32_t sparsity,
                              bool build_seed_filter) {
  if (seed_kmer_size <= 0 ||
      seed_kmer_size > static_cast<int32_t>(TKmer::max_size())) {
    throw std::invalid_argument("local mapper seed length must be 1..32");
  }
  if (sparsity <= 0) {
    throw std::invalid_argument("local mapper sparsity must be positive");
  }
  seed_kmer_size_ = seed_kmer_size;
  index_sparsity_ = sparsity;
  ContigReader reader(contig_file);
  reader.SetMinLen(min_len)->SetDiscardFlag(contig_flag::kLoop);
  auto sizes = reader.GetNumContigsAndBases();
  refseq_.Clear();
  refseq_.ReserveSequences(sizes.first);
  refseq_.ReserveBases(sizes.second);
  bool contig_reverse = false;
  reader.ReadAll(&refseq_, contig_reverse);

  size_t sz = refseq_.seq_count();
  size_t n_kmers = 0;

#pragma omp parallel for reduction(+ : n_kmers)
  for (size_t i = 0; i < sz; ++i) {
    n_kmers +=
        (refseq_.GetSeqView(i).length() - seed_kmer_size + sparsity) / sparsity;
  }

  // Keep outer routing to a shift and give every worker at least one
  // independently buildable shard.  Non-power-of-two shard counts save some
  // phmap slack at particular table sizes, but add a multiply-range operation
  // to every one of the billions of mapping probes.  The power-of-two choice
  // is topology-independent and won the end-to-end bandwidth benchmark.
  const size_t desired_shards = std::max<size_t>(
      1u, std::min<size_t>(static_cast<size_t>(omp_get_max_threads()),
                           std::max<size_t>(1u, n_kmers)));
  size_t num_shards = 1u;
  while (num_shards < desired_shards) {
    if (num_shards > std::numeric_limits<size_t>::max() / 2u) {
      throw std::length_error("local mapper shard count overflow");
    }
    num_shards <<= 1u;
  }

  index_shards_power_of_two_ =
      (num_shards & (num_shards - 1u)) == 0u;
  unsigned shard_bits = 0u;
  if (index_shards_power_of_two_) {
    size_t value = num_shards;
    while (value > 1u) {
      value >>= 1u;
      ++shard_bits;
    }
  }
  index_.clear();
  index_.resize(num_shards);
  index_shard_shift_ = index_shards_power_of_two_
                           ? static_cast<unsigned>(sizeof(size_t) * 8u) -
                                 shard_bits
                           : static_cast<unsigned>(sizeof(size_t) * 8u);

  const size_t reserve_per_shard =
      (n_kmers + num_shards - 1u) / num_shards;
#pragma omp parallel for schedule(static)
  for (int64_t shard = 0; shard < static_cast<int64_t>(num_shards); ++shard) {
    index_[shard].clear();
    index_[shard].reserve(reserve_per_shard);
    AdviseHugePages(index_[shard].storage_address(),
                    index_[shard].storage_bytes());
  }
  std::vector<std::mutex> shard_locks(num_shards);
  const unsigned seed_padding = 64u - (seed_kmer_size_ << 1u);

#pragma omp parallel
  {
    // Keep a small thread-local batch per shard.  A flush takes one lock for
    // dozens of inserts, while independent shards proceed concurrently.
    std::vector<std::vector<std::pair<TSeedKey, uint64_t>>> pending(
        num_shards);
    const size_t kReservePerPendingShard = 64u;
    const size_t kPendingLimit = num_shards * kReservePerPendingShard;
    for (auto &items : pending) {
      items.reserve(kReservePerPendingShard);
    }
    size_t num_pending = 0;
    auto flush_pending = [&]() {
      for (size_t shard = 0; shard < num_shards; ++shard) {
        auto &items = pending[shard];
        if (items.empty()) continue;
        std::lock_guard<std::mutex> lock(shard_locks[shard]);
        auto &index = index_[shard];
        for (const auto &item : items) {
          auto result = index.emplace(item.first, item.second);
          if (!result.second) result.first->second |= uint64_t{1} << 63u;
        }
        items.clear();
      }
      num_pending = 0;
    };

#pragma omp for
    for (size_t i = 0; i < sz; ++i) {
      const auto contig_view = refseq_.GetSeqView(i);
      const auto contig_address = contig_view.raw_address();
      for (int j = 0, len = contig_view.length(); j + seed_kmer_size <= len;
           j += sparsity) {
        const uint64_t forward =
            GetWord64(contig_address.first, contig_address.second, j,
                      seed_kmer_size_, false);
        uint64_t reverse = kmlib::bit::ReverseComplement<2>(forward);
        if (seed_padding != 0u) reverse <<= seed_padding;
        const bool is_reverse = reverse < forward;
        const TSeedKey seed_key = is_reverse ? reverse : forward;
        const size_t hash = IndexHash(seed_key);
        pending[IndexShard(hash)].emplace_back(
            seed_key,
            EncodeContigOffset(contig_view.id(), j, is_reverse));
        if (++num_pending == kPendingLimit) flush_pending();
      }
    }

    if (num_pending != 0) flush_pending();
  }

  size_t index_size = 0;
  for (const auto &shard : index_) {
    index_size += shard.size();
  }
  if (build_seed_filter) {
    BuildSeedFilter(index_size);
  } else {
    seed_filter_.clear();
    seed_filter_replicas_.clear();
    seed_filter_mask_ = 0u;
  }
  xinfo("Number of contigs: {}, index size: {} in {} shards\n",
        refseq_.seq_count(), index_size, index_.size());
}

std::vector<HashMapper::TSeedKey> HashMapper::CollectUniqueEndpointSeeds(
    int32_t endpoint_range) const {
  const uint64_t range = static_cast<uint64_t>(std::max(endpoint_range, 0));
  std::vector<std::vector<TSeedKey>> shard_keys(index_.size());
#pragma omp parallel for schedule(dynamic, 1)
  for (int64_t shard_id = 0;
       shard_id < static_cast<int64_t>(index_.size()); ++shard_id) {
    auto &keys = shard_keys[static_cast<size_t>(shard_id)];
    keys.reserve(index_[static_cast<size_t>(shard_id)].size() / 4u + 1u);
    for (const auto &entry : index_[static_cast<size_t>(shard_id)]) {
      if ((entry.second >> 63u) != 0u) continue;
      uint32_t contig_id = 0;
      uint32_t contig_offset = 0;
      uint8_t strand = 0;
      DecodeContigOffset(entry.second, contig_id, contig_offset, strand);
      const uint64_t contig_length = refseq_.GetSeqView(contig_id).length();
      if (static_cast<uint64_t>(contig_offset) >= range &&
          static_cast<uint64_t>(contig_offset) + range < contig_length) {
        continue;
      }
      keys.push_back(entry.first);
    }
  }

  size_t total = 0;
  for (const auto &keys : shard_keys) {
    if (keys.size() > std::numeric_limits<size_t>::max() - total) {
      throw std::length_error("local endpoint seed set is too large");
    }
    total += keys.size();
  }
  std::vector<TSeedKey> result;
  result.reserve(total);
  for (auto &keys : shard_keys) {
    result.insert(result.end(), keys.begin(), keys.end());
  }
  return result;
}

void HashMapper::BuildSeedFilter(size_t index_size) {
  if (index_size == 0) {
    seed_filter_.clear();
    seed_filter_replicas_.clear();
    seed_filter_mask_ = 0;
    return;
  }

  size_t bits_per_key = 8u;
  if (const char *configured =
          std::getenv("MEGAHIT_LOCAL_SEED_FILTER_BITS_PER_KEY")) {
    const unsigned long value = std::strtoul(configured, nullptr, 10);
    if (value > 0u && value <= 64u) bits_per_key = value;
  }
  if (index_size >
      (std::numeric_limits<size_t>::max() - 63u) / bits_per_key) {
    throw std::length_error("local-mapping seed filter is too large");
  }
  const size_t target_words =
      DivCeiling(index_size * bits_per_key, size_t{64});
  size_t num_words = 1;
  while (num_words < target_words) {
    if (num_words > std::numeric_limits<size_t>::max() / 2u) {
      throw std::length_error("local-mapping seed filter is too large");
    }
    num_words *= 2u;
  }
  seed_filter_.assign(num_words, 0);
  seed_filter_replicas_.clear();
  seed_filter_mask_ = num_words - 1u;

  // The filter is a large shared random-access array. std::vector::assign()
  // otherwise first-touches every page on the calling thread's NUMA domain,
  // forcing workers on all other domains through the interconnect. Install a
  // topology-derived interleave policy, discard the already-zeroed interior,
  // and let the parallel OR pass below perform the real first touch. These
  // calls are advisory/no-ops on a single-domain or unsupported platform.
  const size_t filter_bytes = seed_filter_.size() * sizeof(seed_filter_[0]);
  const size_t numa_domains = GetNumaTopology().domain_count();
  const bool replica_bytes_fit =
      numa_domains > 1u &&
      numa_domains - 1u <=
          std::numeric_limits<size_t>::max() / filter_bytes &&
      index_size <=
          std::numeric_limits<size_t>::max() /
              (sizeof(TSeedKey) + sizeof(uint64_t)) &&
      (numa_domains - 1u) * filter_bytes <=
          index_size * (sizeof(TSeedKey) + sizeof(uint64_t));
  const bool replicate_per_domain =
      omp_get_max_threads() > 1 && replica_bytes_fit;
  if (replicate_per_domain) {
    BindMemoryPagesToNumaDomain(seed_filter_.data(), filter_bytes, 0u);
    DiscardMemoryPages(seed_filter_.data(), filter_bytes);
  } else if (omp_get_max_threads() > 1 &&
             InterleaveMemoryPages(seed_filter_.data(), filter_bytes)) {
    DiscardMemoryPages(seed_filter_.data(), filter_bytes);
  }
  AdviseHugePages(seed_filter_.data(), filter_bytes);
  // Shards are independent for iteration, but Bloom words are shared.  A
  // relaxed atomic OR is sufficient during construction and avoids a serial
  // pass over a potentially hundred-million-entry index.
  uint64_t usable_entries = 0;
  uint64_t repetitive_entries = 0;
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : usable_entries, repetitive_entries)
  for (int64_t shard_id = 0;
       shard_id < static_cast<int64_t>(index_.size()); ++shard_id) {
    for (const auto &entry : index_[shard_id]) {
      // Repetitive seeds are unconditionally rejected by process_seed below.
      // Omitting them here turns that known rejection into a single Bloom
      // load instead of a random phmap lookup followed by the same rejection.
      // Usable unique seeds are unchanged, so false negatives remain
      // impossible for every candidate that can affect mapping semantics.
      if ((entry.second >> 63u) != 0u) {
        ++repetitive_entries;
        continue;
      }
      ++usable_entries;
      const size_t hash = IndexHash(entry.first);
      __atomic_fetch_or(&seed_filter_[hash & seed_filter_mask_],
                        SeedFilterBits(hash), __ATOMIC_RELAXED);
    }
  }

  if (replicate_per_domain) {
    seed_filter_replicas_.resize(numa_domains - 1u);
    for (size_t domain = 1; domain < numa_domains; ++domain) {
      auto &replica = seed_filter_replicas_[domain - 1u];
      replica.assign(seed_filter_.size(), 0u);
      BindMemoryPagesToNumaDomain(replica.data(), filter_bytes,
                                  static_cast<unsigned>(domain));
      DiscardMemoryPages(replica.data(), filter_bytes);
      AdviseHugePages(replica.data(), filter_bytes);
    }
#pragma omp parallel for schedule(static)
    for (int64_t replica_id = 0;
         replica_id <
         static_cast<int64_t>(seed_filter_replicas_.size()); ++replica_id) {
      std::copy(seed_filter_.begin(), seed_filter_.end(),
                seed_filter_replicas_[replica_id].begin());
    }
  }

  xinfo("Local seed membership filter: {} bytes for {} usable and {} "
        "repetitive indexed k-mers ({} NUMA-local copies)\n",
        filter_bytes, usable_entries, repetitive_entries,
        seed_filter_replicas_.size() + 1u);
}

const uint64_t *HashMapper::LocalSeedFilter() const {
  if (seed_filter_.empty()) return nullptr;
  if (seed_filter_replicas_.empty()) return seed_filter_.data();
  const unsigned domain = CurrentNumaDomain();
  if (domain == 0u) return seed_filter_.data();
  const size_t replica = static_cast<size_t>(domain - 1u);
  return replica < seed_filter_replicas_.size()
             ? seed_filter_replicas_[replica].data()
             : seed_filter_.data();
}

void HashMapper::BuildEndpointSeedFilter(int32_t endpoint_range) {
  endpoint_range_ = std::max(endpoint_range, 0);
  endpoint_seed_filter_.clear();
  endpoint_seed_filter_replicas_.clear();
  endpoint_seed_filter_mask_ = 0;
  if (endpoint_range_ == 0 || index_.empty()) return;

  uint64_t endpoint_entry_upper_bound = 0;
#pragma omp parallel for reduction(+ : endpoint_entry_upper_bound)
  for (int64_t contig_id = 0;
       contig_id < static_cast<int64_t>(refseq_.seq_count()); ++contig_id) {
    const uint64_t length = refseq_.GetSeqView(contig_id).length();
    if (length < static_cast<uint64_t>(seed_kmer_size_)) continue;
    const uint64_t positions =
        (length - seed_kmer_size_) / index_sparsity_ + 1u;
    const uint64_t left_end = std::min<uint64_t>(
        positions,
        DivCeiling(static_cast<uint64_t>(endpoint_range_),
                   static_cast<uint64_t>(index_sparsity_)));
    const uint64_t right_threshold =
        length > static_cast<uint64_t>(endpoint_range_)
            ? length - endpoint_range_
            : 0u;
    const uint64_t right_begin = std::min<uint64_t>(
        positions,
        DivCeiling(right_threshold,
                   static_cast<uint64_t>(index_sparsity_)));
    endpoint_entry_upper_bound +=
        left_end + (positions - right_begin) -
        (left_end > right_begin ? left_end - right_begin : 0u);
  }
  if (endpoint_entry_upper_bound == 0u) return;

  // A read-level gate sees many seeds, so false positives are verified against
  // the exact full index before admitting a read. This controls extra work and
  // never changes mapping correctness.
  size_t bits_per_endpoint_key = 8u;
  if (const char *configured =
          std::getenv("MEGAHIT_LOCAL_ENDPOINT_FILTER_BITS_PER_KEY")) {
    const unsigned long value = std::strtoul(configured, nullptr, 10);
    if (value > 0u && value <= 64u) bits_per_endpoint_key = value;
  }
  if (endpoint_entry_upper_bound >
      (std::numeric_limits<size_t>::max() - 63u) /
          bits_per_endpoint_key) {
    throw std::length_error("local endpoint seed filter is too large");
  }
  const size_t target_words = DivCeiling(
      static_cast<size_t>(endpoint_entry_upper_bound) * bits_per_endpoint_key,
      size_t{64});
  size_t num_words = 1u;
  while (num_words < target_words) {
    if (num_words > std::numeric_limits<size_t>::max() / 2u) {
      throw std::length_error("local endpoint seed filter is too large");
    }
    num_words <<= 1u;
  }
  endpoint_seed_filter_.assign(num_words, 0u);
  endpoint_seed_filter_mask_ = num_words - 1u;

  const size_t filter_bytes = num_words * sizeof(uint64_t);
  const size_t numa_domains = GetNumaTopology().domain_count();
  const bool replicate_per_domain =
      omp_get_max_threads() > 1 && numa_domains > 1u &&
      numa_domains - 1u <=
          std::numeric_limits<size_t>::max() / filter_bytes;
  if (replicate_per_domain) {
    BindMemoryPagesToNumaDomain(endpoint_seed_filter_.data(), filter_bytes,
                                0u);
    DiscardMemoryPages(endpoint_seed_filter_.data(), filter_bytes);
  } else if (omp_get_max_threads() > 1 &&
             InterleaveMemoryPages(endpoint_seed_filter_.data(),
                                   filter_bytes)) {
    DiscardMemoryPages(endpoint_seed_filter_.data(), filter_bytes);
  }
  AdviseHugePages(endpoint_seed_filter_.data(), filter_bytes);

  uint64_t endpoint_entries = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : endpoint_entries)
  for (int64_t shard_id = 0;
       shard_id < static_cast<int64_t>(index_.size()); ++shard_id) {
    for (const auto &entry : index_[shard_id]) {
      if ((entry.second >> 63u) != 0u) continue;
      uint32_t contig_id, contig_offset;
      uint8_t strand;
      DecodeContigOffset(entry.second, contig_id, contig_offset, strand);
      const uint64_t contig_length = refseq_.GetSeqView(contig_id).length();
      if (contig_offset >= static_cast<uint32_t>(endpoint_range_) &&
          static_cast<uint64_t>(contig_offset) + endpoint_range_ <
              contig_length) {
        continue;
      }
      ++endpoint_entries;
      const size_t hash = IndexHash(entry.first);
      __atomic_fetch_or(
          &endpoint_seed_filter_[hash & endpoint_seed_filter_mask_],
          SeedFilterBits(hash), __ATOMIC_RELAXED);
    }
  }

  if (replicate_per_domain) {
    endpoint_seed_filter_replicas_.resize(numa_domains - 1u);
    for (size_t domain = 1; domain < numa_domains; ++domain) {
      auto &replica = endpoint_seed_filter_replicas_[domain - 1u];
      replica.assign(num_words, 0u);
      BindMemoryPagesToNumaDomain(replica.data(), filter_bytes,
                                  static_cast<unsigned>(domain));
      DiscardMemoryPages(replica.data(), filter_bytes);
      AdviseHugePages(replica.data(), filter_bytes);
    }
#pragma omp parallel for schedule(static)
    for (int64_t replica_id = 0;
         replica_id <
         static_cast<int64_t>(endpoint_seed_filter_replicas_.size());
         ++replica_id) {
      std::copy(endpoint_seed_filter_.begin(), endpoint_seed_filter_.end(),
                endpoint_seed_filter_replicas_[replica_id].begin());
    }
  }
  xinfo("Local endpoint seed filter: {} bytes for {} unique indexed "
        "k-mers, range {}, {} NUMA-local copies\n",
        filter_bytes, endpoint_entries, endpoint_range_,
        endpoint_seed_filter_replicas_.size() + 1u);
}

const uint64_t *HashMapper::LocalEndpointSeedFilter() const {
  if (endpoint_seed_filter_.empty()) return nullptr;
  if (endpoint_seed_filter_replicas_.empty()) {
    return endpoint_seed_filter_.data();
  }
  const unsigned domain = CurrentNumaDomain();
  if (domain == 0u) return endpoint_seed_filter_.data();
  const size_t replica = static_cast<size_t>(domain - 1u);
  return replica < endpoint_seed_filter_replicas_.size()
             ? endpoint_seed_filter_replicas_[replica].data()
             : endpoint_seed_filter_.data();
}

bool HashMapper::MayMapToEndpoint(
    const SeqPackage::SeqView &seq_view,
    EndpointSeedWitness *witness) const {
  const auto address = seq_view.raw_address();
  return MayMapToEndpointRaw(address.first, address.second,
                             seq_view.length(), witness);
}

bool HashMapper::MayMapToEndpoint(const uint32_t *packed_words,
                                  unsigned length,
                                  EndpointSeedWitness *witness) const {
  return MayMapToEndpointRaw(packed_words, 0u, length, witness);
}

bool HashMapper::MayMapToEndpointRaw(const uint32_t *query_words,
                                     unsigned query_shift,
                                     unsigned length,
                                     EndpointSeedWitness *witness) const {
  if (witness != nullptr) witness->valid = false;
  if (endpoint_range_ <= 0) return true;
  const uint64_t *const filter = LocalEndpointSeedFilter();
  if (filter == nullptr || length < static_cast<unsigned>(seed_kmer_size_) ||
      length < 50u) {
    return false;
  }
  TKmer initial_seed(query_words, query_shift, seed_kmer_size_);
  uint64_t seed_forward = SeedKey(initial_seed);
  const unsigned seed_padding = 64u - (seed_kmer_size_ << 1u);
  const uint64_t seed_mask =
      seed_padding == 0u ? UINT64_MAX : UINT64_MAX << seed_padding;
  uint64_t seed_reverse = kmlib::bit::ReverseComplement<2>(seed_forward);
  if (seed_padding != 0u) seed_reverse <<= seed_padding;

  struct EndpointProbe {
    TSeedKey key;
    size_t hash;
    int32_t end_position;
    uint8_t query_strand;
  };
  constexpr unsigned kProbeBatch = 32u;
  std::array<EndpointProbe, kProbeBatch> probes;
  std::array<uint8_t, kProbeBatch> positives;
  unsigned num_probes = 0u;

  auto consume = [&]() {
    unsigned num_positive = 0u;
    for (unsigned probe_id = 0; probe_id < num_probes; ++probe_id) {
      EndpointProbe &probe = probes[probe_id];
      const uint64_t bits = SeedFilterBits(probe.hash);
      if ((filter[probe.hash & endpoint_seed_filter_mask_] & bits) != bits) {
        continue;
      }
      index_[IndexShard(probe.hash)].prefetch_hash(probe.hash);
      positives[num_positive++] = static_cast<uint8_t>(probe_id);
    }
    for (unsigned positive_id = 0; positive_id < num_positive;
         ++positive_id) {
      const EndpointProbe &probe = probes[positives[positive_id]];
      const auto &shard = index_[IndexShard(probe.hash)];
      const auto found = shard.find(probe.key, probe.hash);
      if (found == shard.end() || (found->second >> 63u) != 0u) continue;
      uint32_t contig_id, contig_offset;
      uint8_t strand;
      DecodeContigOffset(found->second, contig_id, contig_offset, strand);
      const uint64_t contig_length = refseq_.GetSeqView(contig_id).length();
      if (contig_offset < static_cast<uint32_t>(endpoint_range_) ||
          static_cast<uint64_t>(contig_offset) + endpoint_range_ >=
              contig_length) {
        if (witness != nullptr) {
          witness->index_value = found->second;
          witness->end_position = probe.end_position;
          witness->query_strand = probe.query_strand;
          witness->valid = true;
        }
        return true;
      }
    }
    num_probes = 0u;
    return false;
  };

  for (unsigned end_position = seed_kmer_size_ - 1u;
       end_position < length; ++end_position) {
    if (end_position >= static_cast<unsigned>(seed_kmer_size_)) {
      const unsigned absolute = query_shift + end_position;
      const uint8_t base = static_cast<uint8_t>(
          (query_words[absolute / SeqPackage::kBasesPerWord] >>
           SeqPackage::TVector::bit_shift(
               absolute % SeqPackage::kBasesPerWord)) &
          3u);
      seed_forward =
          (seed_forward << 2u) | (uint64_t(base) << seed_padding);
      seed_reverse =
          ((seed_reverse >> 2u) | (uint64_t(3u - base) << 62u)) & seed_mask;
    }
    EndpointProbe &probe = probes[num_probes++];
    probe.query_strand = seed_forward <= seed_reverse ? 0u : 1u;
    probe.key = probe.query_strand == 0u ? seed_forward : seed_reverse;
    probe.hash = IndexHash(probe.key);
    probe.end_position = static_cast<int32_t>(end_position);
    __builtin_prefetch(
        filter + (probe.hash & endpoint_seed_filter_mask_), 0, 3);
    if (num_probes == kProbeBatch && consume()) return true;
  }
  return num_probes != 0u && consume();
}

int32_t HashMapper::Match(const uint32_t *query_first_word,
                          unsigned query_shift, int query_from, int query_to,
                          size_t contig_id, int ref_from, int ref_to,
                          bool strand) const {
  auto contig_view = refseq_.GetSeqView(contig_id);
  auto ref_ptr_and_offset = contig_view.raw_address();
  const uint32_t *ref_first_word = ref_ptr_and_offset.first;
  int ref_shift = ref_ptr_and_offset.second;

  int match_len = query_to - query_from + 1;
  int threshold = lround(similarity_ * match_len);

  for (int i = query_from; i <= query_to; i += 32) {
    int len = std::min(32, query_to - i + 1);
    uint64_t qw = GetWord64(query_first_word, query_shift, i, len, 0);
    int ref_i = strand == 0 ? ref_from + i - query_from
                            : ref_to - (i + len - 1 - query_from);
    uint64_t rw = GetWord64(ref_first_word, ref_shift, ref_i, len, strand);

    match_len -= Mismatch64(qw, rw);

    if (match_len < threshold) {
      return 0;
    }
  }

  return match_len;
}

MappingRecord HashMapper::TryMap(
    const SeqPackage::SeqView &seq_view,
    const EndpointSeedWitness *witness) const {
  auto address = seq_view.raw_address();
  return TryMapRaw(address.first, address.second, seq_view.length(),
                   seq_view.id(), witness);
}

MappingRecord HashMapper::TryMap(const uint32_t *packed_words,
                                 unsigned length,
                                 uint64_t query_id,
                                 const EndpointSeedWitness *witness) const {
  return TryMapRaw(packed_words, 0, length, query_id, witness);
}

MappingRecord HashMapper::TryMapRaw(const uint32_t *query_words,
                                    unsigned query_shift, unsigned length,
                                    uint64_t query_id,
                                    const EndpointSeedWitness *witness) const {
  MappingRecord bad_record;
  bad_record.valid = false;

  const int len = static_cast<int>(length);
  if (len < seed_kmer_size_ || len < 50)
    return bad_record;  // too short reads not reliable

  // small vector optimization
  static const int kArraySize = 3;
  std::array<MappingRecord, kArraySize> mapping_records;
  int n_mapping_records = 0;
  std::unique_ptr<std::vector<MappingRecord>> v_mapping_records;

  // The mapper seed is at most 32 bases and its historical two-uint32 Kmer
  // representation is exactly one big-endian uint64_t.  Keep forward and
  // reverse-complement seeds in that scalar form: rolling, canonical compare,
  // and conversion to the table key then need no word loop or temporary Kmer.
  TKmer initial_seed(query_words, query_shift, seed_kmer_size_);
  uint64_t seed_forward = SeedKey(initial_seed);
  const unsigned seed_padding = 64u - (seed_kmer_size_ << 1u);
  const uint64_t seed_mask =
      seed_padding == 0u ? UINT64_MAX : UINT64_MAX << seed_padding;
  uint64_t seed_reverse =
      kmlib::bit::ReverseComplement<2>(seed_forward);
  if (seed_padding != 0u) seed_reverse <<= seed_padding;
  const uint64_t *const local_seed_filter = LocalSeedFilter();

  struct CandidateDiagonal {
    const uint32_t *ref_words;
    uint32_t ref_shift;
    int64_t diagonal;
    int32_t contig_length;
    uint32_t indexed_end_phase;
    uint8_t strand;
  };
  std::array<CandidateDiagonal, kArraySize> candidate_diagonals;
  const uint32_t sparsity_mask =
      static_cast<uint32_t>(index_sparsity_ - 1);
  const bool power_of_two_sparsity =
      (static_cast<uint32_t>(index_sparsity_) & sparsity_mask) == 0u;

  auto process_seed_value = [&](uint64_t value, uint8_t query_strand, int i) {
    uint32_t contig_id, contig_offset;
    uint8_t contig_strand;
    DecodeContigOffset(value, contig_id, contig_offset, contig_strand);

    auto contig_view = refseq_.GetSeqView(contig_id);
    assert(contig_offset < contig_view.length());

    uint8_t mapping_strand = contig_strand ^ query_strand;
    int32_t contig_from = mapping_strand == 0
                              ? contig_offset - (i - seed_kmer_size_ + 1)
                              : contig_offset - (len - 1 - i);
    int32_t contig_to = mapping_strand == 0
                            ? contig_offset + seed_kmer_size_ - 1 + len - 1 - i
                            : contig_offset + i;
    contig_from = std::max(contig_from, 0);
    contig_to =
        std::min(static_cast<int32_t>(contig_view.length() - 1), contig_to);

    if (contig_to - contig_from + 1 < len &&
        contig_to - contig_from + 1 < min_mapped_len_) {
      return;  // clipped alignment is considered iff its length >=
      // min_mapped_len_
    }

    int32_t query_from = mapping_strand == 0 ? i - (seed_kmer_size_ - 1) -
                                                   (contig_offset - contig_from)
                                             : i - (contig_to - contig_offset);
    int32_t query_to = mapping_strand == 0 ? i - (seed_kmer_size_ - 1) +
                                                 (contig_to - contig_offset)
                                           : i + (contig_offset - contig_from);

    assert(query_from >= 0 &&
           static_cast<uint32_t>(query_from) < length);
    assert(query_to >= 0 &&
           static_cast<uint32_t>(query_to) < length);

    auto rec = MappingRecord{contig_id,  contig_from,
                             contig_to,  query_id,
                             query_from, query_to,
                             0,          mapping_strand,
                             true};
    auto end = mapping_records.begin() + n_mapping_records;
    if (std::find(mapping_records.begin(), end, rec) == end) {
      if (n_mapping_records < kArraySize) {
        const auto ref_address = contig_view.raw_address();
        const int64_t diagonal =
            mapping_strand == 0
                ? static_cast<int64_t>(contig_from) - query_from
                : static_cast<int64_t>(contig_to) + query_from;
        int64_t indexed_end_phase =
            mapping_strand == 0
                ? static_cast<int64_t>(seed_kmer_size_ - 1) - diagonal
                : diagonal;
        indexed_end_phase %= index_sparsity_;
        if (indexed_end_phase < 0) indexed_end_phase += index_sparsity_;
        mapping_records[n_mapping_records] = rec;
        candidate_diagonals[n_mapping_records] = CandidateDiagonal{
            ref_address.first,
            ref_address.second,
            diagonal,
            static_cast<int32_t>(contig_view.length()),
            static_cast<uint32_t>(indexed_end_phase),
            mapping_strand};
        ++n_mapping_records;
      } else {
        if (v_mapping_records.get() == nullptr) {
          v_mapping_records.reset(new std::vector<MappingRecord>(1, rec));
        } else {
          v_mapping_records->push_back(rec);
        }
      }
    }
  };

  // The endpoint gate already performed an authoritative table lookup.  Seed
  // the exact mapper with that same candidate before replaying the read.  The
  // normal diagonal proof then suppresses every redundant seed on this
  // alignment while all unexplained seeds are still probed in historical
  // order, so alternative mappings and tie semantics are unchanged.
  // Keep an exact A/B path for regression tests.  The optimized path is the
  // default; disabling it only skips reuse of the already-proven endpoint
  // seed and otherwise executes the historical mapper unchanged.
  if (std::getenv("MEGAHIT_DISABLE_ENDPOINT_WITNESS") == nullptr &&
      witness != nullptr && witness->valid) {
    process_seed_value(witness->index_value, witness->query_strand,
                       witness->end_position);
  }

  // The blocked membership filter is deliberately much smaller than the
  // contig index, but at CAMI scale it is still a random 128 MiB working set.
  // A scalar seed loop carries only one outstanding miss.  Form a small
  // software pipeline: compute independent rolling seeds, prefetch their
  // filter words, then consume them in the original order.  This changes no
  // candidate or tie semantics while exposing memory-level parallelism.
  struct SeedProbe {
    TSeedKey key;
    uint64_t forward_key;
    size_t hash;
    int end_position;
    uint8_t strand;
    uint8_t checked_candidate_count;
  };

  auto seed_is_proven_redundant = [&](const SeedProbe &probe) {
    for (int candidate_id = 0; candidate_id < n_mapping_records;
         ++candidate_id) {
      const CandidateDiagonal &candidate =
          candidate_diagonals[candidate_id];
      const uint32_t phase =
          power_of_two_sparsity
              ? static_cast<uint32_t>(probe.end_position) & sparsity_mask
              : static_cast<uint32_t>(probe.end_position) % index_sparsity_;
      if (phase != candidate.indexed_end_phase) continue;

      const int64_t ref_begin =
          candidate.strand == 0
              ? candidate.diagonal + probe.end_position -
                    (seed_kmer_size_ - 1)
              : candidate.diagonal - probe.end_position;
      if (ref_begin < 0 ||
          ref_begin + seed_kmer_size_ > candidate.contig_length) {
        continue;
      }
      const uint64_t expected =
          GetWord64(candidate.ref_words, candidate.ref_shift,
                    static_cast<int>(ref_begin), seed_kmer_size_,
                    candidate.strand != 0u);
      if (expected == probe.forward_key) return true;
    }
    return false;
  };

  constexpr unsigned kSeedProbeBatch = 32;
  std::array<SeedProbe, kSeedProbeBatch> probes;
  std::array<uint8_t, kSeedProbeBatch> positive_probe_ids;
  unsigned num_probes = 0;
  auto consume_probes = [&]() {
    unsigned num_positive = 0;
    for (unsigned probe_id = 0; probe_id < num_probes; ++probe_id) {
      SeedProbe &probe = probes[probe_id];
      if (!SeedMayContain(probe.hash, local_seed_filter)) {
        continue;
      }
      if (seed_is_proven_redundant(probe)) continue;
      probe.checked_candidate_count =
          static_cast<uint8_t>(n_mapping_records);
      positive_probe_ids[num_positive++] = static_cast<uint8_t>(probe_id);
      // The first pass consumes the already-prefetched Bloom words, then
      // launches independent control-byte and slot prefetches into the
      // appropriate shard.  The second pass below retains original seed
      // order, but no longer serializes every positive probe on DRAM.
      index_[IndexShard(probe.hash)].prefetch_hash(probe.hash);
    }
    for (unsigned positive_id = 0; positive_id < num_positive;
         ++positive_id) {
      const SeedProbe &probe = probes[positive_probe_ids[positive_id]];
      // A second diagonal check is needed only when an earlier positive in
      // this same batch discovered a new candidate.  Previously every
      // positive repeated the reference gather even when the candidate set
      // was unchanged.
      if (probe.checked_candidate_count !=
              static_cast<uint8_t>(n_mapping_records) &&
          seed_is_proven_redundant(probe)) {
        continue;
      }
      const auto &index = index_[IndexShard(probe.hash)];
      const auto iter = index.find(probe.key, probe.hash);
      if (iter == index.end() || (iter->second >> 63u) != 0u) continue;
      process_seed_value(iter->second, probe.strand, probe.end_position);
    }
    num_probes = 0;
  };

  for (int i = seed_kmer_size_ - 1; i < len; ++i) {
    if (i >= seed_kmer_size_) {
      const unsigned absolute = query_shift + static_cast<unsigned>(i);
      const uint8_t ch = static_cast<uint8_t>(
          (query_words[absolute / SeqPackage::kBasesPerWord] >>
           SeqPackage::TVector::bit_shift(
               absolute % SeqPackage::kBasesPerWord)) &
          3u);
      seed_forward =
          (seed_forward << 2u) | (uint64_t(ch) << seed_padding);
      seed_reverse =
          ((seed_reverse >> 2u) | (uint64_t(3u - ch) << 62u)) &
          seed_mask;
    }

    const uint8_t query_strand = seed_forward <= seed_reverse ? 0u : 1u;
    SeedProbe &probe = probes[num_probes++];
    probe.key = query_strand == 0u ? seed_forward : seed_reverse;
    probe.forward_key = seed_forward;
    probe.hash = IndexHash(probe.key);
    probe.end_position = i;
    probe.strand = query_strand;
    probe.checked_candidate_count = 0u;
    if (local_seed_filter != nullptr) {
      __builtin_prefetch(
          local_seed_filter + (probe.hash & seed_filter_mask_), 0, 3);
    }
    if (num_probes == kSeedProbeBatch) consume_probes();
  }
  if (num_probes != 0) consume_probes();

  if (n_mapping_records == 0) {
    return bad_record;
  }

  MappingRecord *best = &bad_record;
  int32_t max_match = 0;

#define CHECK_BEST_UNIQ(rec)                                              \
  do {                                                                    \
    int32_t match_bases =                                                 \
        Match(query_words, query_shift, rec.query_from, rec.query_to,     \
              rec.contig_id, rec.contig_from, rec.contig_to, rec.strand); \
    if (match_bases == max_match) {                                       \
      best = &bad_record;                                                 \
    } else if (match_bases > max_match) {                                 \
      max_match = match_bases;                                            \
      int32_t mismatch = rec.query_to - rec.query_from + 1 - match_bases; \
      rec.mismatch = mismatch;                                            \
      best = &rec;                                                        \
    }                                                                     \
  } while (0)

  if (v_mapping_records.get() != nullptr) {
    if (v_mapping_records->size() > 1) {
      std::sort(v_mapping_records->begin(), v_mapping_records->end());
      v_mapping_records->resize(
          std::unique(v_mapping_records->begin(), v_mapping_records->end()) -
          v_mapping_records->begin());
    }

    for (auto &rec : *v_mapping_records) {
      CHECK_BEST_UNIQ(rec);
    }
  }

  for (int i = 0; i < n_mapping_records; ++i) {
    auto &rec = mapping_records[i];
    CHECK_BEST_UNIQ(rec);
  }

#undef CHECK_BEST_UNIQ

  return *best;
}
