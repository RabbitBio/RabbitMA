/**
 * @file hash_graph.cpp
 * @brief
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#include "idba/hash_graph.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <limits>
#include <stdexcept>

#include <omp.h>
#include "bit_operation.h"
#include "idba/contig_builder.h"
#include "idba/contig_info.h"
#include "idba/hash_graph_vertex.h"
#include "idba/kmer.h"
#include "idba/sequence.h"
#include "kmlib/kmsort.h"
#include "utils/histgram.h"

using namespace std;

#include <iostream>

namespace {

inline uint32_t EncodeEndpoint(uint32_t vertex_index, bool is_reverse) {
  if (vertex_index >= (UINT32_MAX - 1u) / 2u) {
    return HashGraphVertex::kUncacheableNeighbor;
  }
  return (vertex_index << 1u) | uint32_t(is_reverse);
}

inline void CacheUniqueNext(HashGraphVertexAdaptor &from,
                            uint32_t neighbor_code) {
  uint32_t &slot = from.next_neighbor();
  // Record transition agreement only.  GetNextVertexAdaptor already checks
  // the final out-degree before consulting this cache, so doing two popcounts
  // for every read occurrence is redundant.  A vertex whose degree becomes
  // non-unique cannot become unique again during HashGraph construction.
  if (neighbor_code == HashGraphVertex::kUncacheableNeighbor) {
    slot = HashGraphVertex::kUncacheableNeighbor;
  } else if (slot == HashGraphVertex::kNoCachedNeighbor) {
    slot = neighbor_code;
  } else if (slot != neighbor_code) {
    slot = HashGraphVertex::kUncacheableNeighbor;
  }
}

inline void CacheTransition(HashGraphVertexAdaptor &previous,
                            uint32_t previous_index, bool previous_reverse,
                            HashGraphVertexAdaptor &current,
                            uint32_t current_index, bool current_reverse) {
  const uint32_t current_code =
      EncodeEndpoint(current_index, current_reverse);
  CacheUniqueNext(previous, current_code);
  HashGraphVertexAdaptor reverse_current = current;
  reverse_current.ReverseComplement();
  const uint32_t reverse_previous_code =
      EncodeEndpoint(previous_index, !previous_reverse);
  CacheUniqueNext(reverse_current, reverse_previous_code);
}

}  // namespace

int64_t HashGraph::InsertKmers(const Sequence &seq) {
  if (seq.size() < kmer_size_) return 0;

  IdbaKmer kmer(kmer_size_), rev_kmer(kmer_size_);
  int length = 0;
  int64_t num_kmers = 0;
  uint32_t previous_index = vertex_table_type::kNull;
  bool previous_reverse = false;
  for (uint64_t i = 0; i < seq.size(); ++i) {
    const uint8_t base = seq[i];
    kmer.ShiftAppend(base);
    rev_kmer.ShiftPreappend(3u - (base & 3u));
    length = (base < 4) ? length + 1 : 0;

    if (length < (int)kmer_size_) {
      previous_index = vertex_table_type::kNull;
      continue;
    }

    const bool is_reverse = rev_kmer < kmer;
    const IdbaKmer &key = is_reverse ? rev_kmer : kmer;
    uint32_t vertex_index;
    HashGraphVertex &vertex =
        vertex_table_.find_or_insert_key(key, &vertex_index);
    vertex.count() += 1;
    HashGraphVertexAdaptor adaptor(&vertex, is_reverse);

    if (length > (int)kmer_size_ && seq[i - kmer_size_] < 4)
      adaptor.in_edges().Add(3 - seq[i - kmer_size_]);
    if (i + 1 < seq.size() && seq[i + 1] < 4)
      adaptor.out_edges().Add(seq[i + 1]);

    if (previous_index != vertex_table_type::kNull) {
      HashGraphVertexAdaptor previous(
          &vertex_table_.value_at(previous_index), previous_reverse);
      CacheTransition(previous, previous_index, previous_reverse, adaptor,
                      vertex_index, is_reverse);
    }
    previous_index = vertex_index;
    previous_reverse = is_reverse;

    ++num_kmers;
  }

  return num_kmers;
}

namespace {

inline uint8_t PackedBaseAt(const uint64_t *words, uint32_t index) {
  return static_cast<uint8_t>((words[index >> 5u] >>
                               ((index & 31u) << 1u)) & 3u);
}

template <unsigned Words>
inline void LoadPackedKmer(const uint64_t *source, uint32_t start,
                           uint32_t kmer_size, uint64_t *output) {
  const uint32_t source_word = start >> 5u;
  const uint32_t shift = (start & 31u) << 1u;
  for (unsigned word = 0; word < Words; ++word) {
    uint64_t value = source[source_word + word] >> shift;
    if (shift != 0u) {
      value |= source[source_word + word + 1u] << (64u - shift);
    }
    output[word] = value;
  }
  const uint32_t tail = kmer_size & 31u;
  if (tail != 0u) {
    output[Words - 1u] &= (uint64_t{1} << (tail << 1u)) - 1u;
  }
}

template <unsigned Words>
inline void LoadPackedKmerPair(const uint64_t *source, uint32_t start,
                               uint32_t kmer_size, uint64_t *forward,
                               uint64_t *reverse) {
  LoadPackedKmer<Words>(source, start, kmer_size, forward);
  const uint32_t tail = kmer_size & 31u;

  for (unsigned word = 0; word < Words; ++word) {
    uint64_t value = forward[word];
    bit_operation::ReverseComplement(value);
    reverse[Words - 1u - word] = value;
  }
  if (tail != 0u) {
    const uint32_t align_shift = (32u - tail) << 1u;
    for (unsigned word = 0; word + 1u < Words; ++word) {
      reverse[word] = (reverse[word] >> align_shift) |
                      (reverse[word + 1u] << (64u - align_shift));
    }
    reverse[Words - 1u] >>= align_shift;
  }
}

template <unsigned Words>
struct CompactKmerOccurrence {
  uint64_t key[Words];
  uint32_t ordinal;
  uint32_t count;
  uint8_t in_edges;
  uint8_t out_edges;
  uint8_t flags;

  enum : uint8_t { kReverse = 1u, kHasNext = 2u };

  static const int n_bytes = Words * sizeof(uint64_t);

  int kth_byte(int byte) const {
    return static_cast<int>(
        (key[byte >> 3u] >> ((byte & 7u) << 3u)) & 0xFFu);
  }

  bool operator<(const CompactKmerOccurrence &other) const {
    for (int word = Words - 1; word >= 0; --word) {
      if (key[word] != other.key[word]) return key[word] < other.key[word];
    }
    return false;
  }

  bool SameKey(const CompactKmerOccurrence &other) const {
    for (unsigned word = 0; word < Words; ++word) {
      if (key[word] != other.key[word]) return false;
    }
    return true;
  }

  void SetKey(const IdbaKmer &kmer) {
    for (unsigned word = 0; word < Words; ++word) key[word] = kmer.word(word);
  }

  IdbaKmer GetKey(uint32_t kmer_size) const {
    IdbaKmer result;
    result.AssignWords(key, kmer_size);
    return result;
  }
};

template <unsigned Words>
struct OccurrenceOrdinalRadixTraits {
  static const int n_bytes = sizeof(uint32_t);
  int kth_byte(const CompactKmerOccurrence<Words> &record, int byte) const {
    return static_cast<int>(
        (record.ordinal >> (byte << 3u)) & UINT32_C(0xFF));
  }
  bool operator()(const CompactKmerOccurrence<Words> &lhs,
                  const CompactKmerOccurrence<Words> &rhs) const {
    return lhs.ordinal < rhs.ordinal;
  }
};

struct CompactHashOrder {
  uint32_t hash;
  uint32_t ordinal;

  static const int n_bytes = sizeof(uint32_t);
  int kth_byte(int byte) const {
    return static_cast<int>((hash >> (byte << 3u)) & UINT32_C(0xFF));
  }
  bool operator<(const CompactHashOrder &other) const {
    return hash < other.hash;
  }
};

struct AggregateNeighborCandidates {
  static constexpr uint64_t kNone = UINT64_MAX;
  static constexpr uint64_t kAmbiguous = UINT64_MAX - 1u;
  uint64_t next[2]{kNone, kNone};

  void Add(bool strand, uint64_t candidate) {
    uint64_t &slot = next[strand ? 1 : 0];
    if (slot == kNone)
      slot = candidate;
    else if (slot != candidate)
      slot = kAmbiguous;
  }
};

template <unsigned Words>
uint32_t CompactLookupHash(const CompactKmerOccurrence<Words> &record) {
  const uint64_t hash =
      XXH3_64bits(record.key, Words * sizeof(uint64_t));
  return static_cast<uint32_t>(hash ^ (hash >> 32u));
}

template <>
uint32_t CompactLookupHash<1>(const CompactKmerOccurrence<1> &record) {
  return static_cast<uint32_t>(record.key[0] ^ (record.key[0] >> 32u));
}

template <unsigned Words>
int64_t InsertPackedKmersBulkImpl(HashGraph *graph,
                                  const HashGraphPackedRead *reads,
                                  size_t num_reads,
                                  HashGraphBulkStats *stats) {
  using Record = CompactKmerOccurrence<Words>;
  const uint32_t kmer_size = graph->kmer_size();
  uint64_t total_occurrences = 0;
  for (size_t read_id = 0; read_id < num_reads; ++read_id) {
    if (reads[read_id].length >= kmer_size) {
      total_occurrences += reads[read_id].length - kmer_size + 1u;
    }
  }
  if (total_occurrences == 0) return 0;
  if (total_occurrences > std::numeric_limits<uint32_t>::max() ||
      total_occurrences > std::numeric_limits<size_t>::max()) {
    throw std::length_error("local bulk k-mer record stream is too large");
  }

  std::vector<Record> records;
  records.resize(static_cast<size_t>(total_occurrences));
  std::vector<CompactHashOrder> hash_order;
  hash_order.resize(static_cast<size_t>(total_occurrences));
  const double generate_begin = omp_get_wtime();
  uint32_t ordinal = 0;
  for (size_t read_id = 0; read_id < num_reads; ++read_id) {
    const uint64_t *forward_words = reads[read_id].forward_words;
    const uint32_t sequence_length = reads[read_id].length;
    if (sequence_length < kmer_size) continue;

    IdbaKmer kmer;
    IdbaKmer rev_kmer;
    const uint32_t num_kmers = sequence_length - kmer_size + 1u;
    kmer.AssignPackedBases(forward_words, 0, kmer_size);
    rev_kmer = kmer;
    rev_kmer.ReverseComplement();
    for (uint32_t start = 0; start < num_kmers; ++start) {
      const bool is_reverse = rev_kmer < kmer;
      const IdbaKmer &key = is_reverse ? rev_kmer : kmer;
      Record &record = records[static_cast<size_t>(ordinal)];
      record.SetKey(key);
      record.ordinal = ordinal;
      record.count = 1;
      record.in_edges = 0;
      record.out_edges = 0;
      record.flags = static_cast<uint8_t>(is_reverse ? Record::kReverse : 0u);

      if (start != 0) {
        const uint8_t edge = static_cast<uint8_t>(
            1u << (3u - PackedBaseAt(forward_words, start - 1u)));
        if (is_reverse)
          record.out_edges |= edge;
        else
          record.in_edges |= edge;
      }
      if (start + kmer_size < sequence_length) {
        const uint8_t edge = static_cast<uint8_t>(
            1u << PackedBaseAt(forward_words, start + kmer_size));
        if (is_reverse)
          record.in_edges |= edge;
        else
          record.out_edges |= edge;
      }
      if (start + 1u < num_kmers) record.flags |= Record::kHasNext;
      hash_order[static_cast<size_t>(ordinal)] =
          {CompactLookupHash(record), ordinal};
      ++ordinal;

      if (start + 1u < num_kmers) {
        const uint8_t next_base =
            PackedBaseAt(forward_words, start + kmer_size);
        kmer.ShiftAppend(next_base);
        rev_kmer.ShiftPreappend(3u - next_base);
      }
    }
  }
  if (stats != nullptr) stats->generate_seconds += omp_get_wtime() - generate_begin;

  double phase_begin = omp_get_wtime();
  kmlib::kmsort(hash_order.begin(), hash_order.end());
  if (stats != nullptr) stats->key_sort_seconds += omp_get_wtime() - phase_begin;

  phase_begin = omp_get_wtime();
  std::vector<Record> aggregates;
  aggregates.reserve(records.size());
  std::vector<uint32_t> occurrence_group(records.size());
  for (size_t hash_begin = 0; hash_begin < hash_order.size();) {
    size_t hash_end = hash_begin + 1u;
    while (hash_end < hash_order.size() &&
           hash_order[hash_end].hash == hash_order[hash_begin].hash) {
      ++hash_end;
    }

    Record aggregate =
        records[static_cast<size_t>(hash_order[hash_begin].ordinal)];
    bool one_exact_key = true;
    for (size_t i = hash_begin + 1u; i < hash_end; ++i) {
      const Record &record =
          records[static_cast<size_t>(hash_order[i].ordinal)];
      if (!aggregate.SameKey(record)) {
        one_exact_key = false;
        break;
      }
      aggregate.count += record.count;
      aggregate.in_edges |= record.in_edges;
      aggregate.out_edges |= record.out_edges;
      aggregate.ordinal = std::min(aggregate.ordinal, record.ordinal);
    }

    if (one_exact_key) {
      const uint32_t aggregate_id =
          static_cast<uint32_t>(aggregates.size());
      for (size_t i = hash_begin; i < hash_end; ++i) {
        occurrence_group[hash_order[i].ordinal] = aggregate_id;
      }
      aggregates.push_back(aggregate);
    } else {
      // XXH3 collisions are rare, but correctness must not depend on that.
      // Sort only the colliding hash run by the complete active key.
      std::sort(hash_order.begin() + hash_begin,
                hash_order.begin() + hash_end,
                [&records](const CompactHashOrder &lhs,
                           const CompactHashOrder &rhs) {
                  return records[static_cast<size_t>(lhs.ordinal)] <
                         records[static_cast<size_t>(rhs.ordinal)];
                });
      for (size_t key_begin = hash_begin; key_begin < hash_end;) {
        size_t key_end = key_begin + 1u;
        const Record &key_record =
            records[static_cast<size_t>(hash_order[key_begin].ordinal)];
        while (key_end < hash_end && key_record.SameKey(
                                        records[static_cast<size_t>(
                                            hash_order[key_end].ordinal)])) {
          ++key_end;
        }

        Record collision_aggregate = key_record;
        for (size_t i = key_begin + 1u; i < key_end; ++i) {
          const Record &record =
              records[static_cast<size_t>(hash_order[i].ordinal)];
          collision_aggregate.count += record.count;
          collision_aggregate.in_edges |= record.in_edges;
          collision_aggregate.out_edges |= record.out_edges;
          collision_aggregate.ordinal =
              std::min(collision_aggregate.ordinal, record.ordinal);
        }
        const uint32_t aggregate_id =
            static_cast<uint32_t>(aggregates.size());
        for (size_t i = key_begin; i < key_end; ++i) {
          occurrence_group[hash_order[i].ordinal] = aggregate_id;
        }
        aggregates.push_back(collision_aggregate);
        key_begin = key_end;
      }
    }
    hash_begin = hash_end;
  }

  std::vector<AggregateNeighborCandidates> neighbor_candidates(
      aggregates.size());
  // Reduction creates aggregates in hash order.  Recover historical
  // first-occurrence order during the occurrence scan that adjacency needs
  // anyway, instead of materializing and radix-sorting another U-entry
  // permutation.  aggregate.ordinal is the exact first occurrence, so this
  // adds only one predictable comparison per occurrence.
  std::vector<uint32_t> first_order;
  first_order.reserve(aggregates.size());
  for (uint32_t occurrence = 0; occurrence < records.size(); ++occurrence) {
    const Record &record = records[occurrence];
    const uint32_t from_group = occurrence_group[occurrence];
    if (aggregates[from_group].ordinal == occurrence) {
      first_order.push_back(from_group);
    }
    if ((record.flags & Record::kHasNext) == 0) continue;
    const uint32_t next_occurrence = occurrence + 1u;
    const uint32_t to_group = occurrence_group[next_occurrence];
    const bool from_reverse = (record.flags & Record::kReverse) != 0;
    const bool to_reverse =
        (records[next_occurrence].flags & Record::kReverse) != 0;
    neighbor_candidates[from_group].Add(
        from_reverse, (uint64_t(to_group) << 1u) | uint64_t(to_reverse));
    neighbor_candidates[to_group].Add(
        !to_reverse,
        (uint64_t(from_group) << 1u) | uint64_t(!from_reverse));
  }
  records.swap(aggregates);
  const size_t aggregate_count = records.size();
  std::vector<Record>().swap(aggregates);
  std::vector<CompactHashOrder>().swap(hash_order);
  std::vector<uint32_t>().swap(occurrence_group);
  if (stats != nullptr) stats->reduce_seconds += omp_get_wtime() - phase_begin;

  phase_begin = omp_get_wtime();
  std::vector<uint32_t> vertex_indices(records.size());
  for (const uint32_t aggregate : first_order) {
    const Record &record = records[aggregate];
    graph->InsertAggregate(record.GetKey(kmer_size), record.count,
                           record.in_edges, record.out_edges,
                           &vertex_indices[aggregate]);
  }
  graph->FinishBulkOccurrences(
      !first_order.empty() &&
          records[first_order.back()].ordinal + 1u < total_occurrences);
  for (uint32_t aggregate = 0; aggregate < records.size(); ++aggregate) {
    for (unsigned strand = 0; strand < 2; ++strand) {
      const uint64_t candidate = neighbor_candidates[aggregate].next[strand];
      if (candidate == AggregateNeighborCandidates::kNone) continue;
      if (candidate == AggregateNeighborCandidates::kAmbiguous) {
        // Occurrence-wise insertion permanently invalidates the transition
        // cache as soon as a second distinct successor is observed.  The
        // reduced graph has the same edges either way, but retaining the
        // historical cache state is required both for the exact traversal
        // fast path and for tie-equivalent validation.  An out-of-range
        // neighbor code is deliberately converted to kUncacheableNeighbor by
        // SetAggregateNeighbor().
        graph->SetAggregateNeighbor(vertex_indices[aggregate], strand != 0,
                                    UINT32_MAX, false);
        continue;
      }
      const uint32_t neighbor_group = static_cast<uint32_t>(candidate >> 1u);
      graph->SetAggregateNeighbor(vertex_indices[aggregate], strand != 0,
                                  vertex_indices[neighbor_group],
                                  (candidate & 1u) != 0);
    }
  }
  if (stats != nullptr) {
    stats->replay_seconds += omp_get_wtime() - phase_begin;
    stats->occurrences += total_occurrences;
    stats->aggregates += aggregate_count;
  }
  return static_cast<int64_t>(total_occurrences);
}

}  // namespace

int64_t HashGraph::InsertPackedKmersBulk(const HashGraphPackedRead *reads,
                                         size_t num_reads,
                                         HashGraphBulkStats *stats) {
  const uint32_t words = (kmer_size_ + 31u) >> 5u;
  switch (words) {
    case 1: return InsertPackedKmersBulkImpl<1>(this, reads, num_reads, stats);
    case 2: return InsertPackedKmersBulkImpl<2>(this, reads, num_reads, stats);
    case 3: return InsertPackedKmersBulkImpl<3>(this, reads, num_reads, stats);
    case 4: return InsertPackedKmersBulkImpl<4>(this, reads, num_reads, stats);
    case 5: return InsertPackedKmersBulkImpl<5>(this, reads, num_reads, stats);
    case 6: return InsertPackedKmersBulkImpl<6>(this, reads, num_reads, stats);
    case 7: return InsertPackedKmersBulkImpl<7>(this, reads, num_reads, stats);
    case 8: return InsertPackedKmersBulkImpl<8>(this, reads, num_reads, stats);
    default:
      throw std::length_error("unsupported local bulk k-mer width");
  }
}

void HashGraph::InsertAggregate(const IdbaKmer &key, uint32_t count,
                                uint8_t in_edges, uint8_t out_edges,
                                uint32_t *vertex_index) {
  HashGraphVertex &vertex =
      vertex_table_.find_or_insert_key(key, vertex_index);
  vertex.count() += static_cast<int32_t>(count);
  vertex.in_edges() =
      static_cast<uint8_t>(vertex.in_edges()) | in_edges;
  vertex.out_edges() =
      static_cast<uint8_t>(vertex.out_edges()) | out_edges;
  if (vertex.out_edges().size() > 1) {
    vertex.next_neighbor(false) = HashGraphVertex::kUncacheableNeighbor;
  }
  if (vertex.in_edges().size() > 1) {
    vertex.next_neighbor(true) = HashGraphVertex::kUncacheableNeighbor;
  }
}

void HashGraph::AccumulateResolvedOccurrence(uint32_t vertex_index,
                                             uint8_t in_edges,
                                             uint8_t out_edges) {
  // Raw insertion observes a pending resize at the beginning of every
  // occurrence, including a duplicate reached through the transition cache.
  // Vertex IDs do not move, so the already resolved index remains valid.
  vertex_table_.prepare_cached_occurrence();
  HashGraphVertex &vertex = vertex_table_.value_at(vertex_index);
  ++vertex.count();
  vertex.in_edges() =
      static_cast<uint8_t>(vertex.in_edges()) | in_edges;
  vertex.out_edges() =
      static_cast<uint8_t>(vertex.out_edges()) | out_edges;
  if (vertex.out_edges().size() > 1) {
    vertex.next_neighbor(false) = HashGraphVertex::kUncacheableNeighbor;
  }
  if (vertex.in_edges().size() > 1) {
    vertex.next_neighbor(true) = HashGraphVertex::kUncacheableNeighbor;
  }
}

bool HashGraph::ReplayResolvedTransition(uint32_t from_code, uint8_t base,
                                         uint32_t *to_code) {
  const uint32_t from_index = from_code >> 1u;
  const bool from_reverse = (from_code & 1u) != 0u;
  HashGraphVertexAdaptor from(&vertex_table_.value_at(from_index),
                              from_reverse);
  base &= 3u;

  if (!from.out_edges()[base]) {
    const bool already_had_out_edge = !from.out_edges().empty();
    from.out_edges().Add(base);
    if (already_had_out_edge) {
      from.next_neighbor() = HashGraphVertex::kUncacheableNeighbor;
    }
    return false;
  }

  uint32_t resolved = from.next_neighbor();
  if (resolved == HashGraphVertex::kUncacheableNeighbor) {
    // A branch-table entry is an optimization hint owned by the ordinary
    // packed-reader path. Cross-k replay deliberately accepts only the much
    // stronger single-neighbor certificate; forked topology returns to the
    // exact prefix/suffix signature lookup below.
    return false;
  } else if (resolved >= HashGraphVertex::kUncacheableNeighbor) {
    return false;
  }

  vertex_table_.prepare_cached_occurrence();
  ++vertex_table_.value_at(resolved >> 1u).count();
  *to_code = resolved;
  return true;
}

void HashGraph::SetAggregateNeighbor(uint32_t vertex_index, bool strand,
                                     uint32_t neighbor_index,
                                     bool neighbor_strand) {
  HashGraphVertexAdaptor from(&vertex_table_.value_at(vertex_index), strand);
  if (from.out_edges().size() != 1) {
    from.next_neighbor() = HashGraphVertex::kUncacheableNeighbor;
    return;
  }
  from.next_neighbor() = EncodeEndpoint(neighbor_index, neighbor_strand);
}

void HashGraph::ObserveTransition(uint32_t from_code, uint32_t to_code) {
  const uint32_t from_index = from_code >> 1u;
  const uint32_t to_index = to_code >> 1u;
  HashGraphVertexAdaptor from(&vertex_table_.value_at(from_index),
                              (from_code & 1u) != 0u);
  HashGraphVertexAdaptor to(&vertex_table_.value_at(to_index),
                            (to_code & 1u) != 0u);
  CacheTransition(from, from_index, (from_code & 1u) != 0u, to, to_index,
                  (to_code & 1u) != 0u);
  if (from.next_neighbor() == HashGraphVertex::kUncacheableNeighbor) {
    RememberBranchTransition(from_index, (from_code & 1u) != 0u,
                             to.last_base(), to_code);
  }
  HashGraphVertexAdaptor reverse_to = to;
  reverse_to.ReverseComplement();
  if (reverse_to.next_neighbor() ==
      HashGraphVertex::kUncacheableNeighbor) {
    HashGraphVertexAdaptor reverse_from = from;
    reverse_from.ReverseComplement();
    RememberBranchTransition(to_index, (to_code & 1u) == 0u,
                             reverse_from.last_base(), from_code ^ 1u);
  }
}

void HashGraph::FinalizeAggregateTransitions() {
  for (uint32_t index = 0; index < vertex_table_.size(); ++index) {
    HashGraphVertex &vertex = vertex_table_.value_at(index);
    if (vertex.out_edges().size() > 1u) {
      vertex.next_neighbor(false) = HashGraphVertex::kUncacheableNeighbor;
    }
    if (vertex.in_edges().size() > 1u) {
      vertex.next_neighbor(true) = HashGraphVertex::kUncacheableNeighbor;
    }
  }
}

bool HashGraph::SameGraphState(const HashGraph &other,
                               std::string *difference) const {
  if (kmer_size_ != other.kmer_size_ ||
      vertex_table_.size() != other.vertex_table_.size()) {
    if (difference != nullptr) {
      *difference = "size " + std::to_string(vertex_table_.size()) +
                    " vs " + std::to_string(other.vertex_table_.size());
    }
    return false;
  }
  for (uint32_t i = 0; i < vertex_table_.size(); ++i) {
    const HashGraphVertex &vertex = vertex_table_.value_at(i);
    const HashGraphVertex *other_vertex =
        other.vertex_table_.find_value(vertex.key());
    if (other_vertex == nullptr || vertex.count() != other_vertex->count() ||
        uint8_t(vertex.in_edges()) != uint8_t(other_vertex->in_edges()) ||
        uint8_t(vertex.out_edges()) != uint8_t(other_vertex->out_edges())) {
      if (difference != nullptr) {
        *difference = "vertex " + std::to_string(i) + " count/edges " +
                      std::to_string(vertex.count()) + "/" +
                      std::to_string(uint8_t(vertex.in_edges())) + "/" +
                      std::to_string(uint8_t(vertex.out_edges()));
      }
      return false;
    }
  }
  return true;
}

bool HashGraph::SameTraversalState(const HashGraph &other,
                                   std::string *difference) const {
  if (!SameGraphState(other, difference)) return false;
  if (vertex_table_.bucket_count() != other.vertex_table_.bucket_count()) {
    if (difference != nullptr) {
      *difference = "buckets " +
                    std::to_string(vertex_table_.bucket_count()) + " vs " +
                    std::to_string(other.vertex_table_.bucket_count());
    }
    return false;
  }
  if (vertex_table_.rehash_pending() !=
      other.vertex_table_.rehash_pending()) {
    if (difference != nullptr) {
      *difference = std::string("pending rehash ") +
                    (vertex_table_.rehash_pending() ? "1" : "0") +
                    " vs " +
                    (other.vertex_table_.rehash_pending() ? "1" : "0");
    }
    return false;
  }
  std::vector<uint32_t> left;
  std::vector<uint32_t> right;
  vertex_table_.traversal_indices(&left);
  other.vertex_table_.traversal_indices(&right);
  for (size_t i = 0; i < left.size(); ++i) {
    const HashGraphVertex &lhs = vertex_table_.value_at(left[i]);
    const HashGraphVertex &rhs = other.vertex_table_.value_at(right[i]);
    if (!lhs.EqualsKey(rhs.key())) {
      if (difference != nullptr) {
        *difference = "traversal position " + std::to_string(i) +
                      " vertex IDs " + std::to_string(left[i]) + " vs " +
                      std::to_string(right[i]);
      }
      return false;
    }
  }
  return true;
}

bool HashGraph::SameTransitionCache(const HashGraph &other,
                                    std::string *difference) const {
  if (vertex_table_.size() != other.vertex_table_.size()) return false;
  for (uint32_t i = 0; i < vertex_table_.size(); ++i) {
    const HashGraphVertex &lhs = vertex_table_.value_at(i);
    const HashGraphVertex &rhs = other.vertex_table_.value_at(i);
    if (!lhs.EqualsKey(rhs.key())) {
      if (difference != nullptr) {
        *difference = "vertex insertion order at " + std::to_string(i);
      }
      return false;
    }
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      const uint8_t oriented_out = static_cast<uint8_t>(
          strand == 0u ? lhs.out_edges() : lhs.in_edges());
      if (lhs.next_neighbor(strand != 0u) !=
          rhs.next_neighbor(strand != 0u)) {
        if (difference != nullptr) {
          *difference = "cache vertex " + std::to_string(i) + " strand " +
                        std::to_string(strand) + " current " +
                        std::to_string(lhs.next_neighbor(strand != 0u)) +
                        " expected " +
                        std::to_string(rhs.next_neighbor(strand != 0u)) +
                        " edges " +
                        std::to_string(oriented_out);
        }
        return false;
      }
    }
  }
  return true;
}

template <unsigned Words, bool EmitGroups>
int64_t HashGraph::InsertPackedKmersImpl(const uint64_t *forward_words,
                                         uint32_t sequence_length,
                                         bool increment_count,
                                         uint32_t *forward_groups,
                                         uint32_t *reverse_groups) {
  if (sequence_length < kmer_size_) return 0;

  // Active-word key state. The common repeated-path case obtains the next
  // canonical vertex directly from the transition cache and never needs to
  // materialize this window. Branches/new paths load it directly from the
  // packed read, avoiding rolling work on every cache hit.
  uint64_t kmer[Words + 1u] = {};
  uint64_t rev_kmer[Words + 1u] = {};

  uint32_t previous_index = vertex_table_type::kNull;
  bool previous_reverse = false;
  bool force_exact_lookup = false;
  uint64_t branch_transition_hits = 0;
  uint64_t branch_transition_misses = 0;
  constexpr bool use_branch_transition_cache = true;
  static const bool use_cached_path_replay =
      std::getenv("MEGAHIT_DISABLE_LOCAL_CACHED_REPLAY") == nullptr;
  const uint32_t num_kmers = sequence_length - kmer_size_ + 1u;

  for (uint32_t start = 0; start < num_kmers; ++start) {
    uint32_t vertex_index = vertex_table_type::kNull;
    bool is_reverse = false;
    bool followed_cached_transition = false;

    // Once an oriented vertex has exactly one outgoing transition, the next
    // canonical vertex is already known.  Following that compact ID avoids a
    // full XXH3 plus bucket-chain probe for every repeated read occurrence.
    // The edge bit for this transition was installed on the previous loop
    // iteration, so a newly observed different base has already made the
    // degree non-unique and necessarily takes the exact lookup fallback.
    if (!force_exact_lookup &&
        previous_index != vertex_table_type::kNull) {
      HashGraphVertexAdaptor previous(
          &vertex_table_.value_at(previous_index), previous_reverse);
      const uint32_t cached = previous.next_neighbor();
      if (cached < HashGraphVertex::kUncacheableNeighbor) {
        // Preserve the historical rehash boundary even though this duplicate
        // occurrence no longer needs to probe the table.
        vertex_table_.prepare_cached_occurrence();
        vertex_index = cached >> 1u;
        is_reverse = (cached & 1u) != 0u;
        followed_cached_transition = true;
      } else if (use_branch_transition_cache &&
                 cached == HashGraphVertex::kUncacheableNeighbor) {
        const uint8_t transition_base =
            PackedBaseAt(forward_words, start + kmer_size_ - 1u);
        uint32_t branch_neighbor = HashGraphVertex::kNoCachedNeighbor;
        if (FindBranchTransition(previous_index, previous_reverse,
                                 transition_base, &branch_neighbor)) {
          vertex_table_.prepare_cached_occurrence();
          vertex_index = branch_neighbor >> 1u;
          is_reverse = (branch_neighbor & 1u) != 0u;
          followed_cached_transition = true;
          ++branch_transition_hits;
        } else {
          ++branch_transition_misses;
        }
      }
    }
    force_exact_lookup = false;

    if (vertex_index == vertex_table_type::kNull) {
      LoadPackedKmerPair<Words>(forward_words, start, kmer_size_, kmer,
                                rev_kmer);
      for (int word = static_cast<int>(Words) - 1; word >= 0; --word) {
        if (rev_kmer[word] != kmer[word]) {
          is_reverse = rev_kmer[word] < kmer[word];
          break;
        }
      }
      const uint64_t *key = is_reverse ? rev_kmer : kmer;
      vertex_table_.find_or_insert_packed<Words>(
          key, kmer_size_, packed_hash_fixed_accumulator_,
                                                  &vertex_index);
    }
    HashGraphVertex &vertex = vertex_table_.value_at(vertex_index);
    if (increment_count) vertex.count() += 1;
    HashGraphVertexAdaptor adaptor(&vertex, is_reverse);

    if (EmitGroups) {
      const uint32_t code = (vertex_index << 1u) | uint32_t(is_reverse);
      if (forward_groups != nullptr) forward_groups[start] = code;
      if (reverse_groups != nullptr) {
        reverse_groups[num_kmers - 1u - start] = code ^ 1u;
      }
    }

    if (start != 0 && !followed_cached_transition) {
      adaptor.in_edges().Add(3u - PackedBaseAt(forward_words, start - 1u));
    }
    if (start + kmer_size_ < sequence_length) {
      const uint8_t next_base =
          PackedBaseAt(forward_words, start + kmer_size_);
      if (!adaptor.out_edges()[next_base]) {
        const bool already_had_out_edge = !adaptor.out_edges().empty();
        adaptor.out_edges().Add(next_base);
        // A cached transition describes the formerly unique edge. Invalidate
        // it at the exact mutation that creates a second edge, instead of
        // leaving a stale-valid window until the next occurrence reaches
        // CacheTransition(). This makes cache validity itself a sufficient
        // hot-loop test while preserving the final graph state.
        if (already_had_out_edge) {
          adaptor.next_neighbor() = HashGraphVertex::kUncacheableNeighbor;
        }
      }
    }

    if (previous_index != vertex_table_type::kNull &&
        !followed_cached_transition) {
      HashGraphVertexAdaptor previous(
          &vertex_table_.value_at(previous_index), previous_reverse);
      CacheTransition(previous, previous_index, previous_reverse, adaptor,
                      vertex_index, is_reverse);
      if (use_branch_transition_cache &&
          previous.next_neighbor() ==
              HashGraphVertex::kUncacheableNeighbor) {
        RememberBranchTransition(previous_index, previous_reverse,
                                 adaptor.last_base(),
                                 EncodeEndpoint(vertex_index, is_reverse));
      }
      HashGraphVertexAdaptor reverse_current = adaptor;
      reverse_current.ReverseComplement();
      if (use_branch_transition_cache &&
          reverse_current.next_neighbor() ==
              HashGraphVertex::kUncacheableNeighbor) {
        HashGraphVertexAdaptor reverse_previous = previous;
        reverse_previous.ReverseComplement();
        RememberBranchTransition(
            vertex_index, !is_reverse, reverse_previous.last_base(),
            EncodeEndpoint(previous_index, !previous_reverse));
      }
    }
    previous_index = vertex_index;
    previous_reverse = is_reverse;

    // Once a read enters topology that earlier occurrences have already
    // resolved, replay the longest certified path in one compact loop.  A
    // cached transition proves both the exact next canonical vertex and the
    // corresponding oriented edge.  Counts are therefore the only semantic
    // state that still changes.  The first missing edge/cache immediately
    // returns to the historical lookup path, which installs any new graph
    // state in precisely the original occurrence order.
    if (!EmitGroups && use_cached_path_replay) {
      // The compact loop contains no insertion.  Apply a pending resize once
      // at the first following occurrence rather than retesting the same
      // false flag for every certified vertex in the run.
      if (start + 1u < num_kmers) {
        vertex_table_.prepare_cached_occurrence();
      }
      while (start + 1u < num_kmers) {
        HashGraphVertex &previous_vertex =
            vertex_table_.value_at(previous_index);
        HashGraphVertexAdaptor previous(&previous_vertex, previous_reverse);
        const uint8_t transition_base =
            PackedBaseAt(forward_words, start + kmer_size_);

        if (!previous.out_edges()[transition_base]) {
          const bool already_had_out_edge = !previous.out_edges().empty();
          previous.out_edges().Add(transition_base);
          if (already_had_out_edge) {
            previous.next_neighbor() =
                HashGraphVertex::kUncacheableNeighbor;
          }
          break;
        }

        uint32_t next_code = previous.next_neighbor();
        if (next_code == HashGraphVertex::kUncacheableNeighbor) {
          if (!use_branch_transition_cache ||
              !FindBranchTransition(previous_index, previous_reverse,
                                    transition_base, &next_code)) {
            ++branch_transition_misses;
            force_exact_lookup = true;
            break;
          }
          ++branch_transition_hits;
        } else if (next_code >= HashGraphVertex::kUncacheableNeighbor) {
          break;
        }

        previous_index = next_code >> 1u;
        previous_reverse = (next_code & 1u) != 0u;
        if (increment_count) {
          ++vertex_table_.value_at(previous_index).count();
        }
        ++start;
      }
    }
  }

  branch_transition_hits_ += branch_transition_hits;
  branch_transition_misses_ += branch_transition_misses;

  return num_kmers;
}

int64_t HashGraph::InsertPackedKmers(const uint64_t *forward_words,
                                     uint32_t sequence_length) {
  switch ((kmer_size_ + 31u) >> 5u) {
    case 1:
      return InsertPackedKmersImpl<1, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 2:
      return InsertPackedKmersImpl<2, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 3:
      return InsertPackedKmersImpl<3, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 4:
      return InsertPackedKmersImpl<4, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 5:
      return InsertPackedKmersImpl<5, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 6:
      return InsertPackedKmersImpl<6, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 7:
      return InsertPackedKmersImpl<7, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    case 8:
      return InsertPackedKmersImpl<8, false>(
          forward_words, sequence_length, true, nullptr, nullptr);
    default:
      throw std::length_error("unsupported local packed k-mer width");
  }
}

int64_t HashGraph::InsertPackedKmersIndexed(
    const uint64_t *forward_words, uint32_t sequence_length,
    uint32_t *forward_groups, uint32_t *reverse_groups) {
  switch ((kmer_size_ + 31u) >> 5u) {
    case 1:
      return InsertPackedKmersImpl<1, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 2:
      return InsertPackedKmersImpl<2, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 3:
      return InsertPackedKmersImpl<3, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 4:
      return InsertPackedKmersImpl<4, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 5:
      return InsertPackedKmersImpl<5, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 6:
      return InsertPackedKmersImpl<6, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 7:
      return InsertPackedKmersImpl<7, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    case 8:
      return InsertPackedKmersImpl<8, true>(
          forward_words, sequence_length, true, forward_groups,
          reverse_groups);
    default:
      throw std::length_error("unsupported local packed k-mer width");
  }
}

int64_t HashGraph::InsertPackedUncountKmers(const uint64_t *forward_words,
                                            uint32_t sequence_length) {
  switch ((kmer_size_ + 31u) >> 5u) {
    case 1:
      return InsertPackedKmersImpl<1, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 2:
      return InsertPackedKmersImpl<2, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 3:
      return InsertPackedKmersImpl<3, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 4:
      return InsertPackedKmersImpl<4, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 5:
      return InsertPackedKmersImpl<5, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 6:
      return InsertPackedKmersImpl<6, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 7:
      return InsertPackedKmersImpl<7, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    case 8:
      return InsertPackedKmersImpl<8, false>(
          forward_words, sequence_length, false, nullptr, nullptr);
    default:
      throw std::length_error("unsupported local packed k-mer width");
  }
}

int64_t HashGraph::InsertUncountKmers(const Sequence &seq) {
  if (seq.size() < kmer_size_) return 0;

  IdbaKmer kmer(kmer_size_), rev_kmer(kmer_size_);
  int length = 0;
  int64_t num_kmers = 0;
  uint32_t previous_index = vertex_table_type::kNull;
  bool previous_reverse = false;
  for (uint64_t i = 0; i < seq.size(); ++i) {
    const uint8_t base = seq[i];
    kmer.ShiftAppend(base);
    rev_kmer.ShiftPreappend(3u - (base & 3u));
    length = (base < 4) ? length + 1 : 0;

    if (length < (int)kmer_size_) {
      previous_index = vertex_table_type::kNull;
      continue;
    }

    const bool is_reverse = rev_kmer < kmer;
    const IdbaKmer &key = is_reverse ? rev_kmer : kmer;

    uint32_t vertex_index;
    HashGraphVertex &vertex =
        vertex_table_.find_or_insert_key(key, &vertex_index);
    HashGraphVertexAdaptor adaptor(&vertex, is_reverse);

    if (length > (int)kmer_size_ && seq[i - kmer_size_] < 4)
      adaptor.in_edges().Add(3 - seq[i - kmer_size_]);
    if (i + 1 < seq.size() && seq[i + 1] < 4)
      adaptor.out_edges().Add(seq[i + 1]);

    if (previous_index != vertex_table_type::kNull) {
      HashGraphVertexAdaptor previous(
          &vertex_table_.value_at(previous_index), previous_reverse);
      CacheTransition(previous, previous_index, previous_reverse, adaptor,
                      vertex_index, is_reverse);
    }
    previous_index = vertex_index;
    previous_reverse = is_reverse;

    ++num_kmers;
  }

  return num_kmers;
}

int64_t HashGraph::Assemble(std::vector<ContigGraphVertex> &unitigs) {
  unitigs.clear();
  unitigs.reserve(vertex_table_.size());
  assembled_endpoint_codes_.clear();
  AssembleFunc func(this, &unitigs, &assembled_endpoint_codes_);
  vertex_table_.for_each(func);
  // The local assembler consumes these contigs and then discards the graph.
  // Clearing every status here is therefore an unobservable full-table pass.
  return unitigs.size();
}

bool HashGraph::BuildUnitigAdjacency(
    std::vector<ContigGraphVertex> &unitigs,
    std::vector<uint32_t> &neighbor_codes, uint64_t *num_edges) {
  // Odd k is the normal local-assembly path and cannot contain a DNA
  // reverse-complement palindrome. Keep the historical endpoint-map path for
  // even k, where one physical endpoint may represent both orientations.
  if ((kmer_size_ & 1u) == 0u ||
      unitigs.size() >= (uint64_t{1} << 30u) ||
      vertex_table_.size() >= (uint64_t{1} << 30u) ||
      assembled_endpoint_codes_.size() != unitigs.size() * 4u) {
    return false;
  }

  constexpr uint32_t kUnitigMarker = uint32_t{1} << 31u;
  const size_t required = unitigs.size() * 8u;
  if (neighbor_codes.size() < required) neighbor_codes.resize(required);
  uint64_t total_degree = 0;

  // First resolve outgoing HashGraph neighbors while the original transition
  // cache is intact. Simple terminals use their direct cached index; only a
  // branching/ambiguous terminal falls back to a k-mer table lookup.
  for (uint32_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      const uint32_t terminal_code =
          assembled_endpoint_codes_[size_t(unitig_id) * 4u + 2u + strand];
      HashGraphVertexAdaptor terminal(
          &vertex_table_.value_at(terminal_code >> 1u),
          (terminal_code & 1u) != 0u);
      const uint32_t cached = terminal.next_neighbor();
      const bool use_cached = terminal.out_edges().size() == 1 &&
                              cached < HashGraphVertex::kUncacheableNeighbor;
      uint8_t edges = static_cast<uint8_t>(terminal.out_edges());
      while (edges != 0u) {
        const uint32_t base = static_cast<uint32_t>(__builtin_ctz(edges));
        edges &= static_cast<uint8_t>(edges - 1u);
        uint32_t next_hash_code;
        if (use_cached) {
          next_hash_code = cached;
        } else {
          IdbaKmer next_kmer = terminal.kmer();
          next_kmer.ShiftAppend(base);
          HashGraphVertexAdaptor next = FindVertexAdaptor(next_kmer);
          if (next.is_null()) return false;
          next_hash_code =
              (vertex_table_.index_of(next.vertex()) << 1u) |
              uint32_t(next.is_reverse());
        }
        neighbor_codes[size_t(unitig_id) * 8u + size_t(strand) * 4u + base] =
            next_hash_code;
      }
      total_degree += terminal.out_edges().size();
    }
  }

  // Unitigs no longer consult HashGraph after this point, so reuse the two
  // transition slots as a compact endpoint-to-unitig map.
  for (uint32_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      const uint32_t endpoint_code =
          assembled_endpoint_codes_[size_t(unitig_id) * 4u + strand];
      const uint32_t code = (unitig_id << 1u) | strand;
      HashGraphVertexAdaptor endpoint(
          &vertex_table_.value_at(endpoint_code >> 1u),
          (endpoint_code & 1u) != 0u);
      uint32_t &slot = endpoint.next_neighbor();
      if (slot < HashGraphVertex::kUncacheableNeighbor &&
          (slot & kUnitigMarker) != 0u &&
          (slot & ~kUnitigMarker) != code) {
        return false;
      }
      slot = kUnitigMarker | code;
    }
  }

  // Translate the temporary oriented HashGraph codes to oriented unitig IDs.
  for (uint32_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
    for (uint32_t strand = 0; strand < 2u; ++strand) {
      ContigGraphVertexAdaptor current(&unitigs[unitig_id], strand != 0);
      uint32_t sole_neighbor = ContigGraphVertex::kNoCachedNeighbor;
      uint8_t edges = static_cast<uint8_t>(current.out_edges());
      while (edges != 0u) {
        const uint32_t base = static_cast<uint32_t>(__builtin_ctz(edges));
        edges &= static_cast<uint8_t>(edges - 1u);
        const size_t slot = size_t(unitig_id) * 8u +
                            size_t(strand) * 4u + base;
        const uint32_t next_hash_code = neighbor_codes[slot];
        uint32_t next_code;
        HashGraphVertexAdaptor next(
            &vertex_table_.value_at(next_hash_code >> 1u),
            (next_hash_code & 1u) != 0u);
        const uint32_t marked = next.next_neighbor();
        if (marked >= HashGraphVertex::kUncacheableNeighbor ||
            (marked & kUnitigMarker) == 0u) {
          return false;
        }
        next_code = marked & ~kUnitigMarker;
        neighbor_codes[slot] = next_code;
        sole_neighbor = next_code;
      }
      current.next_neighbor() =
          current.out_edges().size() == 1
              ? sole_neighbor
              : ContigGraphVertex::kUncacheableNeighbor;
    }
  }

  *num_edges = total_degree / 2u;
  return true;
}

void HashGraph::AssembleFunc::operator()(HashGraphVertex &vertex) {
  const bool odd_k = (hash_graph_->kmer_size() & 1u) != 0u;
  if (odd_k) {
    if (!vertex.status().Lock(0)) return;
    std::vector<uint8_t> &forward_arm =
        hash_graph_->assemble_arm_bases_[0];
    std::vector<uint8_t> &reverse_arm =
        hash_graph_->assemble_arm_bases_[1];
    forward_arm.clear();
    reverse_arm.clear();
    uint32_t unitig_kmer_count = static_cast<uint32_t>(vertex.count());
    HashGraphVertexAdaptor terminals[2] = {
        HashGraphVertexAdaptor(&vertex, false),
        HashGraphVertexAdaptor(&vertex, true)};
    HashGraphVertexAdaptor loop_begin(&vertex);

    // This is the exact historical walk and lock order. Only the physical
    // sequence construction is deferred until both arms are known.
    for (int strand = 0; strand < 2; ++strand) {
      HashGraphVertexAdaptor current(&vertex, strand != 0);
      while (true) {
        HashGraphVertexAdaptor next;
        if (!hash_graph_->GetNextVertexAdaptor(current, next)) break;
        if (hash_graph_->IsPalindromeLoop(current, next)) break;
        if (next == loop_begin) return;
        if (!next.status().LockPreempt(0)) return;

        hash_graph_->assemble_arm_bases_[strand].push_back(next.last_base());
        unitig_kmer_count += static_cast<uint32_t>(next.count());
        current = next;
      }
      terminals[strand] = current;
      loop_begin = current;
      loop_begin.ReverseComplement();
    }

    ContigInfo contig_info;
    const uint32_t kmer_size = hash_graph_->kmer_size();
    // Position-wise counts were already proven unused by local cleaning and
    // output. Preserve exactly the aggregate/count and oriented boundary
    // edges that ContigBuilder would have produced.
    contig_info.in_edges() = terminals[1].out_edges();
    contig_info.out_edges() = terminals[0].out_edges();
    contig_info.set_kmer_size(kmer_size);
    contig_info.set_kmer_count(unitig_kmer_count);
    unitigs_->emplace_back();
    unitigs_->back().take_contig_info(contig_info);

    HashGraphVertexAdaptor begin_forward = terminals[1];
    HashGraphVertexAdaptor begin_reverse = terminals[0];
    begin_forward.ReverseComplement();
    begin_reverse.ReverseComplement();

    Sequence contig;
    contig.resize(kmer_size + forward_arm.size() + reverse_arm.size());
    uint32_t output = 0;
    // Historical construction finishes by reverse-complementing the
    // strand-1 walk, hence its bases appear reversed and complemented.
    for (auto it = reverse_arm.rbegin(); it != reverse_arm.rend(); ++it) {
      contig[output++] = 3u - *it;
    }
    const IdbaKmer seed = vertex.kmer();
    for (uint32_t i = 0; i < kmer_size; ++i) {
      contig[output++] = seed[i];
    }
    for (uint8_t base : forward_arm) contig[output++] = base;
    unitigs_->back().take_contig(contig);

    const auto encode = [this](const HashGraphVertexAdaptor &endpoint) {
      return (hash_graph_->vertex_table_.index_of(endpoint.vertex()) << 1u) |
             uint32_t(endpoint.is_reverse());
    };
    endpoint_codes_->push_back(encode(begin_forward));
    endpoint_codes_->push_back(encode(begin_reverse));
    endpoint_codes_->push_back(encode(terminals[0]));
    endpoint_codes_->push_back(encode(terminals[1]));
    return;
  }

  if (!vertex.status().Lock(0)) return;

  ContigBuilder contig_builder(false);
  contig_builder.Append(HashGraphVertexAdaptor(&vertex));
  HashGraphVertexAdaptor terminals[2] = {
      HashGraphVertexAdaptor(&vertex, false),
      HashGraphVertexAdaptor(&vertex, true)};

  if ((hash_graph_->kmer_size() & 1u) != 0u ||
      !vertex.kmer().IsPalindrome()) {
    HashGraphVertexAdaptor loop_begin(&vertex);
    for (int strand = 0; strand < 2; ++strand) {
      HashGraphVertexAdaptor current(&vertex, strand);
      IdbaKmer begin_kmer;
      if ((hash_graph_->kmer_size() & 1u) == 0u) {
        begin_kmer =
            contig_builder.contig().GetIdbaKmer(0, hash_graph_->kmer_size());
      }

      while (true) {
        HashGraphVertexAdaptor next;
        if (!hash_graph_->GetNextVertexAdaptor(current, next)) break;

        if (hash_graph_->IsPalindromeLoop(current, next)) break;

        // For odd k an oriented k-mer has a unique (vertex,strand) identity;
        // no DNA k-mer can be its own reverse complement.  Comparing that
        // identity is exactly equivalent to constructing/reverse-complementing
        // a 72-byte IdbaKmer on every chain step.
        if (((hash_graph_->kmer_size() & 1u) != 0u)
                ? next == loop_begin
                : hash_graph_->IsLoop(begin_kmer, next))
          return;

        if (!next.status().LockPreempt(0)) return;

        contig_builder.Append(next);
        current = next;
      }

      terminals[strand] = current;

      loop_begin = current;
      loop_begin.ReverseComplement();
      contig_builder.ReverseComplement();
    }
  }

  Sequence contig;
  ContigInfo contig_info;
  contig_builder.Release(contig, contig_info);
  unitigs_->emplace_back();
  unitigs_->back().take_contig(contig);
  unitigs_->back().take_contig_info(contig_info);
  HashGraphVertexAdaptor begin_forward = terminals[1];
  HashGraphVertexAdaptor begin_reverse = terminals[0];
  begin_forward.ReverseComplement();
  begin_reverse.ReverseComplement();
  const auto encode = [this](const HashGraphVertexAdaptor &endpoint) {
    return (hash_graph_->vertex_table_.index_of(endpoint.vertex()) << 1u) |
           uint32_t(endpoint.is_reverse());
  };
  endpoint_codes_->push_back(encode(begin_forward));
  endpoint_codes_->push_back(encode(begin_reverse));
  endpoint_codes_->push_back(encode(terminals[0]));
  endpoint_codes_->push_back(encode(terminals[1]));
}
