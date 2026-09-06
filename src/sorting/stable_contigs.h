#ifndef MEGAHIT_SORTING_STABLE_CONTIGS_H_
#define MEGAHIT_SORTING_STABLE_CONTIGS_H_

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>
#include "sequence/sequence_package.h"

namespace stable_contigs {

struct DropSpan {
  uint32_t input, begin, end, owner;
  bool operator<(const DropSpan &b) const {
    return input != b.input ? input < b.input : begin < b.begin;
  }
};

struct Plan {
  std::vector<uint8_t> carry;
  std::vector<uint8_t> remove;
  std::vector<DropSpan> spans;
  uint64_t candidate_bases{0};
  uint64_t carried_bases{0};
  uint64_t removed_windows{0};
  uint64_t carried_contigs{0};
  uint64_t redundant_sequences{0};
};

inline uint64_t Hash(uint64_t key) {
  key ^= key >> 30; key *= UINT64_C(0xbf58476d1ce4e5b9);
  key ^= key >> 27; key *= UINT64_C(0x94d049bb133111eb);
  return key ^ (key >> 31);
}

struct Seed {
  uint64_t key;
  uint32_t owner;
  uint32_t position;
  bool operator<(const Seed &other) const { return key < other.key; }
};

template <class Visitor>
void VisitSeeds(const SeqPackage::SeqView &sequence, const Visitor &visit) {
  const uint64_t mask = (uint64_t{1} << 62u) - 1u;
  uint64_t forward = 0, reverse = 0;
  for (uint32_t i = 0; i < sequence.length(); ++i) {
    const uint64_t base = sequence.base_at(i);
    forward = ((forward << 2u) | base) & mask;
    reverse = (reverse >> 2u) | ((base ^ 3u) << 60u);
    if (i >= 30u) visit(std::min(forward, reverse), i - 30u, reverse < forward);
  }
}

inline bool Contains(const SeqPackage::SeqView &outer,
                     const SeqPackage::SeqView &inner, const Seed &seed,
                     uint32_t query_position, bool query_reverse) {
  if (inner.length() > outer.length()) return false;
  const bool reverse = (seed.position & 1u) != query_reverse;
  const int64_t start = int64_t(seed.position >> 1u) -
      (reverse ? int64_t(inner.length() - 31u - query_position)
               : int64_t(query_position));
  if (start < 0 || uint64_t(start) + inner.length() > outer.length()) return false;
  for (uint32_t i = 0; i < inner.length(); ++i) {
    const unsigned base = reverse ? inner.base_at(inner.length() - 1u - i) ^ 3u
                                  : inner.base_at(i);
    if (outer.base_at(static_cast<uint32_t>(start) + i) != base) return false;
  }
  return true;
}

inline bool SharesNode(const SeqPackage::SeqView &outer,
                       const SeqPackage::SeqView &inner, const Seed &seed,
                       uint32_t query_position, bool query_reverse, unsigned k) {
  const bool reverse = (seed.position & 1u) != query_reverse;
  const uint32_t a = seed.position >> 1u;
  const uint32_t b = reverse ? inner.length() - 31u - query_position : query_position;
  const auto base = [&](uint32_t i) {
    return reverse ? inner.base_at(inner.length() - 1u - i) ^ 3u : inner.base_at(i);
  };
  unsigned matched = 31u;
  uint32_t left = 0;
  while (matched < k && left < a && left < b &&
         outer.base_at(a - left - 1u) == base(b - left - 1u)) {
    ++left; ++matched;
  }
  uint32_t right = 31u;
  while (matched < k && right < outer.length() - a && right < inner.length() - b &&
         outer.base_at(a + right) == base(b + right)) {
    ++right; ++matched;
  }
  return matched >= k;
}

struct Match {
  uint32_t owner{std::numeric_limits<uint32_t>::max()};
  uint32_t outer_begin{0}, outer_end{0}, inner_begin{0}, inner_end{0};
  bool reverse{false}, safe{false};
};

inline Match ExtendMatch(const SeqPackage::SeqView &outer,
                         const SeqPackage::SeqView &inner, const Seed &seed,
                         uint32_t position, bool reverse) {
  Match match;
  match.owner = seed.owner;
  match.reverse = (seed.position & 1u) != reverse;
  match.outer_begin = seed.position >> 1u;
  match.inner_begin = match.reverse ? inner.length() - 31u - position : position;
  match.outer_end = match.outer_begin + 31u;
  match.inner_end = match.inner_begin + 31u;
  const auto base = [&](uint32_t i) {
    return match.reverse ? inner.base_at(inner.length() - 1u - i) ^ 3u : inner.base_at(i);
  };
  while (match.outer_begin && match.inner_begin &&
         outer.base_at(match.outer_begin - 1u) == base(match.inner_begin - 1u)) {
    --match.outer_begin; --match.inner_begin;
  }
  while (match.outer_end < outer.length() && match.inner_end < inner.length() &&
         outer.base_at(match.outer_end) == base(match.inner_end)) {
    ++match.outer_end; ++match.inner_end;
  }
  return match;
}

// An overlap of k bases contains k-30 possible 31-mer starts. Sampling one
// start every k-30 bases in a candidate therefore cannot miss any graph
// connection. Bloom filters only reject absent seeds; every possible match
// gets a full 62-bit lookup. Uncertain matches disqualify the candidate.
// A complete contained input with no higher multiplicity can be removed
// together with its owner, since seq2sdbg reduces duplicates by maximum.
template <class Multiplicity>
Plan Build(const SeqPackage &sequences, const std::vector<Multiplicity> &multi,
           uint64_t primary_begin, const std::vector<unsigned> &flags,
           unsigned k, unsigned num_threads, bool interior_only = false,
           bool exact_overlap = false, unsigned end_margin = 0) {
  Plan result;
  const uint32_t none = std::numeric_limits<uint32_t>::max();
  const size_t size = sequences.seq_count();
  const unsigned margin = end_margin ? end_margin : 2u * k;
  if (flags.empty() || k < 31u || size >= none || multi.size() != size || primary_begin > size ||
      flags.size() > size - primary_begin) return result;
  result.carry.assign(size, 0);
  result.remove.assign(size, 0);
  std::vector<uint8_t> candidate(size, 0);
  std::unique_ptr<std::atomic<uint8_t>[]> affected(new std::atomic<uint8_t>[size]);
  for (size_t i = 0; i < size; ++i) affected[i].store(0, std::memory_order_relaxed);
  for (size_t i = 0; i < flags.size(); ++i) {
    const size_t id = primary_begin + i;
    const auto sequence = sequences.GetSeqView(id);
    const bool eligible_flag = interior_only ? !(flags[i] & 2u) : flags[i] == 1u;
    const uint64_t minimum_length = interior_only ? uint64_t{2} * margin + 4u * k : 3u * k;
    if (eligible_flag && uint64_t(sequence.length()) >= minimum_length &&
        sequence.length() <= (none >> 1u) && multi[id] != 0) {
      candidate[id] = 1;
      result.candidate_bases += sequence.length();
    }
  }
  if (result.candidate_bases == 0) return result;
  const unsigned workers = std::max(1u, num_threads);
  const unsigned stride = k - 30u;
  constexpr unsigned shards = 256;
  std::vector<std::vector<Seed>> boxes(static_cast<size_t>(workers) * shards);
#pragma omp parallel for schedule(dynamic, 64) num_threads(workers)
  for (int64_t i = 0; i < static_cast<int64_t>(flags.size()); ++i) {
    const uint32_t id = static_cast<uint32_t>(primary_begin + i);
    if (!candidate[id]) continue;
    const auto sequence = sequences.GetSeqView(id);
    VisitSeeds(sequence, [&](uint64_t key, uint32_t position, bool reverse) {
      if (interior_only &&
          (position < margin || uint64_t(position) + 31u + margin > sequence.length())) return;
      if (position % stride == 0u) {
        boxes[static_cast<size_t>(omp_get_thread_num()) * shards + (Hash(key) & 255u)]
            .push_back(Seed{key, id, (position << 1u) | uint32_t(reverse)});
      }
    });
  }
  std::array<std::vector<Seed>, shards> index;
#pragma omp parallel for num_threads(workers)
  for (unsigned shard = 0; shard < shards; ++shard) {
    size_t count = 0;
    for (unsigned worker = 0; worker < workers; ++worker) count += boxes[worker * shards + shard].size();
    auto &bucket = index[shard];
    bucket.reserve(count);
    for (unsigned worker = 0; worker < workers; ++worker) {
      auto &box = boxes[worker * shards + shard];
      bucket.insert(bucket.end(), box.begin(), box.end());
      std::vector<Seed>().swap(box);
    }
    std::sort(bucket.begin(), bucket.end());
    for (size_t i = 0; i < bucket.size();) {
      size_t end = i + 1u;
      while (end < bucket.size() && bucket[end].key == bucket[i].key) ++end;
      if (end - i > (exact_overlap ? 8u : 1u)) {
        for (size_t j = i; j < end; ++j)
          affected[bucket[j].owner].store(1, std::memory_order_relaxed);
      }
      i = end;
    }
  }
  uint64_t seed_count = 0;
  for (const auto &bucket : index) seed_count += bucket.size();
  size_t filter_words = 1;
  while (filter_words < (seed_count + 7u) / 8u) filter_words <<= 1u;
  std::vector<uint64_t> filter(filter_words, 0);
  const auto bits = [](uint64_t hash) {
    return (uint64_t{1} << ((hash >> 32u) & 63u)) |
           (uint64_t{1} << ((hash >> 40u) & 63u)) |
           (uint64_t{1} << ((hash >> 48u) & 63u)) |
           (uint64_t{1} << ((hash >> 56u) & 63u));
  };
#pragma omp parallel for num_threads(workers)
  for (unsigned shard = 0; shard < shards; ++shard) {
    for (const Seed &seed : index[shard]) {
      const uint64_t hash = Hash(seed.key);
      __atomic_fetch_or(&filter[hash & (filter_words - 1u)], bits(hash), __ATOMIC_RELAXED);
    }
  }
  std::vector<uint32_t> contained(size, none);
  std::vector<std::vector<DropSpan>> thread_spans(workers);
#pragma omp parallel for schedule(dynamic, 64) num_threads(workers)
  for (int64_t i = 0; i < static_cast<int64_t>(size); ++i) {
    const auto sequence = sequences.GetSeqView(i);
    if (sequence.length() <= k) continue;
    bool attempted_containment = false;
    uint32_t owner = none;
    Match last_match;
    VisitSeeds(sequence, [&](uint64_t key, uint32_t position, bool reverse) {
      const uint64_t hash = Hash(key), mask = bits(hash);
      if ((filter[hash & (filter_words - 1u)] & mask) != mask) return;
      const auto &bucket = index[hash & 255u];
      auto found = std::lower_bound(bucket.begin(), bucket.end(), Seed{key, 0, 0});
      if (found == bucket.end() || found->key != key) return;
      if (exact_overlap && bucket.end() - found > 8 && found[8].key == key) return;
      do {
        if (found->owner == static_cast<uint32_t>(i) &&
            found->position == ((position << 1u) | uint32_t(reverse))) continue;
        if (exact_overlap && !SharesNode(sequences.GetSeqView(found->owner),
                                         sequence, *found, position, reverse, k)) continue;
        if (found->owner == static_cast<uint32_t>(i)) {
          affected[i].store(1, std::memory_order_relaxed);
          continue;
        }
        // An extra input may continue through either retained endpoint.
        // Remove its exact shared edge interval along with the parent span,
        // provided it introduces neither an internal branch nor higher depth.
        if (interior_only && exact_overlap && !candidate[i]) {
          const auto outer = sequences.GetSeqView(found->owner);
          const bool query_rc = (found->position & 1u) != reverse;
          const uint32_t q = query_rc ? sequence.length() - 31u - position : position;
          const uint32_t p = found->position >> 1u;
          if (last_match.owner != found->owner || last_match.reverse != query_rc ||
              q < last_match.inner_begin || q + 31u > last_match.inner_end ||
              int64_t(p) - q != int64_t(last_match.outer_begin) - last_match.inner_begin) {
            last_match = ExtendMatch(outer, sequence, *found, position, reverse);
            const uint32_t begin = margin + 1u;
            const uint32_t end = outer.length() - margin - k - 1u;
            const uint32_t right_node = last_match.outer_end - k;
            const bool left_branch = last_match.inner_begin && last_match.outer_begin &&
                last_match.outer_begin >= begin && last_match.outer_begin <= end;
            const bool right_branch = last_match.inner_end < sequence.length() &&
                last_match.outer_end < outer.length() && right_node >= begin && right_node <= end;
            const bool covered_interior = last_match.outer_begin < end && right_node > begin;
            last_match.safe = !left_branch && !right_branch &&
                (!covered_interior || multi[i] <= multi[found->owner]);
            if (last_match.safe && covered_interior) {
              uint32_t from = last_match.inner_begin +
                  std::max(last_match.outer_begin, begin) - last_match.outer_begin;
              uint32_t to = last_match.inner_begin +
                  std::min(right_node, end) - last_match.outer_begin;
              if (query_rc) {
                const uint32_t reversed_from = sequence.length() - k - to;
                to = sequence.length() - k - from;
                from = reversed_from;
              }
              thread_spans[omp_get_thread_num()].push_back(
                  {static_cast<uint32_t>(i), from, to, found->owner});
            }
          }
          if (!last_match.safe) affected[found->owner].store(1, std::memory_order_relaxed);
          continue;
        }
        if (!candidate[i] && !attempted_containment) {
          attempted_containment = true;
          if (multi[i] <= multi[found->owner] &&
              Contains(sequences.GetSeqView(found->owner), sequence, *found, position, reverse)) {
            owner = found->owner;
          }
        }
        if (found->owner != owner) affected[found->owner].store(1, std::memory_order_relaxed);
      } while (exact_overlap && ++found != bucket.end() && found->key == key);
    });
    contained[i] = owner;
  }
  for (size_t i = 0; i < size; ++i) {
    if (candidate[i] && !affected[i].load(std::memory_order_relaxed)) {
      result.carry[i] = result.remove[i] = 1;
      ++result.carried_contigs;
      result.carried_bases += sequences.GetSeqView(i).length();
      if (interior_only) result.spans.push_back(
          {static_cast<uint32_t>(i), margin + 1u,
           sequences.GetSeqView(i).length() - margin - k - 1u,
           static_cast<uint32_t>(i)});
    }
  }
  for (size_t i = 0; i < size; ++i) {
    if (contained[i] != none && result.carry[contained[i]]) {
      result.remove[i] = 1;
      ++result.redundant_sequences;
    }
    if (result.remove[i]) {
      result.removed_windows += sequences.GetSeqView(i).length() - k;
      if (interior_only && result.carry[i]) result.removed_windows -= 2u * (margin + 1u);
    }
  }
  if (interior_only) {
    for (const auto &spans : thread_spans) {
      for (const auto &span : spans) {
        if (result.carry[span.owner]) result.spans.push_back(span);
      }
    }
    std::sort(result.spans.begin(), result.spans.end());
    size_t output = 0;
    for (const auto &span : result.spans) {
      if (output && result.spans[output - 1].input == span.input &&
          result.spans[output - 1].end >= span.begin) {
        result.spans[output - 1].end = std::max(result.spans[output - 1].end, span.end);
      } else {
        result.spans[output++] = span;
      }
    }
    result.spans.resize(output);
  }
  return result;
}

// SeqToSdbg stores reversed (not complemented) input. OutputContigs chooses
// the strand whose first edge has the smaller BOSS ID, i.e. the smaller
// colexicographic k-base source. Certified paths have distinct source nodes.
inline std::string OutputSequence(const SeqPackage::SeqView &sequence, unsigned k) {
  const unsigned length = sequence.length();
  bool reverse_complement = false;
  for (unsigned i = 0; i < k; ++i) {
    const unsigned forward = sequence.base_at(length - k + i);
    const unsigned reverse = sequence.base_at(k - 1u - i) ^ 3u;
    if (forward != reverse) {
      reverse_complement = forward > reverse;
      break;
    }
  }
  std::string output(length, 'A');
  for (unsigned i = 0; i < length; ++i) {
    output[i] = "ACGT"[reverse_complement ? sequence.base_at(i) ^ 3u
                                         : sequence.base_at(length - 1u - i)];
  }
  return output;
}

}  // namespace stable_contigs
#endif
