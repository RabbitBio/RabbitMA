#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "iterate/contig_flank_index.h"

namespace {
void Require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

template <class KmerType>
uint64_t Prefix(const KmerType &kmer, unsigned k) {
  uint64_t result = 0;
  for (unsigned i = 0; i < std::min(19u, k); ++i) {
    result = (result << 2u) | kmer.GetBase(i);
  }
  return result;
}

template <class KmerType>
void Check(unsigned k) {
  using Index = ContigFlankIndex<KmerType>;
  using Bloom = typename Index::BlockedBloomFilter;
  using PrefixFilter = typename Index::OrientedPrefixFilter;
  std::vector<typename Index::HashShard> shards(64);
  std::vector<KmerType> queries;
  std::mt19937_64 random(987 + k);
  const size_t inserted = 32768;
  for (size_t i = 0; i < 2 * inserted; ++i) {
    KmerType key;
    for (unsigned pos = 0; pos < k; ++pos) {
      // Shared prefixes deliberately make different shards update the same
      // Bloom words; fully random keys also cover less contended words.
      const unsigned base = i % 4 == 0 && pos < 19 ? 0u : random() & 3u;
      key.ShiftAppend(base, k);
    }
    KmerType reverse = key;
    reverse.ReverseComplement(k);
    if (reverse < key) key = reverse;
    if (i < inserted) {
      shards[i % shards.size()].primary.emplace(
          key, typename Index::StoredFlankInfo{});
    }
    queries.push_back(key);
  }

  Bloom reference;
  PrefixFilter prefix_reference;
  omp_set_num_threads(1);
  reference.Build(shards, inserted);
  prefix_reference.Build(shards, inserted, k);
  for (int threads : {2, 8, 32}) {
    omp_set_num_threads(threads);
    Bloom actual;
    PrefixFilter prefix_actual;
    for (unsigned repeat = 0; repeat < 3; ++repeat) {
      actual.Build(shards, inserted);
      prefix_actual.Build(shards, inserted, k);
      Require(actual.byte_size() == reference.byte_size(), "Bloom size changed");
      Require(prefix_actual.byte_size() == prefix_reference.byte_size(),
              "prefix filter size changed");
      for (size_t i = 0; i < queries.size(); ++i) {
        const auto &key = queries[i];
        const size_t hash = Index::FilterHash(key);
        Require(actual.MayContain(hash) == reference.MayContain(hash),
                "parallel membership differs from serial");
        if (i < inserted) Require(actual.MayContain(hash), "Bloom false negative");
        KmerType oriented = key;
        for (int strand = 0; strand < 2; ++strand) {
          const uint64_t prefix = Prefix(oriented, k);
          Require(prefix_actual.MayContain(prefix) ==
                      prefix_reference.MayContain(prefix),
                  "parallel oriented membership differs from serial");
          if (i < inserted) {
            Require(prefix_actual.MayContain(prefix), "prefix false negative");
          }
          oriented.ReverseComplement(k);
        }
      }
      // Rebuilding an empty filter must discard the preceding graph's bits.
      actual.Build({}, 0);
      prefix_actual.Build({}, 0, k);
      Require(actual.byte_size() == 0 && !actual.MayContain(0),
              "empty Bloom retained bits");
      Require(prefix_actual.byte_size() == 0 && prefix_actual.MayContain(0),
              "empty prefix filter must allow exact lookup");
    }
  }
}
}  // namespace

int main() {
  unsetenv("MEGAHIT_DISABLE_PREFIX_FLANK_FILTER");
  omp_set_dynamic(0);
  Check<Kmer<1, uint64_t>>(18);
  Check<Kmer<3, uint32_t>>(40);
  Check<Kmer<2, uint64_t>>(64);
  Check<Kmer<4, uint64_t>>(120);
  Check<Kmer<5, uint64_t>>(140);
  std::cout << "Parallel flank filters preserve serial membership across "
               "key widths, orientations, collisions and rebuilds\n";
}
