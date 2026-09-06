#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include "localasm/hash_mapper.h"
#include "sequence/io/read_anchor_positions.h"

static std::vector<uint32_t> Pack(const std::string &dna) {
  std::vector<uint32_t> words((dna.size() + 15u) / 16u, 0u);
  for (size_t i = 0; i < dna.size(); ++i) {
    unsigned base = std::string("ACGT").find(dna[i]);
    words[i / 16u] |= base << ((15u - i % 16u) * 2u);
  }
  return words;
}
static std::string Reverse(std::string s) {
  std::reverse(s.begin(), s.end());
  for (char &c : s) c = std::string("TGCA")[std::string("ACGT").find(c)];
  return s;
}
int main(int argc, char **argv) {
  if (argc != 2) return 2;
  omp_set_num_threads(4);
  std::mt19937 rng(1931);
  auto dna = [&](unsigned length) {
    std::string s(length, 'A');
    for (char &c : s) c = "ACGT"[rng() % 4u];
    return s;
  };
  std::vector<std::string> refs;
  for (unsigned i = 0; i < 30u; ++i) refs.push_back(dna(2048));
  refs.push_back(std::string(90, 'A') + dna(1000) + std::string(90, 'T'));
  refs.push_back(refs[0]);
  refs.push_back(Reverse(refs[1]));
  refs.push_back(refs[2].substr(0, 300) + dna(900) + refs[3].substr(0, 300));
  refs.push_back(dna(200) + std::string(28, 'A') + dna(300));
  std::ofstream out(argv[1]);
  size_t total = 0;
  for (size_t i = 0; i < refs.size(); ++i) {
    out << ">k39_" << i << " flag=0 multi=5 len=" << refs[i].size()
        << "\n" << refs[i] << "\n";
    total += refs[i].size();
  }
  out.close();
  std::ofstream info(std::string(argv[1]) + ".info");
  info << refs.size() << ' ' << total << '\n';
  info.close();
  size_t checks = 0;
  for (unsigned sparsity : {1u, 3u, 8u, 16u}) {
    HashMapper mapper;
    mapper.LoadAndBuild(argv[1], 0, 31, sparsity);
    mapper.SetMappingThreshold(50, 0.95);
    mapper.BuildEndpointSeedFilter(348);
    if (!mapper.BuildEndpointMinimizerGate()) return 3;
    for (unsigned trial = 0; trial < 50000u; ++trial) {
      const auto &ref = refs[rng() % refs.size()];
      unsigned length = 20u + rng() % 493u;
      unsigned start = rng() % (ref.size() - length + 1u);
      std::string s = ref.substr(start, length);
      if (trial % 4u == 0u) s = dna(length);
      if (trial % 2u) s = Reverse(s);
      for (unsigned error = 0; error < trial % 7u; ++error) s[rng() % length] = "ACGT"[rng() % 4u];
      auto words = Pack(s);
      uint8_t bitmap[64] = {};
      std::vector<uint32_t> pos;
      std::vector<uint64_t> hash, key;
      ForEachPackedReadAnchor(words.data(), length, 19, 31, &pos, &hash, &key,
          [&](uint64_t, uint32_t p) { bitmap[p / 8u] |= 1u << (p % 8u); });
      HashMapper::EndpointSeedWitness a, b;
      bool x = mapper.MayMapToEndpoint(words.data(), length, &a);
      bool y = mapper.MayMapToEndpointAnchors(words.data(), length, bitmap, &b);
      if (x != y || (x && (a.index_value != b.index_value || a.end_position != b.end_position || a.query_strand != b.query_strand))) {
        std::fprintf(stderr, "Mismatch sparsity %u trial %u length %u gate %d/%d end %d/%d\n%s\n", sparsity, trial, length, x, y, a.end_position, b.end_position, s.c_str());
        return 1;
      }
      ++checks;
    }
  }
  std::printf("Exact earliest endpoint witness: %zu checks passed\n", checks);
  return 0;
}
