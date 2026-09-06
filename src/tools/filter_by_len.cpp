/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "sequence/io/contig/packed_contig.h"
#include "sequence/io/kseq.h"
#include "utils/histgram.h"
#include "utils/utils.h"

#ifndef KSEQ_INITED
#define KSEQ_INITED
// Keep the zlib-specific parser types distinct from the mgz parser in
// fastx_reader.h when the executable is optimized across translation units.
namespace {
KSEQ_INIT(gzFile, gzread)
}
#endif

namespace {

void FilterFastx(gzFile input, unsigned min_len, Histgram<long long> *hist) {
  kseq_t *seq = kseq_init(input);
  while (kseq_read(seq) >= 0) {
    if (seq->seq.l >= min_len) {
      hist->insert(seq->seq.l);
      pprintf(">{s} {s}\n{s}\n", seq->name.s,
              seq->comment.s == nullptr ? "" : seq->comment.s, seq->seq.s);
    }
  }
  kseq_destroy(seq);
}

void FilterFastxFile(const char *path, unsigned min_len,
                     Histgram<long long> *hist) {
  gzFile input = gzopen(path, "r");
  if (input == nullptr) {
    xfatal("Cannot open contig file {s}\n", path);
  }
  FilterFastx(input, min_len, hist);
  gzclose(input);
}

bool FilterPackedFile(const std::string &fasta_path, unsigned min_len,
                      Histgram<long long> *hist) {
  const std::string packed_path = packed_contig::Path(fasta_path);
  std::ifstream input(packed_path, std::ios::binary | std::ios::in);
  if (!input.is_open()) {
    return false;
  }

  packed_contig::FileHeader file_header{};
  input.read(reinterpret_cast<char *>(&file_header), sizeof(file_header));
  if (!input || !packed_contig::IsValid(file_header)) {
    xfatal("Invalid packed contig file: {s}\n", packed_path.c_str());
  }

  std::vector<uint32_t> packed_words;
  std::string sequence;
  static const char bases[] = {'A', 'C', 'G', 'T'};
  while (true) {
    packed_contig::RecordHeader header{};
    input.read(reinterpret_cast<char *>(&header), sizeof(header));
    if (input.gcount() == 0 && input.eof()) {
      break;
    }
    if (!input) {
      xfatal("Truncated packed contig header: {s}\n", packed_path.c_str());
    }

    const size_t words = (static_cast<size_t>(header.length) + 15u) / 16u;
    packed_words.resize(words);
    input.read(reinterpret_cast<char *>(packed_words.data()),
               words * sizeof(uint32_t));
    if (!input) {
      xfatal("Truncated packed contig sequence: {s}\n", packed_path.c_str());
    }
    if (header.length < min_len) {
      continue;
    }

    sequence.resize(header.length);
    for (uint32_t i = 0; i < header.length; ++i) {
      sequence[i] = bases[(packed_words[i >> 4u] >>
                           ((15u - (i & 15u)) << 1u)) &
                          3u];
    }
    hist->insert(header.length);
    pprintf(">k{}_{} flag={} multi={.4} len={}\n{s}\n", header.kmer_size,
            header.id, header.flag, header.multiplicity, header.length,
            sequence.c_str());
  }
  return true;
}

void PrintStats(Histgram<long long> &hist) {
  const long long total_bases = hist.sum();
  pfprintf(
      stderr,
      "{} contigs, total {} bp, min {} bp, max {} bp, avg {} bp, N50 {} bp\n",
      (int)hist.size(), total_bases, hist.minimum(), hist.maximum(),
      int(hist.mean() + 0.5), hist.Nx(total_bases * 0.5));
}

}  // namespace

int main_filter_by_len(int argc, char **argv) {
  if (argc < 2) {
    pfprintf(stderr,
             "Usage: cat contigs.fa | {s} <min_len>\n"
             "   or: {s} <min_len> <contigs.fa>...\n",
             argv[0], argv[0]);
    exit(1);
  }

  unsigned min_len = atoi(argv[1]);
  Histgram<long long> hist;

  if (argc == 2) {
    gzFile input = gzdopen(fileno(stdin), "r");
    FilterFastx(input, min_len, &hist);
    gzclose(input);
  } else {
    for (int i = 2; i < argc; ++i) {
      if (!FilterPackedFile(argv[i], min_len, &hist)) {
        FilterFastxFile(argv[i], min_len, &hist);
      }
    }
  }

  PrintStats(hist);
  return 0;
}
