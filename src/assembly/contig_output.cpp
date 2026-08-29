//
// Created by vout on 11/21/18.
//

#include "contig_output.h"
#include <cassert>
#include "definitions.h"
#include "unitig_graph.h"

namespace {

inline char Complement(char c) {
  if (c >= 0 && c < 4) {
    return 3 - c;
  }
  switch (c) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      assert(false);
  }
  return 0;
}

inline void ReverseComplement(std::string &s) {
  int i, j;
  for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
    std::swap(s[i], s[j]);
    s[i] = Complement(s[i]);
    s[j] = Complement(s[j]);
  }
  if (i == j) {
    s[i] = Complement(s[i]);
  }
}

size_t MinimalRotation(const std::string &s) {
  if (s.empty()) return 0;
  size_t i = 0;
  size_t j = 1;
  size_t matched = 0;
  while (i < s.size() && j < s.size() && matched < s.size()) {
    const char lhs = s[(i + matched) % s.size()];
    const char rhs = s[(j + matched) % s.size()];
    if (lhs == rhs) {
      ++matched;
      continue;
    }
    if (lhs > rhs) {
      i += matched + 1;
      if (i <= j) i = j + 1;
    } else {
      j += matched + 1;
      if (j <= i) j = i + 1;
    }
    matched = 0;
  }
  return std::min(i, j);
}

bool RotationLess(const std::string &lhs, size_t lhs_offset,
                  const std::string &rhs, size_t rhs_offset) {
  assert(lhs.size() == rhs.size());
  for (size_t i = 0; i < lhs.size(); ++i) {
    const char a = lhs[(lhs_offset + i) % lhs.size()];
    const char b = rhs[(rhs_offset + i) % rhs.size()];
    if (a != b) return a < b;
  }
  return false;
}

// A circular unitig has no biological first base, but its old linearized
// representation inherited one from a lock winner.  That made otherwise
// identical multi-thread runs emit different strings and could feed a
// different rotation into the next k round.  Choose the minimal rotation over
// both strands and then restore the k-base circular overlap.
void CanonicalizeLoop(std::string &s, unsigned kmer_k) {
  if (s.size() <= kmer_k) return;
  const size_t core_size = s.size() - kmer_k;
  std::string forward = s.substr(0, core_size);
  std::string reverse = forward;
  ReverseComplement(reverse);
  const size_t forward_offset = MinimalRotation(forward);
  const size_t reverse_offset = MinimalRotation(reverse);
  const bool use_reverse =
      RotationLess(reverse, reverse_offset, forward, forward_offset);
  const std::string &core = use_reverse ? reverse : forward;
  const size_t offset = use_reverse ? reverse_offset : forward_offset;
  s.resize(core_size + kmer_k);
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = core[(offset + i) % core_size];
  }
}

void FoldPalindrome(std::string &s, unsigned kmer_k, bool is_loop) {
  if (is_loop) {
    for (unsigned i = 1; i + kmer_k <= s.length(); ++i) {
      std::string rc = s.substr(i, kmer_k);
      ReverseComplement(rc);
      if (rc == s.substr(i - 1, kmer_k)) {
        assert(i <= s.length() / 2);
        s = s.substr(i, s.length() / 2);
        break;
      }
    }
  } else {
    int num_edges = s.length() - kmer_k;
    assert(num_edges % 2 == 1);
    s.resize((num_edges - 1) / 2 + kmer_k + 1);
  }
}

}  // namespace

void OutputContigs(UnitigGraph &graph, ContigWriter *contig_writer,
                   ContigWriter *final_contig_writer, bool change_only,
                   uint32_t min_standalone) {
  assert(!(change_only && final_contig_writer != nullptr));  // if output
                                                             // changed contigs,
                                                             // must not output
                                                             // final contigs

#pragma omp parallel for
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(graph.active_id(i));
    if (change_only && !adapter.IsChanged()) {
      continue;
    }
    double multi = change_only ? 1
                               : std::min(static_cast<double>(kMaxMul),
                                          adapter.GetAvgDepth());
    std::string ascii_contig = graph.VertexToDNAString(adapter);

    if (adapter.IsLoop()) {
      int flag = contig_flag::kLoop | contig_flag::kStandalone;
      auto writer = contig_writer;

      if (adapter.IsPalindrome()) {
        FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        flag = contig_flag::kStandalone;
      } else if (!change_only) {
        // change_only output is k*.addi.fa, which is consumed as sequence
        // evidence by the next k round.  Unlike the main contig input,
        // seq2sdbg does not extend circular addi records from the previous k
        // to the new k.  Rotating such a record here therefore changes the
        // linear windows visible at the circular cut and can create or remove
        // ordinary edges in the next graph.  Preserve the deterministic
        // legacy-owner cut for this semantic intermediate; canonicalization
        // remains appropriate for regular/user-visible contig output.
        CanonicalizeLoop(ascii_contig, graph.k());
      }

      if (final_contig_writer != nullptr) {
        if (ascii_contig.length() < min_standalone) {
          continue;
        } else {
          writer = final_contig_writer;
        }
      }
      writer->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    } else {
      auto out_file = contig_writer;
      int flag = 0;

      if (adapter.IsStandalone() ||
          (graph.InDegree(adapter) == 0 && graph.OutDegree(adapter) == 0)) {
        if (adapter.IsPalindrome()) {
          FoldPalindrome(ascii_contig, graph.k(), adapter.IsLoop());
        }
        flag = contig_flag::kStandalone;
        if (final_contig_writer != nullptr) {
          if (ascii_contig.length() < min_standalone) {
            continue;
          } else {
            out_file = final_contig_writer;
          }
        }
      }
      out_file->WriteContig(ascii_contig, graph.k(), i, flag, multi);
    }
  }

  // Drain one partial block per producer thread while this phase is still
  // covered by the output timer. Full blocks are emitted on demand.
  contig_writer->Flush();
  if (final_contig_writer != nullptr && final_contig_writer != contig_writer) {
    final_contig_writer->Flush();
  }
}
