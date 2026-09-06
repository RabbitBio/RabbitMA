//
// Created by vout on 4/28/19.
//

#ifndef MEGAHIT_CONTIG_READER_H
#define MEGAHIT_CONTIG_READER_H

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include "sequence/io/fastx_reader.h"
#include "sequence/io/contig/packed_contig.h"

class ContigReader : public FastxReader {
 public:
  explicit ContigReader(const std::string &file_name)
      : FastxReader(file_name), file_name_(file_name) {
    packed_input_.open(packed_contig::Path(file_name_),
                       std::ios::binary | std::ios::in);
    if (packed_input_.is_open()) {
      packed_contig::FileHeader header{};
      packed_input_.read(reinterpret_cast<char *>(&header), sizeof(header));
      if (!packed_input_ || !packed_contig::IsValid(header)) {
        xfatal("Invalid packed contig file: {s}\n",
               packed_contig::Path(file_name_).c_str());
      }
      packed_mode_ = true;
    }
  }
  ContigReader *SetMinLen(unsigned min_len) {
    min_len_ = min_len;
    return this;
  }
  ContigReader *SetExtendLoop(unsigned k_from, unsigned k_to) {
    k_from_ = k_from;
    k_to_ = k_to;
    return this;
  }
  ContigReader *SetDiscardFlag(unsigned flag) {
    discard_flag_ = flag;
    return this;
  }
  ContigReader *SetFlagVector(std::vector<unsigned> *flags) {
    flags_ = flags;
    return this;
  }

  std::pair<int64_t, int64_t> GetNumContigsAndBases() const {
    std::ifstream info_fs(file_name_ + ".info");
    int64_t num_contigs, num_bases;
    info_fs >> num_contigs >> num_bases;
    if (!info_fs) {
      xfatal("Invalid format of contig info file: {s}.info",
             file_name_.c_str());
    }

    // Older ContigWriter versions used `length + is_loop ? 28 : 0` without
    // parentheses.  Due to operator precedence that recorded exactly 28
    // bases per non-empty contig, severely under-reserving the compact
    // sequence vector and forcing repeated whole-vector reallocations in the
    // next k-mer round.  Keep the on-disk metadata byte-for-byte compatible,
    // but recognize that historical signature and use the regular FASTA file
    // size as a conservative, O(1) reserve hint (headers make it an upper
    // bound for MEGAHIT's unwrapped FASTA output).
    const char *disable_file_size_reserve =
        std::getenv("MEGAHIT_DISABLE_FILESIZE_CONTIG_RESERVE");
    if ((disable_file_size_reserve == nullptr ||
         std::strcmp(disable_file_size_reserve, "1") != 0) &&
        num_contigs >= 0 &&
        num_contigs <= std::numeric_limits<int64_t>::max() / 28 &&
        num_bases == num_contigs * 28) {
      struct stat file_stat {};
      if (::stat(file_name_.c_str(), &file_stat) == 0 &&
          S_ISREG(file_stat.st_mode) && file_stat.st_size > 0 &&
          static_cast<uint64_t>(file_stat.st_size) <=
              static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
        num_bases = std::max(num_bases,
                             static_cast<int64_t>(file_stat.st_size));
      }
    }
    return {num_contigs, num_bases};
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override {
    return ReadWithMultiplicity<float>(pkg, nullptr, max_num, max_num_bases,
                                       reverse);
  }

  template <typename TMul>
  int64_t ReadAllWithMultiplicity(SeqPackage *pkg, std::vector<TMul> *mul,
                                  bool reverse) {
    return ReadWithMultiplicity(pkg, mul, kMaxNumSeq, kMaxNumBases, reverse);
  }

  template <typename TMul>
  int64_t ReadWithMultiplicity(SeqPackage *pkg, std::vector<TMul> *mul,
                               int64_t max_num, int64_t max_num_bases,
                               bool reverse) {
    if (packed_mode_) {
      return ReadPackedWithMultiplicity(pkg, mul, max_num, max_num_bases,
                                        reverse);
    }
    bool extend_loop = k_from_ < k_to_ && !(discard_flag_ & contig_flag::kLoop);

    int64_t num_bases = 0;
    for (int64_t ri = 0; ri < max_num; ++ri) {
      auto record = ReadNext();
      if (record) {
        if (record->seq.l < min_len_) {
          --ri;
          continue;
        }
        // comment = "flag=x multi=xx.xxxx"
        unsigned flag = record->comment.s[5] - '0';
        if (discard_flag_ & flag) {
          --ri;
          continue;
        }

        if (extend_loop && (flag & contig_flag::kLoop)) {
          if (record->seq.l < k_to_ + 1U) {
            continue;
          }
          std::string ss(record->seq.s);
          for (unsigned i = k_from_; i < k_to_; ++i) {
            ss.push_back(ss[i]);
          }

          if (reverse) {
            pkg->AppendReversedStringSequence(ss.c_str(), ss.length());
          } else {
            pkg->AppendStringSequence(ss.c_str(), ss.length());
          }
        } else {
          if (reverse) {
            pkg->AppendReversedStringSequence(record->seq.s, record->seq.l);
          } else {
            pkg->AppendStringSequence(record->seq.s, record->seq.l);
          }
        }

        if (mul) {
          mul->push_back(GetMultiplicity<TMul>(record->comment.s));
        }
        if (flags_ != nullptr) flags_->push_back(flag);

        num_bases += record->seq.l;
        if (num_bases >= max_num_bases) {
          return ri + 1;
        }
      } else {
        return ri;
      }
    }
    return max_num;
  }

 private:
  template <typename TMul>
  static TMul GetMultiplicity(const char *fastx_comment) {
    auto m = atof(fastx_comment + 13);
    if (std::is_integral<TMul>::value) {
      return m + .5;
    } else {
      return m;
    }
  }

  template <typename TMul>
  int64_t ReadPackedWithMultiplicity(SeqPackage *pkg,
                                     std::vector<TMul> *mul,
                                     int64_t max_num,
                                     int64_t max_num_bases, bool reverse) {
    const bool extend_loop =
        k_from_ < k_to_ && !(discard_flag_ & contig_flag::kLoop);
    int64_t accepted = 0;
    int64_t num_bases = 0;
    while (accepted < max_num) {
      packed_contig::RecordHeader header{};
      packed_input_.read(reinterpret_cast<char *>(&header), sizeof(header));
      if (packed_input_.gcount() == 0 && packed_input_.eof()) return accepted;
      if (!packed_input_) {
        xfatal("Truncated packed contig header: {s}\n", file_name_.c_str());
      }
      const size_t words = (static_cast<size_t>(header.length) + 15u) / 16u;
      packed_words_.resize(words);
      packed_input_.read(reinterpret_cast<char *>(packed_words_.data()),
                         words * sizeof(uint32_t));
      if (!packed_input_) {
        xfatal("Truncated packed contig sequence: {s}\n",
               file_name_.c_str());
      }
      if (header.length < min_len_ ||
          (discard_flag_ & static_cast<unsigned>(header.flag)) != 0u) {
        continue;
      }

      if (extend_loop &&
          (static_cast<unsigned>(header.flag) & contig_flag::kLoop) != 0u) {
        if (header.length < k_to_ + 1u) continue;
        packed_ascii_.resize(header.length);
        static const char bases[] = {'A', 'C', 'G', 'T'};
        for (uint32_t i = 0; i < header.length; ++i) {
          packed_ascii_[i] = bases[(packed_words_[i >> 4u] >>
                                    ((15u - (i & 15u)) << 1u)) &
                                   3u];
        }
        for (unsigned i = k_from_; i < k_to_; ++i) {
          packed_ascii_.push_back(packed_ascii_[i]);
        }
        if (reverse) {
          pkg->AppendReversedStringSequence(packed_ascii_.data(),
                                            packed_ascii_.size());
        } else {
          pkg->AppendStringSequence(packed_ascii_.data(),
                                    packed_ascii_.size());
        }
      } else if (reverse) {
        pkg->AppendReversedCompactSequence(packed_words_.data(),
                                           header.length);
      } else {
        pkg->AppendCompactSequence(packed_words_.data(), header.length);
      }
      if (mul != nullptr) {
        if (std::is_integral<TMul>::value) {
          mul->push_back(static_cast<TMul>(header.multiplicity + .5f));
        } else {
          mul->push_back(static_cast<TMul>(header.multiplicity));
        }
      }
      ++accepted;
      if (flags_ != nullptr) flags_->push_back(static_cast<unsigned>(header.flag));
      num_bases += header.length;
      if (num_bases >= max_num_bases) return accepted;
    }
    return accepted;
  }

 private:
  unsigned min_len_{0};
  unsigned k_from_{0}, k_to_{0};
  unsigned discard_flag_{0};
  std::vector<unsigned> *flags_{nullptr};
  std::string file_name_;
  bool packed_mode_{false};
  std::ifstream packed_input_;
  std::vector<uint32_t> packed_words_;
  std::string packed_ascii_;
};

#endif  // MEGAHIT_CONTIG_READER_H
