/**
 * @file sequence.h
 * @brief Sequence Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#ifndef __SEQUENCE_SEQUENCE_H_

#define __SEQUENCE_SEQUENCE_H_

#include <stdint.h>

#include <algorithm>
#include <istream>
#include <ostream>
#include <string>

#include "idba/kmer.h"

/**
 * @brief It represents a DNA sequence ({A, C, G, T, N}) as a digit sequence
 * ({0, 1, 2, 3, 4}).
 */
class Sequence {
 public:
  friend std::istream &operator>>(std::istream &stream, Sequence &seq);
  friend std::ostream &operator<<(std::ostream &stream, const Sequence &seq);

  Sequence() {}
  Sequence(const Sequence &seq, int offset = 0,
           size_t length = std::string::npos) {
    Assign(seq, offset, length);
  }
  explicit Sequence(const std::string &seq, int offset = 0,
                    size_t length = std::string::npos) {
    Assign(seq, offset, length);
  }
  Sequence(uint32_t num, uint8_t ch) { Assign(num, ch); }
  explicit Sequence(const IdbaKmer &kmer) { Assign(kmer); }

  ~Sequence() {}

  const Sequence &operator=(const Sequence &seq) {
    Assign(seq);
    return *this;
  }
  const Sequence &operator=(const std::string &seq) {
    Assign(seq);
    return *this;
  }
  const Sequence &operator=(const IdbaKmer &kmer) {
    Assign(kmer);
    return *this;
  }

  const Sequence &operator+=(const Sequence &seq) {
    Append(seq);
    return *this;
  }
  const Sequence &operator+=(uint8_t ch) {
    Append(ch);
    return *this;
  }

  bool operator==(const Sequence &seq) const { return bases_ == seq.bases_; }
  bool operator!=(const Sequence &seq) const { return bases_ != seq.bases_; }
  bool operator<(const Sequence &seq) const { return bases_ < seq.bases_; }
  bool operator>(const Sequence &seq) const { return bases_ > seq.bases_; }

  const Sequence &Assign(const Sequence &seq, int offset = 0,
                         size_t length = std::string::npos) {
    if (&seq != this) bases_.assign(seq.bases_, offset, length);
    return *this;
  }
  const Sequence &Assign(const std::string &s, int offset = 0,
                         size_t length = std::string::npos) {
    bases_.assign(s, offset, length);
    Encode();
    return *this;
  }
  const Sequence &Assign(uint32_t num, uint8_t ch) {
    bases_.assign(num, ch);
    return *this;
  }
  const Sequence &Assign(const IdbaKmer &kmer);

  const Sequence &Append(const Sequence &seq, int offset = 0,
                         size_t length = std::string::npos) {
    bases_.append(seq.bases_, offset, length);
    return *this;
  }
  const Sequence &Append(const std::string &seq, int offset = 0,
                         size_t length = std::string::npos) {
    bases_.append(seq, offset, length);
    return *this;
  }
  const Sequence &Append(uint32_t num, uint8_t ch) {
    bases_.append(num, ch);
    return *this;
  }
  const Sequence &Append(uint8_t ch) {
    bases_.append(1, ch);
    return *this;
  }

  // Append a slice of the reverse-complemented view without first
  // materializing that complete temporary Sequence.  Local graph cleaning
  // repeatedly concatenates short oriented unitigs, so resize once and fill
  // the destination directly.
  const Sequence &AppendReverseComplement(
      const Sequence &seq, size_t offset = 0,
      size_t length = std::string::npos) {
    if (&seq == this) {
      const Sequence copy(seq);
      return AppendReverseComplement(copy, offset, length);
    }
    if (offset > seq.bases_.size()) offset = seq.bases_.size();
    const size_t available = seq.bases_.size() - offset;
    const size_t count = std::min(length, available);
    const size_t old_size = bases_.size();
    bases_.resize(old_size + count);
    for (size_t i = 0; i < count; ++i) {
      const uint8_t base = static_cast<uint8_t>(
          seq.bases_[seq.bases_.size() - 1u - offset - i]);
      bases_[old_size + i] = static_cast<char>(base < 4u ? 3u - base : base);
    }
    return *this;
  }

  const Sequence &AssignReverseComplement(const Sequence &seq) {
    if (&seq == this) return ReverseComplement();
    bases_.clear();
    return AppendReverseComplement(seq);
  }

  const Sequence &ReverseComplement();
  bool IsValid() const;
  bool IsPalindrome() const;

  IdbaKmer GetIdbaKmer(uint32_t offset, uint32_t kmer_size) const;

  uint8_t &operator[](unsigned index) { return (uint8_t &)bases_[index]; }
  const uint8_t &operator[](unsigned index) const {
    return (uint8_t &)bases_[index];
  }
  uint8_t get_base(uint32_t index) const { return (uint8_t)bases_[index]; }
  const uint8_t *encoded_data() const {
    return reinterpret_cast<const uint8_t *>(bases_.data());
  }
  void set_base(uint32_t index, uint8_t ch) { bases_[index] = ch; }

  void swap(Sequence &seq) {
    if (this != &seq) bases_.swap(seq.bases_);
  }

  uint32_t size() const { return bases_.size(); }
  void resize(int new_size) { bases_.resize(new_size); }
  void reserve(size_t capacity) { bases_.reserve(capacity); }
  bool empty() const { return bases_.size() == 0; }

  void clear() { bases_.clear(); }
  std::string str() {
    Sequence tmp = *this;
    tmp.Decode();
    return tmp.bases_;
  }

 protected:
  void Encode();
  void Decode();

  bool IsValid(char ch) const {
    return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 0 ||
           ch == 1 || ch == 2 || ch == 3;
  }

 private:
  std::string bases_;
};

namespace std {
template <>
inline void swap(Sequence &seq1, Sequence &seq2) {
  seq1.swap(seq2);
}
}  // namespace std

std::ostream &WriteFasta(std::ostream &os, const Sequence &seq,
                         const std::string &comment);

#endif
