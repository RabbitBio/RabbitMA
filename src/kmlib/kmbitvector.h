//
// Created by vout on 28/2/2018.
//

#ifndef KMLIB_ATOMIC_BIT_VECTOR_H
#define KMLIB_ATOMIC_BIT_VECTOR_H

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <memory>
#include <new>
#include <thread>
#include <type_traits>

#include <omp.h>

namespace kmlib {
/*!
 * @brief Atomic bit vector: a class that represent a vector of "bits".
 * @details Update of each bit is threads safe via set and get.
 * It can also be used as a vector of bit locks via try_lock, lock and unlock
 */
template <typename WordType = unsigned long>
class AtomicBitVector {
 private:
  // std::atomic is neither copyable nor movable, but the wrapper gives each
  // word an explicit construction point.  Large arrays are allocated as raw
  // virtual memory and placement-constructed by the OpenMP team so physical
  // pages follow first touch across NUMA nodes instead of all landing on the
  // thread that happened to call the vector constructor.
  template <typename T>
  struct AtomicWrapper {
    std::atomic<T> v;
    explicit AtomicWrapper(T a = T()) : v(a) {}
    AtomicWrapper(const AtomicWrapper &rhs)
        : v(rhs.v.load(std::memory_order_relaxed)) {}
    AtomicWrapper &operator=(const AtomicWrapper &rhs) {
      v.store(rhs.v.load(std::memory_order_relaxed),
              std::memory_order_relaxed);
      return *this;
    }
  };

 public:
  using word_type = WordType;
  using size_type = size_t;

 public:
  /*!
   * @brief Constructor
   * @param size the size (number of bits) of the bit vector
   */
  explicit AtomicBitVector(size_type size = 0)
      : size_(size),
        word_count_((size + kBitsPerWord - 1) / kBitsPerWord) {
    AllocateAndZero();
  }
  /*!
   * @brief Construct a bit vector from iterators of words
   * @tparam WordIterator iterator to access words
   * @param first the iterator pointing to the first word
   * @param last the iterator pointing to the last word
   */
  template <typename WordIterator>
  explicit AtomicBitVector(WordIterator first, WordIterator last)
      : size_(static_cast<size_type>(last - first) * kBitsPerWord),
        word_count_(static_cast<size_type>(last - first)) {
    AllocateRaw();
    const int64_t count = static_cast<int64_t>(word_count_);
#pragma omp parallel for schedule(static) if (word_count_ * sizeof(word_type) >= kParallelFirstTouchMinBytes)
    for (int64_t i = 0; i < count; ++i) {
      new (data_array_ + i)
          AtomicWrapper<word_type>(static_cast<word_type>(*(first + i)));
    }
  }
  /*!
   * @brief the move constructor
   */
  AtomicBitVector(AtomicBitVector &&rhs)
      : size_(rhs.size_),
        word_count_(rhs.word_count_),
        data_array_(rhs.data_array_) {
    rhs.size_ = 0;
    rhs.word_count_ = 0;
    rhs.data_array_ = nullptr;
  }
  /*!
   * @brief the move operator
   */
  AtomicBitVector &operator=(AtomicBitVector &&rhs) noexcept {
    if (this == &rhs) {
      return *this;
    }
    Release();
    size_ = rhs.size_;
    word_count_ = rhs.word_count_;
    data_array_ = rhs.data_array_;
    rhs.size_ = 0;
    rhs.word_count_ = 0;
    rhs.data_array_ = nullptr;
    return *this;
  }
  AtomicBitVector(const AtomicBitVector &) = delete;
  AtomicBitVector &operator=(const AtomicBitVector &) = delete;
  ~AtomicBitVector() { Release(); }

  /*!
   * @return the size of the bit vector
   */
  size_type size() const { return size_; }

  size_type word_count() const { return word_count_; }

  word_type load_word(size_type word_index) const {
    assert(word_index < word_count_);
    return data_array_[word_index].v.load(std::memory_order_relaxed);
  }

  // Pointer-chasing kernels can discover the next bit address one stage
  // before the atomic update.  Pulling the containing word into an exclusive
  // cache state here overlaps that ownership request with other independent
  // paths without exposing the internal atomic storage representation.
  void prefetch_for_write(size_type i) {
    assert(i < size_);
    __builtin_prefetch(data_array_ + i / kBitsPerWord, 1, 1);
  }

  // Bulk builders can accumulate a complete word in a register and perform
  // one store instead of up to kBitsPerWord atomic read-modify-writes.  The
  // caller must have exclusive ownership of the word.
  void store_word(size_type word_index, word_type value) {
    assert(word_index < word_count_);
    data_array_[word_index].v.store(value, std::memory_order_relaxed);
  }

  static size_type bits_per_word() { return kBitsPerWord; }

  /*!
   * @brief set the i-th bit to 1
   * @param i the index of the bit to be set to 1
   */
  void set(size_type i) {
    word_type mask = word_type(1) << (i % kBitsPerWord);
    data_array_[i / kBitsPerWord].v.fetch_or(mask, std::memory_order_relaxed);
  }

  /*!
   * @brief set the i-th bit to 0
   * @param i the index of the bit to be set to 0
   */
  void unset(size_type i) {
    word_type mask = ~(word_type(1) << (i % kBitsPerWord));
    data_array_[i / kBitsPerWord].v.fetch_and(mask, std::memory_order_relaxed);
  }

  // Atomically clear a bit and report whether this call performed the 1->0
  // transition.  Frontier builders use the transition itself as an exact,
  // allocation-free deduplication event.
  bool try_unset(size_type i) {
    const word_type bit = word_type(1) << (i % kBitsPerWord);
    const word_type old = data_array_[i / kBitsPerWord].v.fetch_and(
        ~bit, std::memory_order_relaxed);
    return (old & bit) != 0;
  }

  /*!
   * @param i the index of the bit
   * @return value of the i-th bit
   */
  bool at(size_type i) const {
    return !!(data_array_[i / kBitsPerWord].v.load(std::memory_order_relaxed) &
              (word_type(1) << i % kBitsPerWord));
  }

  /*!
   * @param i the index of the bit
   * @return whether the i-th bit has been locked successfully
   */
  bool try_lock(size_type i) {
    word_type mask = word_type(1) << (i % kBitsPerWord);
    word_type old_val = data_array_[i / kBitsPerWord].v.fetch_or(
        mask, std::memory_order_acquire);
    return !(old_val & mask);
  }

  /*!
   * @brief lock the i-th bit
   * @param i the bit to lock
   */
  void lock(size_type i) {
    unsigned retry = 0;
    while (!try_lock(i)) {
      if (++retry > 64) {
        retry = 0;
        std::this_thread::yield();
      }
    }
    assert(at(i));
  }

  /*!
   * @brief unlock the i-th bit
   * @param i the index of the bits
   */
  void unlock(size_type i) {
    auto mask = word_type(1) << (i % kBitsPerWord);
    auto old_val = data_array_[i / kBitsPerWord].v.fetch_and(
        ~mask, std::memory_order_release);
    assert(old_val & mask);
    (void)(old_val);
  }

  /*!
   * @brief reset the size of the bit vector and clear all bits
   * @param size the new size of the bit vector
   */
  void reset(size_type size) {
    if (size == size_) {
      reset();
      return;
    }
    AtomicBitVector replacement(size);
    swap(replacement);
  }

  void reset() {
    const int64_t count = static_cast<int64_t>(word_count_);
#pragma omp parallel for schedule(static) if (word_count_ * sizeof(word_type) >= kParallelFirstTouchMinBytes)
    for (int64_t i = 0; i < count; ++i) {
      data_array_[i].v.store(0, std::memory_order_relaxed);
    }
  }

  /*!
   * @brief swap with another bit vector
   * @param rhs the target to swap
   */
  void swap(AtomicBitVector &rhs) {
    std::swap(size_, rhs.size_);
    std::swap(word_count_, rhs.word_count_);
    std::swap(data_array_, rhs.data_array_);
  }

 private:
  void AllocateRaw() {
    if (word_count_ == 0) {
      data_array_ = nullptr;
      return;
    }
    data_array_ = static_cast<AtomicWrapper<word_type> *>(
        ::operator new[](word_count_ * sizeof(AtomicWrapper<word_type>)));
  }

  void AllocateAndZero() {
    AllocateRaw();
    const int64_t count = static_cast<int64_t>(word_count_);
#pragma omp parallel for schedule(static) if (word_count_ * sizeof(word_type) >= kParallelFirstTouchMinBytes)
    for (int64_t i = 0; i < count; ++i) {
      new (data_array_ + i) AtomicWrapper<word_type>(0);
    }
  }

  void Release() noexcept {
    static_assert(
        std::is_trivially_destructible<AtomicWrapper<word_type>>::value,
        "Atomic bit-vector words must remain trivially destructible");
    ::operator delete[](data_array_);
    data_array_ = nullptr;
    word_count_ = 0;
    size_ = 0;
  }

  static const unsigned kBitsPerByte = 8;
  static const unsigned kBitsPerWord = sizeof(word_type) * kBitsPerByte;
  static const size_t kParallelFirstTouchMinBytes = size_t{1} << 20u;
  size_type size_;
  size_type word_count_;
  AtomicWrapper<word_type> *data_array_;
  static_assert(sizeof(AtomicWrapper<word_type>) == sizeof(word_type), "");
};

}  // namespace kmlib

using AtomicBitVector = kmlib::AtomicBitVector<>;

#endif  // KMLIB_ATOMIC_BIT_VECTOR_H
