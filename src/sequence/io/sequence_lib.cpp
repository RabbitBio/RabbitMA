//
// Created by vout on 6/29/19.
//

#include "sequence_lib.h"
#include "async_sequence_reader.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <omp.h>

#include "kmlib/kmbit.h"
#include "sequence/io/read_chunk_index.h"
#include "utils/startup_affinity.h"

namespace {

struct InputLibrary {
  std::string metadata;
  std::string type;
  std::string file_name1;
  std::string file_name2;
};

bool ReadLibraryPath(std::istream &input, std::string *path) {
  input >> std::ws;
  if (!input.good()) {
    return false;
  }

  if (input.peek() != '"') {
    return static_cast<bool>(input >> *path);
  }

  input.get();
  path->clear();
  bool escaped = false;
  char ch = '\0';
  while (input.get(ch)) {
    if (escaped) {
      path->push_back(ch);
      escaped = false;
    } else if (ch == '\\') {
      escaped = true;
    } else if (ch == '"') {
      return true;
    } else {
      path->push_back(ch);
    }
  }
  return false;
}

struct BuiltLibrary {
  int64_t num_reads{0};
  int64_t num_bases{0};
  unsigned max_read_len{0};
  uint64_t num_words{0};
  std::vector<PackedReadChunk> chunks;
  std::string part_path;
  std::string error;
};

uint64_t RegularFileSize(const std::string &path) {
  struct stat st;
  if (path == "-" || stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode) ||
      st.st_size <= 0) {
    return 0;
  }
  return static_cast<uint64_t>(st.st_size);
}

uint64_t InputWorkBytes(const InputLibrary &lib) {
  return RegularFileSize(lib.file_name1) + RegularFileSize(lib.file_name2);
}

unsigned InputStreamCount(const InputLibrary &lib) {
  return lib.type == "pe" ? 2u : 1u;
}

/**
 * Global budget for private decompression and parsing pools.
 *
 * OpenMP library workers are the FASTQ parser/packer population.  Decoder
 * helper threads are leased in addition to those workers, and every lease is
 * bounded by ceil(total auxiliary slots / active libraries).  A library may
 * divide its lease between gzip decoding and ordered FASTQ packing.  This lets
 * one large input use nearly all requested CPUs without multiplying `-t` by
 * the number of files when several libraries are active.
 */
class AuxiliaryThreadBudget {
 public:
  AuxiliaryThreadBudget(unsigned total, unsigned active_libraries)
      : available_(total), total_(total), max_per_library_(
            active_libraries == 0
                ? 0
                : DivCeiling<unsigned>(total, active_libraries)) {}

  unsigned Acquire(unsigned requested) {
    unsigned available = available_.load(std::memory_order_relaxed);
    while (available != 0 && requested != 0) {
      const unsigned grant =
          std::min(requested, std::min(available, max_per_library_));
      if (available_.compare_exchange_weak(
              available, available - grant, std::memory_order_acq_rel,
              std::memory_order_relaxed)) {
        return grant;
      }
    }
    return 0;
  }

  void Release(unsigned count) {
    if (count != 0) {
      available_.fetch_add(count, std::memory_order_release);
    }
  }

  unsigned total() const { return total_; }
  unsigned max_per_library() const { return max_per_library_; }

 private:
  std::atomic<unsigned> available_;
  unsigned total_;
  unsigned max_per_library_;
};

class AuxiliaryThreadLease {
 public:
  AuxiliaryThreadLease(AuxiliaryThreadBudget *budget, unsigned requested)
      : budget_(budget),
        count_(budget == nullptr ? 0 : budget->Acquire(requested)) {}
  ~AuxiliaryThreadLease() {
    if (budget_ != nullptr) {
      budget_->Release(count_);
    }
  }
  unsigned count() const { return count_; }

 private:
  AuxiliaryThreadBudget *budget_;
  unsigned count_;
};

struct LibraryIoPlan {
  unsigned file1_threads{0};
  unsigned file2_threads{0};
};

LibraryIoPlan SplitDecoderThreads(const InputLibrary &lib,
                                  unsigned decoder_threads) {
  LibraryIoPlan plan;
  if (decoder_threads == 0) {
    return plan;
  }
  if (lib.type != "pe") {
    plan.file1_threads = decoder_threads;
    return plan;
  }

  const uint64_t size1 = RegularFileSize(lib.file_name1);
  const uint64_t size2 = RegularFileSize(lib.file_name2);
  if (size1 == 0 && size2 != 0) {
    plan.file2_threads = decoder_threads;
    return plan;
  }
  if (size2 == 0 && size1 != 0) {
    plan.file1_threads = decoder_threads;
    return plan;
  }
  if (decoder_threads == 1) {
    if (size2 > size1) {
      plan.file2_threads = 1;
    } else {
      plan.file1_threads = 1;
    }
    return plan;
  }

  // Keep both mate streams live, then divide the remaining work in
  // proportion to compressed bytes.  This is input-work scheduling, not a
  // machine- or dataset-specific threshold.
  const long double total_size =
      static_cast<long double>(size1) + static_cast<long double>(size2);
  unsigned file1_threads = decoder_threads / 2;
  if (total_size > 0) {
    file1_threads = static_cast<unsigned>(
        static_cast<long double>(decoder_threads) * size1 / total_size +
        0.5L);
  }
  file1_threads = std::max(1u, std::min(decoder_threads - 1, file1_threads));
  plan.file1_threads = file1_threads;
  plan.file2_threads = decoder_threads - file1_threads;
  return plan;
}

/**
 * Encode FASTX records directly into MEGAHIT's binary-library format.
 *
 * The old parallel build path first appended every base to a SeqPackage and
 * then traversed that package again to realign each record into `[len][words]`.
 * RabbitFX's useful lesson here is to hand bounded chunks directly from the
 * parser to the consumer.  This writer does exactly one base pass and keeps a
 * reusable output chunk; memory is bounded by one chunk (or one exceptionally
 * long read), independent of input size.
 */
class PackedBinaryWriter {
 public:
  explicit PackedBinaryWriter(std::ostream *output)
      : output_(output), buffer_(8u << 20u), used_(0) {
    for (int i = 0; i < 10; ++i) {
      dna_map_[static_cast<unsigned char>("ACGTNacgtn"[i])] =
          static_cast<unsigned char>("0123201232"[i] - '0');
    }
  }

  void Append(const char *sequence, unsigned length, BuiltLibrary *result) {
    int begin = 0;
    int end = length;
    FastxReader::TrimN(sequence, length, &begin, &end);

    uint32_t packed_length = static_cast<uint32_t>(end - begin);
    const char *packed_sequence = sequence + begin;
    // SequencePackage represents an empty post-trim read as one 'A'.  Match
    // that historical format and its metadata exactly.
    if (packed_length == 0) {
      packed_length = 1;
      packed_sequence = nullptr;
    }

    const size_t num_words =
        DivCeiling<size_t>(packed_length, SeqPackage::kBasesPerWord);
    const size_t record_bytes = sizeof(uint32_t) * (num_words + 1);
    Ensure(record_bytes);
    std::memcpy(buffer_.data() + used_, &packed_length, sizeof(packed_length));
    used_ += sizeof(packed_length);

    for (size_t word_id = 0; word_id < num_words; ++word_id) {
      uint32_t word = 0;
      const unsigned offset = word_id * SeqPackage::kBasesPerWord;
      const unsigned take = std::min<unsigned>(
          SeqPackage::kBasesPerWord, packed_length - offset);
      if (packed_sequence != nullptr) {
        for (unsigned j = 0; j < take; ++j) {
          word |= static_cast<uint32_t>(
                      dna_map_[static_cast<unsigned char>(
                          packed_sequence[offset + j])])
                  << SeqPackage::TVector::bit_shift(j);
        }
      }
      std::memcpy(buffer_.data() + used_, &word, sizeof(word));
      used_ += sizeof(word);
    }

    ++result->num_reads;
    result->num_bases += packed_length;
    result->max_read_len = std::max(result->max_read_len, packed_length);
    total_words_ += num_words + 1u;
    chunk_bases_ += packed_length;
    chunk_max_read_len_ = std::max(chunk_max_read_len_, packed_length);
    if (total_words_ - chunk_word_begin_ >= kChunkWords) {
      FinishChunk(result);
    }
  }

  /**
   * Commit a batch that was encoded by a worker thread.
   *
   * `packed_lengths` stays in original FASTX order.  Updating the global
   * counters and chunk boundaries here, rather than in workers, preserves the
   * exact historical `.bin` layout and read-index partitioning even when
   * batches finish out of order.
   */
  void AppendEncodedBatch(const std::vector<char> &encoded,
                          const std::vector<uint32_t> &packed_lengths,
                          BuiltLibrary *result) {
    Flush();
    if (!encoded.empty()) {
      output_->write(encoded.data(), encoded.size());
      if (!*output_) {
        throw std::runtime_error("failed writing temporary binary library");
      }
    }
    for (uint32_t packed_length : packed_lengths) {
      const size_t num_words =
          DivCeiling<size_t>(packed_length, SeqPackage::kBasesPerWord);
      ++result->num_reads;
      result->num_bases += packed_length;
      result->max_read_len =
          std::max(result->max_read_len, packed_length);
      total_words_ += num_words + 1u;
      chunk_bases_ += packed_length;
      chunk_max_read_len_ =
          std::max(chunk_max_read_len_, packed_length);
      if (total_words_ - chunk_word_begin_ >= kChunkWords) {
        FinishChunk(result);
      }
    }
  }

  void Finish(BuiltLibrary *result) {
    FinishChunk(result);
    result->num_words = total_words_;
    Flush();
    if (!*output_) {
      throw std::runtime_error("failed writing temporary binary library");
    }
  }

 private:
  static constexpr uint64_t kChunkWords =
      (uint64_t{1} << 20u) / sizeof(uint32_t);

  void FinishChunk(BuiltLibrary *result) {
    if (chunk_read_begin_ == static_cast<uint64_t>(result->num_reads)) {
      return;
    }
    result->chunks.push_back(PackedReadChunk{
        chunk_word_begin_, total_words_, chunk_read_begin_,
        static_cast<uint64_t>(result->num_reads), chunk_bases_,
        chunk_max_read_len_});
    chunk_word_begin_ = total_words_;
    chunk_read_begin_ = static_cast<uint64_t>(result->num_reads);
    chunk_bases_ = 0;
    chunk_max_read_len_ = 0;
  }

  void Ensure(size_t bytes) {
    if (bytes <= buffer_.size() - used_) {
      return;
    }
    Flush();
    if (bytes > buffer_.size()) {
      buffer_.resize(bytes);
    }
  }

  void Flush() {
    if (used_ != 0) {
      output_->write(buffer_.data(), used_);
      used_ = 0;
    }
  }

  std::ostream *output_;
  std::vector<char> buffer_;
  size_t used_;
  uint64_t total_words_{0};
  uint64_t chunk_word_begin_{0};
  uint64_t chunk_read_begin_{0};
  uint64_t chunk_bases_{0};
  unsigned chunk_max_read_len_{0};
  unsigned char dna_map_[256]{};
};

/**
 * Signals a valid FASTX input whose layout is not suitable for the parallel
 * FASTQ path.  The caller discards the partial part file and replays the
 * library through kseq, so unusual FASTA/FASTQ layouts retain their existing
 * semantics instead of becoming a new compatibility restriction.
 */
class UnsupportedParallelFastq : public std::runtime_error {
 public:
  explicit UnsupportedParallelFastq(const std::string &message)
      : std::runtime_error(message) {}
};

constexpr size_t kFastqBatchBytes = size_t{8} << 20u;
constexpr size_t kFastqReadBytes = size_t{4} << 20u;
constexpr size_t kFastqPipelineByteBudget = size_t{512} << 20u;
constexpr unsigned kMaxFastqParserWorkers = static_cast<unsigned>(
    kFastqPipelineByteBudget / kFastqBatchBytes);

size_t FastqLineLength(const std::vector<char> &data, size_t begin,
                       size_t end) {
  if (end > begin && data[end - 1] == '\r') {
    --end;
  }
  return end - begin;
}

/**
 * Split an ordered decompressed stream only at complete FASTQ records.
 *
 * This small state machine supports wrapped sequence and quality lines; it
 * does not assume four physical lines per record.  It scans only line
 * boundaries on the producer thread.  The more expensive sequence packing is
 * performed later by independent workers.  If the input is FASTA or uses a
 * layout outside kseq-compatible FASTQ, the caller can safely restart through
 * the legacy parser because output is staged in a per-library part file.
 */
class FastqRecordBatchReader {
 public:
  FastqRecordBatchReader(const std::string &path, unsigned gzip_threads)
      : stream_(mgz_open(path, false, gzip_threads)) {
    if (stream_ == nullptr) {
      throw std::runtime_error("cannot open FASTQ input " + path);
    }
    data_.reserve(kFastqBatchBytes + kFastqReadBytes);
  }

  ~FastqRecordBatchReader() { mgz_close(stream_); }

  bool Next(std::vector<char> *batch, uint64_t *record_count = nullptr,
            uint64_t record_limit = 0) {
    batch->clear();
    if (record_count != nullptr) {
      *record_count = 0;
    }
    if (finished_) {
      return false;
    }

    while (true) {
      const size_t emit_end = ScanCompleteLines(record_limit);
      if (emit_end != 0) {
        Emit(emit_end, batch, record_count);
        return true;
      }

      if (eof_) {
        if (line_begin_ < data_.size()) {
          const size_t begin = line_begin_;
          const size_t end = data_.size();
          line_begin_ = end;
          ProcessLine(begin, FastqLineLength(data_, begin, end), end);
        }
        if (state_ != State::kHeader) {
          throw UnsupportedParallelFastq(
              "truncated or non-FASTQ record in parallel parser");
        }
        if (records_seen_ == 0) {
          // A batch may end exactly on the previous record and leave only
          // trailing empty lines for the final refill.
          data_.clear();
          finished_ = true;
          return false;
        }
        if (!data_.empty()) {
          // All bytes have now been validated.  Include harmless trailing
          // empty lines; workers ignore them just as kseq does.
          Emit(data_.size(), batch, record_count);
          finished_ = true;
          return !batch->empty();
        }
        finished_ = true;
        return false;
      }

      ReadMore();
    }
  }

 private:
  enum class State { kHeader, kSequence, kQuality };

  size_t ScanCompleteLines(uint64_t record_limit) {
    while (line_begin_ < data_.size()) {
      const void *newline = std::memchr(data_.data() + line_begin_, '\n',
                                        data_.size() - line_begin_);
      if (newline == nullptr) {
        break;
      }
      const size_t line_end =
          static_cast<const char *>(newline) - data_.data();
      const size_t after = line_end + 1;
      const size_t begin = line_begin_;
      line_begin_ = after;
      ProcessLine(begin, FastqLineLength(data_, begin, line_end), after);
      if ((record_limit != 0 && records_seen_ >= record_limit) ||
          (record_limit == 0 && last_record_end_ >= kFastqBatchBytes)) {
        return last_record_end_;
      }
    }
    return 0;
  }

  void ProcessLine(size_t begin, size_t length, size_t after) {
    switch (state_) {
      case State::kHeader:
        if (length == 0) {
          return;
        }
        if (data_[begin] != '@') {
          throw UnsupportedParallelFastq(
              "input is not line-oriented FASTQ");
        }
        sequence_length_ = 0;
        quality_length_ = 0;
        state_ = State::kSequence;
        return;

      case State::kSequence:
        if (length != 0 && data_[begin] == '+') {
          if (sequence_length_ == 0) {
            FinishRecord(after);
          } else {
            state_ = State::kQuality;
          }
          return;
        }
        if (sequence_length_ >
            std::numeric_limits<size_t>::max() - length) {
          throw UnsupportedParallelFastq("FASTQ sequence is too long");
        }
        sequence_length_ += length;
        return;

      case State::kQuality:
        if (quality_length_ >
            std::numeric_limits<size_t>::max() - length) {
          throw UnsupportedParallelFastq("FASTQ quality is too long");
        }
        quality_length_ += length;
        if (quality_length_ > sequence_length_) {
          throw UnsupportedParallelFastq(
              "FASTQ sequence and quality lengths differ");
        }
        if (quality_length_ == sequence_length_) {
          FinishRecord(after);
        }
        return;
    }
  }

  void FinishRecord(size_t after) {
    state_ = State::kHeader;
    last_record_end_ = after;
    ++records_seen_;
  }

  void ReadMore() {
    const size_t old_size = data_.size();
    data_.resize(old_size + kFastqReadBytes);
    const int got = mgz_read(stream_, data_.data() + old_size,
                             static_cast<unsigned>(kFastqReadBytes));
    if (got < 0) {
      const std::string message =
          stream_->error.empty() ? "FASTQ decompression failed"
                                 : stream_->error;
      data_.resize(old_size);
      throw std::runtime_error(message);
    }
    if (got == 0) {
      data_.resize(old_size);
      eof_ = true;
      return;
    }
    data_.resize(old_size + static_cast<size_t>(got));
  }

  void Emit(size_t end, std::vector<char> *batch, uint64_t *record_count) {
    if (records_seen_ == 0) {
      throw UnsupportedParallelFastq("FASTQ batch contains no records");
    }
    if (record_count != nullptr) {
      *record_count = records_seen_;
    }
    batch->swap(data_);
    data_.assign(batch->begin() + end, batch->end());
    batch->resize(end);
    data_.reserve(kFastqBatchBytes + kFastqReadBytes);
    line_begin_ = 0;
    last_record_end_ = 0;
    records_seen_ = 0;
    // Emission is possible only immediately after a complete record.
    state_ = State::kHeader;
    sequence_length_ = 0;
    quality_length_ = 0;
  }

  mgzFile stream_{nullptr};
  std::vector<char> data_;
  size_t line_begin_{0};
  size_t last_record_end_{0};
  size_t sequence_length_{0};
  size_t quality_length_{0};
  uint64_t records_seen_{0};
  State state_{State::kHeader};
  bool eof_{false};
  bool finished_{false};
};

struct RawFastqBatch {
  uint64_t ordinal{0};
  std::vector<char> data;
  std::vector<char> mate_data;
  bool paired{false};
};

struct EncodedFastqBatch {
  uint64_t ordinal{0};
  std::vector<char> data;
  std::vector<uint32_t> packed_lengths;
  std::exception_ptr error;
};

bool NextFastqLine(const std::vector<char> &data, size_t *position,
                   const char **line, size_t *length) {
  if (*position >= data.size()) {
    return false;
  }
  const size_t begin = *position;
  const void *newline =
      std::memchr(data.data() + begin, '\n', data.size() - begin);
  size_t end = newline == nullptr
                   ? data.size()
                   : static_cast<const char *>(newline) - data.data();
  *position = newline == nullptr ? data.size() : end + 1;
  if (end > begin && data[end - 1] == '\r') {
    --end;
  }
  *line = data.data() + begin;
  *length = end - begin;
  return true;
}

void AppendPackedSequence(const char *sequence, size_t length,
                          const unsigned char *dna_map,
                          EncodedFastqBatch *output) {
  if (length > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw UnsupportedParallelFastq("FASTQ sequence exceeds MEGAHIT limits");
  }
  int begin = 0;
  int end = static_cast<int>(length);
  FastxReader::TrimN(sequence, static_cast<unsigned>(length), &begin, &end);

  uint32_t packed_length = static_cast<uint32_t>(end - begin);
  const char *packed_sequence = sequence + begin;
  if (packed_length == 0) {
    packed_length = 1;
    packed_sequence = nullptr;
  }

  const size_t num_words =
      DivCeiling<size_t>(packed_length, SeqPackage::kBasesPerWord);
  const size_t old_size = output->data.size();
  output->data.resize(old_size + sizeof(uint32_t) * (num_words + 1));
  char *destination = output->data.data() + old_size;
  std::memcpy(destination, &packed_length, sizeof(packed_length));
  destination += sizeof(packed_length);

  for (size_t word_id = 0; word_id < num_words; ++word_id) {
    uint32_t word = 0;
    const unsigned offset = word_id * SeqPackage::kBasesPerWord;
    const unsigned take = std::min<unsigned>(
        SeqPackage::kBasesPerWord, packed_length - offset);
    if (packed_sequence != nullptr) {
      for (unsigned j = 0; j < take; ++j) {
        word |= static_cast<uint32_t>(
                    dna_map[static_cast<unsigned char>(
                        packed_sequence[offset + j])])
                << SeqPackage::TVector::bit_shift(j);
      }
    }
    std::memcpy(destination, &word, sizeof(word));
    destination += sizeof(word);
  }
  output->packed_lengths.push_back(packed_length);
}

EncodedFastqBatch EncodeFastqData(const std::vector<char> &input_data,
                                  uint64_t ordinal) {
  EncodedFastqBatch output;
  output.ordinal = ordinal;
  output.data.reserve(input_data.size() / 4u);
  output.packed_lengths.reserve(input_data.size() / 256u);

  unsigned char dna_map[256]{};
  for (int i = 0; i < 10; ++i) {
    dna_map[static_cast<unsigned char>("ACGTNacgtn"[i])] =
        static_cast<unsigned char>("0123201232"[i] - '0');
  }

  size_t position = 0;
  std::string joined_sequence;
  while (position < input_data.size()) {
    const char *header = nullptr;
    size_t header_length = 0;
    if (!NextFastqLine(input_data, &position, &header, &header_length)) {
      break;
    }
    if (header_length == 0) {
      continue;
    }
    if (header[0] != '@') {
      throw UnsupportedParallelFastq(
          "FASTQ worker did not find a record header");
    }

    const char *single_sequence = nullptr;
    size_t single_length = 0;
    size_t sequence_length = 0;
    unsigned sequence_lines = 0;
    joined_sequence.clear();

    while (true) {
      const char *line = nullptr;
      size_t length = 0;
      if (!NextFastqLine(input_data, &position, &line, &length)) {
        throw UnsupportedParallelFastq("FASTQ record has no plus line");
      }
      if (length != 0 && line[0] == '+') {
        break;
      }
      if (sequence_length >
          std::numeric_limits<size_t>::max() - length) {
        throw UnsupportedParallelFastq("FASTQ sequence is too long");
      }
      if (sequence_lines == 0) {
        single_sequence = line;
        single_length = length;
      } else {
        if (sequence_lines == 1) {
          joined_sequence.assign(single_sequence, single_length);
        }
        joined_sequence.append(line, length);
      }
      ++sequence_lines;
      sequence_length += length;
    }

    size_t quality_length = 0;
    while (quality_length < sequence_length) {
      const char *quality = nullptr;
      size_t length = 0;
      if (!NextFastqLine(input_data, &position, &quality, &length)) {
        throw UnsupportedParallelFastq(
            "FASTQ record has a truncated quality field");
      }
      if (quality_length >
          std::numeric_limits<size_t>::max() - length) {
        throw UnsupportedParallelFastq("FASTQ quality is too long");
      }
      quality_length += length;
      if (quality_length > sequence_length) {
        throw UnsupportedParallelFastq(
            "FASTQ sequence and quality lengths differ");
      }
    }

    const char *sequence = nullptr;
    if (sequence_lines <= 1) {
      sequence = single_sequence == nullptr ? "" : single_sequence;
    } else {
      sequence = joined_sequence.data();
    }
    AppendPackedSequence(sequence, sequence_length, dna_map, &output);
  }
  return output;
}

size_t PackedRecordBytes(uint32_t packed_length) {
  return sizeof(uint32_t) *
         (DivCeiling<size_t>(packed_length, SeqPackage::kBasesPerWord) + 1u);
}

EncodedFastqBatch EncodeFastqBatch(RawFastqBatch input) {
  EncodedFastqBatch left = EncodeFastqData(input.data, input.ordinal);
  if (!input.paired) {
    return left;
  }

  EncodedFastqBatch right = EncodeFastqData(input.mate_data, input.ordinal);
  if (left.packed_lengths.size() != right.packed_lengths.size()) {
    throw UnsupportedParallelFastq(
        "paired FASTQ batches contain different record counts");
  }

  EncodedFastqBatch interleaved;
  interleaved.ordinal = input.ordinal;
  interleaved.data.reserve(left.data.size() + right.data.size());
  interleaved.packed_lengths.reserve(left.packed_lengths.size() * 2u);
  size_t left_offset = 0;
  size_t right_offset = 0;
  for (size_t i = 0; i < left.packed_lengths.size(); ++i) {
    const uint32_t left_length = left.packed_lengths[i];
    const uint32_t right_length = right.packed_lengths[i];
    const size_t left_bytes = PackedRecordBytes(left_length);
    const size_t right_bytes = PackedRecordBytes(right_length);
    const size_t destination = interleaved.data.size();
    interleaved.data.resize(destination + left_bytes + right_bytes);
    std::memcpy(interleaved.data.data() + destination,
                left.data.data() + left_offset, left_bytes);
    std::memcpy(interleaved.data.data() + destination + left_bytes,
                right.data.data() + right_offset, right_bytes);
    interleaved.packed_lengths.push_back(left_length);
    interleaved.packed_lengths.push_back(right_length);
    left_offset += left_bytes;
    right_offset += right_bytes;
  }
  if (left_offset != left.data.size() || right_offset != right.data.size()) {
    throw std::runtime_error("invalid packed paired FASTQ batch");
  }
  return interleaved;
}

/**
 * Bounded producer/worker/ordered-commit pipeline.
 *
 * The caller is both the decompressed-byte producer and ordered writer.  It
 * never permits more than a fixed byte budget of raw batches in flight.
 * Worker completion order is deliberately decoupled from semantic order by
 * `ordinal`, which makes the generated binary library byte-for-byte stable.
 */
class ParallelFastqPacker {
 public:
  ParallelFastqPacker(unsigned worker_count, bool paired,
                      PackedBinaryWriter *writer, BuiltLibrary *result)
      : worker_count_(std::max(1u, worker_count)),
        max_in_flight_(std::max<unsigned>(
            2u, std::min<unsigned>(
                    static_cast<unsigned>(
                        kFastqPipelineByteBudget /
                        (kFastqBatchBytes * (paired ? 2u : 1u))),
                    std::max(worker_count_, worker_count_ * 2u)))),
        writer_(writer),
        result_(result) {
    workers_.reserve(worker_count_);
    for (unsigned i = 0; i < worker_count_; ++i) {
      workers_.emplace_back(&ParallelFastqPacker::WorkerLoop, this);
    }
  }

  ~ParallelFastqPacker() { StopAndJoin(); }

  void Run(FastqRecordBatchReader *reader) {
    uint64_t next_submit = 0;
    uint64_t next_commit = 0;
    unsigned in_flight = 0;
    size_t bytes_in_flight = 0;
    std::deque<size_t> submitted_bytes;
    try {
      std::vector<char> bytes;
      while (reader->Next(&bytes)) {
        const size_t batch_bytes = bytes.size();
        while (in_flight >= max_in_flight_ ||
               (in_flight != 0 &&
                (bytes_in_flight >= kFastqPipelineByteBudget ||
                 batch_bytes >
                     kFastqPipelineByteBudget - bytes_in_flight))) {
          Commit(next_commit++);
          --in_flight;
          bytes_in_flight -= submitted_bytes.front();
          submitted_bytes.pop_front();
        }
        RawFastqBatch batch;
        batch.ordinal = next_submit++;
        batch.data.swap(bytes);
        Submit(std::move(batch));
        ++in_flight;
        bytes_in_flight += batch_bytes;
        submitted_bytes.push_back(batch_bytes);
      }
      while (in_flight != 0) {
        Commit(next_commit++);
        --in_flight;
        bytes_in_flight -= submitted_bytes.front();
        submitted_bytes.pop_front();
      }
      StopAndJoin();
    } catch (...) {
      StopAndJoin();
      throw;
    }
  }

  void RunPaired(FastqRecordBatchReader *left_reader,
                 FastqRecordBatchReader *right_reader) {
    uint64_t next_submit = 0;
    uint64_t next_commit = 0;
    unsigned in_flight = 0;
    size_t bytes_in_flight = 0;
    std::deque<size_t> submitted_bytes;
    try {
      while (true) {
        while (in_flight >= max_in_flight_) {
          Commit(next_commit++);
          --in_flight;
          bytes_in_flight -= submitted_bytes.front();
          submitted_bytes.pop_front();
        }

        RawFastqBatch batch;
        uint64_t left_records = 0;
        const bool have_left =
            left_reader->Next(&batch.data, &left_records);
        if (!have_left) {
          uint64_t right_records = 0;
          std::vector<char> extra_right;
          if (right_reader->Next(&extra_right, &right_records, 1)) {
            throw UnsupportedParallelFastq(
                "paired FASTQ files contain different record counts");
          }
          break;
        }

        uint64_t right_records = 0;
        const bool have_right = right_reader->Next(
            &batch.mate_data, &right_records, left_records);
        if (!have_right || left_records != right_records) {
          throw UnsupportedParallelFastq(
              "paired FASTQ files contain different record counts");
        }
        const size_t batch_bytes =
            batch.data.size() + batch.mate_data.size();
        while (in_flight != 0 &&
               (bytes_in_flight >= kFastqPipelineByteBudget ||
                batch_bytes >
                    kFastqPipelineByteBudget - bytes_in_flight)) {
          Commit(next_commit++);
          --in_flight;
          bytes_in_flight -= submitted_bytes.front();
          submitted_bytes.pop_front();
        }
        batch.ordinal = next_submit++;
        batch.paired = true;
        Submit(std::move(batch));
        ++in_flight;
        bytes_in_flight += batch_bytes;
        submitted_bytes.push_back(batch_bytes);
      }

      while (in_flight != 0) {
        Commit(next_commit++);
        --in_flight;
        bytes_in_flight -= submitted_bytes.front();
        submitted_bytes.pop_front();
      }
      StopAndJoin();
    } catch (...) {
      StopAndJoin();
      throw;
    }
  }

 private:
  void Submit(RawFastqBatch batch) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      tasks_.emplace_back(std::move(batch));
    }
    task_ready_.notify_one();
  }

  void Commit(uint64_t ordinal) {
    EncodedFastqBatch batch;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      result_ready_.wait(lock, [&] {
        return completed_.find(ordinal) != completed_.end();
      });
      auto found = completed_.find(ordinal);
      batch = std::move(found->second);
      completed_.erase(found);
    }
    if (batch.error) {
      std::rethrow_exception(batch.error);
    }
    writer_->AppendEncodedBatch(batch.data, batch.packed_lengths, result_);
  }

  void WorkerLoop() {
    // Workers may be created by an OpenMP library task.  In that case the
    // child would otherwise inherit one OpenMP place and all parser workers
    // could silently contend for the same CPU.
    ResetThreadAffinityToStartupMask();
    while (true) {
      RawFastqBatch task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        task_ready_.wait(lock, [&] { return stopping_ || !tasks_.empty(); });
        if (tasks_.empty()) {
          if (stopping_) {
            return;
          }
          continue;
        }
        task = std::move(tasks_.front());
        tasks_.pop_front();
      }

      EncodedFastqBatch output;
      const uint64_t ordinal = task.ordinal;
      try {
        output = EncodeFastqBatch(std::move(task));
      } catch (...) {
        output.ordinal = ordinal;
        output.error = std::current_exception();
      }
      {
        std::lock_guard<std::mutex> lock(mutex_);
        completed_.emplace(output.ordinal, std::move(output));
      }
      result_ready_.notify_all();
    }
  }

  void StopAndJoin() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (joined_) {
        return;
      }
      stopping_ = true;
      tasks_.clear();
    }
    task_ready_.notify_all();
    for (std::thread &worker : workers_) {
      if (worker.joinable()) {
        worker.join();
      }
    }
    joined_ = true;
  }

  unsigned worker_count_;
  unsigned max_in_flight_;
  PackedBinaryWriter *writer_;
  BuiltLibrary *result_;
  std::deque<RawFastqBatch> tasks_;
  std::map<uint64_t, EncodedFastqBatch> completed_;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable task_ready_;
  std::condition_variable result_ready_;
  bool stopping_{false};
  bool joined_{false};
};

struct FastqPipelinePlan {
  unsigned gzip_threads{0};
  unsigned parser_workers{0};
};

FastqPipelinePlan PlanFastqPipeline(unsigned auxiliary_threads,
                                    bool parallel_gzip,
                                    unsigned input_streams) {
  FastqPipelinePlan plan;
  input_streams = std::max(1u, input_streams);
  if (auxiliary_threads <= input_streams) {
    plan.gzip_threads = auxiliary_threads;
    return plan;
  }

  // The in-flight queue, not the number of input streams, enforces the 512 MiB
  // raw-byte cap.  Keep the compute split independent of SE/PE layout; paired
  // batches simply admit half as many simultaneous tasks.
  const unsigned max_parser_workers = kMaxFastqParserWorkers;
  const unsigned parallel_decoder_minimum = input_streams * 2u;
  if (parallel_gzip &&
      auxiliary_threads > parallel_decoder_minimum) {
    plan.parser_workers =
        std::min(max_parser_workers, auxiliary_threads / 2u);
    if (auxiliary_threads - plan.parser_workers <
        parallel_decoder_minimum) {
      plan.parser_workers =
          auxiliary_threads - parallel_decoder_minimum;
    }
    plan.gzip_threads = auxiliary_threads - plan.parser_workers;
  } else {
    // zlib has one bounded producer per input stream.  Give the remaining
    // fair-share slots to parsing/packing rather than leaving them idle.
    plan.gzip_threads = input_streams;
    plan.parser_workers = std::min(
        max_parser_workers, auxiliary_threads - input_streams);
  }
  return plan;
}

void PackLibrary(const InputLibrary &lib, std::ostream *output,
                 AuxiliaryThreadBudget *auxiliary_budget,
                 bool intra_file_parallel,
                 bool allow_parallel_fastq,
                 BuiltLibrary *result) {
  // Independent libraries already expose coarse-grained parallelism.  Keep
  // every gzip stream and packed-output chunk bounded instead of materializing
  // either a whole input or an intermediate SeqPackage.
  // Parallel gzip recovery does extra block-finding and window work.  It is a
  // win when fine-grained per-file parallelism dominates, but it is wasteful
  // when independent libraries already expose the wider scheduling level.
  // In the latter case request one efficient zlib producer per gzip stream;
  // this is a topology choice and does not depend on file contents or a
  // machine-specific threshold.
  const bool parallel_fastq_candidate =
      allow_parallel_fastq && intra_file_parallel &&
      lib.file_name1 != "-" &&
      (lib.type != "pe" || lib.file_name2 != "-");
  const unsigned requested_auxiliary =
      parallel_fastq_candidate || intra_file_parallel
          ? std::numeric_limits<unsigned>::max()
          : (lib.type == "pe" ? 2u : 1u);
  AuxiliaryThreadLease auxiliary_lease(auxiliary_budget,
                                       requested_auxiliary);

  if (parallel_fastq_candidate) {
    const unsigned input_streams = InputStreamCount(lib);
    const FastqPipelinePlan pipeline = PlanFastqPipeline(
        auxiliary_lease.count(), intra_file_parallel, input_streams);
    if (pipeline.parser_workers != 0) {
      xinfo("Parallel FASTQ library {s}: {} decoder and {} parser/packer "
            "workers across {} input stream(s)\n",
            lib.metadata.c_str(), pipeline.gzip_threads,
            pipeline.parser_workers, input_streams);
      PackedBinaryWriter writer(output);
      ParallelFastqPacker packer(pipeline.parser_workers, lib.type == "pe",
                                 &writer, result);
      const LibraryIoPlan io_plan =
          SplitDecoderThreads(lib, pipeline.gzip_threads);
      FastqRecordBatchReader reader1(lib.file_name1,
                                     io_plan.file1_threads);
      if (lib.type == "pe") {
        FastqRecordBatchReader reader2(lib.file_name2,
                                       io_plan.file2_threads);
        packer.RunPaired(&reader1, &reader2);
      } else {
        packer.Run(&reader1);
      }
      writer.Finish(result);
      if ((lib.type == "pe" || lib.type == "interleaved") &&
          result->num_reads % 2 != 0) {
        throw std::runtime_error(
            "paired library has an odd number of reads: " + lib.metadata);
      }
      return;
    }
  }

  const LibraryIoPlan io_plan =
      SplitDecoderThreads(lib, auxiliary_lease.count());
  PackedBinaryWriter writer(output);
  if (lib.type == "pe") {
    FastxReader mate1(lib.file_name1, false, io_plan.file1_threads);
    FastxReader mate2(lib.file_name2, false, io_plan.file2_threads);
    while (true) {
      kseq_t *read1 = mate1.ReadNext();
      kseq_t *read2 = mate2.ReadNext();
      if (read1 == nullptr || read2 == nullptr) {
        break;
      }
      writer.Append(read1->seq.s, read1->seq.l, result);
      writer.Append(read2->seq.s, read2->seq.l, result);
    }
  } else {
    FastxReader reader(lib.file_name1, false, io_plan.file1_threads);
    while (kseq_t *read = reader.ReadNext()) {
      writer.Append(read->seq.s, read->seq.l, result);
    }
  }
  writer.Finish(result);

  if ((lib.type == "pe" || lib.type == "interleaved") &&
      result->num_reads % 2 != 0) {
    throw std::runtime_error("paired library has an odd number of reads: " +
                             lib.metadata);
  }
}

void BuildOneLibrary(const InputLibrary &lib, const std::string &part_path,
                     AuxiliaryThreadBudget *auxiliary_budget,
                     bool intra_file_parallel,
                     BuiltLibrary *result) {
  result->part_path = part_path;
  {
    std::ofstream part_file(part_path,
                            std::ofstream::binary | std::ofstream::out |
                                std::ofstream::trunc);
    if (!part_file.is_open()) {
      throw std::runtime_error("cannot create temporary binary library " +
                               part_path);
    }
    try {
      PackLibrary(lib, &part_file, auxiliary_budget, intra_file_parallel, true,
                  result);
      return;
    } catch (const UnsupportedParallelFastq &e) {
      xwarn("Parallel FASTQ parser fallback for {s}: {s}\n",
            lib.file_name1.c_str(), e.what());
    }
  }

  // The fast path writes only to this private part file.  Reopening with
  // truncation makes fallback exact even when an unusual layout is detected
  // late in the input.
  *result = BuiltLibrary{};
  result->part_path = part_path;
  std::ofstream part_file(part_path,
                          std::ofstream::binary | std::ofstream::out |
                              std::ofstream::trunc);
  if (!part_file.is_open()) {
    throw std::runtime_error("cannot recreate temporary binary library " +
                             part_path);
  }
  PackLibrary(lib, &part_file, auxiliary_budget, intra_file_parallel, false,
              result);
}

void AppendFile(const std::string &path, std::ostream *output,
                std::vector<char> *buffer) {
  std::ifstream input(path, std::ifstream::binary | std::ifstream::in);
  if (!input.is_open()) {
    throw std::runtime_error("cannot reopen temporary binary library " + path);
  }
  while (input) {
    input.read(buffer->data(), buffer->size());
    const std::streamsize size = input.gcount();
    if (size != 0) {
      output->write(buffer->data(), size);
    }
  }
  if (!input.eof() || !*output) {
    throw std::runtime_error("failed merging temporary binary library " +
                             path);
  }
}

}  // namespace

void SequenceLibCollection::Build(const std::string &lib_file,
                                  const std::string &out_prefix,
                                  unsigned num_threads) {
  std::ifstream lib_config(lib_file);

  if (!lib_config.is_open()) {
    xfatal("File to open read_lib file: {}\n", lib_file.c_str());
  }

  std::vector<InputLibrary> input_libs;
  std::string metadata;
  while (std::getline(lib_config, metadata)) {
    InputLibrary lib;
    lib.metadata = metadata;
    if (!(lib_config >> lib.type)) {
      xfatal("Missing read library type after: {s}\n", metadata.c_str());
    }
    bool paths_valid = true;
    if (lib.type == "pe") {
      paths_valid = ReadLibraryPath(lib_config, &lib.file_name1) &&
                    ReadLibraryPath(lib_config, &lib.file_name2);
    } else if (lib.type == "se" || lib.type == "interleaved") {
      paths_valid = ReadLibraryPath(lib_config, &lib.file_name1);
    } else {
      xerr("Cannot identify read library type {}\n", lib.type.c_str());
      xfatal("Valid types: pe, se, interleaved\n");
    }
    if (!paths_valid) {
      xfatal("Malformed read library entry: {s}\n", metadata.c_str());
    }
    input_libs.emplace_back(std::move(lib));
    std::getline(lib_config, metadata);  // eliminate the "\n"
  }

  if (input_libs.empty()) {
    xfatal("No read libraries found in {s}\n", lib_file.c_str());
  }

  const unsigned worker_count = std::max(
      1u, std::min(num_threads, static_cast<unsigned>(input_libs.size())));
  const unsigned auxiliary_thread_count =
      num_threads > worker_count ? num_threads - worker_count : 0;
  AuxiliaryThreadBudget auxiliary_budget(auxiliary_thread_count, worker_count);
  uint64_t active_input_streams = 0;
  for (const InputLibrary &lib : input_libs) {
    active_input_streams += InputStreamCount(lib);
  }
  // Prefer the lower-work zlib path when independent files already expose
  // more parallelism than each file could profitably receive internally.
  // This topology-only comparison scales with `-t`, library count, and paired
  // inputs without a machine- or dataset-specific file-count threshold.
  const bool fine_parallelism_dominates =
      active_input_streams != 0 &&
      auxiliary_thread_count / active_input_streams > active_input_streams;
  const bool intra_file_parallel = fine_parallelism_dominates;
  xinfo("Building {} read libraries with {} library workers and {} "
        "globally bounded decoder/parser slots across {} input streams; "
        "backend: {s}\n",
        input_libs.size(), worker_count, auxiliary_budget.total(),
        active_input_streams,
        intra_file_parallel ? "parallel intra-file gzip"
                            : "file-parallel bounded zlib");
  std::vector<BuiltLibrary> built(input_libs.size());

  std::vector<size_t> task_order(input_libs.size());
  std::vector<uint64_t> task_work_bytes(input_libs.size());
  std::iota(task_order.begin(), task_order.end(), size_t{0});
  for (size_t i = 0; i < input_libs.size(); ++i) {
    task_work_bytes[i] = InputWorkBytes(input_libs[i]);
  }
  std::stable_sort(task_order.begin(), task_order.end(),
                   [&](size_t lhs, size_t rhs) {
                     return task_work_bytes[lhs] > task_work_bytes[rhs];
                   });

  if (worker_count == 1 && input_libs.size() == 1) {
    // A private part file makes the speculative parallel FASTQ parser safely
    // restartable.  Successful output is renamed in place, so the normal path
    // still performs no second read/copy pass.
    const std::string part_path = out_prefix + ".bin.part.0";
    try {
      BuildOneLibrary(input_libs[0], part_path, &auxiliary_budget,
                      intra_file_parallel, &built[0]);
      if (std::rename(part_path.c_str(), (out_prefix + ".bin").c_str()) != 0) {
        throw std::runtime_error("cannot publish temporary binary library " +
                                 part_path);
      }
    } catch (const std::exception &e) {
      std::remove(part_path.c_str());
      xfatal("Failed to build read library: {s}\n", e.what());
    }
  } else if (worker_count == 1) {
    // A one-thread, multi-library invocation has no auxiliary slots for the
    // parallel pipeline.  Preserve the direct serial writer and avoid part
    // file merge traffic.
    std::ofstream bin_file(out_prefix + ".bin",
                           std::ofstream::binary | std::ofstream::out);
    if (!bin_file.is_open()) {
      xfatal("Cannot create binary read library {s}.bin\n",
             out_prefix.c_str());
    }
    try {
      for (size_t i = 0; i < input_libs.size(); ++i) {
        PackLibrary(input_libs[i], &bin_file, &auxiliary_budget,
                    intra_file_parallel, false, &built[i]);
      }
    } catch (const std::exception &e) {
      bin_file.close();
      xfatal("Failed to build read library: {s}\n", e.what());
    }
    bin_file.close();
  } else {
#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (int64_t task = 0; task < static_cast<int64_t>(task_order.size());
         ++task) {
      const size_t i = task_order[task];
      const std::string part_path =
          out_prefix + ".bin.part." + std::to_string(i);
      try {
        BuildOneLibrary(input_libs[i], part_path, &auxiliary_budget,
                        intra_file_parallel, &built[i]);
      } catch (const std::exception &e) {
        built[i].error = e.what();
      }
    }

    for (const auto &result : built) {
      if (!result.error.empty()) {
        for (const auto &part : built) {
          if (!part.part_path.empty()) {
            std::remove(part.part_path.c_str());
          }
        }
        xfatal("Failed to build read library: {s}\n", result.error.c_str());
      }
    }

    std::ofstream bin_file(out_prefix + ".bin",
                           std::ofstream::binary | std::ofstream::out);
    if (!bin_file.is_open()) {
      xfatal("Cannot create binary read library {s}.bin\n",
             out_prefix.c_str());
    }
    std::vector<char> copy_buffer(8u << 20u);
    try {
      for (const auto &result : built) {
        AppendFile(result.part_path, &bin_file, &copy_buffer);
        std::remove(result.part_path.c_str());
      }
    } catch (const std::exception &e) {
      for (const auto &result : built) {
        std::remove(result.part_path.c_str());
      }
      xfatal("Failed to merge read libraries: {s}\n", e.what());
    }
    bin_file.close();
  }

  int64_t total_reads = 0;
  int64_t total_bases = 0;
  unsigned global_max_read_len = 0;
  std::vector<SequenceLib> libs;
  for (size_t i = 0; i < input_libs.size(); ++i) {
    const int64_t begin_index = total_reads;
    total_reads += built[i].num_reads;
    total_bases += built[i].num_bases;
    global_max_read_len =
        std::max(global_max_read_len, built[i].max_read_len);
    xinfo("Lib {} ({s}): {s}, {} reads, {} max length\n", i,
          input_libs[i].metadata.c_str(), input_libs[i].type.c_str(),
          built[i].num_reads, built[i].max_read_len);
    libs.emplace_back(nullptr, begin_index, total_reads,
                      built[i].max_read_len, input_libs[i].type != "se",
                      input_libs[i].metadata);
  }

  std::ofstream lib_info_file(out_prefix + ".lib_info");
  lib_info_file << total_bases << ' ' << total_reads << '\n';

  for (auto &lib : libs) {
    lib.DumpMetadata(lib_info_file);
  }
  lib_info_file.close();

  std::vector<PackedReadChunk> merged_chunks;
  uint64_t word_base = 0;
  uint64_t read_base = 0;
  bool have_complete_chunks = true;
  for (const BuiltLibrary &result : built) {
    if (result.num_reads != 0 && result.chunks.empty()) {
      have_complete_chunks = false;
    }
    for (const PackedReadChunk &local : result.chunks) {
      merged_chunks.push_back(PackedReadChunk{
          word_base + local.word_begin, word_base + local.word_end,
          read_base + local.read_begin, read_base + local.read_end,
          local.num_bases, local.max_read_len});
    }
    word_base += result.num_words;
    read_base += static_cast<uint64_t>(result.num_reads);
  }
  if (have_complete_chunks && !merged_chunks.empty() &&
      !PublishPackedReadChunkIndex(
          out_prefix + ".bin", static_cast<uint64_t>(total_reads),
          static_cast<uint64_t>(total_bases), global_max_read_len,
          merged_chunks)) {
    xwarn("Could not publish packed-read chunk index for {s}.bin\n",
          out_prefix.c_str());
  }
}

namespace {

/**
 * Map and bulk-load a binary read library in bounded parallel chunks.
 *
 * The serial reader decodes ~33 M records one at a time (bit-repacking each
 * read, and bit-reversing it for the sorting stages).  The previous bulk path
 * removed that CPU bottleneck but simultaneously materialized the complete
 * file, an 8/12-byte per-read locator/length index, and the final packed
 * package.  At hundreds of millions of reads those temporary objects were
 * much larger than the useful 2-bit data.
 *
 * A private file mapping keeps I/O in the kernel page cache.  SequencePackage
 * validates and indexes it in one linear pass, then repacks byte-bounded
 * chunks in parallel and drops each completed mapped range.  Resident source
 * memory is therefore bounded by active chunks rather than input size.
 */
bool TryBulkBinaryLoad(const std::string &bin_path, int64_t total_bases,
                       int64_t num_reads, unsigned max_read_len,
                       SeqPackage *pkg, bool reverse_seq) {
  using TWord = SeqPackage::TWord;
  if (num_reads <= 0 || total_bases <= 0 || max_read_len == 0) {
    return false;
  }

  const int fd = open(bin_path.c_str(), O_RDONLY);
  if (fd < 0) {
    return false;
  }
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size <= 0 ||
      st.st_size % static_cast<off_t>(sizeof(TWord)) != 0) {
    close(fd);
    return false;
  }
  const size_t total_bytes = static_cast<size_t>(st.st_size);
  const size_t total_words = total_bytes / sizeof(TWord);
  void *mapping =
      mmap(nullptr, total_bytes, PROT_READ, MAP_PRIVATE, fd, 0);
  if (mapping == MAP_FAILED) {
    close(fd);
    return false;
  }
#if defined(MADV_SEQUENTIAL)
  madvise(mapping, total_bytes, MADV_SEQUENTIAL);
#endif
  const bool loaded = pkg->AssignMappedBinaryRecords(
      static_cast<const TWord *>(mapping), total_words,
      static_cast<size_t>(num_reads), static_cast<uint64_t>(total_bases),
      max_read_len, reverse_seq, omp_get_max_threads());
  munmap(mapping, total_bytes);
  close(fd);
  return loaded;
}

}  // namespace

void SequenceLibCollection::ReadMetadata(SeqPackage *data_holder) {
  std::ifstream lib_info_file(path_ + ".lib_info");
  int64_t total_bases, num_reads;
  bool is_paired;
  std::string metadata;
  libs_.clear();
  lib_info_file >> total_bases >> num_reads;
  std::getline(lib_info_file, metadata);  // eliminate the "\n"

  while (std::getline(lib_info_file, metadata)) {
    int64_t start, end;
    int max_read_len;
    lib_info_file >> start >> end >> max_read_len >> is_paired;
    libs_.emplace_back(data_holder, start, end, max_read_len, is_paired,
                       metadata);
    std::getline(lib_info_file, metadata);  // eliminate the "\n"
  }
}

void SequenceLibCollection::Read(SeqPackage *pkg, bool reverse_seq) {
  ReadMetadata(pkg);
  const SizeInfo size_info = GetSizeInfo();

  pkg->Clear();
  if (TryBulkBinaryLoad(path_ + ".bin", size_info.num_bases,
                        size_info.num_reads, size_info.max_read_len, pkg,
                        reverse_seq)) {
    xinfo("After reading, sizeof seq_package: {}\n", pkg->size_in_byte());
    return;
  }

  pkg->ReserveSequences(size_info.num_reads);
  pkg->ReserveBases(size_info.num_bases);
  BinaryReader reader(path_ + ".bin");

  xinfo("Before reading, sizeof seq_package: {}\n", pkg->size_in_byte());
  reader.ReadAll(pkg, reverse_seq);
  xinfo("After reading, sizeof seq_package: {}\n", pkg->size_in_byte());
}

std::pair<int64_t, int64_t> SequenceLibCollection::GetSize() const {
  std::ifstream lib_info_file(path_ + ".lib_info");
  int64_t total_bases, num_reads;
  lib_info_file >> total_bases >> num_reads;
  return {total_bases, num_reads};
}

SequenceLibCollection::SizeInfo SequenceLibCollection::GetSizeInfo() const {
  std::ifstream lib_info_file(path_ + ".lib_info");
  SizeInfo info{0, 0, 0};
  std::string line;
  lib_info_file >> info.num_bases >> info.num_reads;
  std::getline(lib_info_file, line);
  while (std::getline(lib_info_file, line)) {
    int64_t start, end;
    unsigned max_read_len;
    bool paired;
    lib_info_file >> start >> end >> max_read_len >> paired;
    info.max_read_len = std::max(info.max_read_len, max_read_len);
    std::getline(lib_info_file, line);
  }
  return info;
}
