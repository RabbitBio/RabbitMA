//
// Created by vout on 6/24/19.
//

#ifndef MEGAHIT_CONTIG_WRITER_H
#define MEGAHIT_CONTIG_WRITER_H

#include <definitions.h>
#include <omp.h>
#include <atomic>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include "utils/utils.h"

class ContigWriter {
 public:
  explicit ContigWriter(const std::string file_name) : file_name_(file_name) {
    file_ = xfopen(file_name.c_str(), "w");
    const char *disable_buffering =
        std::getenv("MEGAHIT_DISABLE_BUFFERED_CONTIG_OUTPUT");
    buffered_output_ =
        disable_buffering == nullptr || std::strcmp(disable_buffering, "1") != 0;
    if (buffered_output_) {
      num_thread_buffers_ = static_cast<size_t>(omp_get_max_threads());
      if (num_thread_buffers_ == 0) {
        num_thread_buffers_ = 1;
      }
      thread_buffers_.reset(new ThreadBuffer[num_thread_buffers_]);
    }
  }

  ~ContigWriter() {
    int64_t n_contigs = n_contigs_.load(std::memory_order_relaxed);
    int64_t n_bases = n_bases_.load(std::memory_order_relaxed);
    if (buffered_output_) {
      Flush();
      for (size_t i = 0; i < num_thread_buffers_; ++i) {
        n_contigs += thread_buffers_[i].n_contigs;
        n_bases += thread_buffers_[i].n_bases;
      }
    }
    std::FILE *info_file = xfopen((file_name_ + ".info").c_str(), "w");
    pfprintf(info_file, "{} {}\n", n_contigs, n_bases);
    fclose(info_file);
    fclose(file_);
  }

  // Call only after all producer threads have left their parallel region.
  // This drains the final partial block from every thread while the output
  // stage is still being timed, without adding a lock to each record.
  void Flush() {
    if (buffered_output_) {
      for (size_t i = 0; i < num_thread_buffers_; ++i) {
        FlushThreadBuffer(thread_buffers_[i]);
      }
    }
    if (UNLIKELY(std::fflush(file_) != 0)) {
      xfatal("Cannot flush contigs to {s}\n", file_name_.c_str());
    }
  }

  void WriteContig(const std::string &ascii_contig, unsigned k_size,
                   long long id, int flag, double multi) {
    if (UNLIKELY(!buffered_output_)) {
      pfprintf(file_, ">k{}_{} flag={} multi={.4} len={}\n{s}\n", k_size, id,
               flag, multi, ascii_contig.length(), ascii_contig.c_str());
      n_contigs_.fetch_add(1, std::memory_order_relaxed);
      n_bases_.fetch_add(
          ascii_contig.length() + (flag & contig_flag::kLoop) ? 28 : 0,
          std::memory_order_relaxed);
      return;
    }

    ThreadBuffer &thread_buffer = GetThreadBuffer();
    char header[kHeaderBufferSize];
    int header_size = std::snprintf(
        header, sizeof(header),
        ">k%u_%lld flag=%d multi=%.4f len=%zu\n", k_size, id, flag, multi,
        ascii_contig.length());
    AppendHeader(thread_buffer.data, header, sizeof(header), header_size,
                 k_size, id, flag, multi, ascii_contig.length());
    thread_buffer.data.append(ascii_contig);
    thread_buffer.data.push_back('\n');
    ++thread_buffer.n_contigs;
    // Preserve the historical .info accounting expression byte-for-byte.
    thread_buffer.n_bases +=
        ascii_contig.length() + (flag & contig_flag::kLoop) ? 28 : 0;
    if (UNLIKELY(thread_buffer.data.size() >= kFlushSize)) {
      FlushThreadBuffer(thread_buffer);
    }
  }

  void WriteLocalContig(const std::string &ascii_contig,
                        int64_t origin_contig_id, int strand,
                        int64_t contig_id) {
    if (UNLIKELY(!buffered_output_)) {
      pfprintf(file_, ">lc_{}_strand_{}_id_{} flag=0 multi=1\n{s}\n",
               origin_contig_id, strand, contig_id, ascii_contig.c_str());
      n_contigs_.fetch_add(1, std::memory_order_relaxed);
      n_bases_.fetch_add(ascii_contig.length(), std::memory_order_relaxed);
      return;
    }

    ThreadBuffer &thread_buffer = GetThreadBuffer();
    char header[kHeaderBufferSize];
    int header_size = std::snprintf(
        header, sizeof(header),
        ">lc_%lld_strand_%d_id_%lld flag=0 multi=1\n",
        static_cast<long long>(origin_contig_id), strand,
        static_cast<long long>(contig_id));
    AppendLocalHeader(thread_buffer.data, header, sizeof(header), header_size,
                      origin_contig_id, strand, contig_id);
    thread_buffer.data.append(ascii_contig);
    thread_buffer.data.push_back('\n');
    ++thread_buffer.n_contigs;
    thread_buffer.n_bases += ascii_contig.length();
    if (UNLIKELY(thread_buffer.data.size() >= kFlushSize)) {
      FlushThreadBuffer(thread_buffer);
    }
  }

 private:
  static const size_t kFlushSize = 1u << 20u;
  static const size_t kHeaderBufferSize = 192;

  struct ThreadBuffer {
    std::string data;
    int64_t n_contigs{0};
    int64_t n_bases{0};
  };

  ThreadBuffer &GetThreadBuffer() {
    size_t thread_id = 0;
    if (omp_in_parallel()) {
      thread_id = static_cast<size_t>(omp_get_thread_num());
    }
    if (UNLIKELY(thread_id >= num_thread_buffers_)) {
      xfatal("Contig writer has {} thread buffers, but thread {} attempted "
             "to write\n",
             num_thread_buffers_, thread_id);
    }
    ThreadBuffer &thread_buffer = thread_buffers_[thread_id];
    if (UNLIKELY(thread_buffer.data.capacity() < kFlushSize)) {
      thread_buffer.data.reserve(kFlushSize);
    }
    return thread_buffer;
  }

  void FlushThreadBuffer(ThreadBuffer &thread_buffer) {
    if (thread_buffer.data.empty()) {
      return;
    }
    const size_t bytes = thread_buffer.data.size();
    if (UNLIKELY(std::fwrite(thread_buffer.data.data(), 1, bytes, file_) !=
                 bytes)) {
      xfatal("Cannot write contigs to {s}\n", file_name_.c_str());
    }
    thread_buffer.data.clear();
  }

  static void AppendHeader(std::string &output, const char *stack_header,
                           size_t stack_capacity, int header_size,
                           unsigned k_size, long long id, int flag,
                           double multi, size_t contig_length) {
    if (UNLIKELY(header_size < 0)) {
      xfatal("Cannot format contig header\n");
    }
    if (LIKELY(static_cast<size_t>(header_size) < stack_capacity)) {
      output.append(stack_header, static_cast<size_t>(header_size));
      return;
    }
    std::unique_ptr<char[]> large_header(new char[header_size + 1]);
    std::snprintf(large_header.get(), static_cast<size_t>(header_size) + 1,
                  ">k%u_%lld flag=%d multi=%.4f len=%zu\n", k_size, id, flag,
                  multi, contig_length);
    output.append(large_header.get(), static_cast<size_t>(header_size));
  }

  static void AppendLocalHeader(std::string &output, const char *stack_header,
                                size_t stack_capacity, int header_size,
                                int64_t origin_contig_id, int strand,
                                int64_t contig_id) {
    if (UNLIKELY(header_size < 0)) {
      xfatal("Cannot format local contig header\n");
    }
    if (LIKELY(static_cast<size_t>(header_size) < stack_capacity)) {
      output.append(stack_header, static_cast<size_t>(header_size));
      return;
    }
    std::unique_ptr<char[]> large_header(new char[header_size + 1]);
    std::snprintf(large_header.get(), static_cast<size_t>(header_size) + 1,
                  ">lc_%lld_strand_%d_id_%lld flag=0 multi=1\n",
                  static_cast<long long>(origin_contig_id), strand,
                  static_cast<long long>(contig_id));
    output.append(large_header.get(), static_cast<size_t>(header_size));
  }

  std::string file_name_;
  std::FILE *file_;
  bool buffered_output_{false};
  size_t num_thread_buffers_{0};
  std::unique_ptr<ThreadBuffer[]> thread_buffers_;
  std::atomic<int64_t> n_contigs_{0};
  std::atomic<int64_t> n_bases_{0};
};

#endif  // MEGAHIT_CONTIG_WRITER_H
