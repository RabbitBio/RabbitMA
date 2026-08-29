//
// Created by vout on 4/28/19.
//

#include "fastx_reader.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cassert>
#include <array>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

#ifdef MEGAHIT_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif

#ifdef MEGAHIT_HAVE_RAPIDGZIP
#include "sequence/io/rapidgzip_adapter.h"
#endif

#include "utils/startup_affinity.h"

class AsyncGzipReader {
 private:
  static const size_t kChunkBytes = 256u << 10u;
  static const unsigned kChunkCount = 4;
  static const unsigned kNoChunk = kChunkCount;

  struct Chunk {
    std::unique_ptr<unsigned char[]> data;
    size_t size{0};
  };

 public:
  explicit AsyncGzipReader(gzFile gz) : gz_(gz) {
    // zlib's own input buffer and the output ring are deliberately bounded;
    // the latter is large enough to amortize synchronization but remains
    // tiny compared with a read library.
    gzbuffer(gz_, static_cast<unsigned>(kChunkBytes));
    for (unsigned i = 0; i < kChunkCount; ++i) {
      chunks_[i].data.reset(new unsigned char[kChunkBytes]);
      free_.push_back(i);
    }
    worker_ = std::thread(&AsyncGzipReader::Run, this);
  }

  ~AsyncGzipReader() { Stop(); }

  int Read(void *dst, unsigned len) {
    unsigned char *out = static_cast<unsigned char *>(dst);
    size_t copied = 0;
    std::unique_lock<std::mutex> lock(mutex_);
    while (copied < len) {
      if (current_ == kNoChunk) {
        ready_cv_.wait(lock, [&] { return !ready_.empty() || done_; });
        if (ready_.empty()) {
          return copied != 0 ? static_cast<int>(copied) : (failed_ ? -1 : 0);
        }
        current_ = ready_.front();
        ready_.pop_front();
        current_pos_ = 0;
      }

      Chunk &chunk = chunks_[current_];
      const size_t take =
          std::min<size_t>(len - copied, chunk.size - current_pos_);
      std::memcpy(out + copied, chunk.data.get() + current_pos_, take);
      copied += take;
      current_pos_ += take;
      if (current_pos_ == chunk.size) {
        chunk.size = 0;
        free_.push_back(current_);
        current_ = kNoChunk;
        free_cv_.notify_one();
      }
    }
    return static_cast<int>(copied);
  }

 private:
  void Run() {
    ResetThreadAffinityToStartupMask();
    while (true) {
      unsigned chunk_id;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        free_cv_.wait(lock, [&] { return stop_ || !free_.empty(); });
        if (stop_) {
          return;
        }
        chunk_id = free_.front();
        free_.pop_front();
      }

      const int got = gzread(gz_, chunks_[chunk_id].data.get(),
                             static_cast<unsigned>(kChunkBytes));
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (got > 0) {
          chunks_[chunk_id].size = static_cast<size_t>(got);
          ready_.push_back(chunk_id);
        } else {
          free_.push_back(chunk_id);
          failed_ = got < 0;
          done_ = true;
        }
      }
      ready_cv_.notify_one();
      if (got <= 0) {
        return;
      }
    }
  }

  void Stop() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    free_cv_.notify_all();
    ready_cv_.notify_all();
    if (worker_.joinable()) {
      worker_.join();
    }
  }

  gzFile gz_{};
  std::array<Chunk, kChunkCount> chunks_;
  std::deque<unsigned> free_;
  std::deque<unsigned> ready_;
  std::mutex mutex_;
  std::condition_variable free_cv_;
  std::condition_variable ready_cv_;
  std::thread worker_;
  unsigned current_{kNoChunk};
  size_t current_pos_{0};
  bool done_{false};
  bool failed_{false};
  bool stop_{false};
};

MgzStream::~MgzStream() = default;

namespace {

constexpr size_t kRapidGzipChunkBytes = size_t{4} << 20u;

bool IsRegularGzipFile(const std::string &path) {
  if (path == "-") {
    return false;
  }
  const int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    return false;
  }
  unsigned char magic[2]{};
  const ssize_t got = pread(fd, magic, sizeof(magic), 0);
  struct stat st;
  const bool regular = fstat(fd, &st) == 0 && S_ISREG(st.st_mode);
  close(fd);
  return regular && got == static_cast<ssize_t>(sizeof(magic)) &&
         magic[0] == 0x1f && magic[1] == 0x8b;
}

#ifdef MEGAHIT_HAVE_LIBDEFLATE

/**
 * Decompress a whole (possibly multi-member / BGZF) gzip file into memory.
 * Returns false on any surprise so the caller can fall back to zlib
 * streaming, which also transparently handles non-gzip files.
 */
bool LoadWholeGzip(const char *path, std::vector<char> *out) {
  const int fd = open(path, O_RDONLY);
  if (fd < 0) {
    return false;
  }
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size < 18) {
    close(fd);
    return false;
  }
  const size_t csize = static_cast<size_t>(st.st_size);
  std::unique_ptr<char[]> cbuf;
  try {
    cbuf.reset(new char[csize]);
  } catch (const std::bad_alloc &) {
    close(fd);
    return false;
  }
  size_t got_total = 0;
  while (got_total < csize) {
    const ssize_t got =
        pread(fd, cbuf.get() + got_total, csize - got_total, got_total);
    if (got <= 0) {
      close(fd);
      return false;
    }
    got_total += got;
  }
  close(fd);

  const unsigned char *u = reinterpret_cast<const unsigned char *>(cbuf.get());
  if (u[0] != 0x1f || u[1] != 0x8b) {
    return false;  // not gzip; zlib fallback reads it as plain text
  }

  libdeflate_decompressor *d = libdeflate_alloc_decompressor();
  if (d == nullptr) {
    return false;
  }
  try {
    out->resize(csize * 4 + (size_t{16} << 20u));
  } catch (const std::bad_alloc &) {
    libdeflate_free_decompressor(d);
    return false;
  }
  size_t in_pos = 0;
  size_t out_pos = 0;
  bool ok = true;
  while (in_pos < csize) {
    if (u[in_pos] != 0x1f) {
      // Tolerate zero padding after the last member; anything else fails.
      bool all_zero = true;
      for (size_t i = in_pos; i < csize; ++i) {
        if (cbuf[i] != 0) {
          all_zero = false;
          break;
        }
      }
      ok = all_zero;
      break;
    }
    size_t actual_in = 0;
    size_t actual_out = 0;
    const auto result = libdeflate_gzip_decompress_ex(
        d, cbuf.get() + in_pos, csize - in_pos, out->data() + out_pos,
        out->size() - out_pos, &actual_in, &actual_out);
    if (result == LIBDEFLATE_INSUFFICIENT_SPACE) {
      try {
        out->resize(out->size() * 2 + (size_t{64} << 20u));
      } catch (const std::bad_alloc &) {
        ok = false;
        break;
      }
      continue;
    }
    if (result != LIBDEFLATE_SUCCESS || actual_in == 0) {
      ok = false;
      break;
    }
    in_pos += actual_in;
    out_pos += actual_out;
  }
  libdeflate_free_decompressor(d);
  if (!ok) {
    out->clear();
    return false;
  }
  out->resize(out_pos);
  return true;
}

#endif  // MEGAHIT_HAVE_LIBDEFLATE

}  // namespace

mgzFile mgz_open(const std::string &file_name, bool allow_whole_gzip,
                 unsigned gzip_threads) {
  auto *s = new MgzStream();
#ifdef MEGAHIT_HAVE_LIBDEFLATE
  if (allow_whole_gzip && file_name != "-" &&
      std::getenv("MEGAHIT_DISABLE_LIBDEFLATE") == nullptr &&
      LoadWholeGzip(file_name.c_str(), &s->mem)) {
    return s;
  }
#endif
  s->mem.clear();
  s->mem.shrink_to_fit();
#ifdef MEGAHIT_HAVE_RAPIDGZIP
  if (gzip_threads >= 2 &&
      std::getenv("MEGAHIT_DISABLE_RAPIDGZIP") == nullptr &&
      IsRegularGzipFile(file_name)) {
    // rapidgzip creates its private worker pool lazily.  If this reader is
    // opened by an OpenMP worker, child threads would otherwise inherit that
    // worker's single-place affinity and silently serialize on one CPU.
    ResetThreadAffinityToStartupMask();
    s->rapid = MegahitRapidGzipOpen(file_name, gzip_threads,
                                    kRapidGzipChunkBytes, &s->error);
    if (s->rapid != nullptr) {
      s->error.clear();
      return s;
    }
    // Opening the optimized backend is deliberately fail-soft.  The zlib
    // fallback below preserves support for unusual gzip variants and older
    // platforms; actual decode errors after opening remain fatal.
    s->error.clear();
  }
#endif
  s->gz = file_name == "-" ? gzdopen(fileno(stdin), "r")
                           : gzopen(file_name.c_str(), "r");
  if (s->gz == nullptr) {
    delete s;
    return nullptr;
  }
  if (gzip_threads >= 1) {
    try {
      s->async.reset(new AsyncGzipReader(s->gz));
    } catch (const std::exception &) {
      // Bounded overlap is optional; retain the exact synchronous stream if
      // its small ring or worker thread cannot be created.
      s->async.reset();
    }
  }
  return s;
}

int mgz_read(mgzFile f, void *buf, unsigned len) {
#ifdef MEGAHIT_HAVE_RAPIDGZIP
  if (f->rapid != nullptr) {
    return MegahitRapidGzipRead(f->rapid, buf, len, &f->error);
  }
#endif
  if (f->async) {
    return f->async->Read(buf, len);
  }
  if (f->gz) {
    const int got = gzread(f->gz, buf, len);
    if (got < 0) {
      int error_number = Z_OK;
      const char *message = gzerror(f->gz, &error_number);
      f->error = message == nullptr ? "gzip decompression failed" : message;
    }
    return got;
  }
  const size_t n = std::min<size_t>(len, f->mem.size() - f->pos);
  std::memcpy(buf, f->mem.data() + f->pos, n);
  f->pos += n;
  return static_cast<int>(n);
}

void mgz_close(mgzFile f) {
  if (f == nullptr) {
    return;
  }
  f->async.reset();
#ifdef MEGAHIT_HAVE_RAPIDGZIP
  MegahitRapidGzipClose(f->rapid);
  f->rapid = nullptr;
#endif
  if (f->gz) {
    gzclose(f->gz);
  }
  delete f;
}

FastxReader::FastxReader(const std::string &file_name,
                         bool allow_whole_gzip, unsigned gzip_threads) {
  fp_ = mgz_open(file_name, allow_whole_gzip, gzip_threads);
  if (fp_ == nullptr) {
    throw std::invalid_argument("Cannot open file " + file_name);
  }
  kseq_reader_ = kseq_init(fp_);
  assert(kseq_reader_ != nullptr);
}

FastxReader::~FastxReader() {
  if (kseq_reader_) {
    kseq_destroy(kseq_reader_);
  }
  if (fp_) {
    mgz_close(fp_);
  }
}

int64_t FastxReader::Read(SeqPackage *pkg, int64_t max_num,
                          int64_t max_num_bases, bool reverse) {
  int64_t num_bases = 0;
  for (int64_t i = 0; i < max_num; ++i) {
    auto record = ReadNext();
    if (record) {
      int b = 0, e = record->seq.l;
      if (trim_n_) {
        TrimN(record->seq.s, record->seq.l, &b, &e);
      }

      if (reverse) {
        pkg->AppendReversedStringSequence(record->seq.s + b, e - b);
      } else {
        pkg->AppendStringSequence(record->seq.s + b, e - b);
      }

      num_bases += e - b;
      if (num_bases >= max_num_bases && i % 2 == 1) {
        return i + 1;
      }
    } else {
      return i;
    }
  }
  return max_num;
}

void FastxReader::TrimN(const char *s, int len, int *out_bpos, int *out_epos) {
  *out_bpos = *out_epos = len;
  int i;
  for (i = 0; i < len; ++i) {
    if (s[i] == 'N' || s[i] == 'n') {
      if (*out_bpos < len) {
        break;
      }
    } else {
      if (*out_bpos == len) {
        *out_bpos = i;
      }
    }
  }
  *out_epos = i;
}
