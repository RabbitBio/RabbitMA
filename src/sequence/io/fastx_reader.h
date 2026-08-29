//
// Created by vout on 4/28/19.
//

#ifndef MEGAHIT_FASTX_READER_H
#define MEGAHIT_FASTX_READER_H

#include <zlib.h>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "base_reader.h"
#include "definitions.h"
#include "kseq.h"

/**
 * Input stream for FASTX files: whole-file libdeflate where explicitly
 * allowed, bounded rapidgzip for intra-file parallelism, or a zlib stream
 * fallback for stdin, plain files and unsupported gzip variants.
 */
class AsyncGzipReader;
struct MegahitRapidGzipHandle;

struct MgzStream {
  gzFile gz{nullptr};
  MegahitRapidGzipHandle *rapid{nullptr};
  std::vector<char> mem;
  size_t pos{0};
  std::unique_ptr<AsyncGzipReader> async;
  std::string error;
  ~MgzStream();
};
typedef MgzStream *mgzFile;

mgzFile mgz_open(const std::string &file_name,
                 bool allow_whole_gzip = true,
                 unsigned gzip_threads = 0);
void mgz_close(mgzFile f);
int mgz_read(mgzFile f, void *buf, unsigned len);

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(mgzFile, mgz_read)
#endif

class FastxReader : public BaseSequenceReader {
 public:
  explicit FastxReader(const std::string &file_name,
                       bool allow_whole_gzip = true,
                       unsigned gzip_threads = 0);
  virtual ~FastxReader();
  virtual int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
                       bool reverse);
  static void TrimN(const char *s, int len, int *out_bpos, int *out_epos);

  kseq_t *ReadNext() {
    const int status = kseq_reader_ ? kseq_read(kseq_reader_) : -1;
    if (status >= 0) {
      return kseq_reader_;
    }
    if (status == -3 && fp_ != nullptr && !fp_->error.empty()) {
      throw std::runtime_error(fp_->error);
    }
    return nullptr;
  }

  mgzFile fp_{};
  kseq_t *kseq_reader_{};
  bool trim_n_{true};
};

#endif  // MEGAHIT_FASTX_READER_H
