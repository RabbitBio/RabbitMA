//
// Created by vout on 5/11/19.
//

#ifndef MEGAHIT_PAIR_END_FASTX_READER_H
#define MEGAHIT_PAIR_END_FASTX_READER_H

#include <memory>
#include <thread>
#include "fastx_reader.h"
#include "utils/startup_affinity.h"

class PairedFastxReader : public BaseSequenceReader {
 public:
  PairedFastxReader(const std::string &file1, const std::string &file2,
                    bool allow_whole_gzip = true,
                    unsigned async_streams = 0) {
    if (!allow_whole_gzip) {
      // Library-level workers parse mate 1 while a bounded producer inflates
      // mate 2.  This exposes both gzip streams without materializing either
      // file and still preserves paired record order in the consumer.
      readers_[0].reset(new FastxReader(file1, false, async_streams >= 2));
      readers_[1].reset(new FastxReader(file2, false, async_streams >= 1));
      return;
    }
    // Each constructor may decompress a whole gzip file; do the two mates
    // concurrently.
    std::exception_ptr err;
    std::thread t([&] {
      ResetThreadAffinityToStartupMask();
      try {
        readers_[1].reset(new FastxReader(file2, true));
      } catch (...) {
        err = std::current_exception();
      }
    });
    try {
      readers_[0].reset(new FastxReader(file1, true));
    } catch (...) {
      t.join();
      throw;
    }
    t.join();
    if (err) {
      std::rethrow_exception(err);
    }
  }

  int64_t Read(SeqPackage *pkg, int64_t max_num, int64_t max_num_bases,
               bool reverse) override;

 private:
  std::unique_ptr<FastxReader> readers_[2];
  bool trim_n_{true};
};

#endif  // MEGAHIT_PAIR_END_FASTX_READER_H
