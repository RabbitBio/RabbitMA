#include "sequence/io/rapidgzip_adapter.h"

#include <algorithm>
#include <climits>
#include <exception>
#include <memory>
#include <new>
#include <string>

#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>

struct MegahitRapidGzipHandle {
  using Reader = rapidgzip::ParallelGzipReader<rapidgzip::ChunkData>;

  explicit MegahitRapidGzipHandle(std::unique_ptr<Reader> input)
      : reader(std::move(input)) {}

  std::unique_ptr<Reader> reader;
};

MegahitRapidGzipHandle *MegahitRapidGzipOpen(const std::string &path,
                                             unsigned decoder_threads,
                                             size_t chunk_bytes,
                                             std::string *error) {
  if (error != nullptr) {
    error->clear();
  }
  try {
    const size_t parallelization = std::max(2u, decoder_threads);
    auto file_reader =
        std::make_unique<rapidgzip::StandardFileReader>(path);
    auto reader = std::make_unique<MegahitRapidGzipHandle::Reader>(
        std::move(file_reader), parallelization, chunk_bytes);

    // MEGAHIT consumes every stream exactly once and never seeks.  Retaining
    // the random-access window index would only increase peak memory with the
    // compressed input size.
    reader->setKeepIndex(false);
    reader->setCRC32Enabled(true);
    reader->setStatisticsEnabled(false);
    reader->setShowProfileOnDestruction(false);
    return new MegahitRapidGzipHandle(std::move(reader));
  } catch (const std::exception &e) {
    if (error != nullptr) {
      *error = e.what();
    }
  } catch (...) {
    if (error != nullptr) {
      *error = "unknown rapidgzip initialization failure";
    }
  }
  return nullptr;
}

int MegahitRapidGzipRead(MegahitRapidGzipHandle *handle, void *buffer,
                         unsigned bytes, std::string *error) {
  if (handle == nullptr || handle->reader == nullptr) {
    if (error != nullptr) {
      *error = "invalid rapidgzip stream";
    }
    return -1;
  }
  try {
    const size_t got = handle->reader->read(
        static_cast<char *>(buffer), static_cast<size_t>(bytes));
    return static_cast<int>(std::min<size_t>(got, INT_MAX));
  } catch (const std::exception &e) {
    if (error != nullptr) {
      *error = e.what();
    }
  } catch (...) {
    if (error != nullptr) {
      *error = "unknown rapidgzip decompression failure";
    }
  }
  return -1;
}

void MegahitRapidGzipClose(MegahitRapidGzipHandle *handle) { delete handle; }
