#ifndef MEGAHIT_RAPIDGZIP_ADAPTER_H
#define MEGAHIT_RAPIDGZIP_ADAPTER_H

#include <cstddef>
#include <string>

// Opaque boundary between MEGAHIT's C++11 sources and rapidgzip's C++17
// implementation.  All exceptions are caught inside the adapter.
struct MegahitRapidGzipHandle;

MegahitRapidGzipHandle *MegahitRapidGzipOpen(const std::string &path,
                                             unsigned decoder_threads,
                                             size_t chunk_bytes,
                                             std::string *error);

int MegahitRapidGzipRead(MegahitRapidGzipHandle *handle, void *buffer,
                         unsigned bytes, std::string *error);

void MegahitRapidGzipClose(MegahitRapidGzipHandle *handle);

#endif  // MEGAHIT_RAPIDGZIP_ADAPTER_H
