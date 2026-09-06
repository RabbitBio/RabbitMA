#ifndef MEGAHIT_UTILS_STREAM_STORE_H
#define MEGAHIT_UTILS_STREAM_STORE_H

#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(__SSE2__) && !defined(MEGAHIT_DISABLE_STREAM_STORE_INTRINSICS)
#include <emmintrin.h>
#endif

inline bool StreamStoresAvailable() {
#if defined(__SSE2__) && !defined(MEGAHIT_DISABLE_STREAM_STORE_INTRINSICS)
  return true;
#else
  return false;
#endif
}

// Copy complete cache lines without allocating destination lines in cache.
// Partial boundary lines use ordinary stores, so adjacent worker ranges
// never mix a whole-line streaming store with another worker's bytes.
// Source and destination must not overlap, as with memcpy. A caller must
// issue StreamStoreFence before publishing the copied range to other threads.
inline void StreamStoreCopy(void *destination, const void *source, size_t bytes) {
#if defined(__SSE2__) && !defined(MEGAHIT_DISABLE_STREAM_STORE_INTRINSICS)
  auto *out = static_cast<uint8_t *>(destination);
  const auto *in = static_cast<const uint8_t *>(source);
  size_t prefix = (64u - (reinterpret_cast<uintptr_t>(out) & 63u)) & 63u;
  if (prefix > bytes) prefix = bytes;
  if (prefix != 0u) {
    std::memcpy(out, in, prefix);
    out += prefix;
    in += prefix;
    bytes -= prefix;
  }
  while (bytes >= 64u) {
    const __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in));
    const __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in + 16u));
    const __m128i c = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in + 32u));
    const __m128i d = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in + 48u));
    _mm_stream_si128(reinterpret_cast<__m128i *>(out), a);
    _mm_stream_si128(reinterpret_cast<__m128i *>(out + 16u), b);
    _mm_stream_si128(reinterpret_cast<__m128i *>(out + 32u), c);
    _mm_stream_si128(reinterpret_cast<__m128i *>(out + 48u), d);
    in += 64u;
    out += 64u;
    bytes -= 64u;
  }
  if (bytes != 0u) std::memcpy(out, in, bytes);
#else
  std::memcpy(destination, source, bytes);
#endif
}

inline void StreamStoreFence() {
#if defined(__SSE2__) && !defined(MEGAHIT_DISABLE_STREAM_STORE_INTRINSICS)
  _mm_sfence();
#endif
}

#endif
