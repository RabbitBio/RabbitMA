#ifndef MEGAHIT_KMSORT_SELECTOR_H
#define MEGAHIT_KMSORT_SELECTOR_H

#include <cstdint>
#include <functional>

std::function<void(uint32_t*, int64_t)> SelectSortingFunc(int words_per_substr,
                                                          int extra_words,
                                                          int ignored_low_bytes = 0,
                                                          int ignored_high_bytes = 2);

#endif
