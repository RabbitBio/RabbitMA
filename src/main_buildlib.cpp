#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"

#include <cerrno>
#include <cstdlib>

void DisplayHelp(const char *program) {
  pfprintf(stderr,
           "Usage {s} <read_lib_file> <out_prefix> [num_threads] "
           "[anchor_len window_len]\n",
           program);
}

int main_build_lib(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  if (argc < 3) {
    DisplayHelp(argv[0]);
    exit(1);
  }
  unsigned num_threads = 1;
  if (argc >= 4) {
    char *end = nullptr;
    errno = 0;
    unsigned long value = std::strtoul(argv[3], &end, 10);
    if (errno != 0 || end == argv[3] || *end != '\0' || value == 0) {
      xfatal("Invalid buildlib thread count: {s}\n", argv[3]);
    }
    num_threads = static_cast<unsigned>(value);
  }
  unsigned anchor_len = 0;
  unsigned window_len = 0;
  if (argc == 5 || argc > 6) {
    xfatal("buildlib anchor_len and window_len must be specified together\n");
  }
  if (argc >= 6) {
    char *anchor_end = nullptr;
    char *window_end = nullptr;
    errno = 0;
    const unsigned long parsed_anchor =
        std::strtoul(argv[4], &anchor_end, 10);
    const unsigned long parsed_window =
        std::strtoul(argv[5], &window_end, 10);
    if (errno != 0 || anchor_end == argv[4] || *anchor_end != '\0' ||
        window_end == argv[5] || *window_end != '\0' ||
        parsed_anchor == 0 || parsed_anchor > 31 ||
        parsed_window < parsed_anchor) {
      xfatal("Invalid buildlib anchor/window geometry: {s}/{s}\n", argv[4],
             argv[5]);
    }
    anchor_len = static_cast<unsigned>(parsed_anchor);
    window_len = static_cast<unsigned>(parsed_window);
  }
  SequenceLibCollection::Build(argv[1], argv[2], num_threads, anchor_len,
                               window_len);

  return 0;
}
