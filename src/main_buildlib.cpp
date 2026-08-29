#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"

#include <cerrno>
#include <cstdlib>

void DisplayHelp(const char *program) {
  pfprintf(stderr,
           "Usage {s} <read_lib_file> <out_prefix> [num_threads]\n",
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
  SequenceLibCollection::Build(argv[1], argv[2], num_threads);

  return 0;
}
