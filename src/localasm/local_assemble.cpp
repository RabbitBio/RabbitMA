#include "local_assemble.h"

#include <algorithm>
#include <cassert>
#include <fcntl.h>
#include <iostream>
#include <limits>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <omp.h>
#include "idba/contig_graph.h"
#include "idba/hash_graph.h"
#include "idba/sequence.h"
#include "kmlib/kmbit.h"

#include "hash_mapper.h"
#include "mapping_result_collector.h"
#include "sequence/io/contig/contig_reader.h"
#include "sequence/io/contig/contig_writer.h"
#include "sequence/io/read_chunk_index.h"
#include "sequence/io/sequence_lib.h"
#include "utils/histgram.h"
#include "utils/utils.h"

namespace {

static const int kMaxLocalRange = 650;
using TInsertSize = std::pair<double, double>;

struct PackedReadRecord {
  const uint32_t *words;
  uint32_t length;
};

// Read-only access to buildlib's native [length][packed words] stream.  The
// mapping is virtual; callers discard completed chunk pages so resident input
// is bounded by the active library/chunks instead of the complete read set.
class MappedReadFile {
 public:
  MappedReadFile() = default;
  ~MappedReadFile() { Close(); }
  MappedReadFile(const MappedReadFile &) = delete;
  MappedReadFile &operator=(const MappedReadFile &) = delete;

  bool Open(const std::string &lib_prefix,
            const SequenceLibCollection::SizeInfo &expected) {
    Close();
    path_ = lib_prefix + ".bin";
    if (!LoadPackedReadChunkIndex(path_, &index_) ||
        index_.num_reads != static_cast<uint64_t>(expected.num_reads) ||
        index_.num_bases != static_cast<uint64_t>(expected.num_bases) ||
        index_.max_read_len != expected.max_read_len) {
      return false;
    }

    fd_ = open(path_.c_str(), O_RDONLY | O_CLOEXEC);
    struct stat status;
    if (fd_ < 0 || fstat(fd_, &status) != 0 || status.st_size <= 0 ||
        status.st_size % static_cast<off_t>(sizeof(uint32_t)) != 0) {
      Close();
      return false;
    }
    bytes_ = static_cast<size_t>(status.st_size);
    void *address = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (address == MAP_FAILED) {
      words_ = nullptr;
      Close();
      return false;
    }
    words_ = static_cast<const uint32_t *>(address);
    return true;
  }

  void Close() {
    if (words_ != nullptr) {
      munmap(const_cast<uint32_t *>(words_), bytes_);
      words_ = nullptr;
    }
    if (fd_ >= 0) {
      close(fd_);
      fd_ = -1;
    }
    bytes_ = 0;
    index_ = PackedReadChunkIndex();
  }

  const PackedReadChunkIndex &index() const { return index_; }
  const uint32_t *words() const { return words_; }
  size_t total_words() const { return bytes_ / sizeof(uint32_t); }

  const uint32_t *LocateRead(const PackedReadChunk &chunk,
                             uint64_t read_id) const {
    assert(read_id >= chunk.read_begin && read_id < chunk.read_end);
    const uint32_t *cursor = words_ + chunk.word_begin;
    for (uint64_t id = chunk.read_begin; id < read_id; ++id) {
      const uint32_t len = *cursor++;
      cursor += DivCeiling(static_cast<size_t>(len),
                           SeqPackage::kBasesPerWord);
    }
    return cursor;
  }

  static PackedReadRecord Next(const uint32_t **cursor) {
    const uint32_t len = *(*cursor)++;
    const uint32_t *words = *cursor;
    *cursor += DivCeiling(static_cast<size_t>(len),
                          SeqPackage::kBasesPerWord);
    return {words, len};
  }

  void DropChunk(const PackedReadChunk &chunk) const {
    DiscardMemoryPages(
        const_cast<uint32_t *>(words_ + chunk.word_begin),
        (chunk.word_end - chunk.word_begin) * sizeof(uint32_t));
  }

 private:
  std::string path_;
  int fd_{-1};
  size_t bytes_{0};
  const uint32_t *words_{nullptr};
  PackedReadChunkIndex index_;
};

void LaunchIDBA(const std::vector<Sequence> &reads, size_t num_reads,
                const Sequence &contig_end,
                std::deque<Sequence> &out_contigs,
                std::deque<ContigInfo> &out_contig_infos, uint32_t mink,
                uint32_t maxk, uint32_t step, ContigGraph &contig_graph) {
  int local_range = contig_end.size();
  HashGraph hash_graph;
  out_contigs.clear();
  out_contig_infos.clear();

  uint32_t max_read_len = 0;

  for (size_t read_id = 0; read_id < num_reads; ++read_id) {
    max_read_len = std::max(max_read_len, reads[read_id].size());
  }

  for (uint32_t kmer_size = mink; kmer_size <= std::min(maxk, max_read_len);
       kmer_size += step) {
    hash_graph.clear();
    hash_graph.set_kmer_size(kmer_size);

    for (size_t read_id = 0; read_id < num_reads; ++read_id) {
      const Sequence &read = reads[read_id];
      if (read.size() < kmer_size) continue;
      hash_graph.InsertKmers(read);
    }

    double mean = hash_graph.coverage_percentile(
        1 - 1.0 * local_range / hash_graph.num_vertices());
    double threshold = mean;

    hash_graph.InsertKmers(contig_end);

    for (const auto &out_contig : out_contigs)
      hash_graph.InsertUncountKmers(out_contig);

    hash_graph.Assemble(out_contigs, out_contig_infos);

    contig_graph.clear();
    contig_graph.set_kmer_size(kmer_size);
    contig_graph.Initialize(out_contigs, out_contig_infos);
    contig_graph.RemoveDeadEnd(kmer_size * 2);

    contig_graph.RemoveBubble();
    contig_graph.IterateCoverage(kmer_size * 2, 1, threshold);

    contig_graph.Assemble(out_contigs, out_contig_infos);

    if (out_contigs.size() == 1) {
      break;
    }
  }
}

std::vector<TInsertSize> EstimateInsertSize(
    const HashMapper &mapper, const SequenceLibCollection &lib_collection) {
  std::vector<TInsertSize> insert_sizes(lib_collection.size());
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto lib = lib_collection.GetLib(lib_id);

    if (!lib.IsPaired()) {
      continue;
    }

    Histgram<int> insert_hist;
    std::vector<Histgram<int, NullMutex>> thread_insert_hist(
        omp_get_max_threads());
    const size_t min_hist_size_for_estimation = 1u << 18;
    size_t processed_reads = 0;

    while (insert_hist.size() < min_hist_size_for_estimation &&
           processed_reads < lib.seq_count()) {
      size_t start_read_id = processed_reads;
      processed_reads = std::min(min_hist_size_for_estimation + start_read_id,
                                 lib.seq_count());

      for (auto &local_hist : thread_insert_hist) {
        local_hist.clear();
      }

#pragma omp parallel for
      for (size_t i = start_read_id; i < processed_reads; i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        auto rec1 = mapper.TryMap(seq1);
        auto rec2 = mapper.TryMap(seq2);
        if (rec1.valid && rec2.valid) {
          if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
            int insert_size;

            if (rec1.strand == 0) {
              insert_size = rec2.contig_to + seq2.length() - rec2.query_to -
                            (rec1.contig_from - rec1.query_from);
            } else {
              insert_size = rec1.contig_to + seq1.length() - rec1.query_to -
                            (rec2.contig_from - rec2.query_from);
            }

            if (insert_size >= (int)seq1.length() &&
                insert_size >= (int)seq2.length()) {
              thread_insert_hist[omp_get_thread_num()].insert(insert_size);
            }
          }
        }
      }

      for (const auto &local_hist : thread_insert_hist) {
        insert_hist.MergeFrom(local_hist);
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes[lib_id] = TInsertSize(insert_hist.mean(), insert_hist.sd());

    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id,
          insert_sizes[lib_id].first, insert_sizes[lib_id].second);
  }

  return insert_sizes;
}

std::vector<TInsertSize> EstimateInsertSizeMapped(
    const HashMapper &mapper, const SequenceLibCollection &lib_collection,
    const MappedReadFile &reads) {
  const size_t target_samples = size_t{1} << 18u;
  const int num_threads = std::max(1, omp_get_max_threads());
  const auto &chunks = reads.index().chunks;
  std::vector<TInsertSize> insert_sizes(lib_collection.size());

  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    const auto &lib = lib_collection.GetLib(lib_id);
    if (!lib.IsPaired()) continue;

    Histgram<int> insert_hist;
    std::vector<Histgram<int, NullMutex>> thread_hist(num_threads);
    uint64_t processed_reads = 0;
    while (insert_hist.size() < target_samples &&
           processed_reads < lib.seq_count()) {
      const uint64_t local_begin = processed_reads;
      processed_reads = std::min<uint64_t>(
          processed_reads + target_samples, lib.seq_count());
      const uint64_t range_begin = lib.global_begin() + local_begin;
      const uint64_t range_end = lib.global_begin() + processed_reads;
      const size_t num_pairs = static_cast<size_t>(
          (range_end - range_begin) / 2u);
      std::vector<uint64_t> pair_offsets(num_pairs,
                                         std::numeric_limits<uint64_t>::max());
      std::vector<uint8_t> touched(chunks.size(), 0u);

#pragma omp parallel for schedule(static)
      for (int64_t chunk_id = 0;
           chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
        const PackedReadChunk &chunk = chunks[chunk_id];
        uint64_t first = std::max(chunk.read_begin, range_begin);
        const uint64_t end = std::min(chunk.read_end, range_end);
        if (((first - lib.global_begin()) & 1u) != 0u) ++first;
        if (first >= end) continue;
        touched[chunk_id] = 1u;
        const uint32_t *cursor = reads.LocateRead(chunk, first);
        for (uint64_t read_id = first; read_id < end; read_id += 2u) {
          pair_offsets[(read_id - range_begin) / 2u] =
              static_cast<uint64_t>(cursor - reads.words());
          MappedReadFile::Next(&cursor);
          MappedReadFile::Next(&cursor);
        }
      }

      for (auto &hist : thread_hist) hist.clear();
#pragma omp parallel for schedule(static)
      for (int64_t pair_id = 0;
           pair_id < static_cast<int64_t>(pair_offsets.size()); ++pair_id) {
        const uint64_t offset = pair_offsets[pair_id];
        if (offset == std::numeric_limits<uint64_t>::max()) continue;
        const uint64_t read_id = range_begin + uint64_t(pair_id) * 2u;
        const uint32_t *cursor = reads.words() + offset;
        const PackedReadRecord seq1 = MappedReadFile::Next(&cursor);
        const PackedReadRecord seq2 = MappedReadFile::Next(&cursor);
        const MappingRecord rec1 =
            mapper.TryMap(seq1.words, seq1.length, read_id);
        const MappingRecord rec2 =
            mapper.TryMap(seq2.words, seq2.length, read_id + 1u);
        if (rec1.valid && rec2.valid && rec1.contig_id == rec2.contig_id &&
            rec1.strand != rec2.strand) {
          int insert_size;
          if (rec1.strand == 0) {
            insert_size = rec2.contig_to + seq2.length - rec2.query_to -
                          (rec1.contig_from - rec1.query_from);
          } else {
            insert_size = rec1.contig_to + seq1.length - rec1.query_to -
                          (rec2.contig_from - rec2.query_from);
          }
          if (insert_size >= static_cast<int>(seq1.length) &&
              insert_size >= static_cast<int>(seq2.length)) {
            thread_hist[omp_get_thread_num()].insert(insert_size);
          }
        }
      }
      for (const auto &hist : thread_hist) insert_hist.MergeFrom(hist);

#pragma omp parallel for schedule(static)
      for (int64_t chunk_id = 0;
           chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
        if (touched[chunk_id] != 0u) reads.DropChunk(chunks[chunk_id]);
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes[lib_id] =
        TInsertSize(insert_hist.mean(), insert_hist.sd());
    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id,
          insert_sizes[lib_id].first, insert_sizes[lib_id].second);
  }
  return insert_sizes;
}

int32_t LocalRange(const SequenceLib &lib, const TInsertSize &insert_size) {
  int32_t local_range = lib.GetMaxLength() - 1;

  if (lib.IsPaired() && insert_size.first >= lib.GetMaxLength()) {
    local_range = std::min(2 * insert_size.first,
                           insert_size.first + 3 * insert_size.second);
  }

  if (local_range > kMaxLocalRange) {
    local_range = kMaxLocalRange;
  }

  return local_range;
}

int32_t GetMaxLocalRange(const SequenceLibCollection &lib_collection,
                         const std::vector<TInsertSize> &insert_sizes) {
  int32_t max_local_range = 0;
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    max_local_range =
        std::max(max_local_range, LocalRange(lib, insert_sizes[lib_id]));
  }
  return max_local_range;
}

void MapToContigs(const HashMapper &mapper,
                  const SequenceLibCollection &lib_collection,
                  const std::vector<TInsertSize> &insert_sizes,
                  MappingResultCollector *collector) {
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    int32_t local_range = LocalRange(lib, insert_sizes[lib_id]);
    bool is_paired = lib.IsPaired();

    size_t num_added = 0, num_mapped = 0;

    if (is_paired) {
#pragma omp parallel for reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        auto rec1 = mapper.TryMap(seq1);
        auto rec2 = mapper.TryMap(seq2);

        if (rec1.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec1.contig_id).length();
          num_added += collector->AddSingle(rec1, contig_len, seq1.length(),
                                            local_range);
          num_added += collector->AddMate(rec1, rec2, contig_len, seq2.id(),
                                          local_range);
        }

        if (rec2.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec2.contig_id).length();
          num_added += collector->AddSingle(rec2, contig_len, seq2.length(),
                                            local_range);
          num_added += collector->AddMate(rec2, rec1, contig_len, seq1.id(),
                                          local_range);
        }
      }
    } else {
#pragma omp parallel reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); ++i) {
        auto seq = lib.GetSequenceView(i);
        auto rec = mapper.TryMap(seq);

        if (rec.valid) {
          ++num_mapped;
          num_added += collector->AddSingle(
              rec, mapper.refseq().GetSeqView(rec.contig_id).length(),
              seq.length(), local_range);
        }
      }
    }

    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads to local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
  }
}

void MapToContigsMapped(const HashMapper &mapper,
                        const SequenceLibCollection &lib_collection,
                        const std::vector<TInsertSize> &insert_sizes,
                        const MappedReadFile &reads,
                        MappingResultCollector *collector) {
  const auto &chunks = reads.index().chunks;
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    const auto &lib = lib_collection.GetLib(lib_id);
    const int32_t local_range = LocalRange(lib, insert_sizes[lib_id]);
    const uint64_t lib_begin = lib.global_begin();
    const uint64_t lib_end = lib.global_end();
    uint64_t num_added = 0;
    uint64_t num_mapped = 0;

#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : num_added, num_mapped)
    for (int64_t chunk_id = 0;
         chunk_id < static_cast<int64_t>(chunks.size()); ++chunk_id) {
      const PackedReadChunk &chunk = chunks[chunk_id];
      uint64_t first = std::max(chunk.read_begin, lib_begin);
      const uint64_t end = std::min(chunk.read_end, lib_end);
      if (first >= end) continue;

      if (lib.IsPaired()) {
        if (((first - lib_begin) & 1u) != 0u) ++first;
        if (first >= end) {
          reads.DropChunk(chunk);
          continue;
        }
        const uint32_t *cursor = reads.LocateRead(chunk, first);
        for (uint64_t read_id = first; read_id < end; read_id += 2u) {
          const PackedReadRecord seq1 = MappedReadFile::Next(&cursor);
          const PackedReadRecord seq2 = MappedReadFile::Next(&cursor);
          const MappingRecord rec1 =
              mapper.TryMap(seq1.words, seq1.length, read_id);
          const MappingRecord rec2 =
              mapper.TryMap(seq2.words, seq2.length, read_id + 1u);

          if (rec1.valid) {
            ++num_mapped;
            const int32_t contig_len =
                mapper.refseq().GetSeqView(rec1.contig_id).length();
            num_added += collector->AddSingle(
                rec1, contig_len, seq1.length, local_range);
            num_added += collector->AddMate(
                rec1, rec2, contig_len, read_id + 1u, local_range);
          }
          if (rec2.valid) {
            ++num_mapped;
            const int32_t contig_len =
                mapper.refseq().GetSeqView(rec2.contig_id).length();
            num_added += collector->AddSingle(
                rec2, contig_len, seq2.length, local_range);
            num_added += collector->AddMate(
                rec2, rec1, contig_len, read_id, local_range);
          }
        }
      } else {
        const uint32_t *cursor = reads.LocateRead(chunk, first);
        for (uint64_t read_id = first; read_id < end; ++read_id) {
          const PackedReadRecord seq = MappedReadFile::Next(&cursor);
          const MappingRecord rec =
              mapper.TryMap(seq.words, seq.length, read_id);
          if (rec.valid) {
            ++num_mapped;
            num_added += collector->AddSingle(
                rec, mapper.refseq().GetSeqView(rec.contig_id).length(),
                seq.length, local_range);
          }
        }
      }
      // A chunk is owned by exactly one loop iteration.  Drop its file-backed
      // pages as soon as all records for this library have been consumed,
      // rather than retaining an entire large library until the OpenMP
      // barrier.  This bounds resident input by active chunks independently
      // of how users divide the same reads into files or libraries.
      reads.DropChunk(chunk);
    }

    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads to local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
  }
}

void AssembleAndOutput(const HashMapper &mapper, const SeqPackage &read_pkg,
                       MappingResultCollector &result_collector,
                       const std::string &output_file,
                       const int32_t local_range,
                       const LocalAsmOption &opt) {
  const size_t min_num_reads = read_pkg.max_length() > 0 ?
      local_range / read_pkg.max_length(): 1;
  xinfo("Minimum number of reads to do local assembly: {}\n", min_num_reads);

  Sequence contig_end;
  ContigGraph contig_graph;
  // The vector is private to an OpenMP worker and survives across its dynamic
  // endpoint tasks.  Keeping inactive Sequence objects avoids destroying and
  // reallocating millions of small read strings when the next endpoint starts.
  std::vector<Sequence> reads;
  std::deque<Sequence> out_contigs;
  std::deque<ContigInfo> out_contig_infos;

  ContigWriter local_contig_writer(output_file);

  struct LocalAssemblyTask {
    uint64_t contig_id;
    uint64_t mapping_count;
    uint8_t strand;
  };
  std::vector<LocalAssemblyTask> tasks;
  tasks.reserve(mapper.refseq().seq_count());
  for (uint64_t cid = 0; cid < mapper.refseq().seq_count(); ++cid) {
    for (uint8_t strand = 0; strand < 2; ++strand) {
      const size_t mapping_count =
          result_collector.GetMappingResults(cid, strand).size();
      if (mapping_count > min_num_reads) {
        tasks.push_back({cid, mapping_count, strand});
      }
    }
  }
  // Each contig end is independent.  Scheduling the largest read sets first
  // prevents a handful of expensive endpoints from forming a long serial
  // tail after the graph has otherwise drained, and splitting the two ends
  // doubles useful task granularity without altering an endpoint's inputs.
  std::sort(tasks.begin(), tasks.end(),
            [](const LocalAssemblyTask &lhs, const LocalAssemblyTask &rhs) {
              if (lhs.mapping_count != rhs.mapping_count) {
                return lhs.mapping_count > rhs.mapping_count;
              }
              if (lhs.contig_id != rhs.contig_id) {
                return lhs.contig_id < rhs.contig_id;
              }
              return lhs.strand < rhs.strand;
            });
  xinfo("Scheduling {} independent local-assembly endpoint tasks\n",
        tasks.size());

#pragma omp parallel for private(contig_graph, contig_end, reads, \
                                 out_contigs, out_contig_infos)        \
    schedule(dynamic, 1)
  for (int64_t task_id = 0; task_id < static_cast<int64_t>(tasks.size());
       ++task_id) {
    const uint64_t cid = tasks[task_id].contig_id;
    const uint8_t strand = tasks[task_id].strand;
    auto contig_view = mapper.refseq().GetSeqView(cid);
    int cl = contig_view.length();

    auto mapping_rslts = result_collector.GetMappingResults(cid, strand);

    // collect local reads, convert them into Sequence
    size_t num_local_reads = 0;
    uint64_t last_mapping_pos = -1;
    int pos_count = 0;

    for (const auto &encoded_rslt : mapping_rslts) {
      uint64_t pos = MappingResultCollector::GetContigAbsPos(encoded_rslt);
      pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
      last_mapping_pos = pos;

      if (pos_count <= 3) {
        auto read_view = read_pkg.GetSeqView(
            MappingResultCollector::GetReadId(encoded_rslt));

        if (num_local_reads == reads.size()) {
          reads.emplace_back();
        }
        Sequence &read = reads[num_local_reads++];
        const unsigned read_length = read_view.length();
        read.resize(read_length);
        for (unsigned ri = 0; ri < read_length; ++ri) {
          read[ri] = read_view.base_at(ri);
        }
      }
    }

    if (strand == 0) {
      const int end = std::min(local_range, cl);
      contig_end.resize(end);
      for (int j = 0; j < end; ++j) {
        contig_end[j] = contig_view.base_at(j);
      }
    } else {
      const int begin = std::max(0, cl - local_range);
      contig_end.resize(cl - begin);
      for (int j = begin; j < cl; ++j) {
        contig_end[j - begin] = contig_view.base_at(j);
      }
    }

    out_contigs.clear();
    LaunchIDBA(reads, num_local_reads, contig_end, out_contigs,
               out_contig_infos, opt.kmin, opt.kmax, opt.step, contig_graph);

    for (uint64_t j = 0; j < out_contigs.size(); ++j) {
      if (out_contigs[j].size() > opt.min_contig_len &&
          out_contigs[j].size() > opt.kmax) {
        auto str = out_contigs[j].str();
        local_contig_writer.WriteLocalContig(str, cid, strand, j);
      }
    }
  }
}

}  // namespace

void RunLocalAssembly(const LocalAsmOption &opt) {
  SimpleTimer timer;
  timer.reset();
  timer.start();
  HashMapper mapper;
  mapper.LoadAndBuild(opt.contig_file, opt.min_contig_len, opt.seed_kmer,
                      opt.sparsity);
  mapper.SetMappingThreshold(opt.min_mapping_len, opt.similarity);
  timer.stop();
  xinfo("Hash mapper construction time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  SequenceLibCollection lib_collection;
  SeqPackage read_pkg;
  lib_collection.SetPath(opt.lib_file_prefix);
  const SequenceLibCollection::SizeInfo read_size =
      lib_collection.GetSizeInfo();
  MappedReadFile mapped_reads;
  const bool use_mapped_stream =
      mapped_reads.Open(opt.lib_file_prefix, read_size);
  if (use_mapped_stream) {
    lib_collection.ReadMetadata();
    xinfo("Using {} indexed binary-read chunks; full read package is not "
          "materialized\n",
          mapped_reads.index().chunks.size());
  } else {
    lib_collection.Read(&read_pkg);
  }
  timer.stop();
  xinfo("Read lib time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  const auto insert_sizes =
      use_mapped_stream
          ? EstimateInsertSizeMapped(mapper, lib_collection, mapped_reads)
          : EstimateInsertSize(mapper, lib_collection);
  timer.stop();
  xinfo("Insert size estimation time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();

  MappingResultCollector collector(mapper.refseq().seq_count());
  if (use_mapped_stream) {
    MapToContigsMapped(mapper, lib_collection, insert_sizes, mapped_reads,
                       &collector);
  } else {
    MapToContigs(mapper, lib_collection, insert_sizes, &collector);
  }
  timer.stop();
  xinfo("Mapping time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  collector.Finalize();
  timer.stop();
  xinfo("Mapping result collation time elapsed: {}, retained {} mappings\n",
        timer.elapsed(), collector.size());

  // Mapping is the only consumer of the large contig seed index.  Release it
  // before constructing the compact read package and per-endpoint graphs.
  mapper.ReleaseIndex();

  timer.reset();
  timer.start();
  const size_t original_read_count = static_cast<size_t>(read_size.num_reads);
  const size_t original_read_bytes =
      use_mapped_stream
          ? DivCeiling(static_cast<size_t>(read_size.num_bases),
                       SeqPackage::kBasesPerWord) *
                sizeof(SeqPackage::TWord)
          : read_pkg.size_in_byte();
  std::vector<uint64_t> selected_read_ids =
      collector.CompactReadIds(original_read_count);
  SeqPackage selected_read_pkg;
  if (use_mapped_stream) {
    if (!selected_read_pkg.AssignSelectedMappedBinaryRecords(
            mapped_reads.words(), mapped_reads.total_words(),
            original_read_count, selected_read_ids,
            mapped_reads.index().chunks, omp_get_max_threads())) {
      xfatal("Indexed binary-read stream changed during local mapping\n");
    }
    mapped_reads.Close();
  } else {
    selected_read_pkg.ReserveSequences(selected_read_ids.size());
    if (read_pkg.max_length() != 0 &&
        selected_read_ids.size() <=
            std::numeric_limits<size_t>::max() / read_pkg.max_length()) {
      selected_read_pkg.ReserveBases(selected_read_ids.size() *
                                     read_pkg.max_length());
    }
    for (const uint64_t read_id : selected_read_ids) {
      selected_read_pkg.AppendSequenceView(read_pkg.GetSeqView(read_id));
    }
  }
  const unsigned selected_gap_bits =
      selected_read_pkg.CompactReadOnlyLengthIndex();
  const size_t selected_read_bytes = selected_read_pkg.size_in_byte();
  read_pkg.ReleaseStorage();
  read_pkg = std::move(selected_read_pkg);
  std::vector<uint64_t>().swap(selected_read_ids);
  timer.stop();
  xinfo("Compacted referenced reads: {} / {} reads, {} -> {} bytes "
        "(length gap index {} bits), elapsed {}\n",
        read_pkg.seq_count(), original_read_count, original_read_bytes,
        selected_read_bytes, selected_gap_bits, timer.elapsed());

  timer.reset();
  timer.start();
  int32_t max_local_range = GetMaxLocalRange(lib_collection, insert_sizes);
  AssembleAndOutput(mapper, read_pkg, collector, opt.output_file,
                    max_local_range, opt);
  timer.stop();
  xinfo("Local assembly time elapsed: {}\n", timer.elapsed());
}
