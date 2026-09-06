/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <omp.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

#include "assembly/all_algo.h"
#include "assembly/contig_output.h"
#include "assembly/contig_stat.h"
#include "sequence/io/contig/contig_reader.h"
#include "utils/histgram.h"
#include "utils/options_description.h"
#include "utils/utils.h"

using std::string;

namespace {

struct LocalAsmOption {
  string sdbg_name;
  string output_prefix{"out"};
  int num_cpu_threads{0};
  double host_mem{0};

  int local_width{1000};
  int max_tip_len{-1};
  int min_standalone{200};
  double min_depth{-1};
  bool is_final_round{false};
  int bubble_level{2};
  int merge_len{20};
  double merge_similar{0.98};
  int prune_level{2};
  double disconnect_ratio{0.1};
  double low_local_ratio{0.2};
  int cleaning_rounds{5};
  bool output_standalone{false};
  bool careful_bubble{false};

  string contig_file() { return output_prefix + ".contigs.fa"; }
  string standalone_file() { return output_prefix + ".final.contigs.fa"; }
  string addi_contig_file() { return output_prefix + ".addi.fa"; }
  string bubble_file() { return output_prefix + ".bubble_seq.fa"; }
} opt;

bool PreferUnitigFirstInitialTips(const SDBG &dbg, uint64_t max_tip_len) {
  if (max_tip_len == 0 || dbg.size() == 0) return true;

  // Edge-first trimming performs bounded walks from the endpoint frontier at
  // every geometric threshold, then unitig construction traverses the
  // surviving graph again.  Unitig-first traverses the graph once and runs
  // the same mark/delete/contract thresholds on the much smaller unitig
  // graph.  It wins when the conservative endpoint-walk upper bound is small
  // relative to one full edge stream.  On tip-dense early graphs, retain the
  // edge-first path because pruning before compression materially shrinks the
  // expensive raw-unitig construction.
  //
  // Write the comparison as division rather than sentinel_count * length so
  // it remains exact for arbitrarily large graph metadata without overflow.
  const uint64_t quarter_graph = dbg.size() / 4u;
  return dbg.TipSentinelCount() <= quarter_graph / max_tip_len;
}

void ParseAsmOption(int argc, char *argv[]) {
  OptionsDescription desc;

  desc.AddOption("sdbg_name", "s", opt.sdbg_name,
                 "succinct de Bruijn graph name");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("num_cpu_threads", "t", opt.num_cpu_threads,
                 "number of cpu threads");
  desc.AddOption("host_mem", "", opt.host_mem,
                 "available memory budget in bytes");
  desc.AddOption("max_tip_len", "", opt.max_tip_len,
                 "max length for tips to be removed. -1 for 2k");
  desc.AddOption(
      "min_standalone", "", opt.min_standalone,
      "min length of a standalone contig to output to final.contigs.fa");
  desc.AddOption("bubble_level", "", opt.bubble_level, "bubbles level 0-3");
  desc.AddOption("merge_len", "", opt.merge_len,
                 "merge complex bubbles of length <= merge_len * k");
  desc.AddOption("merge_similar", "", opt.merge_similar,
                 "min similarity of complex bubble merging");
  desc.AddOption("prune_level", "", opt.prune_level,
                 "strength of low local depth contig pruning (0-3)");
  desc.AddOption("disconnect_ratio", "", opt.disconnect_ratio,
                 "ratio threshold for disconnecting contigs");
  desc.AddOption("low_local_ratio", "", opt.low_local_ratio,
                 "ratio to define low depth contigs");
  desc.AddOption("cleaning_rounds", "", opt.cleaning_rounds,
                 "number of rounds of graphs cleaning");
  desc.AddOption("min_depth", "", opt.min_depth,
                 "if prune_level >= 2, permanently remove low local coverage "
                 "unitigs under this threshold");
  desc.AddOption("is_final_round", "", opt.is_final_round,
                 "this is the last iteration");
  desc.AddOption("output_standalone", "", opt.output_standalone,
                 "output standalone contigs to *.final.contigs.fa");
  desc.AddOption("careful_bubble", "", opt.careful_bubble,
                 "remove bubble carefully");

  try {
    desc.Parse(argc, argv);
    if (opt.sdbg_name.empty()) {
      throw std::logic_error("no succinct de Bruijn graph name!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " -s sdbg_name -o output_prefix"
              << std::endl;
    std::cerr << "options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

}  // namespace

int main_assemble(int argc, char **argv) {
  AutoMaxRssRecorder recorder;
  ParseAsmOption(argc, argv);

  // Loading consists of independent SDBG shards, so configure the requested
  // OpenMP team before constructing the graph rather than after it is loaded.
  if (opt.num_cpu_threads == 0) {
    opt.num_cpu_threads = omp_get_max_threads();
  }
  omp_set_num_threads(opt.num_cpu_threads);

  SDBG dbg;
  SimpleTimer timer;

  // graph loading
  timer.reset();
  timer.start();
  xinfo("Loading succinct de Bruijn graph: {s}", opt.sdbg_name.c_str());
  dbg.LoadFromFile(opt.sdbg_name.c_str());
  uint64_t topology_memory_budget = 0;
  if (std::isfinite(opt.host_mem) && opt.host_mem > 0 &&
      opt.host_mem <
          static_cast<double>(std::numeric_limits<uint64_t>::max())) {
    topology_memory_budget = static_cast<uint64_t>(opt.host_mem);
  }
  const uint64_t expected_tip_walk_length =
      opt.max_tip_len < 0
          ? uint64_t{2} * dbg.k()
          : static_cast<uint64_t>(std::max(0, opt.max_tip_len));
  if (dbg.BuildTopologyCache(topology_memory_budget,
                             expected_tip_walk_length)) {
    xinfo("SDBG topology cache: {} bytes\n", dbg.TopologyCacheBytes());
  }
  timer.stop();
  xinfoc("Done. Time elapsed: {}\n", timer.elapsed());
  xinfo("Number of Edges: {}; K value: {}\n", dbg.size(), dbg.k());

  SeqPackage isolated;
  std::vector<mul_t> isolated_multi;
  std::vector<uint64_t> isolated_depth_counts;
  const string isolated_path = opt.sdbg_name + ".isolated.fa";
  struct stat isolated_stat {};
  if (::stat(packed_contig::Path(isolated_path).c_str(), &isolated_stat) == 0) {
    ContigReader reader(isolated_path);
    reader.ReadAllWithMultiplicity(&isolated, &isolated_multi, false);
    isolated_depth_counts.assign(kMaxMul + 1u, 0);
    for (size_t i = 0; i < isolated.seq_count(); ++i) {
      const uint64_t length = isolated.GetSeqView(i).length();
      if (length <= dbg.k()) xfatal("Invalid isolated path length\n");
      isolated_depth_counts[isolated_multi[i]] += 2u * (length - dbg.k());
    }
    xinfo("Restoring certified isolated paths: {} contigs / {} bases\n",
          isolated.seq_count(), isolated.base_count());
  }
  auto bridges = std::make_shared<ContigBridges>();
  bridges->Load(opt.sdbg_name, dbg);
  if (!bridges->empty()) {
    if (bridges->margin() < static_cast<unsigned>(std::max(0, opt.cleaning_rounds)))
      xfatal("Bridge endpoint margin is smaller than cleaning-round count\n");
    bridges->AddHistogram(&isolated_depth_counts);
    xinfo("Restoring {} compressed path interiors\n", bridges->size());
  }

  // set cpu threads
  xinfo("Number of CPU threads: {}\n", opt.num_cpu_threads);

  // set tip len
  if (opt.max_tip_len == -1) {
    opt.max_tip_len = dbg.k() * 2;
  }
  // set min depth
  if (opt.min_depth <= 0) {
    opt.min_depth = sdbg_pruning::InferMinDepth(dbg, isolated_depth_counts);
    xinfo("min depth set to {.3}\n", opt.min_depth);
  }

  std::unique_ptr<UnitigGraph> graph_storage;
  const bool force_unitig_first =
      std::getenv("MEGAHIT_FORCE_UNITIG_FIRST_TIPS") != nullptr ||
      std::getenv("MEGAHIT_EXPERIMENTAL_UNITIG_FIRST_TIPS") != nullptr;
  const bool disable_unitig_first =
      std::getenv("MEGAHIT_DISABLE_UNITIG_FIRST_TIPS") != nullptr;
  const bool structurally_sparse_tips = PreferUnitigFirstInitialTips(
      dbg, static_cast<uint64_t>(std::max(0, opt.max_tip_len)));
  const bool unitig_first =
      !bridges->empty() || (!disable_unitig_first && (force_unitig_first || structurally_sparse_tips));
  if (unitig_first) {
    xinfo("Initial-tip plan: unitig-first ({} sentinels, {} edges, max "
          "length {})\n",
          dbg.TipSentinelCount(), dbg.size(), opt.max_tip_len);
    timer.reset();
    timer.start();
    graph_storage.reset(new UnitigGraph(&dbg, true));
    graph_storage->AttachBridges(bridges);
    timer.stop();
    xinfo("raw unitig graph size: {}, time for building: {.3}\n",
          graph_storage->size(), timer.elapsed());

    if (opt.max_tip_len > 0) {
      timer.reset();
      timer.start();
      const uint64_t removed =
          RemoveInitialTips(*graph_storage, opt.max_tip_len);
      timer.stop();
      xinfo("Initial unitig tips removed: {} edge-bases; time elapsed(sec): "
            "{.3}\n",
            removed, timer.elapsed());
    }
    graph_storage->FinalizeInitialTipCompression();
  } else {
    xinfo("Initial-tip plan: edge-frontier ({} sentinels, {} edges, max "
          "length {})\n",
          dbg.TipSentinelCount(), dbg.size(), opt.max_tip_len);
    // Historical edge-level initial pruning remains the reference path.
    if (opt.max_tip_len > 0) {
      timer.reset();
      timer.start();
      sdbg_pruning::RemoveTips(dbg, opt.max_tip_len);
      timer.stop();
      xinfo("Tips removal done! Time elapsed(sec): {.3}\n", timer.elapsed());
    }

    timer.reset();
    timer.start();
    graph_storage.reset(new UnitigGraph(&dbg));
    timer.stop();
  }
  UnitigGraph &graph = *graph_storage;
  xinfo("unitig graph size: {}, time for building: {.3}\n", graph.size(),
        unitig_first ? 0.0 : timer.elapsed());
  CalcAndPrintStat(graph);

  // set up bubble
  ContigWriter bubble_writer(opt.bubble_file(), !opt.is_final_round);
  NaiveBubbleRemover naiver_bubble_remover;
  ComplexBubbleRemover complex_bubble_remover;
  complex_bubble_remover.SetMergeSimilarity(opt.merge_similar)
      .SetMergeLevel(opt.merge_len);
  Histgram<int64_t> bubble_hist;
  if (opt.careful_bubble) {
    naiver_bubble_remover.SetCarefulThreshold(0.2).SetWriter(&bubble_writer);
    complex_bubble_remover.SetCarefulThreshold(0.2).SetWriter(&bubble_writer);
  }

  // graph cleaning
  for (int round = 1; round <= opt.cleaning_rounds; ++round) {
    xinfo("Graph cleaning round {}\n", round);
    bool changed = false;
    if (round > 1) {
      timer.reset();
      timer.start();
      uint32_t num_tips = RemoveTips(graph, opt.max_tip_len);
      changed |= num_tips > 0;
      timer.stop();
      xinfo("Tips removed: {}, time: {.3}\n", num_tips, timer.elapsed());
    }
    // remove bubbles
    if (opt.bubble_level >= 1) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = naiver_bubble_remover.PopBubbles(graph, true);
      timer.stop();
      xinfo("Number of bubbles removed: {}, Time elapsed(sec): {.3}\n",
            num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }
    // remove complex bubbles
    if (opt.bubble_level >= 2) {
      timer.reset();
      timer.start();
      uint32_t num_bubbles = complex_bubble_remover.PopBubbles(graph, true);
      timer.stop();
      xinfo("Number of complex bubbles removed: {}, Time elapsed(sec): {}\n",
            num_bubbles, timer.elapsed());
      changed |= num_bubbles > 0;
    }

    // disconnect
    timer.reset();
    timer.start();
    uint32_t num_disconnected =
        DisconnectWeakLinks(graph, opt.disconnect_ratio);
    timer.stop();
    xinfo("Number unitigs disconnected: {}, time: {.3}\n", num_disconnected,
          timer.elapsed());
    changed |= num_disconnected > 0;

    // excessive pruning
    uint32_t num_excessive_pruned = 0;
    if (opt.prune_level >= 3) {
      timer.reset();
      timer.start();
      num_excessive_pruned = RemoveLowDepth(graph, opt.min_depth);
      num_excessive_pruned += naiver_bubble_remover.PopBubbles(graph, true);
      if (opt.bubble_level >= 2 && opt.merge_len > 0) {
        num_excessive_pruned += complex_bubble_remover.PopBubbles(graph, true);
      }
      timer.stop();
      xinfo("Unitigs removed in (more-)excessive pruning: {}, time: {.3}\n",
            num_excessive_pruned, timer.elapsed());
    } else if (opt.prune_level >= 2) {
      timer.reset();
      timer.start();
      RemoveLocalLowDepth(graph, opt.min_depth, opt.max_tip_len,
                          opt.local_width, std::min(opt.low_local_ratio, 0.1),
                          true, &num_excessive_pruned);
      timer.stop();
      xinfo("Unitigs removed in excessive pruning: {}, time: {.3}\n",
            num_excessive_pruned, timer.elapsed());
    }
    if (!changed) break;
  }

  ContigStat stat = CalcAndPrintStat(graph);

  // output contigs
  ContigWriter contig_writer(opt.contig_file(), !opt.is_final_round);
  ContigWriter standalone_writer(opt.standalone_file());

  if (!(opt.is_final_round &&
        opt.prune_level >=
            1)) {  // otherwise output after local low depth pruning
    timer.reset();
    timer.start();

    OutputContigs(graph, &contig_writer,
                  opt.output_standalone ? &standalone_writer : nullptr, false,
                  opt.min_standalone);
    timer.stop();
    xinfo("Time to output: {}\n", timer.elapsed());
  }

  // remove local low depth & output as contigs
  if (opt.prune_level >= 1) {
    ContigWriter addi_contig_writer(opt.addi_contig_file(),
                                    !opt.is_final_round);

    timer.reset();
    timer.start();
    uint32_t num_removed = IterateLocalLowDepth(
        graph, opt.min_depth, opt.max_tip_len, opt.local_width,
        opt.low_local_ratio, opt.is_final_round);

    uint32_t n_bubbles = 0;
    if (opt.bubble_level >= 2 && opt.merge_len > 0) {
      complex_bubble_remover.SetWriter(nullptr);
      n_bubbles = complex_bubble_remover.PopBubbles(graph, false);
      timer.stop();
    }
    xinfo(
        "Number of local low depth unitigs removed: {}, complex bubbles "
        "removed: {}, time: {}\n",
        num_removed, n_bubbles, timer.elapsed());
    CalcAndPrintStat(graph);

    if (!opt.is_final_round) {
      OutputContigs(graph, &addi_contig_writer, nullptr, true, 0);
    } else {
      OutputContigs(graph, &contig_writer,
                    opt.output_standalone ? &standalone_writer : nullptr, false,
                    opt.min_standalone);
    }

  }

  // These are separate linear connected components with constant coverage.
  // Only initial tip trimming and global low-depth pruning can remove them;
  // bubbles, weak links and local-depth pruning require graph neighbours.
  // Their full edge histogram above still participates in depth inference.
  const double isolated_begin = omp_get_wtime();
  uint64_t isolated_written = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : isolated_written)
  for (int64_t i = 0; i < static_cast<int64_t>(isolated.seq_count()); ++i) {
    const auto sequence = isolated.GetSeqView(i);
    if (opt.max_tip_len > 0 && sequence.length() - dbg.k() <
                                  static_cast<unsigned>(std::max(2, opt.max_tip_len))) continue;
    if (opt.cleaning_rounds > 0 && opt.prune_level >= 3 &&
        isolated_multi[i] < opt.min_depth) continue;
    if (opt.output_standalone && sequence.length() <
                                    static_cast<unsigned>(std::max(0, opt.min_standalone))) continue;
    std::string ascii(sequence.length(), 'A');
    for (unsigned j = 0; j < sequence.length(); ++j) ascii[j] = "ACGT"[sequence.base_at(j)];
    auto &writer = opt.output_standalone ? standalone_writer : contig_writer;
    writer.WriteContig(ascii, dbg.k(), static_cast<int64_t>(graph.size()) + i,
                       contig_flag::kStandalone, isolated_multi[i]);
    ++isolated_written;
  }
  if (isolated.seq_count()) {
    contig_writer.Flush();
    standalone_writer.Flush();
    xinfo("Restored {} isolated contigs in {.4} s\n", isolated_written,
          omp_get_wtime() - isolated_begin);
  }

  return 0;
}
