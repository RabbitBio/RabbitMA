//
// Created by vout on 11/21/18.
//

#ifndef MEGAHIT_TIP_REMOVER_H
#define MEGAHIT_TIP_REMOVER_H

#include <cstdint>

class UnitigGraph;

uint32_t RemoveTips(UnitigGraph &graph, uint32_t max_tip_len);

// Initial SDBG tip pruning has no depth-ratio predicate: at each geometric
// length barrier it removes every short isolated/dead-end path.  Once maximal
// simple paths have already been compressed, the same operation can be
// expressed on UnitigGraph vertices without walking every constituent SDBG
// edge again.  This entry point deliberately remains separate from the later
// depth-aware RemoveTips() used during graph-cleaning rounds.
uint64_t RemoveInitialTips(UnitigGraph &graph, uint32_t max_tip_len);

#endif  // MEGAHIT_TIP_REMOVER_H
