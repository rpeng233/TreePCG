#ifndef INCLUDE_LOW_STRETCH_TREE_H__
#define INCLUDE_LOW_STRETCH_TREE_H__

#include <vector>
#include "binary_heap.h"
#include "common.h"
#include "graph.h"

void ComputeStretch(const TreeR& tree,
                    const EdgeListR& es,
                    std::vector<double>& stretch);

#endif // INCLUDE_LOW_STRETCH_TREE_H__
