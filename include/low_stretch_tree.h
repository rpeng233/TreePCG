#ifndef INCLUDE_LOW_STRETCH_TREE_H__
#define INCLUDE_LOW_STRETCH_TREE_H__

#include <limits>
#include <queue>
#include <vector>
#include "binary_heap.h"
#include "common.h"
#include "graph.h"

template <typename TreeType, typename ArcT>
TreeType DijkstraTree(const Graph2<ArcT>& graph, size_t root) {
  size_t n = graph.n;
  TreeType tree(n);

  std::vector<bool> finished(n, false);
  std::vector<double> dist(n, std::numeric_limits<double>::max());
  std::priority_queue<std::pair<double, size_t> > queue;

  dist[root] = 0;
  queue.push(std::make_pair(0, root));
  size_t count = 0;
  while (count < n) {
    size_t u = queue.top().second;
    queue.pop();
    if (finished[u]) {
      continue;
    }
    finished[u] = true;
    count++;
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const ArcR& arc = graph.arcs[i];
      double tmp = dist[u] + (double) arc.Resistance();
      if (tmp < dist[arc.v]) {
        dist[arc.v] = tmp;
        tree.SetParent(arc.v, u, arc.Resistance());
        queue.push(std::make_pair(-tmp, arc.v));
      }
    }
  }

  return tree;
}

void ComputeStretch(const std::vector<TreeVertexR>& vs,
                    const std::vector<EdgeR>& es,
                    std::vector<double>& stretch);

#endif // INCLUDE_LOW_STRETCH_TREE_H__
