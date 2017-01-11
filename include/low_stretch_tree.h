#ifndef INCLUDE_LOW_STRETCH_TREE_H__
#define INCLUDE_LOW_STRETCH_TREE_H__

#include <limits>
#include <queue>
#include <vector>
#include "common.h"
#include "graph.h"

inline
Tree DijkstraTree(const Graph2& graph, size_t root) {
  size_t n = graph.n;
  Tree tree(n);

  std::vector<bool> finished(n, false);
  std::vector<double> dist(n, std::numeric_limits<double>::max());
  std::priority_queue<std::pair<double, size_t>> queue;

  dist[root] = 0;
  queue.push(std::make_pair(0, root));
  for (size_t i = 0; i < n; i++) {
    size_t u = queue.top().second;
    if (finished[u]) {
      continue;
    }
    queue.pop();
    finished[u] = true;
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const Arc& arc = graph.arcs[i];
      double tmp = dist[u] + arc.resistance;
      if (tmp < dist[arc.v]) {
        dist[arc.v] = tmp;
        tree.vertices[arc.v].parent = u;
        tree.vertices[arc.v].parent_resistance = arc.resistance;
        queue.push(std::make_pair(-tmp, arc.v));
      }
    }
  }

  return tree;
}

#endif // INCLUDE_LOW_STRETCH_TREE_H__
