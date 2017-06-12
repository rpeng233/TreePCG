#ifndef INCLUDE_LOW_STRETCH_TREE_H__
#define INCLUDE_LOW_STRETCH_TREE_H__

#include <vector>
#include "binary_heap.h"
#include "common.h"
#include "disjoint_set.h"
#include "graph.h"

class LCAVertex {
public:
  DisjointSetNode ds_node;
  double dist;
  size_t ancestor;
  std::vector<size_t> children;
  std::vector<size_t> incident_edges;
  bool finished;

  LCAVertex () {
    finished = false;
    dist = 0;
  }
};

void LCA(std::vector<LCAVertex>& tree,
         const std::vector<EdgeR>& es,
         std::vector<double>& strs,
         size_t root);


template <typename TVtxT>
void ComputeStretch(const Tree<TVtxT>& tree,
                    const EdgeListR& es,
                    std::vector<double>& stretch) {
  const std::vector<TVtxT>& vs = tree.vertices;
  std::vector<LCAVertex> lca_tree(vs.size());
  size_t root = 0;
  size_t found_root = 0;

  stretch.clear();
  stretch.resize(es.Size(), -1);

  for (size_t i = 0; i < vs.size(); i++) {
    if (vs[i].parent != i) {
      lca_tree[vs[i].parent].children.push_back(i);
      lca_tree[i].dist = vs[i].parent_resistance;
    } else {
      root = i;
      found_root++;
    }
  }
  assert(found_root == 1);

  for (size_t i = 0; i < es.Size(); i++) {
    lca_tree[es[i].u].incident_edges.push_back(i);
    lca_tree[es[i].v].incident_edges.push_back(i);
  }

  LCA(lca_tree, es.edges, stretch, root);
}

#endif // INCLUDE_LOW_STRETCH_TREE_H__
