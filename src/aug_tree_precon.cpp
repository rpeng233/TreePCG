#include <numeric>
#include <random>
#include "akpw.h"
#include "aug_tree_precon.h"
#include "stretch.h"

void AugTreePrecon(const EdgeListR& es,
                   CholeskySolver& precon,
                   size_t k) {
  EdgeList<EdgeR> tree_es;
  EdgeList<EdgeR> off_tree_es;
  TreeR tree;

  AKPW(es, tree_es);
  AdjacencyArray<ArcR> g(tree_es);
  DijkstraTree(g, es.n / 2, tree);
  g.FreeMemory();

  for (size_t i = 0; i < es.Size(); i++) {
    const EdgeR& e = es[i];
    if (tree[e.u].parent == e.v || tree[e.v].parent == e.u) {
      continue;
    }
    off_tree_es.AddEdge(e);
  }

  vector<double> stretches(off_tree_es.Size());
  vector<size_t> indices(off_tree_es.Size());

  ComputeStretch(tree, off_tree_es, stretches);

  FLOAT total_stretch = std::accumulate(stretches.begin(), stretches.end(), 0);
  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);
  AdjacencyMap aug_tree(tree);

  for (size_t i = 0; i < stretches.size(); i++) {
    FLOAT p = k * stretches[i] / total_stretch;
    if (unif01(rng) < p) {
      aug_tree.AddEdgeR(off_tree_es.edges[i]);
    }
  }

  Cholesky(aug_tree, precon.cholesky_factor);
}

