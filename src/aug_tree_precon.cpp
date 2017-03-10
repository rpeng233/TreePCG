#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <random>
#include "akpw.h"
#include "aug_tree_precon.h"
#include "stretch.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

void AugTreePrecon(const EdgeListR& es,
                   CholeskySolver& precon,
                   size_t k) {
  EdgeList<EdgeR> tree_es;
  EdgeList<EdgeR> off_tree_es;
  TreeR tree;

  AKPW2(es, tree_es);
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
  cerr << "Total stretch: " << total_stretch << "\n";
  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);
  AdjacencyMap aug_tree(tree);

  size_t count = 0;
  for (size_t i = 0; i < stretches.size(); i++) {
    FLOAT p = k * stretches[i] / total_stretch;
    if (unif01(rng) < p) {
      aug_tree.AddEdgeR(off_tree_es.edges[i]);
      count++;
    }
  }

  cerr << "Added " << count << " off-tree edges\n";

  Cholesky(aug_tree, precon.cholesky_factor);
}

void AugTreePrecon2(const EdgeListR& es,
                    CholeskySolver& precon,
                    size_t k) {
  EdgeList<EdgeR> tree_es;
  EdgeList<EdgeR> off_tree_es;
  TreeR tree;

  AKPW2(es, tree_es);
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
  AdjacencyMap aug_tree(tree);

  std::mt19937 rng(std::random_device{}());
  std::discrete_distribution<unsigned>
    sample(stretches.begin(), stretches.end());
  std::uniform_real_distribution<> unif01(0, 1);

  for (size_t i = 0; i < off_tree_es.Size(); i++) {
    off_tree_es[i].resistance *= k * stretches[i] / total_stretch;
  }

  size_t count = 0;
  std::unordered_set<size_t> set;
  for (size_t i = 0; i < k; i++) {
    size_t e = sample(rng);
    aug_tree.AddEdgeR(off_tree_es[e]);
    if (set.find(e) == set.end()) {
      set.insert(e);
      count++;
    }
  }

  cerr << "Added " << count << " off-tree edges\n";

  Cholesky(aug_tree, precon.cholesky_factor);
}

void AugTreePrecon3(const EdgeListR& es,
                    CholeskySolver& precon,
                    size_t k) {
  EdgeList<EdgeR> tree_es;
  EdgeList<EdgeR> off_tree_es;
  TreeR tree;

  AKPW2(es, tree_es);
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

  AdjacencyMap aug_tree(tree);

  struct {
    vector<double> *strs;
    bool operator()(size_t e1, size_t e2) {
      return (*strs)[e1] > (*strs)[e2];
    }
  } comp;
  comp.strs = &stretches;

  vector<size_t> ids(off_tree_es.Size());
  for (size_t i = 0; i < ids.size(); i++) {
    ids[i] = i;
  }

  std::sort(ids.begin(), ids.end(), comp);

  for (size_t i = 0; i < k; i++) {
    aug_tree.AddEdgeR(off_tree_es[ids[i]]);
  }

  cerr << "Added " << k << " off-tree edges\n";

  Cholesky(aug_tree, precon.cholesky_factor);
}
