#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <random>
#include "akpw.h"
#include "aug_tree_precon.h"
#include "graph.h"
#include "stretch.h"

using std::cout;
using std::cerr;
using std::endl;

void AugTreePrecon(const EdgeListR& tree_es,
                   const std::vector<double>& stretches,
                   const EdgeListR& off_tree_es,
                   EdgeListC& precon,
                   size_t k) {
  FLOAT total_stretch = std::accumulate(stretches.begin(), stretches.end(), 0);
  cerr << "Total stretch: " << total_stretch << "\n";
  cerr << "Average stretch: " << total_stretch / off_tree_es.Size() << "\n";

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);

  struct {
    const std::vector<double> *strs;
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

  precon = tree_es;
  size_t count = 0;
  EdgeC tmp;

  for (size_t i = 0; i < k; i++) {
    tmp = off_tree_es[ids[i]];
    total_stretch -= stretches[ids[i]];
    precon.AddEdge(tmp);
    count++;
  }

  cerr << "Added the top " << count << " stretched off-tree edges\n";
  count = 0;

  for (size_t i = k; i < ids.size(); i++) {
    FLOAT p = k * stretches[ids[i]] / total_stretch;
    if (unif01(rng) < p) {
      tmp = off_tree_es[ids[i]];
      precon.AddEdge(tmp);
      count++;
    }
  }

  cerr << "Added " << count << " sampled off-tree edges\n";
}

void AugTreePrecon(const EdgeListR& es,
                   const EdgeListR& tree_es,
                   EdgeListC& precon,
                   size_t k) {
  AdjacencyArray<ArcR> g(tree_es);
  TreeR tree;

  DijkstraTree(g, es.n / 2, tree);
  g.FreeMemory();

  EdgeListR off_tree_es;
  for (size_t i = 0; i < es.Size(); i++) {
    const EdgeR& e = es[i];
    if (tree[e.u].parent == e.v || tree[e.v].parent == e.u) {
      continue;
    }
    off_tree_es.AddEdge(e);
  }

  vector<double> stretches(off_tree_es.Size());
  ComputeStretch(tree, off_tree_es, stretches);

  AugTreePrecon(tree_es, stretches, off_tree_es, precon, k);
}

void AugTreePrecon(const EdgeListR& es,
                   EdgeListC& precon,
                   size_t k) {
  EdgeListR tree_es;

  AKPW(es, tree_es);
  AugTreePrecon(es, tree_es, precon, k);
}

// void AugTreePrecon(const EdgeListR& es,
//                    CholeskySolver& precon,
//                    size_t k) {
//   EdgeList<EdgeR> tree_es;
//   EdgeList<EdgeR> off_tree_es;
//   TreeR tree;
// 
//   AKPW(es, tree_es);
//   AdjacencyArray<ArcR> g(tree_es);
//   DijkstraTree(g, es.n / 2, tree);
//   g.FreeMemory();
// 
//   for (size_t i = 0; i < es.Size(); i++) {
//     const EdgeR& e = es[i];
//     if (tree[e.u].parent == e.v || tree[e.v].parent == e.u) {
//       continue;
//     }
//     off_tree_es.AddEdge(e);
//   }
// 
//   vector<double> stretches(off_tree_es.Size());
//   vector<size_t> indices(off_tree_es.Size());
// 
//   ComputeStretch(tree, off_tree_es, stretches);
// 
//   FLOAT total_stretch = std::accumulate(stretches.begin(), stretches.end(), 0);
//   cerr << "Total stretch: " << total_stretch << "\n";
//   std::mt19937 rng(std::random_device{}());
//   std::uniform_real_distribution<> unif01(0, 1);
//   AdjacencyMap aug_tree(tree);
// 
//   size_t count = 0;
//   for (size_t i = 0; i < stretches.size(); i++) {
//     FLOAT p = k * stretches[i] / total_stretch;
//     if (unif01(rng) < p) {
//       aug_tree.AddEdgeR(off_tree_es.edges[i]);
//       count++;
//     }
//   }
// 
//   cerr << "Added " << count << " off-tree edges\n";
// 
//   Cholesky(aug_tree, precon.cholesky_factor);
// }

void AugTreePrecon2(const EdgeListR& es,
                    EdgeListC& precon,
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
  AdjacencyMap aug_tree(tree);

  std::mt19937 rng(std::random_device{}());
  std::discrete_distribution<unsigned>
    sample(stretches.begin(), stretches.end());
  std::uniform_real_distribution<> unif01(0, 1);

  for (size_t i = 0; i < off_tree_es.Size(); i++) {
    off_tree_es[i].resistance *= k * stretches[i] / total_stretch;
  }

  precon = tree_es;

  size_t count = 0;
  std::unordered_set<size_t> set;
  EdgeC tmp;
  for (size_t i = 0; i < k; i++) {
    size_t e = sample(rng);
    tmp = off_tree_es[e];
    precon.AddEdge(tmp);
    if (set.find(e) == set.end()) {
      set.insert(e);
      count++;
    }
  }

  cerr << "Added " << count << " off-tree edges\n";
}

void AugTreePrecon3(const EdgeListR& es,
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
