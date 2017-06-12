#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include "aug_tree_chain.h"
#include "cholmod_solver.h"
#include "partial_cholesky.h"
#include "stretch.h"

void AugTreeChain(const EdgeListR& tree_es,
                  const EdgeListR& off_tree_es_,
                  std::vector<PreconLevel>& chain,
                  CholmodSolver& base_solver) {
  AdjacencyArray<ArcR> g(tree_es);
  Tree<PCholVertex> tree;
  DijkstraTree(g, tree_es.n / 2, tree);
  g.FreeMemory();

  EdgeListR off_tree_es(off_tree_es_);
  std::vector<double> stretch;
  ComputeStretch(tree, off_tree_es, stretch);
  EdgeListR new_tree_es;
  EdgeListR sampled_es;
  std::vector<double> sampled_stretch;

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);

  size_t n = tree.n;
  std::unordered_map<size_t,size_t> new_id;

  for (size_t i = 0; i < n; i++) {
    tree[i].parent_resistance /= 1 * log(n) * log(n);
  }

  size_t lvl = 0;
  while (n > 5000 && lvl < 3) {
    assert(lvl == chain.size());
    chain.push_back(PreconLevel());
    chain[lvl].es.n = n;

    // construct the current level graph

    chain[lvl].invD.resize(n);
    for (size_t i = 0; i < tree.n; i++) {
      if (tree[i].parent != i) {
        chain[lvl].es.AddEdge(
            EdgeC(i, tree[i].parent, 1 / tree[i].parent_resistance)
        );
        chain[lvl].invD[i] += 1 / tree[i].parent_resistance;
        chain[lvl].invD[tree[i].parent] += 1 / tree[i].parent_resistance;
      }
    }

    EdgeC tmp;
    for (size_t i = 0; i < off_tree_es.Size(); i++) {
      tmp = off_tree_es[i];
      chain[lvl].es.AddEdge(tmp);
      chain[lvl].invD[tmp.u] += tmp.conductance;
      chain[lvl].invD[tmp.v] += tmp.conductance;
    }

    for (size_t i = 0; i < n; i++) {
      assert(chain[lvl].invD[i] != 0);
      chain[lvl].invD[i] = 1 / chain[lvl].invD[i];
    }

    // sample the off-tree edges

    sampled_es.Clear();
    sampled_es.n = n;
    sampled_stretch.clear();
    double total_stretch = std::accumulate(stretch.begin(), stretch.end(), 0);
    for (size_t i = 0; i < off_tree_es.Size(); i++) {
      EdgeR& e = off_tree_es[i];
      double p = (off_tree_es.Size() / 4.0) * stretch[i] / total_stretch;
      if (p < 1) {
        if (unif01(rng) < p) {
          sampled_es.AddEdge(EdgeR(e.u, e.v, e.resistance / p));
          sampled_stretch.push_back(stretch[i]);
          tree[e.u].ref_count++;
          tree[e.v].ref_count++;
        }
      } else {
        sampled_es.AddEdge(EdgeR(e.u, e.v, e.resistance));
        sampled_stretch.push_back(stretch[i]);
        tree[e.u].ref_count++;
        tree[e.v].ref_count++;
      }
    }

    PartialCholesky(tree, sampled_es, chain[lvl].pchol.cholesky_factor);

    chain[lvl].r.resize(n);
    chain[lvl].s.resize(n);
    chain[lvl].d.resize(n);
    chain[lvl].y.resize(n);

    // prepare the next level

    n = 0;
    new_id.clear();
    for (size_t i = 0; i < tree.n; i++) {
      if (tree[i].eliminated) {
        continue;
      }
      new_id[i] = chain[lvl].id.size();
      chain[lvl].id.push_back(i);
      n++;
    }

    chain[lvl].x_small.resize(n);
    chain[lvl].b_small.resize(n);

    new_tree_es.Clear();
    new_tree_es.n = n;
    for (size_t i = 0; i < tree.n; i++) {
      if (tree[i].eliminated || tree[i].parent == i) {
        continue;
      }
      new_tree_es.AddEdge(
          EdgeR(new_id.at(i), new_id.at(tree[i].parent), tree[i].parent_resistance / 1)
      );
    }

    assert(new_tree_es.Size() == n - 1);

    tree.Resize(n);
    for (size_t i = 0; i < n - 1; i++) {
      EdgeR& e = new_tree_es[i];
      tree.SetParentR(e.u, e.v, e.resistance);
    }

    off_tree_es.Clear();
    off_tree_es.n = n;
    stretch.clear();
    for (size_t i = 0; i < sampled_es.Size(); i++) {
      EdgeR& e = sampled_es[i];
      if (e.resistance < 0) {
        continue;
      }
      off_tree_es.AddEdge(
          EdgeR(new_id.at(e.u), new_id.at(e.v), e.resistance)
      );
      stretch.push_back(sampled_stretch[i]);
    }

    chain[lvl].iter = chain[lvl].es.Size() / (off_tree_es.Size() + tree.n - 1);

    std::cout << "level " << lvl << ": "
              << chain[lvl].es.n << " vertices, "
              << chain[lvl].es.Size() << " edges, "
              << "# of recursive calls: " << chain[lvl].iter
              << '\n';

    assert(stretch.size() == off_tree_es.Size());

    lvl++;
  }

  EdgeListC base_case;

  base_case.n = tree.n;
  for (size_t i = 0; i < tree.n; i++) {
    if (tree[i].parent != i) {
      base_case.AddEdge(EdgeC(i, tree[i].parent, 1 / tree[i].parent_resistance));
    }
  }
  EdgeC tmp;
  for (size_t i = 0; i < off_tree_es.Size(); i++) {
    tmp = off_tree_es[i];
    base_case.AddEdge(tmp);
  }

  std::cout << "level " << lvl << ": "
            << base_case.n << " vertices, "
            << base_case.Size() << " edges\n";

  base_solver.Initialize(base_case);
}
