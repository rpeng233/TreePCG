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
                  std::vector<size_t>& new_id,
                  std::vector<size_t>& original_id,
                  CholmodSolver& base_solver) {
  // construct tree structure from list of tree edges

  AdjacencyArray<ArcR> g(tree_es);
  Tree<PCholVertex> tree;
  DijkstraTree(g, tree_es.n / 2, tree);
  g.FreeMemory();

  // initialize off-tree edges and stretches

  EdgeListR off_tree_es(off_tree_es_);
  std::vector<double> stretch;

  ComputeStretch(tree, off_tree_es, stretch);
  for (size_t i = 0; i < off_tree_es.Size(); i++) {
    const EdgeR& e = off_tree_es[i];
    tree[e.u].ref_count++;
    tree[e.v].ref_count++;
  }

  // needed for sampling

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);
  EdgeListR sampled_es;
  std::vector<double> sampled_stretch;

  const size_t N = tree_es.n;
  size_t n = tree.n;

  // scale up the tree

  for (size_t i = 0; i < n; i++) {
    tree[i].parent_resistance /= 1 * log(n) * log(n);
  }


  std::vector<size_t> eorder; // elimination order

  size_t lvl = 0;
  while (n > 1000 && lvl < 20) {
    assert(lvl == chain.size());
    chain.push_back(PreconLevel());
    chain[lvl].es.n = N;

    // save the current level graph

    // chain[lvl].invD.resize(N);
    chain[lvl].n = n;
    for (size_t i = 0; i < tree.n; i++) {
      if (!tree[i].eliminated && tree[i].parent != i) {
        chain[lvl].es.AddEdge(
            EdgeC(i, tree[i].parent, 1 / tree[i].parent_resistance)
        );
        // chain[lvl].invD[i] += 1 / tree[i].parent_resistance;
        // chain[lvl].invD[tree[i].parent] += 1 / tree[i].parent_resistance;
      }
    }

    for (size_t i = 0; i < off_tree_es.Size(); i++) {
      chain[lvl].es.AddEdge(
          EdgeC(off_tree_es[i].u,
                off_tree_es[i].v,
                1 / off_tree_es[i].resistance)
      );
      // chain[lvl].invD[tmp.u] += tmp.conductance;
      // chain[lvl].invD[tmp.v] += tmp.conductance;
    }

    // for (size_t i = 0; i < N; i++) {
    //   if (chain[lvl].invD[i] != 0) {
    //     chain[lvl].invD[i] = 1 / chain[lvl].invD[i];
    //   }
    // }

    // sample the off-tree edges

    sampled_es.Clear();
    sampled_es.n = N;
    sampled_stretch.clear();
    double total_stretch = std::accumulate(stretch.begin(), stretch.end(), 0);
    for (size_t i = 0; i < off_tree_es.Size(); i++) {
      EdgeR& e = off_tree_es[i];
      double p = (off_tree_es.Size() / 4.0) * stretch[i] / total_stretch;
      if (p < 1) {
        if (unif01(rng) < p) {
          sampled_es.AddEdge(EdgeR(e.u, e.v, e.resistance / p));
          sampled_stretch.push_back(stretch[i]);
        } else {
          tree[e.u].ref_count--;
          tree[e.v].ref_count--;
        }
      } else {
        sampled_es.AddEdge(e);
        sampled_stretch.push_back(stretch[i]);
      }
    }

    for (size_t i = 0; i < tree.n; i++) {
      assert(tree[i].ref_count >= 0);
    }

    PartialCholesky(tree, sampled_es, chain[lvl].pchol.cholesky_factor, eorder);

    chain[lvl].r.resize(chain[lvl].n);
    chain[lvl].s.resize(chain[lvl].n);
    chain[lvl].d.resize(chain[lvl].n);
    chain[lvl].y.resize(chain[lvl].n);

    // prepare the next level

    n = 0;
    for (size_t i = 0; i < tree.n; i++) {
      if (!tree[i].eliminated) {
        tree[i].parent_resistance /= 4;
        n++;
      }
    }

    // new_tree_es.Clear();
    // new_tree_es.n = n;
    // for (size_t i = 0; i < tree.n; i++) {
    //   if (tree[i].eliminated || tree[i].parent == i) {
    //     continue;
    //   }
    //   new_tree_es.AddEdge(
    //       EdgeR(new_id.at(i), new_id.at(tree[i].parent), tree[i].parent_resistance / 1)
    //   );
    // }

    // assert(new_tree_es.Size() == n - 1);

    // tree.Resize(n);
    // for (size_t i = 0; i < n - 1; i++) {
    //   EdgeR& e = new_tree_es[i];
    //   tree.SetParentR(e.u, e.v, e.resistance);
    // }

    off_tree_es.Clear();
    off_tree_es.n = N;
    stretch.clear();
    for (size_t i = 0; i < sampled_es.Size(); i++) {
      EdgeR& e = sampled_es[i];
      if (e.resistance < 0) {
        continue;
      }
      off_tree_es.AddEdge(e);
      stretch.push_back(sampled_stretch[i]);
    }

    if (off_tree_es.Size() + n - 1 == 0) {
      break;
    }

    chain[lvl].iter = chain[lvl].es.Size() / (off_tree_es.Size() + n - 1);

    std::cout << "level " << lvl << ": "
              << chain[lvl].n << " vertices, "
              << chain[lvl].es.Size() << " edges, "
              << "# of recursive calls: " << chain[lvl].iter
              << '\n';

    assert(stretch.size() == off_tree_es.Size());

    lvl++;
  }

  n = 0;
  for (size_t i = 0; i < tree.n; i++) {
    if (!tree[i].eliminated) {
      eorder.push_back(i);
      n++;
    }
  }

  assert(eorder.size() == N);

  new_id.resize(N);
  original_id.resize(N);

  for (size_t i = 0; i < N; i++) {
    new_id[eorder[N - i - 1]] = i;
    original_id[i] = eorder[N - i - 1];
  }

  for (size_t lvl = 0; lvl < chain.size(); lvl++) {
    chain[lvl].es.n = chain[lvl].n;
    chain[lvl].invD.resize(chain[lvl].n);
    assert(chain[lvl].r.size() == chain[lvl].n);
    assert(chain[lvl].s.size() == chain[lvl].n);
    assert(chain[lvl].d.size() == chain[lvl].n);
    assert(chain[lvl].y.size() == chain[lvl].n);
    for (size_t i = 0; i < chain[lvl].es.Size(); i++) {
      EdgeC& e = chain[lvl].es[i];
      e.u = new_id[e.u];
      e.v = new_id[e.v];
      assert(e.u < chain[lvl].es.n);
      assert(e.v < chain[lvl].es.n);
      chain[lvl].invD[e.u] += e.conductance;
      chain[lvl].invD[e.v] += e.conductance;
    }

    for (size_t i = 0; i < chain[lvl].invD.size(); i++) {
      chain[lvl].invD[i] = 1 / chain[lvl].invD[i];
    }

    CholeskyFactor& pchol = chain[lvl].pchol.cholesky_factor;
    pchol.n = chain[lvl].n;
    for (size_t i = 0; i < pchol.elims.size(); i++) {
      pchol.elims[i].v = new_id[pchol.elims[i].v];
      assert(pchol.elims[i].v < pchol.n);
    }
    for (size_t i = 0; i < pchol.elim_arcs.size(); i++) {
      ArcC& a = pchol.elim_arcs[i];
      a.v = new_id[a.v];
    }

  }

  EdgeListC base_case;

  base_case.n = n;
  for (size_t i = 0; i < tree.n; i++) {
    if (!tree[i].eliminated && tree[i].parent != i) {
      base_case.AddEdge(
          EdgeC(new_id[i], new_id[tree[i].parent], 1 / tree[i].parent_resistance)
      );
    }
  }

  for (size_t i = 0; i < off_tree_es.Size(); i++) {
    EdgeR& e = off_tree_es[i];
    base_case.AddEdge(
        EdgeC(new_id[e.u], new_id[e.v], 1 / e.resistance)
    );
  }

  std::cout << "level " << lvl << ": "
            << n << " vertices, "
            << base_case.Size() << " edges\n";

  base_solver.Initialize(base_case);
}
