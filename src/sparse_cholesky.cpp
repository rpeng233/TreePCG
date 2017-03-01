#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include "cholesky.h"
#include "sparse_cholesky.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void SparseCholesky(AdjacencyMap& graph,
                    size_t k,
                    CholeskyFactor& cholesky_factor) {
  size_t n = graph.n;
  cholesky_factor.n = n;
  vector<size_t> order(n);

  for (size_t i = 0; i < n; i++) {
    order[i] = i;
  }
  std::random_shuffle(order.begin(), order.end());

  typedef std::map<size_t, FLOAT>::const_iterator IterType;

  vector<std::map<size_t, FLOAT>>& neighbor_map = graph.neighbor_map;
  vector<EliminatedVertex>& elims = cholesky_factor.elims;
  vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;

  std::mt19937 rng(std::random_device{}());
  vector<size_t> neighbors;
  vector<double> weights;

  elims.resize(n);
  for (size_t i = 0; i < n - 1; i++) {
    size_t u = order[i];
    elims[i].v = u;
    elims[i].degree = 0;
    elims[i].first_arc = elim_arcs.size();
    neighbors.clear();
    weights.clear();
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      elims[i].degree += it->second;
      elim_arcs.push_back(ArcC(it->first, it->second));
      neighbors.push_back(it->first);
      weights.push_back(static_cast<double>(it->second));
      neighbor_map[it->first].erase(u);
    }
    std::discrete_distribution<unsigned>
      sample(weights.begin(), weights.end());
    std::uniform_int_distribution<unsigned>
      uniform(0, neighbors.size() - 1);

    for (size_t j = 0; j < neighbors.size() * k; j++) {
      size_t a1 = sample(rng);
      size_t a2 = uniform(rng);
      if (a1 == a2) continue;
      FLOAT w1 = weights[a1] / k;
      FLOAT w2 = weights[a2] / k;
      FLOAT c = w1 * w2 / (w1 + w2);
      size_t v1 = neighbors[a1];
      size_t v2 = neighbors[a2];
      neighbor_map[v1][v2] += c;
      neighbor_map[v2][v1] += c;
    }
  }
  elims[n - 1].first_arc = elim_arcs.size();
}
