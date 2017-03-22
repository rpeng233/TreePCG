#include <cassert>
#include <iostream>
#include <random>
#include <map>
#include <vector>
#include "cholesky.h"
#include "graph.h"
#include "sparse_cholesky.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

struct Upper {
  size_t n;
  std::vector<std::map<size_t, FLOAT>> neighbor_map;

  Upper(const EdgeList<EdgeC>& es) {
    n = es.n;
    neighbor_map.resize(n);
    for (size_t i = 0; i < es.edges.size(); i++) {
      size_t u = std::min(es[i].u, es[i].v);
      size_t v = std::max(es[i].u, es[i].v);
      neighbor_map[u][v] += es[i].conductance;
    }
  }
};

void SparseCholesky(EdgeListC& es,
                    size_t k,
                    CholeskyFactor& cholesky_factor) {
  size_t n = es.n;
  cholesky_factor.n = n;
  vector<size_t> new_id(n);

  for (size_t i = 0; i < n; i++) {
    new_id[i] = i;
  }
  std::random_shuffle(new_id.begin(), new_id.end());

  for (size_t i = 0; i < es.Size(); i++) {
    es[i].u = new_id[es[i].u];
    es[i].v = new_id[es[i].v];
  }

  vector<size_t> old_id(n);

  for (size_t i = 0; i < n; i++) {
    old_id[new_id[i]] = i;
  }

  typedef std::map<size_t, FLOAT>::const_iterator IterType;

  Upper g(es);

  vector<std::map<size_t, FLOAT>>& neighbor_map = g.neighbor_map;
  vector<EliminatedVertex>& elims = cholesky_factor.elims;
  vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;

  std::mt19937 rng(std::random_device{}());
  vector<size_t> neighbors;
  vector<double> weights;

  elims.resize(n);
  for (size_t i = 0; i < n - 1; i++) {
    if (i % 1000 == 0) {
      cout << i << endl;
    }
    elims[i].v = i;
    elims[i].degree = 0;
    elims[i].first_arc = elim_arcs.size();
    neighbors.clear();
    weights.clear();
    for (IterType it = neighbor_map[i].begin();
         it != neighbor_map[i].end();
         ++it) {
      elims[i].degree += it->second;
      elim_arcs.push_back(ArcC(it->first, it->second));
      neighbors.push_back(it->first);
      weights.push_back(static_cast<double>(it->second));
      // neighbor_map[it->first].erase(u);
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
      size_t v1 = std::min(neighbors[a1], neighbors[a2]);
      size_t v2 = std::max(neighbors[a1], neighbors[a2]);
      neighbor_map[v1][v2] += c;
    }
    neighbor_map[i].clear();
  }
  elims[n - 1].first_arc = elim_arcs.size();

  cerr << elim_arcs.size() << '\n';

  for (size_t i = 0; i < elim_arcs.size(); i++) {
    elim_arcs[i].v = old_id[elim_arcs[i].v];
  }

  for (size_t i = 0; i < elims.size(); i++) {
    elims[i].v = old_id[elims[i].v];
  }
}

void SparseCholesky2(EdgeListC& es,
                     size_t k,
                     CholeskyFactor& cholesky_factor) {
  size_t n = es.n;
  cholesky_factor.n = n;
  vector<size_t> new_id(n);

  for (size_t i = 0; i < n; i++) {
    new_id[i] = i;
  }
  std::random_shuffle(new_id.begin(), new_id.end());

  vector<size_t> old_id(n);
  for (size_t i = 0; i < n; i++) {
    old_id[new_id[i]] = i;
  }

  for (size_t i = 0; i < es.Size(); i++) {
    es[i].u = new_id[es[i].u];
    es[i].v = new_id[es[i].v];
  }

  UpperTriangular g(es);

  vector<EliminatedVertex>& elims = cholesky_factor.elims;
  vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;
  vector<FLOAT> spa(n, 0);
  vector<bool> is_nonzero(n, false);
  vector<size_t> indices;

  struct {
    bool operator()(const size_t v, const ArcC& a) {
      return v < a.v;
    }
  } cmp;

  typedef vector<ArcC>::iterator IterType;
  indices.reserve(n);
  elims.resize(n);

  std::mt19937 rng(std::random_device{}());
  vector<vector<size_t>> prevs(n);

  for (size_t current = 0; current < n - 1; current++) {
    elims[current].v = current;
    elims[current].degree = 0;
    elims[current].first_arc = elim_arcs.size();
    /* load the row in the sparse accumulator */
    for (size_t i = g.first_arc[current]; i < g.first_arc[current + 1]; i++) {
      const ArcC& a = g.arcs[i];
      spa[a.v] = a.conductance;
      is_nonzero[a.v] = true;
      indices.push_back(a.v);
    }
    /* go through previously eliminated vertices */
    // for (size_t previous = 0; previous < current; previous++) {
    for (size_t i = 0; i < prevs[current].size(); i++) {
      size_t previous = prevs[current][i];
      const IterType first = elim_arcs.begin() + elims[previous].first_arc;
      const IterType last = elim_arcs.begin() + elims[previous + 1].first_arc;
      const IterType upper = std::upper_bound(first, last, current, cmp);
      if (upper == first || (upper - 1)->v != current) {
        /* not edge between previous and current */
        continue;
      }
      FLOAT w1 = (upper - 1)->conductance / k;
      size_t lol = elims[previous + 1].first_arc - elims[previous].first_arc;
      /* go through each neighbor of previous that came after current */
      for (IterType it = upper; it != last; ++it) {
        FLOAT w2 = it->conductance / k;
        FLOAT c = w1 * w2 / (w1 + w2);
        // FLOAT c = w1 * w2 / elims[previous].degree;
        std::binomial_distribution<unsigned> sample(
            k * lol,
            (1.0 / lol) * (it->conductance + (upper - 1)->conductance) / elims[previous].degree
        );
        size_t count = sample(rng);
        spa[it->v] += c * count;
        if (count != 0 && !is_nonzero[it->v]) {
          is_nonzero[it->v] = true;
          indices.push_back(it->v);
        }
      }
    }

    std::sort(indices.begin(), indices.end());
    for (size_t i = 0; i < indices.size(); i++) {
      elim_arcs.push_back(ArcC(indices[i], spa[indices[i]]));
      elims[current].degree += spa[indices[i]];
      prevs[indices[i]].push_back(current);
      spa[indices[i]] = 0;
      is_nonzero[indices[i]] = false;
    }

    indices.clear();
  }
  elims[n - 1].first_arc = elim_arcs.size();

  cerr << elim_arcs.size() << '\n';

  for (size_t i = 0; i < elim_arcs.size(); i++) {
    elim_arcs[i].v = old_id[elim_arcs[i].v];
  }

  for (size_t i = 0; i < elims.size(); i++) {
    elims[i].v = old_id[elims[i].v];
  }
}
