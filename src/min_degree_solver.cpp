#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>
#include "min_degree_solver.h"

using std::vector;

class DegreeLess {
  const std::vector<std::map<size_t, FLOAT> > *neighbor_map;
public:
  DegreeLess(const std::vector<std::map<size_t, FLOAT> > *neighbor_map) {
    assert(neighbor_map != NULL);
    this->neighbor_map = neighbor_map;
  }

  bool operator() (size_t u, size_t v) const {
    return (*neighbor_map)[u].size() < (*neighbor_map)[v].size();
  }
};

static inline
vector<FLOAT> eliminate_rhs(const vector<EliminatedVertex>& elims,
                            const vector<FLOAT>& rhs_) {
  vector<FLOAT> rhs(rhs_);
  vector<FLOAT> rhs_elims(elims.size());

  for (size_t i = 0; i < elims.size(); i++) {
    rhs_elims[i] = rhs[elims[i].v];
    for (vector<Arc>::const_iterator it = elims[i].neighbors.begin();
         it != elims[i].neighbors.end();
         ++it) {
      rhs[it->v] += rhs[elims[i].v] / (elims[i].degree * it->resistance);
    }
    rhs[elims[i].v] = 0;
  }

  return rhs_elims;
}

static inline
void back_substitution(const vector<EliminatedVertex>& elims,
                       const vector<FLOAT>& rhs_elims,
                       vector<FLOAT>& x) {
  for (size_t i = elims.size(); i > 0; i--) {
    const EliminatedVertex& e = elims[i - 1];
    x[e.v] = rhs_elims[i - 1];
    for (vector<Arc>::const_iterator it = e.neighbors.begin();
         it != e.neighbors.end();
         ++it) {
      x[e.v] += x[it->v] / it->resistance;
    }
    x[e.v] /= e.degree;
  }
}

MinDegreeSolver::MinDegreeSolver(Graph3& graph) {
  n = graph.n;
  vector<size_t> vs(n);
  DegreeLess compare(&graph.neighbor_map);

  for (size_t i = 0; i < n; i++) {
    vs[i] = i;
  }

  typedef std::map<size_t, FLOAT>::const_iterator IterType;

  vector<std::map<size_t, FLOAT> >& neighbor_map = graph.neighbor_map;
  elims.resize(n - 1);

  for (size_t i = 0; i < n - 1; i++) {
    std::vector<size_t>::iterator min_it =
      std::min_element(vs.begin(), vs.end() - i, compare);
    size_t u = *min_it;
    elims[i].v = u;
    *min_it = vs[n - 1 - i];
    FLOAT degree = 0;
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      degree += 1 / it->second;
      elims[i].neighbors.push_back(Arc(it->first, it->second));
    }
    elims[i].degree = degree;
    for (IterType it1 = neighbor_map[u].begin();
         it1 != neighbor_map[u].end();
         ++it1) {
      for (IterType it2 = it1;
           it2 != neighbor_map[u].end();
           ++it2) {
        if (it1 == it2) continue;
        size_t v1 = it1->first;
        size_t v2 = it2->first;
        const FLOAT& r1 = it1->second;
        const FLOAT& r2 = it2->second;
        std::map<size_t, FLOAT>::iterator r12 = neighbor_map[v1].find(v2);
        FLOAT new_r;
        if (r12 == neighbor_map[v1].end()) {
          new_r = degree * r1 * r2;
        } else {
          const FLOAT& r = neighbor_map[v1][v2];
          // assert(neighbor_map[v1][v2] == neighbor_map[v2][v1]);
          new_r = 1 / (1 / r + 1 / (degree * r1 * r2));
        }
        neighbor_map[v1][v2] = new_r;
        neighbor_map[v2][v1] = new_r;
      }
    }
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      neighbor_map[it->first].erase(u);
    }
  }
}

void MinDegreeSolver::solve(const vector<FLOAT>& b, vector<FLOAT>& x) const {
  assert(b.size() == x.size());
  assert(b.size() == n);

  vector<FLOAT> rhs_elims = eliminate_rhs(elims, b);
  back_substitution(elims, rhs_elims, x);

  FLOAT sum = 0;
  for (size_t i = 0; i < n; i++) {
    sum += x[i];
  }
  sum /= n;
  for (size_t i = 0; i < n; i++) {
    x[i] -= sum;
  }
}
