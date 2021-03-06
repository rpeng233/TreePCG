#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>
#include "binary_heap.h"
#include "cholesky.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

class DegreeGreater {
  const std::vector<std::map<size_t, FLOAT> >& neighbor_map;
public:
  DegreeGreater(const std::vector<std::map<size_t, FLOAT> >& neighbor_map_)
  : neighbor_map(neighbor_map_) {
  }

  bool operator() (size_t u, size_t v) const {
    return neighbor_map[u].size() > neighbor_map[v].size();
  }
};

/*
class DegreeLess {
  const std::vector<std::map<size_t, FLOAT> >& neighbor_map;
public:
  DegreeLess(const std::vector<std::map<size_t, FLOAT> >& neighbor_map_)
  : neighbor_map(neighbor_map_) {
  }

  bool operator() (size_t u, size_t v) const {
    return neighbor_map[u].size() < neighbor_map[v].size();
  }
};

CholeskySolver::CholeskySolver(AdjacencyMap& graph, int brute_force) {
  n = graph.n;
  vector<size_t> vs(n);
  DegreeLess compare(graph.neighbor_map);

  for (size_t i = 0; i < n; i++) {
    vs[i] = i;
  }

  typedef std::map<size_t, FLOAT>::const_iterator IterType;

  vector<std::map<size_t, FLOAT> >& neighbor_map = graph.neighbor_map;
  elims.resize(n);

  for (size_t i = 0; i < n - 1; i++) {
    std::vector<size_t>::iterator min_it =
      std::min_element(vs.begin(), vs.end() - i, compare);
    size_t u = *min_it;
    elims[i].v = u;
    *min_it = vs[n - 1 - i];
    FLOAT degree = 0;
    elims[i].first_arc = elim_arcs.size();
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      degree += it->second;
      elim_arcs.push_back(ArcC(it->first, it->second));
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
        const FLOAT& c1 = it1->second;
        const FLOAT& c2 = it2->second;
        FLOAT new_c = c1 * c2 / degree;
        neighbor_map[v1][v2] += new_c;
        neighbor_map[v2][v1] += new_c;
      }
    }
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      neighbor_map[it->first].erase(u);
    }
  }
  elims[n - 1].first_arc = elim_arcs.size();
}
*/

void Cholesky(AdjacencyMap& graph, CholeskyFactor& cholesky_factor) {
  size_t n = graph.n;
  cholesky_factor.n = n;
  DegreeGreater less(graph.neighbor_map);
  BinaryHeap<DegreeGreater> heap(n, &less);

  typedef std::map<size_t, FLOAT>::const_iterator IterType;
  vector<std::map<size_t, FLOAT> >& neighbor_map = graph.neighbor_map;
  vector<EliminatedVertex>& elims = cholesky_factor.elims;
  vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;

  elims.resize(n);

  for (size_t i = 0; i < n - 1; i++) {
    size_t u = heap.Top();
    heap.Pop();
    elims[i].v = u;
    elims[i].degree = 0;
    elims[i].first_arc = elim_arcs.size();
    for (IterType it = neighbor_map[u].begin();
         it != neighbor_map[u].end();
         ++it) {
      elims[i].degree += it->second;
      elim_arcs.push_back(ArcC(it->first, it->second));
      neighbor_map[it->first].erase(u);
      heap.BubbleUp(it->first);
    }
    for (IterType it1 = neighbor_map[u].begin();
         it1 != neighbor_map[u].end();
         ++it1) {
      for (IterType it2 = it1;
           it2 != neighbor_map[u].end();
           ++it2) {
        if (it1 == it2) continue;
        size_t v1 = it1->first;
        size_t v2 = it2->first;
        const FLOAT& c1 = it1->second;
        const FLOAT& c2 = it2->second;
        FLOAT new_c = c1 * c2 / elims[i].degree;
        neighbor_map[v1][v2] += new_c;
        heap.BubbleDown(v1);
        neighbor_map[v2][v1] += new_c;
        heap.BubbleDown(v2);
      }
    }
  }
  elims[n - 1].first_arc = elim_arcs.size();

  // for (size_t i = 0; i < n - 1; i++) {
  //   for (size_t j = elims[i].first_arc; j < elims[i + 1].first_arc; j++) {
  //     elim_arcs[j].conductance /= elims[i].degree;
  //   }
  // }
}

void CholeskySolver::ForwardSubstitution(const vector<FLOAT>& b,
                                         vector<FLOAT>& y) const {
  // vector<FLOAT> b(b_);
  const vector<EliminatedVertex>& elims = cholesky_factor.elims;
  const vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;

  y = b;
  for (size_t i = 0; i < elims.size() - 1; i++) {
    FLOAT tmp = y[elims[i].v] / elims[i].degree;
    for (size_t j = elims[i].first_arc; j < elims[i + 1].first_arc; j++) {
      y[elim_arcs[j].v] += tmp * elim_arcs[j].conductance;
    }
  }
}

void CholeskySolver::BackSubstitution(const vector<FLOAT>& y,
                                       vector<FLOAT>& x) const {
  const vector<EliminatedVertex>& elims = cholesky_factor.elims;
  const vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;

  for (size_t i = elims.size() - 1; i > 0; i--) {
    const EliminatedVertex& e = elims[i - 1];
    // x[e.v] = y[i - 1];
    x[e.v] = y[e.v];
    for (size_t j = e.first_arc; j < elims[i].first_arc; j++) {
      x[e.v] += x[elim_arcs[j].v] * elim_arcs[j].conductance;
    }
    x[e.v] /= e.degree;
  }
}

void CholeskySolver::Solve(const vector<FLOAT>& b, vector<FLOAT>& x) const {
  assert(b.size() == x.size());
  assert(b.size() == cholesky_factor.n);

  vector<FLOAT> y(cholesky_factor.n);
  ForwardSubstitution(b, y);
  BackSubstitution(y, x);

  // FLOAT sum = 0;
  // for (size_t i = 0; i < n; i++) {
  //   sum += x[i];
  // }
  // sum /= n;
  // for (size_t i = 0; i < n; i++) {
  //   x[i] -= sum;
  // }
}
