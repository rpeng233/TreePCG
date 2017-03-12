#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include "cholesky.h"
#include "graph.h"
#include "incomplete_cholesky.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void IncompleteCholesky(EdgeListC& es,
                        size_t k,
                        CholeskyFactor& cholesky_factor) {
  size_t n = es.n;
  cholesky_factor.n = n;
  UpperTriangular g(es);

  struct {
    bool operator()(const size_t v, const ArcC& a) {
      return v < a.v;
    }
  } cmp;

  typedef vector<ArcC>::iterator IterType;
  vector<EliminatedVertex>& elims = cholesky_factor.elims;
  vector<ArcC>& elim_arcs = cholesky_factor.elim_arcs;
  vector<FLOAT> spa(n, 0);
  vector<bool> is_nonzero(n, false);
  vector<size_t> indices;

  indices.reserve(n);
  elims.resize(n);

  for (size_t current = 0; current < n - 1; current++) {
    elims[current].v = current;
    elims[current].degree = 0;
    elims[current].first_arc = elim_arcs.size();

    /* load the current row into the sparse accumulator */
    for (size_t i = g.first_arc[current]; i < g.first_arc[current + 1]; i++) {
      const ArcC& a = g.arcs[i];
      spa[a.v] = a.conductance;
      is_nonzero[a.v] = true;
      indices.push_back(a.v);
    }

    /* go through each previously eliminated vertices */
    for (size_t previous = 0; previous < current; previous++) {
      IterType first = elim_arcs.begin() + elims[previous].first_arc;
      IterType last = elim_arcs.begin() + elims[previous + 1].first_arc;
      IterType upper = std::upper_bound(first, last, current, cmp);

      if (upper == first || (upper - 1)->v != current) {
        /* no edge between previous and current */
        continue;
      }

      FLOAT w1 = (upper - 1)->conductance;
      size_t lol = elims[previous + 1].first_arc - elims[previous].first_arc;

      /* go through each neighbor of previous that came after current */
      for (IterType it = upper; it != last; ++it) {
        if (!is_nonzero[it->v]) {
          /* Don't introduce any fill */
          continue;
        }
        FLOAT w2 = it->conductance;
        spa[it->v] += w1 * w2 / elims[previous].degree;
      }
    }

    std::sort(indices.begin(), indices.end());
    for (size_t i = 0; i < indices.size(); i++) {
      elim_arcs.push_back(ArcC(indices[i], spa[indices[i]]));
      elims[current].degree += spa[indices[i]];
    }

    indices.clear();
    std::fill(spa.begin(), spa.end(), 0);
    std::fill(is_nonzero.begin(), is_nonzero.end(), false);
  }
  elims[n - 1].first_arc = elim_arcs.size();

  cerr << elim_arcs.size() << '\n';
}
