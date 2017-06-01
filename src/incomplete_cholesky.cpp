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
                        CholeskyFactor& cholesky_factor) {
  size_t n = es.n;
  cholesky_factor.n = n;
  AdjacencyArray<ArcC> g(es);

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

  elim_arcs.reserve(es.Size());
  indices.reserve(n);
  elims.resize(n);
  vector<size_t> prevs;

  for (size_t current = 0; current < n - 1; current++) {
    elims[current].v = current;
    elims[current].degree = 0;
    elims[current].first_arc = elim_arcs.size();

    /* load the current row into the sparse accumulator */
    prevs.clear();
    for (size_t i = g.first_arc[current]; i < g.first_arc[current + 1]; i++) {
      const ArcC& a = g.arcs[i];
      if (a.v < current) {
        prevs.push_back(a.v);
        continue;
      }
      spa[a.v] = a.conductance;
      is_nonzero[a.v] = true;
      indices.push_back(a.v);
    }
    std::sort(prevs.begin(), prevs.end());

    /* go through each previously eliminated vertices */
    for (size_t i = 0; i < prevs.size(); i++) {
      size_t previous = prevs[i];
      IterType first = elim_arcs.begin() + elims[previous].first_arc;
      IterType last = elim_arcs.begin() + elims[previous + 1].first_arc;
      IterType upper = std::upper_bound(first, last, current, cmp);
      // IterType upper = first;
      // while (upper->v <= current) {
      //   ++upper;
      // }

      if (upper == first || (upper - 1)->v != current) {
        /* no edge between previous and current */
        continue;
      }

      FLOAT w1 = (upper - 1)->conductance;

      /* go through each neighbor of previous that came after current */
      for (IterType it = upper; it != last; ++it) {
        if (!is_nonzero[it->v]) {
          /* Don't introduce any fill */
          continue;
          // is_nonzero[it->v] = true;
          // indices.push_back(it->v);
        }
        FLOAT w2 = it->conductance;
        spa[it->v] += w1 * w2 / elims[previous].degree;
        // spa[it->v] += w1 * w2 / (w1 + w2);
      }
    }

    std::sort(indices.begin(), indices.end());
    for (size_t i = 0; i < indices.size(); i++) {
      elim_arcs.push_back(ArcC(indices[i], spa[indices[i]]));
      elims[current].degree += spa[indices[i]];
      spa[indices[i]] = 0;
      is_nonzero[indices[i]] = false;
    }

    indices.clear();
  }
  elims[n - 1].first_arc = elim_arcs.size();

  cerr << elim_arcs.size() << '\n';
}

void IncompleteCholesky(EdgeListC& es,
                        double droptol,
                        CholeskyFactor& cholesky_factor) {
  size_t n = es.n;
  cholesky_factor.n = n;
  AdjacencyArray<ArcC> g(es);

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

  elim_arcs.reserve(es.Size());
  indices.reserve(n);
  elims.resize(n);
  vector<vector<size_t>> prevs(n);

  for (size_t current = 0; current < n - 1; current++) {
    elims[current].v = current;
    elims[current].degree = 0;
    elims[current].first_arc = elim_arcs.size();

    /* load the current row into the sparse accumulator */
    // prevs.clear();
    for (size_t i = g.first_arc[current]; i < g.first_arc[current + 1]; i++) {
      const ArcC& a = g.arcs[i];
      if (a.v < current) {
        // prevs.push_back(a.v);
        continue;
      }
      spa[a.v] = a.conductance;
      is_nonzero[a.v] = true;
      indices.push_back(a.v);
    }
    // std::sort(prevs.begin(), prevs.end());

    /* go through each previously eliminated vertices */
    // for (size_t previous = 0; previous < current; previous++) {
    for (size_t i = 0; i < prevs[current].size(); i++) {
      size_t previous = prevs[current][i];
      IterType first = elim_arcs.begin() + elims[previous].first_arc;
      IterType last = elim_arcs.begin() + elims[previous + 1].first_arc;
      IterType upper = std::upper_bound(first, last, current, cmp);
      // IterType upper = first;
      // while (upper->v <= current) {
      //   ++upper;
      // }

      if (upper == first || (upper - 1)->v != current) {
        /* no edge between previous and current */
        continue;
      }

      FLOAT w1 = (upper - 1)->conductance;

      /* go through each neighbor of previous that came after current */
      for (IterType it = upper; it != last; ++it) {
        FLOAT w2 = it->conductance;
        spa[it->v] += w1 * w2 / elims[previous].degree;
        if (!is_nonzero[it->v]) {
          is_nonzero[it->v] = true;
          indices.push_back(it->v);
        }
      }
    }

    std::sort(indices.begin(), indices.end());
    double sum = 0;
    for (size_t i = 0; i < indices.size(); i++) {
      sum += spa[indices[i]];
    }

    for (size_t i = 0; i < indices.size(); i++) {
      if (spa[indices[i]] < droptol * sum) {
        continue;
      }
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
}
