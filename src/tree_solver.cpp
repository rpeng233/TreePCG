#include <iostream>
#include <utility>
#include "tree_solver.h"

using std::vector;
using std::pair;

static inline
pair<size_t, vector<ElimnatedLeaf> > eliminate_leaves(vector<TreeSolverVertex>& vs) {
  vector<ElimnatedLeaf> elims;
  size_t root;
  size_t found_root = 0;

  elims.reserve(vs.size() - 1);
  for (size_t i = 0; i < vs.size(); i++) {
    if (vs[i].parent == i) {
      root = i;
      found_root++;
      continue;
    }
    if (vs[i].degree != 0 || vs[i].eliminated) {
      continue;
    }
    size_t cur = i;
    do {
      size_t p = vs[cur].parent;
      assert(p != cur);
      assert(vs[p].eliminated == false);
      elims.push_back(ElimnatedLeaf(cur, p, vs[cur].parent_resistance));
      vs[p].degree--;
      vs[cur].eliminated = true;
      cur = p;
    } while (vs[cur].degree == 0 && vs[cur].parent != cur);
  }

  assert(found_root == 1);

  return std::make_pair(root, elims);
}

static inline
void ForwardSubstitution(const vector<ElimnatedLeaf>& elims,
                         const vector<FLOAT>& b_,
                         vector<FLOAT> &y) {
  // vector<FLOAT> b(b_);
  y = b_;

  for (size_t i = 0; i < elims.size(); i++) {
    // y[i] = b[elims[i].v];
    y[elims[i].parent] += y[elims[i].v];
  }
}

static inline
void BackSubstitution(const vector<ElimnatedLeaf>& elims,
                      const vector<FLOAT>& y,
                      vector<FLOAT>& x) {
  for (size_t i = elims.size(); i > 0; i--) {
    const ElimnatedLeaf& e = elims[i - 1];
    x[e.v] = y[e.v] * e.resistance + x[e.parent];
  }
}

TreeSolver::TreeSolver(const Tree<TreeSolverVertex>& tree) {
  n = tree.n;
  vector<TreeSolverVertex> vs(tree.vertices);
  pair<size_t, vector<ElimnatedLeaf> > p = eliminate_leaves(vs);
  root = p.first;
  elims = p.second;
  assert(elims.size() == tree.n - 1);
}

void TreeSolver::Solve(const vector<FLOAT>& b, vector<FLOAT>& x) {
  assert(b.size() == x.size());
  assert(b.size() == n);

  vector<FLOAT> y(n);
  ForwardSubstitution(elims, b, y);
  x[root] = 0;
  BackSubstitution(elims, y, x);

  FLOAT sum = 0;
  for (size_t i = 0; i < n; i++) {
    sum += x[i];
  }
  sum /= n;
  for (size_t i = 0; i < n; i++) {
    x[i] -= sum;
  }
}
