#include "TreeSolver.h"
#include <tuple>
#include <utility>

using std::vector;

static std::pair<size_t, vector<ElimnatedLeaf>>
eliminate_leaves(vector<TreeVertex>& vs) {
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
    if (vs[i].children_count != 0 || vs[i].eliminated) {
      continue;
    }
    size_t cur = i;
    do {
      size_t p = vs[cur].parent;
      assert(p != cur);
      assert(vs[p].eliminated == false);
      elims.emplace_back(cur, p, vs[cur].parent_resistance);
      vs[p].children_count--;
      vs[cur].eliminated = true;
      cur = p;
    } while (vs[cur].children_count == 0 && vs[cur].parent != cur);
  }

  assert(found_root == 1);

  return std::make_pair(root, std::move(elims));
}

static vector<FLOAT> eliminate_rhs(const vector<ElimnatedLeaf>& elims,
                                   const vector<FLOAT>& rhs_) {
  vector<FLOAT> rhs(rhs_);
  vector<FLOAT> rhs_elims(elims.size());

  for (size_t i = 0; i < elims.size(); i++) {
    rhs_elims[i] = rhs[elims[i].v];
    rhs[elims[i].parent] += rhs[elims[i].v];
  }

  return rhs_elims;
}

static void back_substitution(const vector<ElimnatedLeaf>& elims,
                              const vector<FLOAT>& rhs_elims,
                              vector<FLOAT>& x) {
  for (size_t i = elims.size(); i > 0; i--) {
    auto& e = elims[i - 1];
    x[e.v] = rhs_elims[i - 1] * e.resistance + x[e.parent];
  }
}

TreeSolver::TreeSolver(const Tree& tree) {
  n = tree.n;
  vector<TreeVertex> vs(tree.vertices);
  std::tie(root, elims) = eliminate_leaves(vs);
  assert(elims.size() == tree.n - 1);
}

void TreeSolver::solve(const vector<FLOAT>& b, vector<FLOAT>& x) {
  assert(b.size() == x.size());
  assert(b.size() == n);

  vector<FLOAT> rhs_elims = eliminate_rhs(elims, b);
  x[root] = 0;
  back_substitution(elims, rhs_elims, x);

  FLOAT sum = 0;
  for (const auto& f : x) {
    sum += f;
  }
  sum /= n;
  for (auto& f : x) {
    f -= sum;
  }
}
