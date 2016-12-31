#ifndef INCLUDE_TREE_SOLVER_H__
#define INCLUDE_TREE_SOLVER_H__

#include <tuple>
#include <vector>
#include "common.h"
#include "graph.h"

class ElimnatedLeaf {
 public:
  size_t v;
  size_t parent;
  FLOAT resistance;

  ElimnatedLeaf(size_t v_, size_t p, FLOAT r) {
    v = v_;
    parent = p;
    resistance = r;
  }
};

class TreeSolver {
 public:
  TreeSolver(const Tree& tree);
  void solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x);

 private:
  size_t n;
  size_t root;
  std::vector<ElimnatedLeaf> elims;
};

#endif  // INCLUDE_TREE_SOLVER_H_
