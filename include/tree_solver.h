#ifndef INCLUDE_TREE_SOLVER_H__
#define INCLUDE_TREE_SOLVER_H__

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

struct TreeSolverVertex {
  size_t parent;
  size_t degree; // right now the degree doesn't count the parent
  FLOAT parent_resistance;
  bool eliminated;

  TreeSolverVertex()
  : parent(0), degree(0), parent_resistance(0.0), eliminated(false)
  { }

  FLOAT ParentResistance() const {
    return parent_resistance;
  }

  FLOAT ParentConductance() const {
    return 1.0 / parent_resistance;
  }

  void SetParentR(size_t p, FLOAT r) {
    parent = p;
    parent_resistance = r;
  }

  void SetParentC(size_t p, FLOAT c) {
    parent = p;
    parent_resistance = 1.0 / c;
  }
};

template <>
inline void Tree<TreeSolverVertex>::SetParentR(size_t v, size_t p, FLOAT r) {
  size_t old_p = vertices[v].parent;
  if (old_p != v) {
    vertices[old_p].degree--;
  }
  vertices[v].parent = p;
  vertices[v].parent_resistance = r;
  vertices[p].degree++;
}


class TreeSolver {
 public:
  TreeSolver(const Tree<TreeSolverVertex>& tree);
  void Solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x);

 private:
  size_t n;
  size_t root;
  std::vector<ElimnatedLeaf> elims;
};

#endif  // INCLUDE_TREE_SOLVER_H_
