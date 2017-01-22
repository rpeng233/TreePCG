#ifndef INCLUDE_AUG_TREE_SOLVER_H
#define INCLUDE_AUG_TREE_SOLVER_H

#include "min_degree_solver.h"

class AugTreeSolver {
private:
  MinDegreeSolver solver;

public:
  /*
  AugTreeSolver(const TreePlusEdgesR& tree) {
    AdjacencyMap g(tree);
    solver = MinDegreeSolver(g, 1);
  }
  */
  void solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x) const {
    solver.solve(b, x);
  }
};

#endif // INCLUDE_AUG_TREE_SOLVER_H
