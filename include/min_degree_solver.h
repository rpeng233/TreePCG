#ifndef INCLUDE_MIN_DEGREE_SOLVER_H_
#define INCLUDE_MIN_DEGREE_SOLVER_H_

#include "common.h"
#include "graph.h"

class EliminatedVertex {
public:
  size_t v;
  FLOAT degree;
  std::vector<Arc> neighbors;
};

class MinDegreeSolver {
public:
  MinDegreeSolver(Graph3& graph);
  void solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x) const;

private:
  size_t n;
  std::vector<EliminatedVertex> elims;
};

#endif // INCLUDE_MIN_DEGREE_SOLVER_H_
