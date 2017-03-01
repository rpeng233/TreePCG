#ifndef INCLUDE_CHOLESKY_SOLVER_H_
#define INCLUDE_CHOLESKY_SOLVER_H_

#include "common.h"
#include "graph.h"

struct EliminatedVertex {
  size_t v;
  size_t first_arc;
  FLOAT degree;
};

class CholeskySolver {
public:
  CholeskySolver() {
    n = 0;
  }

  void Factorize(AdjacencyMap& graph);
  void Solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x) const;

  CholeskySolver(AdjacencyMap& graph) {
    Factorize(graph);
  }
  // CholeskySolver(AdjacencyMap& graph, int brute_force);

private:
  size_t n;
  std::vector<EliminatedVertex> elims;
  std::vector<ArcC> elim_arcs;

  void ForwardSubstitution(const std::vector<FLOAT>& rhs_,
                           std::vector<FLOAT>& y) const;
  void BackSubstitution(const std::vector<FLOAT>& rhs_elims,
                        std::vector<FLOAT>& x) const;
};

#endif  // INCLUDE_CHOLESKY_SOLVER_H_
