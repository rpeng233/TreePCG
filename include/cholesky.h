#ifndef INCLUDE_CHOLESKY_SOLVER_H_
#define INCLUDE_CHOLESKY_SOLVER_H_

#include "common.h"
#include "graph.h"

struct EliminatedVertex {
  size_t v;
  size_t first_arc;
  FLOAT degree;
};

struct CholeskyFactor {
  size_t n;
  std::vector<EliminatedVertex> elims;
  std::vector<ArcC> elim_arcs;
};

void Cholesky(AdjacencyMap& graph, CholeskyFactor& chol);

class CholeskySolver {
public:
  CholeskyFactor cholesky_factor;

  CholeskySolver() { }

  CholeskySolver(AdjacencyMap& graph) {
    Cholesky(graph, cholesky_factor);
  }

  void Factorize(AdjacencyMap& graph) {
    Cholesky(graph, cholesky_factor);
  }

  void Solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x) const;

private:
  void ForwardSubstitution(const std::vector<FLOAT>& b,
                           std::vector<FLOAT>& y) const;
  void BackSubstitution(const std::vector<FLOAT>& y,
                        std::vector<FLOAT>& x) const;
};

#endif  // INCLUDE_CHOLESKY_SOLVER_H_
