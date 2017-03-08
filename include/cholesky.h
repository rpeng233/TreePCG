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

struct UpperTriangular {
  size_t n;
  std::vector<ArcC> arcs;
  std::vector<size_t> first_arc;

  UpperTriangular() {
    n = 0;
  }

  UpperTriangular(const EdgeListC& es) {
    BuildGraph(es);
  }

  void BuildGraph(const EdgeListC& es) {
    n = es.n;
    size_t m = es.edges.size();

    arcs.resize(m);
    first_arc.resize(n + 1);
    std::vector<size_t> degrees(n);

    for (typename std::vector<EdgeC>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      degrees[std::min(it->u, it->v)]++;
    }

    first_arc[0] = 0;
    for (size_t i = 0; i < n - 1; i++) {
      first_arc[i + 1] = degrees[i] + first_arc[i];
    }
    first_arc[n] = m;

    size_t tmp_index;
    for (typename std::vector<EdgeC>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      size_t u = std::min(it->u, it->v);
      size_t v = std::max(it->u, it->v);
      arcs[first_arc[u]].v = v;
      arcs[first_arc[u]].conductance = it->conductance;
      first_arc[u]++;
    }

    first_arc[0] = 0;
    for (size_t i = 0; i < n - 1; i++) {
      first_arc[i + 1] = degrees[i] + first_arc[i];
    }
    first_arc[n] = m;
  }

  void FreeMemory() {
    n = 0;
    arcs = std::vector<ArcC>();
    first_arc = std::vector<size_t>();
  }
};

#endif  // INCLUDE_CHOLESKY_SOLVER_H_
