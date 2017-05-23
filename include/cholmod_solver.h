#ifndef INCLUDE_CHOLMOD_SOLVER_H_
#define INCLUDE_CHOLMOD_SOLVER_H_

#include "cholmod.h"
#include "graph.h"

class CholmodSolver {
private:
  size_t n;
  cholmod_common common_;
  cholmod_common *common;
  cholmod_factor *factor;
  cholmod_dense *b = NULL;
  cholmod_dense *x = NULL;
  cholmod_dense *Y = NULL;
  cholmod_dense *E = NULL;

public:
  CholmodSolver(const EdgeListC& es, cholmod_common *c=nullptr);

  void Solve(const std::vector<FLOAT>& b_, std::vector<FLOAT>& x_);

  ~CholmodSolver() {
    cholmod_free_factor(&factor, common);
    cholmod_free_dense(&b, common);
    cholmod_free_dense(&x, common);
    cholmod_free_dense(&Y, common);
    cholmod_free_dense(&E, common);
    cholmod_finish(common);
  }
};

#endif  // INCLUDE_CHOLMOD_SOLVER_H_
