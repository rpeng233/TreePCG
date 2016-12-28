#ifndef IDENTITYSOLVER_H
#define IDENTITYSOLVER_H

#include <iostream>
#include <vector>
#include "AbstractSolver.h"

class IdentitySolver {
public:
  void solve(
      const std::vector<FLOAT>& b,
      std::vector<FLOAT>& x,
      FLOAT tol=1e-6,
      int maxit=-1
  ) const {
    // x = b;
    for (size_t i = 0; i < x.size(); i++) {
      x[i] = b[i];
    }
  }
};

#endif
