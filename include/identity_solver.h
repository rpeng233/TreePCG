#ifndef INCLUDE_IDENTITY_SOLVER_H_
#define INCLUDE_IDENTITY_SOLVER_H_

#include <iostream>
#include <vector>
#include "common.h"

class IdentitySolver {
 public:
  void Solve(
      const std::vector<FLOAT>& b,
      std::vector<FLOAT>& x,
      FLOAT tol = 1e-6,
      int maxit = -1) const {
    // x = b;
    for (size_t i = 0; i < x.size(); i++) {
      x[i] = b[i];
    }
  }
};

#endif  // INCLUDE_IDENTITY_SOLVER_H_
