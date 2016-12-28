#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <cmath>
#include <iostream>
#include <type_traits>
#include "common.h"
#include "linalg.h"
#include "matrix.h"

template <typename Preconditioner>
class PCGSolver {

public:

  PCGSolver(Preconditioner p, Matrix A_) {
    preconditioner = p;
    A = A_;
  }


  void solve(
      const std::vector<FLOAT>& b,
      std::vector<FLOAT>& x,
      FLOAT tol=1e-6,
      int maxit=-1
  ) const {
    size_t n = A.n;
    vector<FLOAT> r(n);
    vector<FLOAT> q(n);
    vector<FLOAT> d(n, 0);
    vector<FLOAT> s(n, 0);
    FLOAT delta_old;
    FLOAT delta_new;
    FLOAT alpha;
    FLOAT beta;

    // size_t i = 0;
    mv(-1, A, x, 1, b, r);                    // r = b - A * x
    FLOAT res = (r * r);
    preconditioner.solve(r, d, tol, maxit);   // Solve P * d = r
    delta_new = r * d;
    for (;;) {
      mv(1, A, d, 0, q);                      // q = A * d
      alpha = delta_new / (d * q);
      axpy(alpha, d, x);                      // x = alpha * d + x
      mv(-1, A, x, 1, b, r);                  // r = b - A * x
      res = sqrt(r * r);
      if (res < tol) return;
      preconditioner.solve(r, s, tol, maxit); // Solve P * s = r
      delta_old = delta_new;
      delta_new = r * s;
      beta = delta_new / delta_old;
      axpby(1, s, beta, d);                   // d = s + beta * d
    }
  }

private:

  Preconditioner preconditioner;
  Matrix A;
};

#endif
