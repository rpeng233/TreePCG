#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include <cmath>
#include <iostream>
#include "common.h"
#include "linear_algebra.h"
#include "matrix.h"

template <typename MatrixType, typename Preconditioner>
class PCGSolver {
public:

  PCGSolver(const Preconditioner& p, const MatrixType& A_)
    : preconditioner(p), A(A_) { }

  PCGSolver(Preconditioner&& p, MatrixType&& A_)
    : preconditioner(std::move(p)), A(std::move(A_)) {}

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
    preconditioner.solve(r, d);   // Solve P * d = r
    delta_new = r * d;
    for (;;) {
      mv(1, A, d, 0, q, q);                      // q = A * d
      alpha = delta_new / (d * q);
      axpy(alpha, d, x, x);                   // x = alpha * d + x
      mv(-1, A, x, 1, b, r);                  // r = b - A * x
      res = sqrt(r * r);
      if (res < tol) return;
      preconditioner.solve(r, s); // Solve P * s = r
      delta_old = delta_new;
      delta_new = r * s;
      beta = delta_new / delta_old;
      axpy(beta, d, s, d);                    // d = s + beta * d
    }
  }

private:

  Preconditioner preconditioner;
  MatrixType A;
};

#endif
