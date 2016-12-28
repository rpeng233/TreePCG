#include <iostream>
#include <vector>
#include "common.h"
#include "matrix.h"
#include "IdentitySolver.h"
#include "PCGSolver.h"

int main(void) {
  size_t n = 5;
  Matrix A(n, n);
  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n, 0);

  for (size_t i = 0; i < n; i++) {
    A.addNonZero(i, i, i + 1);
    b[i] = 1;
  }

  PCGSolver<IdentitySolver> s(IdentitySolver(), A);

  s.solve(b, x);

  for (const FLOAT& f : x) {
    std::cout << f << '\n';
  }
  return 0;
}
