#include <iostream>
#include <random>
#include <vector>
#include "common.h"
#include "graph.h"
#include "matrix.h"
#include "IdentitySolver.h"
#include "PCGSolver.h"
#include "TreeSolver.h"

int main(void) {
  size_t n = 4096 - 1;
  TreePlusEdges tree(n);

  // an unweighted complete binary tree
  for (size_t i = 0; i * 2 + 2 < n; i++) {
    tree.setParent(i * 2 + 1, i, 1);
    tree.setParent(i * 2 + 2, i, 1);
  }

  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n);

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif(-5, 5);

  // random demand
  FLOAT sum = 0;
  for (size_t i = 0; i < n - 1; i++) {
    FLOAT tmp = unif(rng);
    sum += tmp;
    b[i] = tmp;
  }
  b[n - 1] = -sum;

  TreeSolver s(tree);
  s.solve(b, x);

  // r = b - A * x
  std::vector<FLOAT> r(n);
  mv(-1, tree, x, 1, b, r);

  std::cout << r * r << std::endl;

  return 0;
}

/*
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
*/
