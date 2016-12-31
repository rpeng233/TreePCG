#include <iostream>
#include <random>
#include <vector>
#include "common.h"
#include "graph.h"
#include "matrix.h"
#include "identity_solver.h"
#include "pcg_solver.h"
#include "tree_solver.h"

void bst(std::mt19937& rng) {
  size_t n = 65535;
  Tree tree(n);
  std::uniform_real_distribution<> weight(1, 10);

  // a complete binary tree
  for (size_t i = 0; i * 2 + 2 < n; i++) {
    tree.setParent(i * 2 + 1, i, weight(rng));
    tree.setParent(i * 2 + 2, i, weight(rng));
  }

  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n);

  // random demand
  std::uniform_real_distribution<> demand(-5, 5);
  FLOAT sum = 0;
  for (size_t i = 0; i < n - 1; i++) {
    FLOAT tmp = demand(rng);
    sum += tmp;
    b[i] = tmp;
  }
  b[n - 1] = -sum;

  TreeSolver s(tree);
  s.solve(b, x);

  // r = b - A * x
  std::vector<FLOAT> r(n);
  mv(-1, tree, x, 1, b, r);

  std::cout << "bst\n";
  std::cout << r * r << std::endl;
}

void pcg(std::mt19937& rng) {
  size_t n = 5;
  Matrix A(n, n);
  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n, 0);

  std::uniform_real_distribution<> weight(1, 10);
  for (size_t i = 0; i < n; i++) {
    A.addNonZero(i, i, weight(rng));
    b[i] = 1;
  }

  PCGSolver<Matrix, IdentitySolver> s(IdentitySolver(), A);

  s.solve(b, x);

  std::vector<FLOAT> r(n);

  mv(-1, A, x, 1, b, r);

  std::cout << "pcg\n";
  std::cout << r * r << std::endl;
}

int main(void) {
  std::mt19937 rng(std::random_device{}());

  bst(rng);
  pcg(rng);
}
