#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include "common.h"
#include "graph.h"
#include "identity_solver.h"
#include "low_stretch_tree.h"
#include "matrix.h"
#include "min_degree_solver.h"
#include "pcg_solver.h"
#include "tree_solver.h"

using std::cout;
using std::cerr;
using std::endl;

template <class Distribution, class RandomEngine>
class RNG {
  Distribution& dist;
  RandomEngine& rng;
public:
  RNG(Distribution& dist_, RandomEngine& rng_) : dist(dist_), rng(rng_) { }

  typename Distribution::result_type operator()() {
    return dist(rng);
  }
};

void min_degree(std::mt19937& rng) {
  size_t k = 70;
  size_t n = k * k;

  EdgeListR es;

  torus(k, k, es);
  // line(n, es);
  Graph3 g(es);

  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n);
  std::uniform_real_distribution<> demand(-10, 10);
  FLOAT sum = 0;
  for (size_t i = 0; i < n - 1; i++) {
    FLOAT tmp = demand(rng);
    sum += tmp;
    b[i] = tmp;
  }
  b[n - 1] = -sum;

  cerr << 1 << endl;
  MinDegreeSolver s(g, 1);
  cerr << 2 << endl;
  s.solve(b, x);
  cerr << 3 << endl;

  std::vector<FLOAT> r(n);
  mv(-1, es, x, 1, b, r);

  std::cout << "min_degree\n";
  std::cout << MYSQRT(r * r) << std::endl;
}

void resistance_vs_conductance(std::mt19937& rng) {
  size_t k = 300;
  size_t n = k * k;

  EdgeListR es;
  std::uniform_real_distribution<> unif_1_100(1, 10);
  RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);
  torus(k, k, es, random_resistance);

  Tree t;

  {
    Graph2 g(es);
    t = DijkstraTree(g, 0);
  }


  EdgeListC es2(es);

  TreeSolver preconditioner(t);
  PCGSolver<EdgeListR, TreeSolver> s(es, preconditioner);
  PCGSolver<EdgeListC, TreeSolver> s2(es2, preconditioner);

  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n);
  std::uniform_real_distribution<> demand(-10, 10);
  FLOAT sum = 0;
  for (size_t i = 0; i < n - 1; i++) {
    FLOAT tmp = demand(rng);
    sum += tmp;
    b[i] = tmp;
  }
  b[n - 1] = -sum;

  std::vector<FLOAT> r(n);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> duration;

  std::cout << "resistance_vs_conductance\n";
  for (size_t i = 0; i < 3; i++) {
    for (auto&f : x) {
      f = 0;
    }
    start = std::chrono::system_clock::now();
    s.solve(b, x);
    end = std::chrono::system_clock::now();
    duration = end - start;
    cout << "resistance: " << duration.count() << "s" << endl;

    mv(-1, es, x, 1, b, r);

    std::cout << MYSQRT(r * r) << std::endl;

    for (auto&f : x) {
      f = 0;
    }
    start = std::chrono::system_clock::now();
    s2.solve(b, x);
    end = std::chrono::system_clock::now();
    duration = end - start;
    cout << "conductance: " << duration.count() << "s" << endl;

    mv(-1, es, x, 1, b, r);

    std::cout << MYSQRT(r * r) << std::endl;
  }
}

void bst(std::mt19937& rng) {
  size_t n = 65535;
  Tree tree(n);
  std::uniform_real_distribution<> weight(1, 100);

  // a complete binary tree
  for (size_t i = 0; i * 2 + 2 < n; i++) {
    tree.SetParent(i * 2 + 1, i, weight(rng));
    tree.SetParent(i * 2 + 2, i, weight(rng));
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
  std::cout << MYSQRT(r * r) << std::endl;
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

  PCGSolver<Matrix, IdentitySolver> s(A, IdentitySolver());

  s.solve(b, x);

  std::vector<FLOAT> r(n);

  mv(-1, A, x, 1, b, r);

  std::cout << "pcg\n";
  std::cout << MYSQRT(r * r) << std::endl;
}

int main(void) {
  std::mt19937 rng(std::random_device{}());

  // bst(rng);
  // pcg(rng);
  resistance_vs_conductance(rng);
  // min_degree(rng);

  return 0;
}
