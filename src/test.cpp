#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include "aug_tree_solver.h"
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
using std::vector;

class Timer {
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::string msg;

public:
  Timer() {
    msg = "Timer uninitialized";
    start = std::chrono::steady_clock::now();
  }

  void tic(const std::string& m) {
    msg = m;
    start = std::chrono::steady_clock::now();
  }

  void toc() {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << msg << ": " << duration.count() << "s" << endl;
  }
};

template <typename Distribution, typename RandomEngine>
class RNG {
  Distribution& dist;
  RandomEngine& rng;
public:
  RNG(Distribution& dist_, RandomEngine& rng_) : dist(dist_), rng(rng_) { }

  typename Distribution::result_type operator()() {
    return dist(rng);
  }
};

class StretchGreater {
  const vector<double>& strs;
public:
  StretchGreater(const vector<double>& strs_) : strs(strs_) { }

  bool operator() (size_t i, size_t j) const {
    return strs[i] > strs[j];
  }
};

void aug_tree_pcg(std::mt19937& rng) {
  Timer timer;

  size_t k = 500;
  size_t n = k * k;

  EdgeListR es;

  std::uniform_real_distribution<> unif_1_100(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);

  torus(k, k, es, random_resistance);

  TreeR t;

  {
    Graph2 g(es);
    timer.tic("dijkstra");
    t = DijkstraTree<TreeR>(g, 0);
    timer.toc();
  }

  EdgeListR otes;
  otes.n = n;

  for (size_t i = 0; i < es.edges.size(); i++) {
    const EdgeR& e = es.edges[i];
    if (t.vertices[e.u].parent == e.v || t.vertices[e.v].parent == e.u) {
      continue;
    }
    otes.AddEdge(e.u, e.v, e.resistance);
  }

  vector<double> strs(otes.edges.size());
  Graph3 g(t);

  timer.tic("computing stretch and adding edges");
  ComputeStretch(t.vertices, otes.edges, strs);

  StretchGreater less(strs);
  vector<size_t> indices(otes.edges.size());

  for (size_t i = 0; i < indices.size(); i++) {
    indices[i] = i;
  }

  std::sort(indices.begin(), indices.end(), less);

  for (size_t i = 0; i < 3 * k; i++) {
    g.AddEdgeR(otes.edges[indices[i]]);
  }
  timer.toc();

  timer.tic("eliminating aug tree");
  MinDegreeSolver precon(g);
  timer.toc();

  PCGSolver<EdgeListR, MinDegreeSolver> s(es, precon);

  EdgeListC es2(es);
  std::vector<FLOAT> b(n);
  std::vector<FLOAT> x(n);
  std::uniform_real_distribution<> unif(0, 100);
  FLOAT sum = 0;
  for (auto& f : x) {
    f = unif(rng);
  }
  mv(-1, es, x, 0, b, b);

  for (auto& f : x) {
    f = 0;
  }

  timer.tic("aug tree pcg");
  s.solve(b, x);
  timer.toc();

  std::vector<FLOAT> r(n);
  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) << std::endl;

  TreeSolver ts(t);
  PCGSolver<EdgeListC, TreeSolver> s2(es2, ts);

  for (auto& f : x) {
    f = 0;
  }

  timer.tic("tree pcg");
  s2.solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) << std::endl;
}

void stretch(std::mt19937& rng) {
  size_t k = 4;
  size_t n = k * k;

  EdgeListR es;

  // std::uniform_real_distribution<> unif_1_100(1, 100);
  // RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);

  torus(k, k, es);

  TreeR t;

  {
    Graph2 g(es);
    t = DijkstraTree<TreeR>(g, 0);
  }

  vector<double> strs(es.edges.size());
  ComputeStretch(t.vertices, es.edges, strs);

  for (size_t i = 0; i < t.vertices.size(); i++) {
    cout << i << ' ' << t.vertices[i].parent << endl;
  }

  for (size_t i = 0; i < es.edges.size(); i++) {
    EdgeR& e = es.edges[i];
    cout << e.u << ' ' << e.v << ' ' << strs[i] << endl;
  }
}

void aug_tree(std::mt19937& rng) {
  size_t k = 50;
  size_t n = k * k;

  EdgeListR es;

  std::uniform_real_distribution<> unif_1_100(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);

  torus(k, k, es, random_resistance);

  TreePlusEdgesR t;

  {
    Graph2 g(es);
    t = DijkstraTree<TreePlusEdgesR>(g, 0);
  }

  for (size_t i = 0; i < es.edges.size(); i++) {
    const EdgeR& e = es.edges[i];
    if (t.vertices[e.u].parent == e.v || t.vertices[e.v].parent == e.u) {
      continue;
    }
    t.AddEdge(e.u, e.v, e.resistance);
  }


  AugTreeSolver s(t);

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

  s.solve(b, x);

  std::vector<FLOAT> r(n);
  mv(-1, es, x, 1, b, r);

  std::cout << "aug_tree\n";
  std::cout << MYSQRT(r * r) << std::endl;
}

void min_degree(std::mt19937& rng) {
  size_t k = 100;
  size_t n = k * k;

  EdgeListR es;
  std::uniform_real_distribution<> unif_1_100(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);

  torus(k, k, es, random_resistance);
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

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> duration;

  start = std::chrono::system_clock::now();
  MinDegreeSolver s(g, 1);
  end = std::chrono::system_clock::now();
  duration = end - start;
  cout << duration.count() << "s" << endl;

  s.solve(b, x);

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

  TreeR t;

  {
    Graph2 g(es);
    t = DijkstraTree<TreeR>(g, 0);
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
  TreeR tree(n);
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
  // resistance_vs_conductance(rng);
  // min_degree(rng);
  // aug_tree(rng);
  aug_tree_pcg(rng);
  // stretch(rng);

  return 0;
}
