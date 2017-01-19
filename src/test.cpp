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

  void tic(const std::string& m="") {
    msg = m;
    start = std::chrono::steady_clock::now();
  }

  void toc() {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << msg << duration.count() << "s" << endl;
  }
};

Timer timer;

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

void aug_tree_pcg(const EdgeListR& es, const vector<FLOAT>& b, size_t k) {
  cout << "===== aug-tree PCG =====\n";
  TreeR t;

  EdgeListR recc;
  recursive_c(k, k, recc);

  {
    Graph2 g(recc);
    t = DijkstraTree<TreeR>(g, 0);
  }

  EdgeListR otes;
  otes.n = es.n;

  for (size_t i = 0; i < es.edges.size(); i++) {
    const EdgeR& e = es.edges[i];
    if (t.vertices[e.u].parent == e.v || t.vertices[e.v].parent == e.u) {
      continue;
    }
    otes.AddEdge(e.u, e.v, e.resistance);
  }

  vector<double> strs(otes.edges.size());
  Graph3 g(t);

  timer.tic("computing stretch and adding edges: ");
  ComputeStretch(t.vertices, otes.edges, strs);

  StretchGreater less(strs);
  vector<size_t> indices(otes.edges.size());

  for (size_t i = 0; i < indices.size(); i++) {
    indices[i] = i;
  }

  std::sort(indices.begin(), indices.end(), less);

  for (size_t i = 0; i < 5 * k; i++) {
    g.AddEdgeR(otes.edges[indices[i]]);
  }
  timer.toc();

  timer.tic("eliminating aug tree: ");
  MinDegreeSolver precon(g);
  timer.toc();

  EdgeListC es2(es);
  PCGSolver<EdgeListC, MinDegreeSolver> s(es2, precon);

  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("aug tree pcg: ");
  s.solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;

  return;
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

/*
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
*/

void min_degree(const EdgeListR& es, const vector<FLOAT>& b) {
  cout << "===== min degree =====\n";
  Graph3 g(es);
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("factorizing: ");
  MinDegreeSolver s(g, 1);
  timer.toc();

  timer.tic("solving: ");
  s.solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);
  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}

void resistance_vs_conductance(const EdgeListR& es, const vector<FLOAT>& b) {
  cout << "===== resistance vs conductance =====\n";
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  EdgeListC es2(es);

  PCGSolver<EdgeListR, IdentitySolver> s(es, IdentitySolver());
  PCGSolver<EdgeListC, IdentitySolver> s2(es2, IdentitySolver());

  for (size_t i = 0; i < 2; i++) {
    for (auto&f : x) {
      f = 0;
    }
    timer.tic("resistance: ");
    s.solve(b, x);
    timer.toc();

    mv(-1, es, x, 1, b, r);
    cout << MYSQRT(r * r) / MYSQRT(b * b) << endl;

    for (auto&f : x) {
      f = 0;
    }
    timer.tic("conductance: ");
    s2.solve(b, x);
    timer.toc();

    mv(-1, es, x, 1, b, r);
    cout << MYSQRT(r * r) / MYSQRT(b * b) << endl;
  }
}

void bst(std::mt19937& rng) {
  cout << "===== Solving random complete BST =====\n";
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

void pcg(const EdgeListR& es, const vector<FLOAT>& b) {
  cout << "===== PCG =====\n";
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  PCGSolver<EdgeListR, IdentitySolver> s(es, IdentitySolver());

  timer.tic();
  s.solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}

int main(void) {
  size_t k = 500;
  size_t n = k * k;

  EdgeListR unweighted_torus;
  EdgeListR weighted_torus;

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> uniform(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937>
    random_resistance(uniform, rng);

  torus(k, k, unweighted_torus);
  torus(k, k, weighted_torus, random_resistance);

  std::vector<FLOAT> unweighted_b(n);
  std::vector<FLOAT> weighted_b(n);
  std::vector<FLOAT> x(n);

  for (auto& f : x) {
    f = uniform(rng);
  }
  mv(1, weighted_torus, x, 0, x, weighted_b);

  for (auto& f : x) {
    f = uniform(rng);
  }
  mv(1, unweighted_torus, x, 0, x, unweighted_b);

  // bst(rng);
  // pcg(unweighted_torus, unweighted_b);
  // resistance_vs_conductance(weighted_torus, weighted_b);
  // min_degree(weighted_torus, weighted_b);
  aug_tree_pcg(unweighted_torus, unweighted_b, k);

  return 0;
}
