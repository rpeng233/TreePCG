#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "akpw.h"
#include "aug_tree_precon.h"
#include "cholesky.h"
#include "common.h"
#include "graph.h"
#include "identity_solver.h"
#include "matrix.h"
#include "pcg_solver.h"
#include "sparse_cholesky.h"
#include "stretch.h"
#include "tree_solver.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

class Timer {
  std::chrono::time_point<std::chrono::steady_clock> start;

public:
  Timer() {
    start = std::chrono::steady_clock::now();
  }

  void tic(const std::string& msg="") {
    cout << msg << endl;
    start = std::chrono::steady_clock::now();
  }

  void toc() {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << "Finished in " << duration.count() << "s" << endl;
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

void akpw(const EdgeList<EdgeR>& es) {
  EdgeList<EdgeR> tree;
  AKPW(es, tree);

  for (const auto& e : tree.edges) {
    cout << e.u << ' ' << e.v << ' ' << e.resistance << '\n';
  }
}

/*
void dijkstra(const EdgeList<EdgeR>& es) {
  struct DV {
  public:
    double distance;
    size_t parent;

    void initialize(size_t i) {
      parent = i;
      distance = std::numeric_limits<double>::max();
    }

    void SetPrev(size_t u) {
      parent = u;
    }
  };

  struct CMP {
    vector<DV>& vs;

    CMP(vector<DV>& vs_) : vs(vs_) { }

    bool operator() (size_t u, size_t v) const {
      return vs[u].distance > vs[v].distance;
    }
  };

  AdjacencyArray g(es);
  vector<DV> vs(es.n);
  for (size_t i = 0; i < vs.size(); i++) {
    vs[i].initialize(i);
  }
  vs[10].distance = 0;
  CMP less(vs);
  BinaryHeap<CMP> queue(es.n, &less);

  Dijkstra(g, queue, vs);

  for (size_t i = 0; i < vs.size(); i++) {
    cout << i << ' ' << vs[i].parent << ' ' << vs[i].distance << '\n';
  }
}
*/

void aug_tree_pcg(const EdgeList<EdgeR>& es, const vector<FLOAT>& b, size_t k) {
  cout << "===== aug-tree PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  Tree<TreeVertexR> t;

  timer.tic("Finding low stretch spanning tree... ");
    EdgeList<EdgeR> tree_es;
    // recursive_c(k, k, tree_es);
    AKPW(es, tree_es);
    AdjacencyArray<ArcR> g(tree_es);
    DijkstraTree<Tree<TreeVertexR>>(g, es.n / 2, t);
  timer.toc();

  EdgeList<EdgeR> off_tree_es;
  off_tree_es.n = es.n;

  for (size_t i = 0; i < es.edges.size(); i++) {
    const EdgeR& e = es.edges[i];
    if (t.vertices[e.u].parent == e.v || t.vertices[e.v].parent == e.u) {
      continue;
    }
    off_tree_es.AddEdge(EdgeR(e.u, e.v, e.resistance));
  }

  vector<double> strs(off_tree_es.edges.size());

  timer.tic("Computing stretch and adding edges... ");
  ComputeStretch(t, off_tree_es, strs);

  AdjacencyMap aug_tree(t);
  size_t count = 0;
  FLOAT total_stretch = std::accumulate(strs.begin(), strs.end(), 0);
  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> unif01(0, 1);

  for (size_t i = 0; i < strs.size(); i++) {
    FLOAT p = 5 * k * strs[i] / total_stretch; // add 5k edges in expectation
    if (unif01(rng) < p) {
      // off_tree_es.edges[i].resistance *= 5 * k * strs[i] / total_stretch;
      count++;
      aug_tree.AddEdgeR(off_tree_es.edges[i]);
    }
  }
  // cout << "Added " << count << " edges\n";
  // cout << k << '\n';

  // StretchGreater less(strs);
  // vector<size_t> indices(off_tree_es.edges.size());
  // for (size_t i = 0; i < indices.size(); i++) {
  //   indices[i] = i;
  // }
  // std::sort(indices.begin(), indices.end(), less);
  // for (size_t i = 0; i < 5 * k; i++) {
  //   aug_tree.AddEdgeR(off_tree_es.edges[indices[i]]);
  // }
  timer.toc();

  timer.tic("Factorizing aug tree... ");
  CholeskySolver precon(aug_tree);
  timer.toc();

  EdgeList<EdgeC> es2(es);
  PCGSolver<EdgeList<EdgeC>, CholeskySolver> s(&es2, &precon);

  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("Aug tree pcg... ");
  s.Solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;

  /*
  std::ofstream esf("es.txt");
  esf << es2.n << '\n';
  for (const auto& e : es2.edges) {
    esf << e.u << ' ' << e.v << ' ' << std::setprecision(17) << e.conductance << '\n';
  }
  esf.close();

  std::ofstream bf("b.txt");
  for (const auto& f : b) {
    bf << std::setprecision(17) << f << '\n';
  }
  bf.close();

  std::ofstream xf("x.txt");
  for (const auto& f : x) {
    xf << std::setprecision(17) << f << '\n';
  }
  xf.close();
  */

  return;
}

void aug_tree_pcg2(const EdgeList<EdgeR>& es, const vector<FLOAT>& b, size_t k) {
  cout << "===== aug-tree PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  CholeskySolver precon;

  timer.tic("Constructing preconditioner... ");
  AugTreePrecon3(es, precon, es.n / 10);
  timer.toc();

  EdgeList<EdgeC> es2(es);
  PCGSolver<EdgeList<EdgeC>, CholeskySolver> s(&es2, &precon);

  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("Aug tree pcg... ");
  s.Solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;

  return;
}

void sparse_cholesky(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== sparse cholesky PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  CholeskySolver precon;

  AdjacencyMap g(es);
  timer.tic("Constructing preconditioner... ");
  SparseCholesky(g, log(es.n), precon.cholesky_factor);
  timer.toc();

  EdgeList<EdgeC> es2(es);
  PCGSolver<EdgeList<EdgeC>, CholeskySolver> s(&es2, &precon);

  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("PCG... ");
  s.Solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;

  return;
}

void stretch(std::mt19937& rng) {
  size_t k = 4;
  size_t n = k * k;

  EdgeList<EdgeR> es;

  // std::uniform_real_distribution<> unif_1_100(1, 100);
  // RNG<std::uniform_real_distribution<>, std::mt19937> random_resistance(unif_1_100, rng);

  grid2(k, k, es);

  Tree<TreeVertexR> t;

  {
    AdjacencyArray<ArcR> g(es);
    DijkstraTree<Tree<TreeVertexR>>(g, 0, t);
  }

  vector<double> strs(es.edges.size());
  ComputeStretch(t, es, strs);

  for (size_t i = 0; i < t.vertices.size(); i++) {
    cout << i << ' ' << t.vertices[i].parent << endl;
  }

  for (size_t i = 0; i < es.edges.size(); i++) {
    EdgeR& e = es.edges[i];
    cout << e.u << ' ' << e.v << ' ' << strs[i] << endl;
  }
}

void min_degree(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== min degree =====\n";
  AdjacencyMap g(es);
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  timer.tic("factorizing: ");
  CholeskySolver s(g);
  timer.toc();

  timer.tic("solving: ");
  s.Solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);
  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}

void resistance_vs_conductance(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== resistance vs conductance =====\n";
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  EdgeList<EdgeC> es2(es);

  auto id = IdentitySolver();
  PCGSolver<EdgeList<EdgeR>, IdentitySolver> s(&es, &id);
  PCGSolver<EdgeList<EdgeC>, IdentitySolver> s2(&es2, &id);

  for (size_t i = 0; i < 2; i++) {
    for (auto&f : x) {
      f = 0;
    }
    timer.tic("resistance: ");
    s.Solve(b, x);
    timer.toc();

    mv(-1, es, x, 1, b, r);
    cout << MYSQRT(r * r) / MYSQRT(b * b) << endl;

    for (auto&f : x) {
      f = 0;
    }
    timer.tic("conductance: ");
    s2.Solve(b, x);
    timer.toc();

    mv(-1, es, x, 1, b, r);
    cout << MYSQRT(r * r) / MYSQRT(b * b) << endl;
  }
}

void bst(std::mt19937& rng) {
  cout << "===== Solving random complete BST =====\n";
  size_t n = 65535;
  Tree<TreeSolverVertex> tree(n);
  std::uniform_real_distribution<> weight(1, 100);

  // a complete binary tree
  for (size_t i = 0; i * 2 + 2 < n; i++) {
    tree.SetParentR(i * 2 + 1, i, weight(rng));
    tree.SetParentR(i * 2 + 2, i, weight(rng));
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
  s.Solve(b, x);

  // r = b - A * x
  std::vector<FLOAT> r(n);
  mv(-1, tree, x, 1, b, r);

  std::cout << "bst\n";
  std::cout << MYSQRT(r * r) << std::endl;
}

void pcg(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== PCG =====\n";
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  auto id = IdentitySolver();
  PCGSolver<EdgeList<EdgeR>, IdentitySolver> s(&es, &id);

  timer.tic();
  s.Solve(b, x);
  timer.toc();

  mv(-1, es, x, 1, b, r);

  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}

int main(void) {
  size_t k = 1000;
  size_t n = k * k;

  EdgeList<EdgeR> unweighted_grid;
  EdgeList<EdgeR> weighted_grid;

  struct {
    bool operator() (const EdgeR& e1, const EdgeR& e2) const {
      return e1.resistance < e2.resistance;
    }
  } less;

  std::sort(weighted_grid.edges.begin(), weighted_grid.edges.end(), less);

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> uniform(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937>
    random_resistance(uniform, rng);

  grid2(k, k, unweighted_grid);
  grid2(k, k, weighted_grid, random_resistance);

  std::vector<FLOAT> unweighted_b(n);
  std::vector<FLOAT> weighted_b(n);
  std::vector<FLOAT> x(n);

  for (auto& f : x) {
    f = uniform(rng);
  }
  mv(1, weighted_grid, x, 0, x, weighted_b);

  for (auto& f : x) {
    f = uniform(rng);
  }
  mv(1, unweighted_grid, x, 0, x, unweighted_b);

  // bst(rng);
  // pcg(unweighted_grid, unweighted_b);
  // resistance_vs_conductance(weighted_grid, weighted_b);
  // min_degree(weighted_grid, weighted_b);
  // aug_tree_pcg(weighted_grid, weighted_b, k);
  aug_tree_pcg2(weighted_grid, weighted_b, k);
  // sparse_cholesky(weighted_grid, weighted_b);
  // akpw(unweighted_grid);

  return 0;
}
