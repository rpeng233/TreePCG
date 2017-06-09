#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "akpw.h"
#include "aug_tree_precon.h"
#include "aug_tree_chain.h"
#include "cholesky.h"
#include "cholmod.h"
#include "cholmod_solver.h"
#include "common.h"
#include "cycle_toggling_solver.h"
#include "flow_gradient_solver.h"
#include "graph.h"
#include "graph_gen.h"
#include "identity_solver.h"
#include "incomplete_cholesky.h"
#include "io.h"
#include "matrix.h"
#include "partial_cholesky.h"
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
  std::string open;
  std::string close;

public:
  Timer(const std::string& o=">>>>> ", const std::string& c="<<<<< ")
    : open(o), close(c) {
    start = std::chrono::steady_clock::now();
  }

  void tic(const std::string& msg="") {
    cout << open << msg << endl;
    start = std::chrono::steady_clock::now();
  }

  void toc() {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << close << "Finished in " << duration.count() << "s\n" << endl;
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

  // for (const auto& e : tree.edges) {
  //   cout << e.u << ' ' << e.v << ' ' << e.resistance << '\n';
  // }
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

void aug_tree_pcg(const EdgeList<EdgeR>& es,
                  const EdgeListR& tree_es,
                  const vector<FLOAT>& b,
                  size_t top,
                  size_t k) {
  cout << "===== aug-tree PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  EdgeListC aug_tree;

  timer.tic("constructing augmented tree... ");
  AugTreePrecon(es, tree_es, aug_tree, top, k);
  timer.toc();

  // FILE *f = fopen("precon.mtx", "w");
  // WriteMtx(f, aug_tree);
  // fclose(f);
  // return;

  timer.tic("factorizing... ");
  cholmod_common common;
  cholmod_start(&common);
  common.supernodal = CHOLMOD_SIMPLICIAL;
  CholmodSolver precon(aug_tree, &common);
  timer.toc();

  EdgeList<EdgeC> es2(es);
  PCGSolver<EdgeList<EdgeC>, CholmodSolver> s(&es2, &precon);

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

  EdgeList<EdgeC> es2(es);
  EdgeList<EdgeC> es3(es);

  timer.tic("Constructing preconditioner... ");
  SparseCholesky3(es3, log(es.n) + 1, precon.cholesky_factor);
  timer.toc();

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

void incomplete_cholesky(const EdgeList<EdgeR>& es, const vector<FLOAT>& b, double eps) {
  cout << "===== incomplete cholesky PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  CholeskySolver precon;

  AdjacencyMap g(es);

  EdgeList<EdgeC> es2(es);
  EdgeList<EdgeC> es3(es);

  timer.tic("Constructing preconditioner... ");
  IncompleteCholesky(es3, eps, precon.cholesky_factor);
  timer.toc();

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
  // size_t n = k * k;

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
    cout << i << '\t' << t.vertices[i].parent << endl;
  }

  for (size_t i = 0; i < es.edges.size(); i++) {
    EdgeR& e = es.edges[i];
    cout << e.u << '\t' << e.v << '\t' << strs[i] << endl;
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

void cholmod(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== cholmod =====\n";
  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  cholmod_common common;
  cholmod_start(&common);
  common.supernodal = CHOLMOD_SIMPLICIAL;

  timer.tic("factorizing: ");
  CholmodSolver s(es, &common);
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

void flow_gradient_descent(const EdgeListR& es, const vector<FLOAT>& b) {
  EdgeListR tree_es;
  EdgeListR off_tree_es;
  AKPW(es, tree_es);

  TreeR tree;
  AdjacencyArray<ArcR> g(tree_es);
  DijkstraTree(g, es.n / 2, tree);
  g.FreeMemory();

  for (size_t i = 0; i < es.Size(); i++) {
    const EdgeR& e = es[i];
    if (tree[e.u].parent == e.v || tree[e.v].parent == e.u) {
      continue;
    }
    off_tree_es.AddEdge(e);
  }

  FlowGradientSolver solver(tree, off_tree_es);

  std::vector<FLOAT> x(es.n);
  std::vector<FLOAT> r(es.n);

  solver.Solve(b, x);

  mv(-1, es, x, 1, b, r);
  std::cout << MYSQRT(r * r) << std::endl;
}

void cycle_toggling(EdgeListR& es, const vector<FLOAT>& b) {
  EdgeListR tree;
  EdgeListR otes;
  AKPW(es, tree);
  for (size_t i = 0; i < tree.Size(); i++) {
    tree[i].resistance /= 1;
  }
  TreeR t;
  AdjacencyArray<ArcR> g(tree);
  DijkstraTree<TreeR>(g, es.n / 2, t);
  for (size_t i = 0; i < es.Size(); i++) {
    EdgeR& e = es[i];
    if (t[e.u].parent == e.v || t[e.v].parent == e.u) {
      e.resistance /= 1;
    }
    otes.AddEdge(e);
  }

  CycleTogglingSolver s(t, otes);

  std::vector<double> x(es.n);
  std::vector<FLOAT> r(es.n);

  s.Solve(b, x);
  // for (auto d : x) {
  //   printf("%5f ", d);
  // }
  // printf("\n");
  // mv(1, es, x, 0, b, r);
  // for (auto d : r) {
  //   printf("%.2f ", d);
  // }
  // printf("\n");
  mv(-1, es, x, 1, b, r);
  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}

void partial_cholesky(EdgeListR& es, const vector<FLOAT>& b) {
  AdjacencyArray<ArcR> g(es);
  Tree<PCholVertex> tree;
  EdgeListR off_tree_es;

  DijkstraTree(g, es.n / 2, tree);
  g.FreeMemory();

  CholeskySolver s;

  for (size_t i = 0; i < es.Size(); i++) {
    const EdgeR& e = es[i];
    if (tree[e.u].parent == e.v || tree[e.v].parent == e.u) {
      continue;
    }
    off_tree_es.AddEdge(e);
    tree[e.u].ref_count++;
    tree[e.v].ref_count++;
  }

  PartialCholesky(tree, off_tree_es, s.cholesky_factor);

  // AdjacencyMap g(es);
  // CholeskySolver s(g);

  std::vector<double> x(es.n);
  std::vector<FLOAT> r(es.n);

  s.Solve(b, x);

  mv(-1, es, x, 1, b, r);
  std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
  std::cout << MYSQRT(r * r) << std::endl;
  std::cout << MYSQRT(b * b) << std::endl;
  std::cout << s.cholesky_factor.n << std::endl;
  std::cout << s.cholesky_factor.elims.size() << std::endl;
  std::cout << s.cholesky_factor.elim_arcs.size() << std::endl;
}

int main(void) {
  size_t k = 1000;

  EdgeList<EdgeR> es;
  EdgeListR tree_es;

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> uniform(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937>
    random_resistance(uniform, rng);

  // grid2(k, k, es);
  // grid2(k, k, es, random_resistance);
  // grid3(k, k, k, es);
  // grid3(k, k, k, es, random_resistance);
  cycle(1000000, es);
  recursive_c(k, k, tree_es);

  // ReadEdgeList(stdin, es);
  // FILE *f = fopen("grid.txt", "w");
  // WriteEdgeList(f, es);
  // fclose(f);
  // f = fopen("grid.mtx", "w");
  // WriteMtx(f, es);
  // fclose(f);
  // return 0;

  size_t n = es.n;

  struct {
    bool operator() (const EdgeR& e1, const EdgeR& e2) const {
      return e1.resistance < e2.resistance;
    }
  } less;

  std::sort(es.edges.begin(), es.edges.end(), less);

  std::vector<FLOAT> random_b(n);
  std::vector<FLOAT> unit_b(n);
  std::vector<FLOAT> x(n);

  unit_b[0] = 100;
  unit_b[n - 1] = -100;

  for (auto& f : x) {
    f = uniform(rng);
  }
  mv(1, tree_es, x, 0, x, random_b);

  // bst(rng);
  // pcg(es, random_b);
  // resistance_vs_conductance(es, random_b);
  // min_degree(es, random_b);
  // aug_tree_pcg(es, tree_es, unit_b, 0, es.Size() / 4);
  // cholmod(es, unit_b);
  // sparse_cholesky(es, random_b);
  // incomplete_cholesky(es, unit_b, 1e-6);
  // akpw(es);
  // flow_gradient_descent(es, unit_b);
  // cycle_toggling(es, unit_b);
  partial_cholesky(es, random_b);

  // WriteMtx(stdout, es);

  return 0;
}
