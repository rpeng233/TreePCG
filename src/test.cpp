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
  Timer(const std::string& o=">>>> ", const std::string& c="<<<< ")
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
    cout << close << "Finished in " << duration.count() << "s" << endl;
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
                  size_t k) {
  cout << "===== aug-tree PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  EdgeListC aug_tree;

  timer.tic("constructing augmented tree... ");
  AugTreePrecon(es, aug_tree, k);
  timer.toc();

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

void incomplete_cholesky(const EdgeList<EdgeR>& es, const vector<FLOAT>& b) {
  cout << "===== incomplete cholesky PCG =====\n";
  cout << "n = " << es.n << ", m = " << es.Size() << endl;
  CholeskySolver precon;

  AdjacencyMap g(es);

  EdgeList<EdgeC> es2(es);
  EdgeList<EdgeC> es3(es);

  timer.tic("Constructing preconditioner... ");
  // SparseCholesky(g, log(es.n) + 1, precon.cholesky_factor);
  IncompleteCholesky(es3, 1e-4, precon.cholesky_factor);
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

// int main(void) {
//   size_t n;
//   size_t m;
//   std::cin >> n >> m;
//   EdgeListR tes;
//   EdgeListR otes;
//   EdgeListR es;
// 
//   tes.n = otes.n = es.n = n;
// 
//   size_t i = 0;
//   while (i < n - 1) {
//     size_t u, v;
//     double r;
//     cin >> u >> v >> r;
//     r /= log(n);
//     es.AddEdge(EdgeR(u, v, r));
//     tes.AddEdge(EdgeR(u, v, r));
//     i++;
//   }
//   while (i < m) {
//     size_t u, v;
//     double r;
//     cin >> u >> v >> r;
//     es.AddEdge(EdgeR(u, v, r));
//     otes.AddEdge(EdgeR(u, v, r));
//     i++;
//   }
// 
//   TreeR t;
//   AdjacencyArray<ArcR> g(tes);
//   DijkstraTree<TreeR>(g, 0, t);
//   CycleTogglingSolver s(t, otes);
// 
//   std::vector<double> x(es.n);
//   std::vector<FLOAT> r(es.n);
//   std::vector<FLOAT> b(es.n);
//   b[0] = -1;
//   b[n - 1] = 1;
// 
//   s.Solve(b, x);
//   mv(-1, es, x, 1, b, r);
//   std::cout << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
// 
//   return 0;
// }

int main(void) {
  size_t k = 3;

  EdgeList<EdgeR> unweighted_grid;
  EdgeList<EdgeR> weighted_grid;
  EdgeList<EdgeR> c;
  EdgeListR rec_c;

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<> uniform(1, 100);
  RNG<std::uniform_real_distribution<>, std::mt19937>
    random_resistance(uniform, rng);

  grid2(k, k, unweighted_grid);
  grid2(k, k, weighted_grid, random_resistance);
  // grid3(k, k, 2, unweighted_grid);
  // grid3(k, k, 2, weighted_grid, random_resistance);

  size_t n = unweighted_grid.n;

  // cycle(n, c);
  // recursive_c(k, k, rec_c);

  struct {
    bool operator() (const EdgeR& e1, const EdgeR& e2) const {
      return e1.resistance < e2.resistance;
    }
  } less;

  std::sort(weighted_grid.edges.begin(), weighted_grid.edges.end(), less);

  std::vector<FLOAT> unweighted_b(n);
  std::vector<FLOAT> weighted_b(n);
  std::vector<FLOAT> unit_b(n);
  std::vector<FLOAT> x(n);

  unit_b[0] = 1;
  unit_b[n - 1] = -1;
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
  // aug_tree_pcg(unweighted_grid, rec_c, unit_b, 250 * sqrt(n));
  // cholmod(unweighted_grid, unit_b);
  // sparse_cholesky(weighted_grid, weighted_b);
  incomplete_cholesky(unweighted_grid, unit_b);
  // akpw(weighted_grid);
  // flow_gradient_descent(unweighted_grid, unit_b);
  // cycle_toggling(unweighted_grid, unit_b);

  IO::WriteMtx(stdout, unweighted_grid);

  return 0;
}
