#include "cholmod.h"
#include "cholmod_solver.h"
#include "graph.h"

CholmodSolver::CholmodSolver(const EdgeListC& es, cholmod_common *c) {
  common = c;
  if (common == NULL) {
    common = &common_;
    cholmod_start(common);
  }

  n = es.n;

  cholmod_triplet *T = cholmod_allocate_triplet(
      n,
      n,
      es.Size() + n,
      1,
      CHOLMOD_REAL,
      common
  );

  std::vector<size_t> col_count(n, 1);
  std::vector<double> degrees(n);
  for (size_t j = 0; j < es.Size(); j++) {
    const EdgeC& e = es[j];
    // col_count[std::min(e.u, e.v)]++;
    col_count[e.u]++;
    col_count[e.v]++;
    degrees[e.u] += e.conductance;
    degrees[e.v] += e.conductance;
  }

  int *Ti = reinterpret_cast<int *>(T->i);
  int *Tj = reinterpret_cast<int *>(T->j);
  double *Tx = reinterpret_cast<double *>(T->x);

  for (size_t i = 0; i < n; i++) {
    Ti[T->nnz] = i;
    Tj[T->nnz] = i;
    Tx[T->nnz] = degrees[i];
    // Tx[T->nnz] = degrees[i] * 1.01;
    if (i == 0) {
      Tx[T->nnz] *=  1.1;
    }
    T->nnz++;
  }

  for (size_t j = 0; j < es.Size(); j++) {
    const EdgeC& e = es[j];
    size_t u = std::min(e.u, e.v);
    size_t v = std::max(e.u, e.v);
    Ti[T->nnz] = u;
    Tj[T->nnz] = v;
    Tx[T->nnz] = -e.conductance;
    T->nnz++;
  }

  cholmod_sparse *A = cholmod_triplet_to_sparse(T, T->nnz, common);
  cholmod_free_triplet(&T, common);

  std::cerr << "Analyzing..." << std::endl;
  factor = cholmod_analyze(A, common);
  std::cerr << "Factorizing..." << std::endl;
  cholmod_factorize(A, factor, common);
  std::cerr << "Done!" << std::endl;

  cholmod_free_sparse(&A, common);

  b = cholmod_zeros(n, 1, CHOLMOD_REAL, common);
}

void CholmodSolver::Solve(const std::vector<FLOAT>& b_,
                          std::vector<FLOAT>& x_) {
  cholmod_dense *b = cholmod_zeros(n, 1, CHOLMOD_REAL, common);
  double* bx = reinterpret_cast<double *>(b->x);

  for (size_t i = 0; i < n; i++) {
    bx[i] = b_[i];
  }

  cholmod_solve2(CHOLMOD_A, factor, b, NULL, &x, NULL, &Y, &E, common);
  // cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b, common);
  double* xx = reinterpret_cast<double *>(x->x);

  for (size_t i = 0; i < n; i++) {
    x_[i] = xx[i];
  }
}
