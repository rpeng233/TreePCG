#include "linear_algebra.h"
#include "precon_chain_solver.h"

inline void PreconChainSolver::PreconSolve(size_t lvl,
                                           const std::vector<double>& b,
                                           std::vector<double>& x) {
  CholeskySolver& pchol = chain[lvl].pchol;
  std::vector<double>& y = chain[lvl].y;

  pchol.ForwardSubstitution(b, y);
  SolveChain(lvl + 1, y, x);
  pchol.BackSubstitution(y, x);
}

void PreconChainSolver::SolveChain(size_t lvl,
                                   const std::vector<double>& b,
                                   std::vector<double>& x) {
  if (lvl >= chain.size()) {
    base_solver.Solve(b, x);
    return;
  }

  EdgeListC& es = chain[lvl].es;
  std::vector<double>& r = chain[lvl].r;
  std::vector<double>& s = chain[lvl].s;
  std::vector<double>& d = chain[lvl].d;

  // double bnorm = MYSQRT(b * b);
  // double res;
  double delta_old;
  double delta_new;
  double alpha;

  size_t n = chain[lvl].n;

  for (size_t i = 0; i < n; i++) {
    x[i] = 0;
  }

  // for (size_t j = 0; j < 0; j++) {
  //   mv(-1, es, x, 1, b, r);
  //   for (size_t i = 0; i < chain[lvl].invD.size(); i++) {
  //     r[i] *= chain[lvl].invD[i];
  //   }
  //   axpy(1, x, r, x, n);
  // }

  mv(-1, es, x, 1, b, r);
  PreconSolve(lvl, r, d);
  delta_new = r * d;
  for (size_t i = 0; i < chain[lvl].iter - 1; i++) {
    mv(1, es, d, 0, r, s);
    alpha = delta_new / dot(d, s, n);
    axpy(alpha, d, x, x, n);
    // mv(-1, es, x, 1, b, r);
    axpy(-alpha, s, r, r, n);
    PreconSolve(lvl, r, s);
    delta_old = delta_new;
    delta_new = dot(r, s, n);
    axpy(delta_new / delta_old, d, s, d, n);
  }

  // for (size_t i = 0; i < chain[lvl].iter; i++) {
  //   mv(-1, es, x, 1, b, r);
  //   PreconSolve(lvl, r, s);
  //   axpy(1, s, x, x, n);
  // }

  // for (size_t j = 0; j < 0; j++) {
  //   mv(-1, es, x, 1, b, r);
  //   for (size_t i = 0; i < chain[lvl].invD.size(); i++) {
  //     r[i] *= chain[lvl].invD[i];
  //   }
  //   axpy(1, x, r, x, n);
  // }

  // mv(-1, es, x, 1, b, r);
  // std::cout << lvl << ' ' << chain[lvl].iter << ' ' << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}
