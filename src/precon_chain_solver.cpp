#include "linear_algebra.h"
#include "precon_chain_solver.h"

inline void PreconChainSolver::PreconSolve(size_t lvl,
                                           const std::vector<double>& b,
                                           std::vector<double>& x) {
  CholeskySolver& pchol = chain[lvl].pchol;
  std::vector<double>& y = chain[lvl].y;
  std::vector<double>& x_small = chain[lvl].x_small;
  std::vector<double>& b_small = chain[lvl].b_small;
  std::vector<size_t>& id = chain[lvl].id;

  pchol.ForwardSubstitution(b, y);
  for (size_t j = 0; j < id.size(); j++) {
    b_small[j] = y[id[j]];
  }
  SolveChain(lvl + 1, b_small, x_small);
  for (size_t j = 0; j < id.size(); j++) {
    x[id[j]] = x_small[j];
  }
  pchol.BackSubstitution(y, x);
}

void PreconChainSolver::SolveChain(size_t lvl,
                                   const std::vector<double>& b,
                                   std::vector<double>& x,
                                   double tol) {

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

  for (size_t j = 0; j < 0; j++) {
    mv(-1, es, x, 1, b, r);
    for (size_t i = 0; i < chain[lvl].invD.size(); i++) {
      r[i] *= chain[lvl].invD[i];
    }
    axpy(1, x, r, x);
  }

  mv(-1, es, x, 1, b, r);
  PreconSolve(lvl, r, d);
  delta_new = r * d;
  for (size_t i = 0; i < chain[lvl].iter - 1; i++) {
    mv(1, es, d, 0, r, s);
    alpha = delta_new / (d * s);
    axpy(alpha, d, x, x);
    // mv(-1, es, x, 1, b, r);
    axpy(-alpha, s, r, r);
    PreconSolve(lvl, r, s);
    delta_old = delta_new;
    delta_new = r * s;
    axpy(delta_new / delta_old, d, s, d);
  }

  // for (size_t i = 0; i < chain[lvl].iter; i++) {
  //   mv(-1, es, x, 1, b, r);
  //   PreconSolve(lvl, r, s);
  //   axpy(1, s, x, x);
  // }


  for (size_t j = 0; j < 2; j++) {
    mv(-1, es, x, 1, b, r);
    for (size_t i = 0; i < chain[lvl].invD.size(); i++) {
      r[i] *= chain[lvl].invD[i];
    }
    axpy(1, x, r, x);
  }

  // mv(-1, es, x, 1, b, r);
  // std::cout << lvl << ' ' << chain[lvl].iter << ' ' << MYSQRT(r * r) / MYSQRT(b * b) << std::endl;
}
