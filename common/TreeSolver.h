#ifndef __TREESOLVER_H__
#define __TREESOLVER_H__

#include <tuple>
#include <vector>
#include "common.h"
#include "graph.h"
#include "pchol.h"

class TreeSolver {
public:
  TreeSolver(const TreePlusEdges& tree) {
    std::tie(std::ignore, pchol) =  partial_cholesky(tree);
    assert(pchol.relabeling.size() == 1);
    // for (auto &e : pchol.elims) {
    //   std::cout << e.v << ' '
    //             << e.nghbr1 << ' '
    //             << e.resistance1 << ' '
    //             << e.rhs << ' '
    //             << e.deg << '\n';
    // }
    // for (auto &it : pchol.relabeling) {
    //   std::cout << it.first << ' ' << it.second << '\n';
    // }
  }

  void solve(
      const std::vector<FLOAT>& b,
      std::vector<FLOAT>& x,
      FLOAT tol=1e-6,
      int maxit=-1
  ) {
    std::vector<FLOAT> partial_x(1);
    std::vector<FLOAT> partial_b = eliminate_rhs(b, pchol);
    assert(partial_b.size() == 1);
    assert(partial_b[0] < 0.001); // sanity check...
    partial_x[0] = 0;
    back_substitution(pchol, partial_x, x);

    FLOAT sum = 0;
    for (const auto& f : x) {
      sum += f;
    }
    sum /= x.size();

    for (auto& f : x) {
      f -= sum;
    }
  }

private:
  PartialCholesky pchol;
};

#endif
