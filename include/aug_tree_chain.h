#ifndef INCLUDE_AUG_TREE_CHAIN_H__
#define INCLUDE_AUG_TREE_CHAIN_H__

#include "cholesky.h"
#include "cholmod_solver.h"
#include "graph.h"

struct PreconLevel {
  EdgeListC es;
  CholeskySolver pchol;
  size_t iter;
  std::vector<size_t> id;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> d;
  std::vector<double> y;
  std::vector<double> x_small;
  std::vector<double> b_small;
  std::vector<double> invD;
};

void AugTreeChain(const EdgeListR& tree_es,
                  const EdgeListR& off_tree_es_,
                  std::vector<PreconLevel>& chain,
                  CholmodSolver& base_solver);

#endif  // INCLUDE_AUG_TREE_CHAIN_H__
