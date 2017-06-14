#ifndef INCLUDE_AUG_TREE_CHAIN_H__
#define INCLUDE_AUG_TREE_CHAIN_H__

#include "cholesky.h"
#include "cholmod_solver.h"
#include "graph.h"

struct PreconLevel {
  EdgeListC es;
  CholeskySolver pchol;
  size_t n;
  size_t iter;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> d;
  std::vector<double> y;
  std::vector<double> invD;
};

void AugTreeChain(const EdgeListR& tree_es,
                  const EdgeListR& off_tree_es_,
                  std::vector<PreconLevel>& chain,
                  std::vector<size_t>& new_id,
                  std::vector<size_t>& original_id,
                  CholmodSolver& base_solver);

#endif  // INCLUDE_AUG_TREE_CHAIN_H__
