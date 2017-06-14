#ifndef INCLUDE_PRECON_CHAIN_SOLVER_H__
#define INCLUDE_PRECON_CHAIN_SOLVER_H__

#include <vector>
#include "graph.h"
#include "aug_tree_chain.h"

class PreconChainSolver {
public:
  PreconChainSolver(const EdgeListR& tree_es, const EdgeListR& off_tree_es) {
    AugTreeChain(tree_es,
                 off_tree_es,
                 chain,
                 new_id,
                 original_id,
                 base_solver);
  }

  void Solve(const std::vector<double>& b,
             std::vector<double>& x) {
    SolveChain(0, b, x);
  }

  void SolveChain(size_t lvl,
                  const std::vector<double>& b,
                  std::vector<double>& x);


  const PreconLevel& GetLevel(size_t lvl) {
    return chain[lvl];
  }

  std::vector<size_t> new_id;
  std::vector<size_t> original_id;
private:
  std::vector<PreconLevel> chain;
  CholmodSolver base_solver;


  inline void PreconSolve(size_t lvl,
                          const std::vector<double>& b,
                          std::vector<double>& x);
};

#endif  // INCLUDE_PRECON_CHAIN_SOLVER_H__
