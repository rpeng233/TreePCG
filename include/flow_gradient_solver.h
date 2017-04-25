#ifndef INCLUDE_FLOW_GRADIENT_SOLVER_H__
#define INCLUDE_FLOW_GRADIENT_SOLVER_H__

#include "disjoint_set.h"
#include "graph.h"

class FlowGradientSolver {
public:
  FlowGradientSolver(const TreeR& tree, const EdgeListR& off_tree_es);
  void Solve(const std::vector<FLOAT>& b, std::vector<FLOAT>& x);

private:
  struct FlowGradientVtx {
    // DisjointSetNode ds_node;
    size_t parent;
    // size_t ancestor;
    double gradient;
    double internal;
    double flow;
    double parent_resistance;
    std::vector<size_t> children;
    std::vector<size_t> incident_edges;
    // bool finished;

    FlowGradientVtx () {
      // finished = false;
    }
  };

  struct FlowGradientEdge {
    size_t u;
    size_t v;
    double resistance;
    double flow;
  };

  size_t root;
  std::vector<FlowGradientVtx> tree;
  std::vector<FlowGradientEdge> off_tree_es;

  void ComputeTreeGradients(std::vector<double>& off_tree_grad, size_t cur);
  void TreeFlow(const std::vector<double>& demand, size_t cur);
  void TreeVoltage(std::vector<double>& tree_vs, size_t cur) const;
  double PrimalEnergy() const ;
  double DualEnergy(const std::vector<double>& b,
                    const std::vector<double>& x) const;

};

#endif  // INCLUDE_FLOW_GRADIENT_SOLVER_H__
