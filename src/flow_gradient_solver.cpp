#include <cstddef>
#include <vector>
#include "flow_gradient_solver.h"
#include "linear_algebra.h"

FlowGradientSolver::FlowGradientSolver(const TreeR& t,
                                       const EdgeListR& off_tree)
: tree(t.n),
  off_tree_es(off_tree.Size()) {
  size_t found_root = 0;
  for (size_t i = 0; i < tree.size(); i++) {
    if (t[i].parent != i) {
      tree[i].parent = t[i].parent;
      tree[i].parent_resistance = t[i].parent_resistance;
      tree[tree[i].parent].children.push_back(i);
    } else {
      root = i;
      found_root++;
    }
  }
  assert(found_root == 1);

  for (size_t i = 0; i < off_tree_es.size(); i++) {
    off_tree_es[i].u = std::min(off_tree[i].u, off_tree[i].v);
    off_tree_es[i].u = std::max(off_tree[i].u, off_tree[i].v);
    off_tree_es[i].resistance = off_tree[i].resistance;
    tree[off_tree_es[i].u].incident_edges.push_back(i);
    tree[off_tree_es[i].v].incident_edges.push_back(i);
  }

  std::cout << tree.size() << ' ' << off_tree_es.size() << '\n';
}

void FlowGradientSolver::ComputeTreeGradients(
    std::vector<double>& off_tree_grad,
    size_t cur
) {
  FlowGradientVtx& u = tree[cur];

  u.gradient = 0;
  for (std::vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    ComputeTreeGradients(off_tree_grad, *it);
    u.gradient += tree[*it].gradient;
  }

  for (std::vector<size_t>::const_iterator it = u.incident_edges.begin();
       it != u.incident_edges.end();
       ++it) {
    FlowGradientEdge& e = off_tree_es[*it];
    double tmp = off_tree_grad[*it] / e.resistance;
    u.gradient += tmp * (e.u == cur ? -1 : 1);
  }
}

void FlowGradientSolver::TreeFlow(const std::vector<double>& demand,
                                  size_t cur) {
  FlowGradientVtx& u = tree[cur];
  u.flow = demand[cur];
  for (std::vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    TreeFlow(demand, *it);
    u.flow += tree[*it].flow;
  }
}

void FlowGradientSolver::TreeVoltage(std::vector<double>& tree_vs,
                                     size_t cur) const {
  const FlowGradientVtx& u = tree[cur];

  if (cur == root) {
    tree_vs[cur] = 0;
  } else {
    tree_vs[cur] = tree_vs[u.parent] + u.flow * u.parent_resistance;
  }

  for (std::vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    TreeVoltage(tree_vs, *it);
  }
}

double FlowGradientSolver::PrimalEnergy() const {
  double energy = 0;
  for (size_t i = 0; i < tree.size(); i++) {
    energy += tree[i].flow * tree[i].flow * tree[i].parent_resistance;
  }

  for (size_t i = 0; i < off_tree_es.size(); i++) {
    const FlowGradientEdge& e = off_tree_es[i];
    energy += e.flow * e.flow * e.resistance;
  }

  return energy;
}

double FlowGradientSolver::DualEnergy(const std::vector<double>& b,
                                      const std::vector<double>& x) const {
  double energy = 0;
  for (size_t i = 0; i < b.size(); i++) {
    energy += b[i] * x[i];
  }
  energy *= 2;

  for (size_t i = 0; i < tree.size(); i++) {
    if (i == root) continue;
    double drop = x[tree[i].parent] - x[i];
    energy -= drop * drop / tree[i].parent_resistance;
  }

  for (size_t i = 0; i < off_tree_es.size(); i++) {
    const FlowGradientEdge& e = off_tree_es[i];
    double drop = x[e.u] - x[e.v];
    energy -= drop * drop / e.resistance;
  }

  return energy;
}

void FlowGradientSolver::Solve(const std::vector<double>& b,
                               std::vector<double>& x) {
  std::vector<double> tree_vs(tree.size());
  std::vector<double> demand(b);

  for (size_t i = 0; i < off_tree_es.size(); i++) {
    off_tree_es[i].flow = 0;
  }

  // compute the feasible tree flow
  TreeFlow(demand, root);

  for (size_t iter = 0; ; iter++) {
    // compute tree voltages
    TreeVoltage(x, root);

    double p = PrimalEnergy();
    double d = DualEnergy(b, x);
    std::cout << iter << '\t'
              << "Primal energy:\t" << p << '\t'
              << "Dual energy:\t" << d << '\t'
              << "Gap:\t" << p - d << '\n';

    // demand = b;
    // for (size_t i = 0; i < off_tree_es.size(); i++) {
    //   FlowGradientEdge& e = off_tree_es[i];
    //   double f = (x[e.u] - x[e.v]) / e.resistance;
    //   demand[e.u] -= f;//e.flow;
    //   demand[e.v] += f;//e.flow;
    // }
    // for (size_t i = 0; i < tree.size(); i++) {
    //   if (i == root) continue;
    //   double f = (x[i] - x[tree[i].parent]) / tree[i].parent_resistance;
    //   demand[i] -= f;//tree[i].flow;
    //   demand[tree[i].parent] += f;//tree[i].flow;
    // }
    // std::cout << "Residual: " << demand * demand << std::endl;

    if (p - d < 0.5 * p) break;
    // if (iter >= 20) break;

    // compute gradient on off-tree edges
    std::vector<double> off_tree_grad(off_tree_es.size());
    double cRc = 0;
    double fRc = 0;
    for (size_t i = 0; i < off_tree_es.size(); i++) {
      FlowGradientEdge& e = off_tree_es[i];
      off_tree_grad[i] = e.flow * e.resistance - (x[e.u] - x[e.v]);
      cRc += off_tree_grad[i] * off_tree_grad[i] * e.resistance;
      fRc += off_tree_grad[i] * e.flow * e.resistance;
    }

    // compute gradient on tree edges
    ComputeTreeGradients(off_tree_grad, root);
    // ComputeOverestmates(off_tree_grad, root);
    for (size_t i = 0; i < tree.size(); i++) {
      // tree[i].gradient -= tree[i].internal * 2;
      tree[i].gradient *= tree[i].parent_resistance;
      fRc += tree[i].flow * tree[i].gradient * tree[i].parent_resistance;
      cRc += tree[i].gradient * tree[i].gradient * tree[i].parent_resistance;
    }

    // update flow on off tree edges
    double step_size = fRc / cRc;
    // std::cout << "step size: " << step_size << std::endl;
    // std::cout << "improvement: " << fRc * fRc / cRc << std::endl;
    for (size_t i = 0; i < off_tree_es.size(); i++) {
      off_tree_es[i].flow -= step_size * off_tree_grad[i];
    }

    demand = b;
    for (size_t i = 0; i < off_tree_es.size(); i++) {
      FlowGradientEdge& e = off_tree_es[i];
      demand[e.u] -= e.flow;
      demand[e.v] += e.flow;
    }

    // fix the demand using a tree flow
    TreeFlow(demand, root);
  }
}
