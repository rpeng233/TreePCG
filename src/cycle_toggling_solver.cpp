#include <cstddef>
#include <cstdio>
#include <random>
#include "cycle_toggling_solver.h"

CycleTogglingSolver::CycleTogglingSolver(const TreeR& t, const EdgeListR& o)
: tree(t.n),
  hld(t.n),
  es(o.Size()),
  stretches(o.Size()) {
  std::vector<HelperNode> helper(t.n);
  size_t found_root = 0;
  for (size_t i = 0; i < tree.size(); i++) {
    tree[i].parent = t[i].parent;
    helper[i].Initialize(i);
    if (t[i].parent != i) {
      tree[i].resistance = t[i].parent_resistance;
      helper[t[i].parent].children.push_back(i);
    } else {
      root = i;
      tree[i].resistance = 0;
      found_root++;
    }
  }

  assert(found_root == 1);
  assert(tree[root].parent == root);
  assert(tree[root].resistance == 0);

  for (size_t i = 0; i < o.Size(); i++) {
    es[i].u = o[i].u;
    es[i].v = o[i].v;
    es[i].resistance = o[i].resistance;
    helper[es[i].u].incident_edges.push_back(i);
    helper[es[i].v].incident_edges.push_back(i);
  }

  HLD(helper, root);
  LCA(helper, root);

  // for (size_t i = 0; i < helper.size(); i++) {
  //   std::cout << "Node " << i << ":"
  //             << " head = " << helper[i].is_head
  //             << " size = " << helper[i].size
  //             << " h = " << helper[i].heavy << std::endl;
  // }
  // for (size_t i = 0; i < hld.size(); i++) {
  //   std::cout << "VNode " << i << ":"
  //             << " p = " << hld[i].parent << ','
  //             << " l = " << hld[i].left_child << ','
  //             << " r = " << hld[i].right_child << ','
  //             << " t = " << (int) hld[i].type << std::endl;
  // }
  // for (size_t i = 0 ; i < es.size(); i++) {
  //   std::cout << "Edge " << i << ":"
  //             << " (" << es[i].u << ", " << es[i].v << "), "
  //             << "lca = " << es[i].lca
  //             << " str = " << stretches[i] << std::endl;
  // }
}

void CycleTogglingSolver::Solve(const std::vector<double>& b,
                                std::vector<double>& x) {
  ComputeTreeFlow(b);
  DecomposeTreeFlow();

  std::mt19937 rng(std::random_device{}());
  std::discrete_distribution<unsigned>
    sample(stretches.begin(), stretches.end());

  for (size_t i = 0; i < 5000000; i++) {
    size_t e = sample(rng);
    Toggle(es[e]);
  }
  Dump();
  tree[root].flow = 0;

  // for (size_t i = 0; i < tree.size(); i++) {
  //   std::cout << "Node " << i << ":"
  //             << " p = " << tree[i].parent << ','
  //             << " r = " << tree[i].resistance << ','
  //             << " f = " << tree[i].flow << std::endl;
  // }

  // for (size_t i = 0; i < es.size(); i++) {
  //   std::cout << "Edge " << i << ":"
  //             << " (" << es[i].u << ", " << es[i].v << "), "
  //             << "f = " << es[i].flow << std::endl;
  // }

  ComputeTreeVoltage(x);
}

void CycleTogglingSolver::HLD(std::vector<HelperNode>& helper, size_t root) {
  DFS(helper, root);
  chain_roots.clear();
  std::vector<size_t> chain;
  for (size_t i = 0; i < helper.size(); i++) {
    if (helper[i].is_head) {
      chain.clear();
      for (size_t j = i; ; j = helper[j].heavy) {
        chain.push_back(j);
        if (helper[j].size == 1) break;
      }
      std::pair<size_t, double> p = BST(chain, 0, chain.size() - 1);
      chain_roots.push_back(p.first);
      hld[p.first].parent = (i == root) ? p.first : tree[i].parent;
      hld[p.first].type = VIRTUAL;
    }
  }
}

double CycleTogglingSolver::Query(size_t v) {
  double rf_sum = 0;
  bool was_left = false;

  for (;;) {
    double resistance_so_far = 0;
    for (;;) {
      if (was_left) {
        rf_sum += resistance_so_far * hld[v].flow_to_left;
      } else {
        rf_sum += hld[v].rf_sum;
        resistance_so_far += hld[v].resistance_to_left;
      }
      was_left = hld[v].type == LEFT;
      if (hld[v].type == VIRTUAL) break;
      v = hld[v].parent;
    }
    if (hld[v].parent == v) break;
    v = hld[v].parent;
  }

  return rf_sum;
}

void CycleTogglingSolver::Update(size_t v, double delta) {
  bool was_left = false;

  for (;;) {
    double rf_delta_so_far = 0;
    for (;;) {
      if (was_left) {
        hld[v].rf_sum += rf_delta_so_far;
      } else {
        hld[v].flow_to_left += delta;
        double rf_delta = hld[v].resistance_to_left * delta;
        hld[v].rf_sum += rf_delta;
        rf_delta_so_far += rf_delta;
      }
      was_left = hld[v].type == LEFT;
      if (hld[v].type == VIRTUAL) break;
      v = hld[v].parent;
    }
    if (hld[v].parent == v) break;
    v = hld[v].parent;
  }
}

void CycleTogglingSolver::Toggle(OffTreeEdge& e) {
  size_t u = std::min(e.u, e.v);
  size_t v = std::max(e.u, e.v);

  double sum_r = e.resistance;
  sum_r += tree[u].resistance_to_root + tree[v].resistance_to_root;
  sum_r -= tree[e.lca].resistance_to_root * 2;

  double sum_rf = e.resistance * e.flow;
  sum_rf += Query(v) - Query(u);

  double delta = sum_rf / sum_r;
  e.flow -= delta;
  Update(u, delta);
  Update(v, -delta);
}

void CycleTogglingSolver::LCA(std::vector<HelperNode>& helper,
                              size_t cur) {
  HelperNode& u = helper[cur];
  u.ds_node.MakeSet();
  u.ancestor = cur;
  tree[cur].resistance_to_root += tree[cur].resistance;
  for (std::vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    tree[*it].resistance_to_root += tree[cur].resistance_to_root;
    LCA(helper, *it);
    u.ds_node.Union(&helper[*it].ds_node);
    HelperNode *tmp = (HelperNode *) (
        (char *) u.ds_node.Find() - offsetof(HelperNode, ds_node)
    );
    tmp->ancestor = cur;
  }
  u.finished = true;
  for (std::vector<size_t>::const_iterator it = u.incident_edges.begin();
       it != u.incident_edges.end();
       ++it) {
    size_t ngbr = es[*it].u ^ es[*it].v ^ cur;
    HelperNode *p = (HelperNode *) (
        (char *) helper[ngbr].ds_node.Find() - offsetof(HelperNode, ds_node)
    );
    if (helper[ngbr].finished) {
      es[*it].lca = p->ancestor;
      double r = tree[p->ancestor].resistance_to_root;
      assert(tree[cur].resistance_to_root >= r);
      assert(tree[ngbr].resistance_to_root >= r);
      stretches[*it] = (
          tree[cur].resistance_to_root + tree[ngbr].resistance_to_root - r * 2
      ) / es[*it].resistance;
    }
  }
}

void CycleTogglingSolver::DFS(std::vector<HelperNode>& helper, size_t cur) {
  preorder.push_back(cur);
  helper[cur].size = 1;
  helper[cur].is_head = true;
  for (size_t i = 0; i < helper[cur].children.size(); i++) {
    size_t v = helper[cur].children[i];
    DFS(helper, v);
    // helper[v].parent = cur;
    helper[cur].size += helper[v].size;
    if (helper[v].size > helper[cur].heavy_size) {
      helper[cur].heavy_size = helper[v].size;
      helper[cur].heavy = v;
    }
  }
  if (helper[cur].heavy != cur) {
    helper[helper[cur].heavy].is_head = false;
  }
}

std::pair<size_t, double>
CycleTogglingSolver::BST(const std::vector<size_t>& chain, size_t l, size_t r) {
  // std::cout << "DFS(" << chain[l] << ", " << chain[r] << ")\n";
  size_t m = (l + r) / 2;
  size_t v = chain[m];
  double total_resistance = tree[v].resistance;

  // preorder.push_back(v);
  hld[v].rf_sum = 0;
  hld[v].flow_to_left = 0;
  hld[v].resistance_to_left = tree[v].resistance;

  if (l < m) {
    std::pair<size_t, double> p = BST(chain, l, m - 1);
    size_t left_child = p.first;
    hld[v].left_child = left_child;
    hld[left_child].parent = v;
    hld[left_child].type = LEFT; 
    hld[v].resistance_to_left += p.second; 
    total_resistance += p.second;
  } else {
    hld[v].left_child = v;
  }

  if (m < r) {
    std::pair<size_t, double> p = BST(chain, m + 1, r);
    size_t right_child = p.first;
    hld[v].right_child = right_child;
    hld[right_child].parent = v;
    hld[right_child].type = RIGHT;
    total_resistance += p.second;
  } else {
    hld[v].right_child = v;
  }

  return std::make_pair(v, total_resistance);
}

void CycleTogglingSolver::DumpChain(size_t v, double flow_from_right) {
  tree[v].flow = hld[v].flow_to_left + flow_from_right;
  if (hld[v].left_child != v) {
    DumpChain(hld[v].left_child, hld[v].flow_to_left + flow_from_right);
  }
  if (hld[v].right_child != v) {
    DumpChain(hld[v].right_child, flow_from_right);
  }
}

void CycleTogglingSolver::Dump() {
  for (size_t i = 0; i < chain_roots.size(); i++) {
    DumpChain(chain_roots[i], 0);
  }
}

void CycleTogglingSolver::ComputeTreeFlow(const std::vector<double>& demand) {
  for (std::vector<size_t>::const_reverse_iterator it = preorder.rbegin();
       it != preorder.rend();
       it++) {
    Node& v = tree[*it];
    v.flow += demand[*it];
    tree[v.parent].flow += v.flow;
  }
}

void CycleTogglingSolver::ComputeTreeVoltage(std::vector<double>& x) {
  for (std::vector<size_t>::const_iterator it = preorder.begin();
       it != preorder.end();
       it++) {
    if (*it == root) {
      // std::cout << "v[" << *it << "] = 0\n";
      x[*it] = 0;
    } else {
      const Node& u = tree[*it];
      x[*it] = x[u.parent] + u.flow * u.resistance;
      // std::cout << *it << ' ' << u.parent << ", f = " << u.flow
      //           << ", r = " << u.resistance << ", v[parent] = " <<x[u.parent] << ", thus v = "
      //           << x[*it] << std::endl;
    }
  }
}

double CycleTogglingSolver::DecomposeChainFlow(size_t v, double flow_from_right) {
  hld[v].flow_to_left = tree[v].flow - flow_from_right;
  hld[v].rf_sum = hld[v].flow_to_left * hld[v].resistance_to_left;
  double total_rf = hld[v].flow_to_left * hld[v].resistance_to_left;
  if (hld[v].left_child != v) {
    double tmp =
      DecomposeChainFlow(hld[v].left_child, tree[v].flow);
    hld[v].rf_sum += tmp;
    total_rf += tmp;
  }
  if (hld[v].right_child != v) {
    total_rf += DecomposeChainFlow(hld[v].right_child, flow_from_right);
  }

  return total_rf;
}

void CycleTogglingSolver::DecomposeTreeFlow() {
  for (size_t i = 0; i < chain_roots.size(); i++) {
    DecomposeChainFlow(chain_roots[i], 0);
  }
}
