#include "cycle_toggling_solver.h"

struct HelperNode {
  DisjointSetNode ds_node;
  size_t ancestor;
  std::vector<size_t> children;
  std::vector<size_t> incident_edges;
  bool finished;
  bool is_head;

  size_t size;
  size_t heavy;
  size_t heavy_size;
};

static void LCA(std::vector<HelperNode>& tree,
                const Vector<OffTreeEdge>& es,
                size_t cur) {
  HelperNode& u = tree[cur];
  u.ds_node.MakeSet();
  u.ancestor = cur;
  for (vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    HelperNode& v = tree[*it];
    v.dist += u.dist;
    LCA(tree, es, *it);
    u.ds_node.Union(&v.ds_node);
    LCAVertex *tmp = (LCAVertex *) (
        (char *) u.ds_node.Find() - offsetof(LCAVertex, ds_node)
    );
    tmp->ancestor = cur;
  }
  u.finished = true;
  for (vector<size_t>::const_iterator it = u.incident_edges.begin();
       it != u.incident_edges.end();
       ++it) {
    LCAVertex& v = tree[es[*it].u ^ es[*it].v ^ id];
    if (v.finished) {
      double d = tree[p->ancestor].dist;
      assert(v.dist >= d);
      assert(u.dist >= d);
      strs[*it] = (v.dist + u.dist - d * 2) / es[*it].resistance;
    }
  }
}

CycleTogglingSolver::CycleTogglingSolver(const TreeR& t, const EdgeListR& o)
: tree(t.n),
  hld(t.n),
  es(o.Size()) {
  std::vector<HelperNode> helper(t.n);
  size_t found_root = 0;
  for (size_t i = 0; i < tree.size(); i++) {
    if (t[i].parent != i) {
      tree[i].parent = t[i].parent;
      tree[i].resistance = t[i].parent_resistance;
      helper[t[i].parent].children.push_back(i);
      helper[i].dist = t[i].parent_resistance;
    } else {
      root = i;
      found_root++;
    }
  }
  assert(found_root == 1);

  for (size_t i = 0; i < o.Size(); i++) {
    es[i].u = o[i].u;
    es[i].v = o[i].v;
    es[i].resistance = o[i].resistance;
    helper[es[i].u].incident_edges.push_back(i);
    helper[es[i].v].incident_edges.push_back(i);
  }

  LCA;
}

void DFS(std::vector<HelperNode>& tree, size_t cur) {
  tree[cur].size = 1;
  for (size_t i = 0; i < tree[cur].children.size(); i++) {
    size_t v = tree[cur].children[i];
    dfs(g, v);
    tree[v].parent = cur;
    tree[cur].size += tree[v].size;
    if (tree[v].size > tree[cur].heavy_size) {
      tree[cur].heavy_size = tree[v].size;
      tree[cur].heavy = v;
    }
  }
  tree[tree[cur].heavy].is_head = false;
}

std::pair<size_t, double>
BST(const std::vector<size_t>& chain, size_t l, size_t r) {
  size_t m = (r + l) / 2;
  size_t root = chain[m];
  double total_resistance = tree[root].resistance;

  preorder.push_back(root);
  hld[root].flow_to_left = 0;
  hld[root].rf_sum = 0;
  hld[root].resistance_to_left = tree[root].resistance;

  if (l < m) {
    std::pair<size_t, double> p = BST(chain, l, m - 1);
    size_t left_child = p.first;
    hld[left_child].parent = root;
    hld[left_child].type = LEFT; 
    hld[root].resistance_to_left += p.second; 
    total_resistance += p.second;
  } else {
    hld[root].left_child = root;
  }

  if (m < r) {
    std::pair<size_t, double> p = BST(chain, m + 1, r);
    size_t right_child = p.first;
    hld[right_child].parent = root;
    hld[right_child].type = RIGHT;
    total_resistance += p.second;
  } else {
    hld[root].right_child = root;
  }

  return std::make_pair(root, total_resistance);
}

void HLD(std::vector<Node>& tree, size_t root){
  dfs(tree, root);
  std::vector<size_t> chain;
  for (size_t i = 0; i < tree.size(); i++) {
    if (tree[i].is_head) {
      chain.clear();
      chain.push_back(i);
      for (size_t j = i; !tree[j].children.empty(); j = tree[j].heavy) {
        chain.push_back(j);
      }
      std::pair<size_t, double> p = BST(chain, 0, chain.size());
      size_t root = p.first;
      chains.push_back(root);
      hld[root].parent = tree[i].parent;
      hld[root].type = VIRTUAL;
    }
  }
}

void Query(size_t v) {
  double rf_sum = 0;
  bool was_left = false;

  while (vtree[v].parent != v) {
    double resistance_so_far = 0;
    for (;;) {
      if (was_left) {
        rf_sum += resistance_so_far * vtree[v].flow_to_left;
      } else {
        rf_sum += vtree[v].rf_sum;
        resistance_so_far += vtree[v].resistance_to_left;
      }
      if (vtree[v].type == LEFT) {
        was_left = true;
      }
      if (vtree[v].type == VIRTUAL) break;
      v = vtree[v].parent;
    }
    v = vtree[v].parent;
  }

  return rf_sum;
}

void Update(size_t v, double delta) {
  while (vtree[v].parent != v) {
    double rf_delta_so_far = 0;
    for (;;) {
      if (was_left) {
        vtree[v].rf_sum += rf_delta_so_far;
      } else {
        vtree[v].flow_to_left += delta;
        double rf_delta = vtree[v].resistance_to_left * delta;
        vtree[v].rf_sum += rf_delta;
        rf_delta_so_far += rf_delta;
      }
      if (vtree[v].type == LEFT) {
        was_left = true;
      }
      if (vtree[v].type == VIRTUAL) break;
      v = vtree[v].parent;
    }

    v = vtree[v].parent;
  }
}

void DumpChain(size_t v, double flow_from_right) {
  tree[v].flow += vtree[v].flow_to_left + flow_from_right;
  DumpChain(vtree[v].left_child, vtree[v].flow_to_left + flow_from_right);
  DumpChain(vtree[v].right_child, flow_from_right);
}

void Dump() {
  for (size_t i = 0; i < chains.size(); i++) {
    DumpChain(chains[i], 0);
  }
}

void Toggle(const EdgeR& e) {
  size_t u = std::min(e.u, e.v);
  size_t v = std::max(e.u, e.v);

  double sum_r = e.resistance;
  sum_r += tree[u].resistance_to_root + tree[v].resistance_to_root
  sum_r -= tree[e.lca].resistance_to_root * 2;

  double sum_rf = e.resistance * e.flow;
  sum_rf += Query(v) - Query(u);

  double delta = sum_rf / sum_r;
  e.flow -= delta;
  Update(u, delta);
  Update(v, -delta);
}
