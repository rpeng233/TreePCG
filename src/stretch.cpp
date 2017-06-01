#include <cstddef>
#include "disjoint_set.h"
#include "graph.h"
#include "stretch.h"

using std::vector;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;

class LCAVertex {
public:
  DisjointSetNode ds_node;
  double dist;
  size_t ancestor;
  std::vector<size_t> children;
  std::vector<size_t> incident_edges;
  bool finished;

  LCAVertex () {
    finished = false;
    dist = 0;
  }
};

struct Frame {
  const size_t id;
  size_t child;

  Frame(size_t i, size_t c) : id(i), child(c) { }
};

static void LCA(vector<LCAVertex>& tree,
                const vector<EdgeR>& es,
                vector<double>& strs,
                size_t root) {
  std::vector<Frame> stack;

  stack.push_back(Frame(root, 0));
  tree[root].ds_node.MakeSet();
  tree[root].ancestor = root;

  while (!stack.empty()) {
    Frame& f = stack.back();
    LCAVertex& u = tree[f.id];
    assert(f.child <= u.children.size());
    if (f.child > 0) {
      LCAVertex& v = tree[u.children[f.child - 1]];
      u.ds_node.Union(&v.ds_node);
      LCAVertex *tmp = (LCAVertex *) (
          (char *) u.ds_node.Find() - offsetof(LCAVertex, ds_node)
      );
      tmp->ancestor = f.id;
    }
    if (f.child < u.children.size()) {
      size_t cid = u.children[f.child];
      f.child++;
      LCAVertex& v = tree[cid];
      v.ds_node.MakeSet();
      v.ancestor = u.children[f.child];
      v.dist += u.dist;
      stack.push_back(Frame(cid, 0));
      continue;
    }
    u.finished = true;
    stack.pop_back();
    for (vector<size_t>::const_iterator it = u.incident_edges.begin();
         it != u.incident_edges.end();
         ++it) {
      LCAVertex& v = tree[es[*it].u ^ es[*it].v ^ f.id];
      assert(es[*it].u == f.id || es[*it].v == f.id);
      if (v.finished) {
        LCAVertex *p = (LCAVertex *) (
            (char *) v.ds_node.Find() - offsetof(LCAVertex, ds_node)
        );
        double d = tree[p->ancestor].dist;
        assert(v.dist >= d);
        assert(u.dist >= d);
        strs[*it] = (v.dist + u.dist - d * 2) / es[*it].resistance;
      }
    }
  }
}

void ComputeStretch(const TreeR& tree,
                    const EdgeListR& es,
                    vector<double>& stretch) {
  const vector<TreeVertexR>& vs = tree.vertices;
  vector<LCAVertex> lca_tree(vs.size());
  size_t root = 0;
  size_t found_root = 0;

  stretch.clear();
  stretch.resize(es.Size(), -1);

  for (size_t i = 0; i < vs.size(); i++) {
    if (vs[i].parent != i) {
      lca_tree[vs[i].parent].children.push_back(i);
      lca_tree[i].dist = vs[i].parent_resistance;
    } else {
      root = i;
      found_root++;
    }
  }
  assert(found_root == 1);

  for (size_t i = 0; i < es.Size(); i++) {
    lca_tree[es[i].u].incident_edges.push_back(i);
    lca_tree[es[i].v].incident_edges.push_back(i);
  }

  LCA(lca_tree, es.edges, stretch, root);
}
