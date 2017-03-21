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

static void LCA(vector<LCAVertex>& tree,
                const vector<EdgeR>& es,
                vector<double>& strs,
                size_t id) {
  LCAVertex& u = tree[id];
  u.ds_node.MakeSet();
  u.ancestor = id;
  for (vector<size_t>::const_iterator it = u.children.begin();
       it != u.children.end();
       ++it) {
    LCAVertex& v = tree[*it];
    v.dist += u.dist;
    LCA(tree, es, strs, *it);
    u.ds_node.Union(&v.ds_node);
    LCAVertex *tmp = (LCAVertex *) (
        (char *) u.ds_node.Find() - offsetof(LCAVertex, ds_node)
    );
    tmp->ancestor = id;
  }
  u.finished = true;
  for (vector<size_t>::const_iterator it = u.incident_edges.begin();
       it != u.incident_edges.end();
       ++it) {
    LCAVertex& v = tree[es[*it].u ^ es[*it].v ^ id];
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
