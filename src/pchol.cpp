#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "pchol.h"

using std::vector;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;

static size_t pchol_leaf(vector<TreeVertex>& vs, vector<PartialElim>& elims)
{
  size_t elim_count = 0;
  for (size_t i = 0; i < vs.size(); i++) {
    if (vs[i].ref_count != 0 || vs[i].eliminated || vs[i].parent == i) {
      continue;
    }
    size_t cur = i;
    do {
      size_t p = vs[cur].parent;
      assert(p != cur);
      assert(vs[p].eliminated == false);
      elims.emplace_back(cur, p, vs[cur].parent_resistance);
      vs[p].ref_count--;
      vs[cur].eliminated = true;
      elim_count++;
      cur = p;
    } while (vs[cur].ref_count == 0 && vs[cur].parent != cur);
  }
  return elim_count;
}


static size_t pchol_tree_path(vector<TreeVertex>& vs,
                              vector<PartialElim>& elims) {
  size_t elim_count = 0;
  for (size_t i = 0; i < vs.size(); i++) {
    size_t j = vs[i].parent;
    size_t k = vs[j].parent;
    // we have i -> j -> k in the tree, we want to eliminated j
    if (vs[i].eliminated
        || vs[j].eliminated
        || vs[j].ref_count != 1 // j has off-tree edges
        || k == j               // j is the root
    ) {
      continue;
    }
    elims.emplace_back(j,
                       i, vs[i].parent_resistance,
                       k, vs[j].parent_resistance);
    vs[i].parent = k;
    vs[i].parent_resistance = vs[i].parent_resistance + vs[j].parent_resistance;
    vs[j].eliminated = true;
    elim_count++;
  }
  return elim_count;
}


static size_t pchol_deg2(vector<TreeVertex>& vs,
                         vector<Edge>& es,
                         vector<PartialElim>& elims) {
  size_t elim_count = 0;
  for (auto& e : es) {
    // tries to eliminate one of the endpoints of e
    if (e.resistance < 0) continue; // e is already eliminated

    if (!vs[e.u].eliminated && vs[e.u].ref_count == 1) {
      // e.v -> e.u -> e.u.parent is a path of length 2
      vs[e.u].eliminated = true;
      vs[e.u].ref_count--;
      elim_count++;
      size_t p = vs[e.u].parent;
      elims.emplace_back(e.u, p, vs[e.u].parent_resistance, e.v, e.resistance);
      FLOAT res1 = vs[e.u].parent_resistance + e.resistance;
      if (vs[e.v].parent == p) {
        // we have a triangle: e.v.parent == e.u.parent
        FLOAT res2 = vs[e.v].parent_resistance;
        vs[e.v].parent_resistance = 1 / (1 / res1 + 1 / res2);
        vs[e.v].ref_count--;
        vs[p].ref_count--;
        e.resistance = -1;
      } else {
        e.resistance = res1;
        e.u = p;
      }
      continue;
    }

    if (!vs[e.v].eliminated && vs[e.v].ref_count == 1) {
      vs[e.v].eliminated = true;
      vs[e.v].ref_count--;
      elim_count++;
      size_t p = vs[e.v].parent;
      elims.emplace_back(e.v, p, vs[e.v].parent_resistance, e.u, e.resistance);
      FLOAT res1 = vs[e.v].parent_resistance + e.resistance;
      if (vs[e.u].parent == p) {
        FLOAT res2 = vs[e.u].parent_resistance;
        vs[e.u].parent_resistance = 1 / (1 / res1 + 1 / res2);
        vs[e.u].ref_count--;
        vs[p].ref_count--;
        e.resistance = -1;
      } else {
        e.resistance = res1;
        e.v = p;
      }
      continue;
    }
  }

  return elim_count;
}

std::pair<TreePlusEdges, PartialCholesky>
partial_cholesky(const TreePlusEdges& tree) {
  // assert(tree.vertices.size() <= rhs_.size());

  PartialCholesky pchol;
  vector<TreeVertex> vs(tree.vertices);
  vector<Edge> es(tree.off_tree_edges);

  size_t elim_count;
  do {
    elim_count = 0;
    elim_count += pchol_leaf(vs, pchol.elims);
    // elim_count += pchol_tree_path(vs, pchol.elims);
    // elim_count += pchol_deg2(vs, es, pchol.elims);
  } while (elim_count > 0);


  size_t next_id = 0;
  for (size_t i = 0; i < vs.size(); i++) {
    if (vs[i].eliminated) continue;
    pchol.relabeling[i] = next_id++;
  }

  TreePlusEdges new_tree(pchol.relabeling.size());
  for (const auto& p : pchol.relabeling) {
    auto id = p.first;
    auto new_id = p.second;
    new_tree.setParent(new_id,
                       pchol.relabeling[vs[id].parent],
                       vs[id].parent_resistance);
  }
  for (auto& e : es) {
    if (e.resistance < 0) continue;
    new_tree.addEdge(pchol.relabeling[e.u],
                     pchol.relabeling[e.v],
                     e.resistance);
  }
  return std::make_pair(
      std::move(new_tree),
      std::move(pchol)
  );
}


void back_substitution(
    const PartialCholesky& pchol,
    const vector<FLOAT>& partial_x,
    vector<FLOAT>& x
) {
  for (const auto& it : pchol.relabeling) {
    x[it.first] = partial_x[it.second];
  }

  for (auto it = pchol.elims.rbegin(); it != pchol.elims.rend(); ++it) {
    if (it->deg == 1) {
      x[it->v] = it->rhs * it->resistance1 + x[it->nghbr1];
    } else if (it->deg == 2) {
      x[it->v] = it->rhs + x[it->nghbr1] + x[it->nghbr2];
      x[it->v] *= (1 / it->resistance1 + 1 / it->resistance2);
    } else {
      cerr << __FILE__ << ":" << __LINE__;
      cerr << ": Invalid value for PartialElim::deg" << endl;
      exit(1);
    }
  }
}


vector<FLOAT> eliminate_rhs(const vector<FLOAT>& rhs_, PartialCholesky& pchol)
{
  vector<FLOAT> rhs(rhs_);
  vector<FLOAT> new_rhs(pchol.relabeling.size());

  for (auto& elim : pchol.elims) {
    if (elim.deg == 1) {
      elim.rhs = rhs[elim.v];
      rhs[elim.nghbr1] += rhs[elim.v];
    } else if (elim.deg == 2) {
      elim.rhs = rhs[elim.v];
      FLOAT w1 = 1 / elim.resistance1;
      FLOAT w2 = 1 / elim.resistance2;
      rhs[elim.nghbr1] += w1 * rhs[elim.v] / (w1 + w2);
      rhs[elim.nghbr2] += w2 * rhs[elim.v] / (w1 + w2);
    } else {
      cerr << __FILE__ << ":" << __LINE__
           << ": Invalid value for PartialElim.deg" << endl;
      exit(1);
    }
  }

  for (const auto& p : pchol.relabeling) {
    new_rhs[p.second] = rhs[p.first];
  }

  return new_rhs;
}
