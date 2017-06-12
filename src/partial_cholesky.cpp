#include "partial_cholesky.h"

void PartialCholesky(Tree<PCholVertex>& tree,
                     EdgeListR& off_tree_es,
                     CholeskyFactor& factor) {
  size_t n = tree.n;
  factor.n = n;
  std::vector<EliminatedVertex>& elims = factor.elims;
  std::vector<ArcC>& elim_arcs = factor.elim_arcs;

  EliminatedVertex ev;

  size_t count;

  do {
    count = 0;

    // eliminate the leaves
    size_t leaf = 0;
    for (size_t i = 0; i < n; i++) {
      if (tree[i].ref_count > 0 || tree[i].eliminated || tree[i].parent == i) {
        continue;
      }
      size_t j = i;
      do {
        size_t p = tree[j].parent;
        assert(p != j);
        assert(tree[p].eliminated == false);

        ev.v = j;
        ev.degree = tree[j].ParentConductance();
        ev.first_arc = elim_arcs.size();
        elims.push_back(ev);
        leaf++;

        elim_arcs.push_back(ArcC(p, tree[j].ParentConductance()));

        tree[j].eliminated = true;
        tree[p].ref_count--;


        j = p;
      } while (tree[j].ref_count == 0 && tree[j].parent != j);
    }

    count += leaf;
    // std::cout << "eliminated " << leaf << " leaves\n";

    // eliminate degree 2 vertices with two tree edges
    size_t tpath = 0;
    for (size_t i = 0; i < n; ) {
      size_t j = tree[i].parent;
      size_t k = tree[j].parent;
      if (i == j ||
          j == k ||
          tree[i].eliminated ||
          tree[j].ref_count != 1) {
        i++;
        continue;
      }

      assert(tree[j].eliminated == false);
      assert(tree[k].eliminated == false);

      // we have i -> j -> k, we want to eliminated j

      ev.v = j;
      ev.degree = tree[i].ParentConductance() + tree[j].ParentConductance();
      ev.first_arc = elim_arcs.size();
      elims.push_back(ev);
      tpath++;

      elim_arcs.push_back(ArcC(i, tree[i].ParentConductance()));
      elim_arcs.push_back(ArcC(k, tree[j].ParentConductance()));

      tree[j].eliminated = true;
      tree[i].parent = k;
      tree[i].parent_resistance += tree[j].ParentResistance();

    }

    count += tpath;
    // std::cout << "eliminated " << tpath << " tree paths\n";

    // eliminate the other degree 2 vertices
    size_t deg2 = 0;
    for (size_t i = 0; i < off_tree_es.Size(); ) {
      EdgeR& e = off_tree_es[i];

      if (e.resistance < 0) {
        i++;
        continue; // e is no longer present
      }

      assert(tree[e.u].eliminated == false);
      assert(tree[e.v].eliminated == false);

      bool flag = false;

      if (tree[e.u].parent == e.v) {
        double c = tree[e.u].ParentConductance();
        c += 1 / e.resistance;
        tree[e.u].parent_resistance = 1 / c;
        e.resistance = -1;
        tree[e.u].ref_count--;
        tree[e.v].ref_count--;
        continue;
      } else if (tree[e.v].parent == e.u) {
        double c = tree[e.v].ParentConductance();
        c += 1 / e.resistance;
        tree[e.v].parent_resistance = 1 / c;
        e.resistance = -1;
        tree[e.u].ref_count--;
        tree[e.v].ref_count--;
        continue;
      }

      // try to eliminate e.u
      if (!tree[e.u].eliminated && tree[e.u].ref_count == 1) {
        size_t p = tree[e.u].parent;

        assert(tree[p].eliminated == false);
        assert(p != e.u);

        ev.v = e.u;
        ev.degree = tree[e.u].ParentConductance() + 1 / e.resistance;
        ev.first_arc = elim_arcs.size();
        elims.push_back(ev);
        deg2++;

        elim_arcs.push_back(ArcC(p, tree[e.u].ParentConductance()));
        elim_arcs.push_back(ArcC(e.v, 1 / e.resistance));

        tree[e.u].eliminated = true;
        e.resistance += tree[e.u].ParentResistance();

        if (tree[e.v].parent == p) {
          double c = tree[e.v].ParentConductance();
          c += 1 / e.resistance;
          tree[e.v].parent_resistance = 1 / c;
          tree[p].ref_count--;
          tree[e.v].ref_count--;
          e.resistance = -1;
          i++;
          continue;
        } else {
          e.u = p;
          flag = true;
        }

      }

      // try the same thing with e.v
      if (!tree[e.v].eliminated && tree[e.v].ref_count == 1) {
        size_t p = tree[e.v].parent;

        assert(tree[p].eliminated == false);
        assert(p != e.v);

        ev.v = e.v;
        ev.degree = tree[e.v].ParentConductance() + 1 / e.resistance;
        ev.first_arc = elim_arcs.size();
        elims.push_back(ev);
        deg2++;

        elim_arcs.push_back(ArcC(p, tree[e.v].ParentConductance()));
        elim_arcs.push_back(ArcC(e.u, 1 / e.resistance));

        tree[e.v].eliminated = true;
        e.resistance += tree[e.v].ParentResistance();

        if (tree[e.u].parent == p) {
          double c = tree[e.u].ParentConductance();
          c += 1 / e.resistance;
          tree[e.u].parent_resistance = 1 / c;
          tree[p].ref_count--;
          tree[e.u].ref_count--;
          e.resistance = -1;
          i++;
          continue;
        } else {
          e.v = p;
          flag = true;
        }
      }

      if (!flag) {
        i++;
      }
    }

    count += deg2;
    // std::cout << "eliminated " << deg2 << " deg 2\n";
  } while (count != 0);

  ev.first_arc = elim_arcs.size();
  elims.push_back(ev);
} 
