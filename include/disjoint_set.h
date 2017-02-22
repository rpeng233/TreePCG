#ifndef INCLUDE_DISJOINT_SET_H
#define INCLUDE_DISJOINT_SET_H

class DisjointSetNode {
private:
  DisjointSetNode *parent;
  int rank;
public:
  DisjointSetNode () {
    parent = this;
    rank = 0;
  }

  void MakeSet() {
    parent = this;
    rank = 0;
  }

  DisjointSetNode *Find() {
    if (parent != this) {
      parent = parent->Find();
    }

    return parent;
  }

  void Union(DisjointSetNode *y) {
    DisjointSetNode *x_root = this->Find();
    DisjointSetNode *y_root = y->Find();

    if (x_root == y_root) return;

    if (x_root->rank < y_root->rank) {
      x_root->parent = y_root;
    } else if (x_root->rank > y_root->rank) {
      y_root->parent = x_root;
    } else {
      y_root->parent = x_root;
      x_root->rank++;
    }
  }
};

#endif // INCLUDE_DISJOINT_SET_H
