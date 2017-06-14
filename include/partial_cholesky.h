#include "common.h"
#include "graph.h"
#include "cholesky.h"

struct PCholVertex {
  size_t parent;
  FLOAT parent_resistance;
  size_t ref_count;
  bool eliminated;


  PCholVertex() :
    parent(0), parent_resistance(0.0), ref_count(0), eliminated(false)
  { }

  void SetParentR(size_t p, FLOAT r) {
    parent = p;
    parent_resistance = r;
  }

  FLOAT ParentResistance() const {
    return parent_resistance;
  }

  FLOAT ParentConductance() const {
    return 1.0 / parent_resistance;
  }
};

template <>
inline void Tree<PCholVertex>::SetParentR(size_t v, size_t p, FLOAT r) {
  vertices[p].ref_count++;
  vertices[v].SetParentR(p, r);
}

void PartialCholesky(Tree<PCholVertex>& tree,
                     EdgeListR& off_tree_es,
                     CholeskyFactor& factor,
                     std::vector<size_t>& eorder);
