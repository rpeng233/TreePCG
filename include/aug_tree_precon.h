#ifndef INCLUDE_AUG_TREE_PRECON_H_
#define INCLUDE_AUG_TREE_PRECON_H_

#include "cholesky.h"
#include "graph.h"

using std::vector;

void AugTreePrecon(const EdgeListR& tree_es,
                   const std::vector<double>& stretches,
                   const EdgeListR& off_tree_es,
                   EdgeListC& precon,
                   size_t k);

void AugTreePrecon(const EdgeListR& es,
                   const EdgeListR& tree_es,
                   EdgeListC& precon,
                   size_t k);

void AugTreePrecon(const EdgeListR& es,
                   EdgeListC& precon,
                   size_t k);

void AugTreePrecon2(const EdgeListR& es,
                   EdgeListC& precon,
                    size_t k);

void AugTreePrecon3(const EdgeListR& es,
                    CholeskySolver& precon,
                    size_t k);

#endif  // INCLUDE_AUG_TREE_PRECON_H_
