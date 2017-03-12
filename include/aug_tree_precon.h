#ifndef INCLUDE_AUG_TREE_PRECON_H_
#define INCLUDE_AUG_TREE_PRECON_H_

#include "cholesky.h"
#include "graph.h"

using std::vector;

void AugTreePrecon(const EdgeListR& es,
                   CholeskySolver& precon,
                   size_t k);
void AugTreePrecon2(const EdgeListR& es,
                    CholeskySolver& precon,
                    size_t k);
void AugTreePrecon3(const EdgeListR& es,
                    CholeskySolver& precon,
                    size_t k);

#endif  // INCLUDE_AUG_TREE_PRECON_H_
