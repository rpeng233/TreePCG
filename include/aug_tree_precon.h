#ifndef INCLUDE_AUG_TREE_PRECON_H_
#define INCLUDE_AUG_TREE_PRECON_H_

#include "cholesky.h"

using std::vector;

void AugTreePrecon(const EdgeListR& es,
                   CholeskySolver& precon,
                   size_t k);

#endif  // INCLUDE_AUG_TREE_PRECON_H_
