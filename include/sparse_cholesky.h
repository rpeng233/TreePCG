#ifndef INCLUDE_SPARSE_CHOLESKY_H
#define INCLUDE_SPARSE_CHOLESKY_H

#include "cholesky.h"
#include "graph.h"

void SparseCholesky(EdgeListC& es,
                    size_t k,
                    CholeskyFactor& cholesky_factor);

void SparseCholesky2(EdgeListC& es,
                     size_t k,
                     CholeskyFactor& cholesky_factor);

#endif  // INCLUDE_SPARSE_CHOLESKY_H
