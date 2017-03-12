#ifndef INCLUDE_INCOMPLETE_CHOLESKY_H_
#define INCLUDE_INCOMPLETE_CHOLESKY_H_

#include "cholesky.h"
#include "graph.h"
void IncompleteCholesky(EdgeListC& es,
                        size_t k,
                        CholeskyFactor& cholesky_factor);

#endif  // INCLUDE_INCOMPLETE_CHOLESKY_H_
