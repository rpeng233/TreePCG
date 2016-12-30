#ifndef __LINALG_H__
#define __LINALG_H__

#include <cassert>
#include <vector>
#include "common.h"
#include "matrix.h"

inline FLOAT operator *(
    const std::vector<FLOAT>& a,
    const std::vector<FLOAT>& b
) {
  assert(a.size() == b.size());
  FLOAT result = 0;

  for (size_t i = 0; i < a.size(); i++) {
    result += a[i] * b[i];
  }

  return result;
}

inline void axpy(
    FLOAT a,
    const std::vector<FLOAT>& x,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  assert(x.size() == y.size());
  assert(x.size() == result.size());

  for (size_t i = 0; i < x.size(); i++) {
    result[i] = a * x[i] + y[i];
  }
}

// inline void axpy(
//     FLOAT a,
//     const std::vector<FLOAT>& x,
//     std::vector<FLOAT>& y
// ) {
//   axpy(a, x, y, y);
// }

inline void mv(
    FLOAT alpha,
    const Matrix& A,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  assert(A.m == x.size());
  assert(x.size() == y.size());
  assert(x.size() == result.size());

  for (FLOAT& f : result) {
    f = 0;
  }

  for (const MatrixElement& e : A.non_zero) {
    result[e.row] += e.value * x[e.column];
  }

  for (size_t i = 0; i < x.size(); i++) {
    result[i] = alpha * result[i] + beta * y[i];
  }
}

// inline void mv(
//     FLOAT alpha,
//     const Matrix& A,
//     const std::vector<FLOAT>& x,
//     FLOAT beta,
//     std::vector<FLOAT>& y
// ) {
//   assert(A.m == x.size());
//   assert(x.size() == y.size());
// 
//   std::vector<FLOAT> tmp(x.size(), 0);
// 
//   for (const MatrixElement& e : A.non_zero) {
//     tmp[e.row] += e.value * x[e.column];
//   }
// 
//   for (size_t i = 0; i < x.size(); i++) {
//     y[i] = alpha * tmp[i] + beta * y[i];
//   }
// }

inline void mv(
    FLOAT alpha,
    const Tree& tree,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  auto n = tree.vertices.size();
  assert(n == x.size());
  assert(n == y.size());
  assert(n == result.size());

  for (auto& f : result) {
    f = 0;
  }

  auto& vs = tree.vertices;
  for (size_t i = 0; i < n; i++) {
    if (vs[i].parent == i) continue;
    FLOAT current = (x[i] - x[vs[i].parent]) / vs[i].parent_resistance;
    result[i] += current;
    result[vs[i].parent] -= current;
  }

  for (size_t i = 0; i < n; i++) {
    result[i] = alpha * result[i] +  beta * y[i];
  }
}

#endif
