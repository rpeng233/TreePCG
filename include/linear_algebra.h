#ifndef INCLUDE_LINEAR_ALGEBRA_H__
#define INCLUDE_LINEAR_ALGEBRA_H__

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

  for (size_t i = 0; i < A.m; i++) {
    result[i] = 0;
  }

  for (std::vector<MatrixElement>::const_iterator it = A.non_zero.begin();
       it != A.non_zero.end();
       ++it) {
    result[it->row] += it->value * x[it->column];
  }

  for (size_t i = 0; i < x.size(); i++) {
    result[i] = alpha * result[i] + beta * y[i];
  }
}

inline void mv(
    FLOAT alpha,
    TreeR& tree,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  size_t n = tree.n;
  assert(n == x.size());
  assert(n == y.size());
  assert(n == result.size());

  for (size_t i = 0; i < n; i++) {
    result[i] = 0;
  }

  std::vector<TreeVertexR>& vs = tree.vertices;
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

inline void mv(
    FLOAT alpha,
    const EdgeListR& es,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  size_t n = es.n;
  assert(n == x.size());
  assert(n == y.size());
  assert(n == result.size());

  for (size_t i = 0; i < n; i++) {
    result[i] = 0;
  }

  const std::vector<EdgeR>& edges = es.edges;
  for (size_t i = 0; i < edges.size(); i++) {
    size_t u = edges[i].u;
    size_t v = edges[i].v;
    const FLOAT& r = edges[i].resistance;
    FLOAT current = (x[u] - x[v]) / r;
    result[u] += current;
    result[v] -= current;
  }

  for (size_t i = 0; i < n; i++) {
    result[i] = alpha * result[i] +  beta * y[i];
  }
}

inline void mv(
    FLOAT alpha,
    const EdgeListC& es,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  size_t n = es.n;
  assert(n == x.size());
  assert(n == y.size());
  assert(n == result.size());

  for (size_t i = 0; i < n; i++) {
    result[i] = 0;
  }

  const std::vector<EdgeC>& edges = es.edges;
  for (size_t i = 0; i < edges.size(); i++) {
    size_t u = edges[i].u;
    size_t v = edges[i].v;
    const FLOAT& c = edges[i].conductance;
    FLOAT current = (x[u] - x[v]) * c;
    result[u] += current;
    result[v] -= current;
  }

  for (size_t i = 0; i < n; i++) {
    result[i] = alpha * result[i] +  beta * y[i];
  }
}

/*
inline void mv(
    FLOAT alpha,
    const Graph& graph,
    const std::vector<FLOAT>& x,
    FLOAT beta,
    const std::vector<FLOAT>& y,
    std::vector<FLOAT>& result
) {
  auto n = graph.n;
  assert(n == x.size());
  assert(n == y.size());
  assert(n == result.size());

  for (auto& f : result) {
    f = 0;
  }

  auto& vs = graph.vertices;
  for (size_t i = 0; i < n; i++) {
    for (const auto& a : vs[i].nghbrs) {
      FLOAT current = (x[i] - x[a.v]) / a.resistance;
      result[i] += current;
    }
  }

  for (size_t i = 0; i < n; i++) {
    result[i] = alpha * result[i] +  beta * y[i];
  }
}
*/

#endif  // INCLUDE_LINEAR_ALGEBRA_H
