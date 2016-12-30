/********************************************************************************
 *
 * Basic Definitions for Matrixrices and Vectors
 *
 * Public API for Vec:
 *     Vec(n):             length-n vector, initialized to 0
 *     .n:                length of vector
 *    [index]:            access the index-th coodinate (0-indexed)
 *     .norm():            norm of vector
 *    Vec*FLOAT:            scale multiplication
 *     Vec+Vec (Vec-Vec):    vector add/subtract
 *     Vec*Vec:            dot multiplication
 *
 * Public API for Matrix:
 *     Matrix(n,m):            n x m sparse matrix initialized to 0
 *     .n (.m):            size of matrix
 *    .valueAddValue(x,y,z):    Add z to value (x,y)
 *     .sortup():            sort the entries for better cache efficiency
 *     .freeMemory():        destroys matrix and frees memory
 *     Matrix*Vec:            Matrix multiplication
 *     Vec*Matrix:            Matrix multiplication, Vec is treated as row vector
 *
 * NOTE:
 *     Matrix behaves like objects in python or javascript. That is, if you do
 *         Matrix A; Mat B=A;
 *     Then A and B will be actually pointing to the same object,
 *     i.e. modifying B will result in A being modified as well!
 *
 *     However, C++ does not have an automatic garbage collection system,
 *     so you have to run .freeMemory() to free memory of matrices that are no longer needed.
 *
 ********************************************************************************/

#include <algorithm>
#include <vector>
using namespace std;

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "common.h"

struct Vec {
  size_t n;
  FLOAT *value;

  Vec() {
    n = 0;
    value = NULL;
  }

  Vec(size_t _n) {
    n = _n;
    value = new FLOAT[n];

#ifdef USE_MPFR
    for (size_t i = 0; i < n; ++i) {
      value[i] = 0.0;
    }
#else
    memset(value, 0, sizeof(FLOAT) * n);
#endif
  }

  Vec(const Vec &o) {
    n = o.n;
    value = new FLOAT[n];

#ifdef USE_MPFR
    for (size_t i = 0; i < n; ++i) {
      value[i] = o.value[i];
    }
#else
    memcpy(value, o.value, n * sizeof(FLOAT));
#endif
  }

  ~Vec() {
    n = 0;
    delete[] value;
    value = NULL;
  }

  FLOAT &operator[](size_t i) const {
#ifndef NO_RANGE_CHECK
    assert(0 <= 0 && i < n);
#endif
    return value[i];
  }

  Vec &operator =(const Vec &o) {
    if (!value || o.n != n) {
      delete[] value;
      n = o.n;
      value = new FLOAT[n];
    }

#ifdef USE_MPFR
    for (size_t i = 0; i < n; ++i) {
      value[i] = o.value[i];
    }
#else
    memcpy(value, o.value, n * sizeof(FLOAT));
#endif
    return (*this);
  }

  FLOAT Norm() const {
    FLOAT sum = 0;
    for (size_t i = 0; i < n; ++i) {
      sum += value[i] * value[i];
    }
    return MYSQRT(sum);
  }
};

struct MatrixElement {
  size_t row;
  size_t column;
  FLOAT value;

  MatrixElement() {
    row = 0;
    column = 0;
    value = 0;
  }

  MatrixElement(size_t _row, size_t _column, FLOAT _value) {
    row = _row;
    column = _column;
    value = _value;
  }

  bool operator <(const MatrixElement &o) const {
    if (this->row != o.row) {
      return this->row < o.row;
    }

    if (this->column != o.column) {
      return this->column < o.column;
    }

    return this->value < o.value;
  }
};


struct Matrix {
  size_t n;
  size_t m;
  vector<MatrixElement> non_zero;

  Matrix() {
    n = 0;
    m = 0;
  }

  Matrix(size_t _n, size_t _m) {
    n = _n;
    m = _m;
  }

  Matrix(const Matrix &o) {
    n = o.n;
    m = o.m;
    non_zero = o.non_zero;
  }

  Matrix& operator =(const Matrix &o) {
    n = o.n;
    m = o.m;
    non_zero = o.non_zero;
    return (*this);
  }

  void addNonZero(size_t row, size_t column, FLOAT value) {
#ifndef NO_RANGE_CHECK
    assert(0 <= row && row < n && 0 <= column && column < m);
#endif
    non_zero.emplace_back(row, column, value);
  }

  void sortAndCombine() {
    if (non_zero.empty()) return;

    sort(non_zero.begin(), non_zero.end());

    size_t new_nnz = 0;
    MatrixElement last(non_zero[0].row, non_zero[0].column, 0);

    for (const auto& nz : non_zero) {
      if (nz.row == last.row && nz.column == last.column) {
        last.value += nz.value;
      } else {
        non_zero[new_nnz] = last;
        new_nnz++;
        last = nz;
      }
    }

    non_zero[new_nnz] = last;
    new_nnz++;
    non_zero.resize(new_nnz);
  }

  Matrix transpose() const {
    Matrix result(m, n);

    for (const auto& nz : non_zero) {
      result.addNonZero(nz.column, nz.row, nz.value);
    }

    result.sortAndCombine();
    return result;
  }
};

/*
FLOAT operator *(const Vec &a, const Vec &b) {
  assert(a.n == b.n);
  FLOAT result = 0;

  for (size_t i = 0; i < a.n; ++i) {
    result += a.value[i] * b.value[i];
  }

  return result;
}

Vec operator *(const Vec &a, FLOAT b) {
  Vec result(a.n);
  for (size_t i = 0; i < a.n; ++i) {
    result[i] = a.value[i] * b;
  }
  return result;
}

Vec operator +(const Vec &a, const Vec &b) {
  assert(a.n == b.n);
  Vec result(a.n);

  for (size_t i = 0; i < a.n; ++i) {
    result[i] = a.value[i] + b.value[i];
  }
  return result;
}

Vec operator -(const Vec &a, const Vec &b) {
  assert(a.n == b.n);
  Vec result(a.n);

  for (size_t i = 0; i < a.n; ++i) {
    result[i] = a.value[i] - b.value[i];
  }
  return result;
}

Vec operator *(const Matrix &a, const Vec &b) {
  assert(a.m == b.n);

  Vec result(a.n);
  for (vector<MatrixElement>::iterator it = a.non_zero.begin();
      it != a.non_zero.end(); ++it) {
    result.value[it->row] += it->value * b.value[it->column];
  }
  return result;
}

Vec operator *(const Vec &a, const Matrix &b) {
  assert(a.n == b.n);
  Vec result(b.m);

  for (vector<MatrixElement>::iterator it = b.non_zero.begin();
      it != b.non_zero.end(); ++it) {
    result.value[it->column] += it->value * a.value[it->row];
  }
  return result;
}
*/

#endif   // __MATRIX_H__
