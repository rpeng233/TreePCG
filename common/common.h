#ifndef __COMMON_H__
#define __COMMON_H__

/*
 * Uncomment the following macro to use __float128 instead of double
 * remember to add -lquadmath to compile options
 */

// #define USE_FLOAT128

/*
 * Uncomment the following macro to use MPFR instead of double
 * remember to add -Wall -ansi -pedantic -lmpfr to compile options
 * overrides USE_FLOAT128
 */

#define USE_MPFR

const int MPFR_PRECISION = 1024;  // bit precision of MPFR

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <sstream>
#ifdef USE_MPFR
#include "mpfrreal.hpp"
#else
#ifdef USE_FLOAT128
extern "C" {
  #include <quadmath.h>
}
#endif
#endif

#ifdef USE_MPFR
  typedef mpfr::real<MPFR_PRECISION> FLOAT;
  typedef mpfr::real<MPFR_PRECISION> FLOAT128;
  #define MYSQRT sqrt
  #define MYABS fabs

  inline double PrintFloat(FLOAT x) {
    std::stringstream s;
    s << std::setprecision(20) << std::scientific << x;
    double ret;
    sscanf(s.str().c_str(), "%lf", &ret);
    return ret;
  }
#else
  #ifdef USE_FLOAT128
    typedef __float128 FLOAT;
    typedef __float128 FLOAT128;
    #define MYSQRT sqrtq
    #define MYABS fabsq
    inline double PrintFloat(FLOAT x) {
      return double(x);
    }
  #else
    typedef double FLOAT;
    typedef __float128 FLOAT128;
    #define MYSQRT sqrt
    #define MYABS fabs
    inline double PrintFloat(FLOAT x) {
      return x;
    }
  #endif
#endif
 
// inline FLOAT Sqr(FLOAT x) {
//   return x * x;
// }
// 
// #define formatf(len, prec, value) setw(len) << setprecision(prec) << fixed << value
// #define formate(len, prec, value) setw(len) << setprecision(prec) << scientific << value
// #define formatd(len, value) setw(len) << fixed << value
// #define formats(len, value) setw(len) << value
// 
#endif
