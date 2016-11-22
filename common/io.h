/********************************************************************************
 * 
 * File IO Interface
 * 
 * Public API:
 *     Mat IO::constructMatrixFromGraph(g):
 *         GraphSP g: the graph
 *         returns a Mat, the Laplacian matrix of graph g
 * 
 *     GraphSP IO::convertLaplacianMatrixToGraphSP(const Mat &A)
 *         Mat A: the Laplacian matrix
 *         returns the graph corresponding to the Laplacian matrix
 *         FIXME: 
 *             Need to check input is indeed Laplacian matrix. 
 *             Currently I'm not doing this bcuz of floating error issues.
 * 
 *     GraphSP IO::readGraph(filename):
 *         string filename: the input file name
 *        Reads input, calls the low stretch tree finder, and returns a GraphSP object.
 *   
 *     GraphSP IO::readGraphSP(filename):
 *         string filename: the input file name
 *         Reads input, returns the GraphSP object.
 *         It assumes that the first n-1 edges form a directed spanning tree rooted at 1.
 * 
 *     Mat IO::readMML(filename):
 *         string filename: the input file name
 *         Reads the Matrix Market format input, and returns the matrix.
 *        It assumes that the input contains only the upper-triangular part of the matrix.
 * 
 *     Mat IO::readMMA(filename):
 *         string filename: the input file name
 *         Reads the Matrix Market format input, and returns the matrix.
 *        It assumes that the input is the adjancy matrix of a graph, and contains only upper-triangular part.
 *         It then constructs the Laplacian matrix corresponding to the graph.
 * 
 *     Mat IO::readMMVec(filename):
 *         string filename: the input file name
 *         Reads the Matrix Market format input, and returns the vector.
 *        It assumes that the input is a dense matrix market format file, with a single column.
 *
 * 
 ********************************************************************************/

#ifndef GENERATORS_IO_H__
#define GENERATORS_IO_H__

#include <utility>
#include <tuple>
#include <string>
#include <vector>
#include "common.h"
#include "matrix.h"
#include "graph.h"
#include "treeFinder.h"

namespace IO {

  FILE * openAsRead(string fileName) {
    FILE *f_in = fopen(fileName.c_str(), "r");
    if (!f_in) {
      printf("File error.\n");
      assert(0);
    }
    return f_in;
  }

  FILE * openAsWrite(string fileName) {
    FILE *f_out = fopen(fileName.c_str(), "w");
    if (!f_out) {
      printf("File error.\n");
      assert(0);
    }
    return f_out;
  }

  Mat constructMatrixFromGraph(const GraphSP g) {
    int n = g.n;
    static FLOAT *s = new FLOAT[MAXN];

#ifdef USE_MPFR
    for (int i = 0; i < n; ++i) {
      s[i] = 0;
    }
#else
    memset(s, 0, sizeof(FLOAT)*(n+1));
#endif

    Mat A(n, n);
    for (int i = 1; i <= n; ++i) {
      for (auto it = g.e[i].begin(); it != g.e[i].end(); ++it) {
        int x = i;
        int y = it->first;
        FLOAT z = it->second;

        x--;
        y--;
        z = 1.0 / z;

        A.entryAddValue(x, y, -z);
        A.entryAddValue(y, x, -z);
        s[x] += z;
        s[y] += z;
      }
    }

    for (auto it = g.o.begin(); it != g.o.end(); ++it) {
      int x = get<0>(*it);
      int y = get<1>(*it);
      FLOAT z = get<2>(*it);

      x--;
      y--;
      z = 1.0/z;
      A.entryAddValue(x, y, -z);
      A.entryAddValue(y, x, -z);

      s[x] += z;
      s[y] += z;
    }

    for (int i = 0; i < n; ++i) {
      A.entryAddValue(i, i, s[i]);
    }
    A.sortup();
    return A;
  }


  GraphSP convertLaplacianMatrixToGraphSP(const Mat &A) {
  // WARNING: won't report error if matrix is not symmetric
    assert(A.n == A.m);
    Graph g(A.n);

    static FLOAT *s = new FLOAT[MAXN];
    for (auto it = A.values.begin(); it != A.values.end(); ++it) {
      int x = it->x;
      int y = it->y;
      FLOAT z = it->z;

      if (x != y) {
        if (z > 0) {
          fprintf(stderr, "Warning: Given matrix is not Laplacian");
          fprintf(stderr, "it has positive off-diagonals!\n");
        }
        s[x] += z;
        z = -1.0/z;
        g.e[x + 1].push_back(make_pair(y + 1, z));
      } else {
        if (z < 0) {
          fprintf(stderr, "Warning: Given matrix is not Laplacian,");
          fprintf(stderr, " it has negative diagonals!\n");
        }
        s[x]+=z;
      }
    }
    //  rep(i,0,A.n-1) if (fabs(s[i])>1e-6)
    //  printf("Warning: Given matrix is not Laplacian,
    //  fpritnf("its row/col sum is not 0!\n");
    //

    GraphSP g2 = TreeFinder::findLowStretchTree(g);
    g.freeMemory();
    return g2;
  }


  GraphSP readGraphSPAsc(string fileName) {
    FILE* f_in = openAsRead(fileName);

    int n, m;
    fscanf(f_in, "%d%d", &n, &m);

    vector< pair<int, FLOAT> > *e = new vector< pair<int, FLOAT> >[n + 1];
    vector< tuple<int, int, FLOAT> > o;
    for (int i = 1; i <= n; ++i) {
      e[i].clear();
    }

    for (int i = 1; i < n; ++i) {
      int x, y;
      double z;
      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      e[x + 1].push_back(make_pair(y + 1, z));
    }

    o.clear();
    for (int i = 1; i <= m - (n - 1); ++i) {
      int x, y;
      double z;
      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      o.push_back(make_tuple(x + 1, y + 1, z));
    }
    fclose(f_in);

    GraphSP g;
    g.n = n;
    g.e = e;
    g.o.swap(o);
    return g;
  }

  GraphSP readGraphAsc(string filename) {
    FILE* f_in = fopen(filename.c_str(), "r");
    if (!f_in) {
      printf("File error.\n");
      assert(0);
    }

    int n, m;
    fscanf(f_in, "%d%d", &n, &m);
    Graph h(n);

    for (int i = 1; i <= n; ++i) {
      h.e[i].clear();
    }

    for (int i = 1; i <= m; ++i) {
      int x, y;
      double z;

      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      h.e[x].push_back(make_pair(y, z));
      h.e[y].push_back(make_pair(x, z));
    }
    fclose(f_in);

    GraphSP g = TreeFinder::findLowStretchTree(h);
    h.freeMemory();
    return g;
  }

  GraphSP readGraphBin(string filename) {
    FILE* f_in = fopen(filename.c_str(), "r");
    if (!f_in) {
      printf("File error.\n");
      assert(0);
    }
    int buffer[3];

    fread(buffer, sizeof(int), 2, f_in);
    int n = buffer[0];
    int m = buffer[1];
    Graph h(n);

    for (int i = 1; i <= n; ++i) {
      h.e[i].clear();
    }

    for (int i = 0; i < m; ++i) {
      fread(buffer, sizeof(int), 3, f_in);

      int x = buffer[0];
      int y = buffer[1];
      int z = buffer[2];

      h.e[x].push_back(make_pair(y, z));
      h.e[y].push_back(make_pair(x, z));
    }
    fclose(f_in);

    GraphSP g = TreeFinder::findLowStretchTree(h);
    h.freeMemory();
    return g;
  }

  // TOFIX Only works on symmetric mm right now
  Mat readMML(string fileName) {
    FILE* f_in = openAsRead(fileName);

    int c;
    char buff[100];

    do {
      fgets(buff, 100, f_in);
      c = getc(f_in);
      ungetc(c, f_in);
    } while (c == '%');

    int n, temp, m;
    fscanf(f_in, "%d%d%d", &n, &temp, &m);

    int x, y;
    double z;

    Mat A = Mat(n, n);
    for (int i = 1; i <= m; ++i) {
      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      x--;
      y--;
      A.entryAddValue(x, y, z);
      if (x != y) {
        A.entryAddValue(y, x, z);
      }
    }
    fclose(f_in);

    A.sortup();
    return A;
  }

  // TOFIX Only works on symmetric mm right now
  Mat readMMA(string fileName) {
    FILE* f_in = openAsRead(fileName);
    int c;
    char buff[100];
    do {
      fgets(buff, 100, f_in);
      c = getc(f_in);
      ungetc(c, f_in);
    } while (c == '%');

    int n, temp, m;
    fscanf(f_in, "%d%d%d", &n, &temp, &m);
    static FLOAT *s = new FLOAT[MAXN];
#ifdef USE_MPFR
    for (int i = 0; i <= n; ++i) {
      s[i] = 0;
    }
#else
    memset(s, 0, sizeof(FLOAT) * (n + 1));
#endif
    int x, y;
    double z;

    Mat A = Mat(n, n);
    for (int i = 1; i <= m; ++i) {
      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      x--;
      y--;

      A.entryAddValue(x, y, -z);
      A.entryAddValue(y, x, -z);
      s[x] += z; s[y] += z;
    }

    fclose(f_in);

    for (int i = 0; i < n; ++i) {
      A.entryAddValue(i, i, s[i]);
    }
    A.sortup();
    return A;
  }

  Vec readMMVec(string filename) {
    FILE* f_in = fopen(filename.c_str(), "r");
    if (!f_in) {
      printf("File error.\n");
      assert(0);
    }

    int c;
    char buff[100];
    do {
      fgets(buff, 100, f_in);
      c = getc(f_in);
      ungetc(c, f_in);
    } while (c == '%');

    int n, m;
    fscanf(f_in, "%d%d", &n, &m);
    assert(m == 1);

    Vec v(n);

    double z;

    for (int i = 0; i < n; ++i) {
      fscanf(f_in, "%lf", &z);
      FLOAT zz = z;
      v[i] = zz;
    }
    fclose(f_in);
    return v;
  }

  void saveMMVec(const Vec &A, string filename) {
    FILE* f_out = fopen(filename.c_str(), "w");
    if (!f_out) {
      printf("Write file error.\n");
      assert(0);
    }

    fprintf(f_out, "%%%%MatrixMarket matrix array real general\n");
    fprintf(f_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");
    fprintf(f_out, "%d %d\n", A.n, 1);
    for (int i =0 ; i < A.n; ++i) {
      fprintf(f_out, "%.16lf\n", printFloat(A[i]));
    }
    fclose(f_out);
  }


  void saveMM(const Mat &A, string filename) {
    A.sortup();
    FILE* f_out = fopen(filename.c_str(), "w");
    if (!f_out) {
      printf("Write file error.\n");
      assert(0);
    }

    fprintf(f_out, "%%%%MatrixMarket matrix coordinate real symmetric\n");
    fprintf(f_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");

    int cnt = 0;
    for (auto it = A.values.begin(); it != A.values.end(); ++it) {
      if (it->x >= it->y) {
        cnt++;
      }
    }

    fprintf(f_out, "%d %d %d\n", A.n, A.m, cnt);

    for (auto it = A.values.begin(); it != A.values.end(); ++it) {
      if (it -> x >= it -> y) {
/*
        char temp[30];
        sprintf(temp, "%0.16lf", printFloat(it -> z));
        int r = strlen(temp);
        while(temp[r - 1] == '0') {
          r--;
        }
        temp[r] = '\0';
*/ 
        fprintf(f_out, "%d %d ", it -> x + 1, it -> y + 1);
        fprintf(f_out, "%0.16lf\n", printFloat(it -> z));
      }
    }
    fclose(f_out);
  }

  void saveMMnonsym(const Mat &A, string filename) {
    A.sortup();
    FILE* f_out = fopen(filename.c_str(), "w");
    if (!f_out) {
      printf("Write file error.\n");
      assert(0);
    }

    fprintf(f_out, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f_out, "%% Generated 01-Feb-2016\n");
    int cnt = 0;
//  richard: why is this not just size?
    for (auto it = A.values.begin(); it != A.values.end(); ++it) {
      cnt++;
    }

    fprintf(f_out, "%d %d %d\n", A.n, A.m, cnt);
    for (auto it = A.values.begin(); it != A.values.end(); ++it) {
      fprintf(f_out, "%d %d ", it -> x + 1, it -> y + 1);
      fprintf(f_out, "%.16lf\n", printFloat(it -> z));
    }
    fclose(f_out);
  }

  // assumes non-symmetric!
  Mat readMMnonsym(string filename) {
    FILE* f_in = fopen(filename.c_str(), "r");
    if (!f_in) {
      printf("File error.\n");
      assert(0);
    }

    int c;
    char buff[100];

    do {
      fgets(buff, 100, f_in);
      c = getc(f_in);
      ungetc(c, f_in);
    } while (c == '%');

    int n, m, nnz;
    fscanf(f_in, "%d%d%d", &n, &m, &nnz);

    int x, y;
    double z;

    Mat A = Mat(n, m);
    for (int i = 0; i < nnz; ++i) {
      fscanf(f_in, "%d%d%lf", &x, &y, &z);
      x--;
      y--;
      A.entryAddValue(x, y, z);
    }
    fclose(f_in);

    A.sortup();
    return A;
  }
}      // namespace IO
#endif  // GENERATORS_IO_H_
