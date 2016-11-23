/********************************************************************************
 * 
 * File IO Interface
 * 
 * We assume that the files are 1-indexed.
 * HOWEVER, all intermediate storages are 0-indexed
 *
 *
 * Note:
 *   matrices store conductances, or weights
 *   graphs store resistances
 *   BUT upon input, everything is in conductances
 *
 *
 *
 * Public API:
 *     Mat IO::ConstructMatrixFromGraph(g):
 *         GraphSP g: the graph
 *         returns a Mat, the Laplacian matrix of graph g
 * 
 *     GraphSP IO::ConvertLaplacianMatrixToGraph(const Mat &A)
 *         Mat A: the Laplacian matrix
 *         returns the graph corresponding to the Laplacian matrix
 *         FIXME: 
 *             Need to check input is indeed Laplacian matrix. 
 *             Currently I'm not doing this due tofloating error issues.
 *             Richard: this is indeed tricky due to floating point issues
 * 
 *     GraphSP IO::ReadGraph(filename):
 *         string filename: the input file name
 *         Reads input, returns a graph object
 *   
 *     Mat IO::ReadMMMat(filename):
 *         string filename: the input file name
 *         Reads the Matrix Market format input, and returns the matrix.
 *        It assumes that the input contains only the upper-triangular part of the matrix.
 * 
 *     Mat IO::ReadMMVec(filename):
 *         string filename: the input file name
 *         Reads the Matrix Market format input, and returns the vector.
 *         It assumes that the input is a dense matrix market format file,
 *         with a single column.
 *
 *
 *
 * FILE FORMATS:
 *
 *  0: Asc
 *  1: binary of ints
 * 
 ********************************************************************************/

#ifndef GENERATORS_IO_H__
#define GENERATORS_IO_H__

#include <utility>
#include <string>
#include <vector>
#include "common.h"
#include "matrix.h"
#include "graph.h"

namespace IO {

  FILE *OpenAsRead(string file_name) {
    FILE *file_in;
    if (ile_name == "__stdin") {
      file_in = stdin;
    } else {
      file_in = fopen(file_name.c_str(), "r");

      if (!file_in) {
        printf("File error.\n");
        assert(0);
      }
    }
    return file_in;
  }

  FILE * OpenAsWrite(string file_name) {
    FILE *file_out;
    if (file_name == "__stdout") {
      file_out = stdout;
    } else {
      file_out = fopen(file_name.c_str(), "w");

      if (!file_out) {
        printf("File error.\n");
        assert(0);
      }
    }
    return file_out;
  }

  void SkipHeader(FILE *file_in) {
    char buff[100], ch;

    do {
      fgets(buff, 100, file_in);
      ch = getc(file_in);
      ungetc(ch, file_in);
    } while (ch == '%');
  }

  inline int ReadInt(FILE *file_in, int format) {
    int result;
    if (format == 0) {
      fscanf(file_in, "%d", &result);
    } else if (format == 1) {
      fread(&result, sizeof(int), 1, file_in);
    }
    return result;
  }

  inline FLOAT ReadFloat(FILE *file_in, int format) {
    double result;
    if (format == 0) {
      fscanf(file_in, "%lf", &result);
    } else if (format == 1) {
      fread(&result, sizeof(double), 1, file_in);
    }
    return FLOAT(result);
  }

  Matrix GraphToMatrix(const Graph graph) {
    Matrix result(graph.n, graph.n);

    for (int u = 0; u < graph.n; ++u) {
      for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
           it != graph.neighbor_list[u].end(); ++it) {
        FLOAT weight = FLOAT(1.0) / it -> resistance;
// off-diagonals
        result.AddNonZero(u, it -> v, -weight);
        result.AddNonZero(it -> v, u, -weight);
// diagonals
        result.AddNonZero(u, u, weight);
        result.AddNonZero(it -> v, it -> v, weight);
      }
    }

    result.SortAndCombine();
    return result;
  }


  Graph MatrixToGraph(const Matrix &matrix) {
  // puts in arcs whose weights are the
  // average of the upper off-diagonals
  //
  // also converts it to resistance
  //

    assert(matrix.n == matrix.m);
    Graph result(matrix.n);

    for (vector<MatrixElement>::iterator it = matrix.non_zero -> begin();
        it != matrix.non_zero -> end(); ++it) {
      if (it -> row != it -> column) {
        (result.neighbor_list[it -> row]).
          push_back(Arc(it -> column, -FLOAT(1.0) / it -> value));
      }
    }

    return result;
  }


  Graph ReadGraph(string file_name, int format) {
// format:
//   0: Asc
//   1: binary of ints

    FILE* file_in = fopen(file_name.c_str(), "r");

    int num_vertices = ReadInt(file_in, format);
    int num_edges = ReadInt(file_in, format);

    Graph result(num_vertices);

    for (int i = 0; i < num_vertices; ++i) {
      result.neighbor_list[i].clear();
    }

    for (int i = 0; i < num_edges; ++i) {
// convert out of 1-indexing
      int row = ReadInt(file_in, format) - 1;
      int column = ReadInt(file_in, format) - 1;

// read float or int depending on format
      FLOAT weight;
      if (format == 0) {
        weight = ReadFloat(file_in, format);
      } else if (format == 1) {
        weight = ReadInt(file_in, format);
      } else {
        assert(0);
      }
// convert to resistances
      FLOAT resistance = FLOAT(1.0) / weight;
      result.neighbor_list[row].push_back(Arc(column, resistance));
      result.neighbor_list[column].push_back(Arc(row, resistance));
    }
    fclose(file_in);

    return result;
  }

  Matrix ReadMMMatrix(string file_name) {
    FILE* file_in = OpenAsRead(file_name);
    SkipHeader(file_in);

    int n = ReadInt(file_in, 0);
    int m = ReadInt(file_in, 0);
    int nnz = ReadInt(file_in, 0);

    Matrix result = Matrix(n, m);
    for (int i = 0; i < nnz; ++i) {
      int row = ReadInt(file_in, 0) - 1;
      int column = ReadInt(file_in, 0) - 1;
      FLOAT value = ReadFloat(file_in, 0);
      result.AddNonZero(row, column, value);
    }
    fclose(file_in);

    result.SortAndCombine();
    return result;
  }

  Vec readMMVec(string file_name) {
    FILE* file_in = fopen(file_name.c_str(), "r");
    SkipHeader(file_in);

    int n = ReadInt(file_in, 0);
    int m = ReadInt(file_in, 0);
    assert(m == 1);

    Vec result(n);

    for (int i = 0; i < n; ++i) {
      result[i] = ReadFloat(file_in, 0);
    }
    fclose(file_in);
    return result;
  }

  void SaveMMVec(const Vec *vec, string file_name) {
    FILE* file_out = OpenAsWrite(file_name);

    fprintf(file_out, "%%%%MatrixMarket matrix array real general\n");
    fprintf(file_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");
    fprintf(file_out, "%d %d\n", vec -> n, 1);
    for (int i =0 ; i < vec -> n; ++i) {
      fprintf(file_out, "%.16lf\n", PrintFloat((*vec)[i]));
    }
    fclose(file_out);
  }

  void SaveMMMatrix(const Matrix matrix, string file_name) {
    matrix.SortAndCombine();
    FILE* file_out = OpenAsWrite(file_name);

    fprintf(file_out, "%%%%MatrixMarket matrix coordinate real symmetric\n");
    fprintf(file_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");

    fprintf(file_out, "%d %d %d\n", matrix.n,
      matrix.m, int((matrix.non_zero) -> size()));
    for (vector<MatrixElement>::iterator it = (matrix.non_zero) -> begin();
        it != (matrix.non_zero) -> end(); ++it) {
      fprintf(file_out, "%d %d ", it -> column + 1, it -> column + 1);
      fprintf(file_out, "%.16lf\n", PrintFloat(it -> value));
    }
    fclose(file_out);
  }

}      // namespace IO
#endif  // GENERATORS_IO_H_
