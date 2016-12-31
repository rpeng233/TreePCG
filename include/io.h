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
#include <sstream>
#include "common.h"
#include "matrix.h"
#include "graph.h"

namespace IO {

  const char ASCII = 0;
  const char BINARY = 1;
  const char DIMACS = 2;

  char ParseFormat(string format) {
    if (format[0] == 'a' || format[0] == 'A') {
      return ASCII;
    } else if (format[0] == 'b' || format[0] == 'B') {
      return BINARY;
    } else if (format[0] == 'd' || format[0] == 'D') {
      return DIMACS;
    } else {
      fprintf(stderr, "UNSUPPORTED FORMAT: %s\n", format.c_str());
      exit(0);
    }
  }

  FILE *OpenAsRead(string file_name) {
    FILE *file_in;
    if (file_name == "__stdin") {
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

  FILE *OpenAsWrite(string file_name) {
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

  void CloseInputFile(FILE *file_in) {
    if(file_in != stdin) {
      fclose(file_in);
    }
  }
 
  void CloseOutputFile(FILE *file_out) {
    if(file_out != stdout) {
      fclose(file_out);
    }
  }

  void SkipHeader(FILE *file_in, char skip) {
    char buff[100], ch;
    
    while(true) {
      ch = getc(file_in);
      ungetc(ch, file_in);
      if (ch != skip) {
        break;
      }
      fgets(buff, 100, file_in);
      fprintf(stderr, "COMMENT: %s", buff);
    }
  }

  void SkipToken(FILE *file_in) {
    char buff[100];
    fscanf(file_in, "%s", buff);
  }

 
  inline int ReadInt(FILE *file_in, int format) {
    int result;
    if (format == ASCII || format == DIMACS) {
      fscanf(file_in, "%d", &result);
    } else if (format == BINARY) {
      fread(&result, sizeof(int), 1, file_in);
    } else {
      assert(0);
    }
    return result;
  }

  inline FLOAT ReadFloat(FILE *file_in, int format) {
    double result;
    if (format == ASCII || format == DIMACS) {
      fscanf(file_in, "%lf", &result);
    } else if (format == BINARY) {
      fread(&result, sizeof(double), 1, file_in);
    } else {
      assert(0);
    }
    return FLOAT(result);
  }

  void WriteFloat(FILE *file_out, int format, FLOAT x) {
    if (format == ASCII || format == DIMACS) {
      fprintf(file_out, "%.16lf", PrintFloat(x));
    } else if (format == BINARY) {
      fwrite(&x, sizeof(double), 1, file_out);
    } else {
      assert(0);
    }
  }

  inline void WriteSpace(FILE *file_out, int format) {
    if (format == ASCII || format == DIMACS) {
      fprintf(file_out, " ");
    }
  }

  inline void WriteNewLine(FILE *file_out, int format) {
    if (format == ASCII || format == DIMACS) {
      fprintf(file_out, "\n");
    }
  }
 
  inline void WriteInt(FILE *file_out, char format, int v) {
    if (format == ASCII || format == DIMACS) {
      fprintf(file_out, "%d", v);
    } else if (format == BINARY) {
      fwrite(&v, sizeof(int), 1, file_out);
    }
  }


  Matrix GraphToMatrix(const Graph graph) {
    Matrix result(graph.n, graph.n);

    for (int u = 0; u < graph.n; ++u) {
      for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
           it != graph.neighbor_list[u].end(); ++it) {
        FLOAT weight = FLOAT(1.0) / it -> resistance;
// off-diagonals
        result.addNonZero(u, it -> v, -weight);
        result.addNonZero(it -> v, u, -weight);
// diagonals
        result.addNonZero(u, u, weight);
        result.addNonZero(it -> v, it -> v, weight);
      }
    }

    result.sortAndCombine();
    return result;
  }


  Graph MatrixToGraph(Matrix &matrix) {
  // puts in arcs whose weights are the
  // average of the upper off-diagonals
  //
  // also converts it to resistance
  //

    assert(matrix.n == matrix.m);
    Graph result(matrix.n);

    for (vector<MatrixElement>::iterator it = matrix.non_zero.begin();
        it != matrix.non_zero.end(); ++it) {
      if (it -> row != it -> column) {
        (result.neighbor_list[it -> row]).
          push_back(Arc(it -> column, -FLOAT(1.0) / it -> value));
      }
    }

    return result;
  }


  Graph ReadGraph(string file_name, char format) {
    FILE* file_in = OpenAsRead(file_name);
    int num_vertices;
    int num_edges;

    if(format == DIMACS) {
      SkipHeader(file_in, 'c');
      char line[1234];
      fgets(line, 100, file_in);
      std::istringstream sin(line);
      string temp1, temp2;
      sin >> temp1 >> temp2 >> num_vertices >> num_edges;
      SkipHeader(file_in, 'n');  // ignore source/sink
    }
    else {
      num_vertices = ReadInt(file_in, format);
      num_edges = ReadInt(file_in, format);
    }

// fprintf(stderr, "V = %d E = %d\n", num_vertices, num_edges); fflush(stderr);

    Graph result(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
      result.neighbor_list[i].clear();
    }

    for (int i = 0; i < num_edges; ++i) {
      if(format == DIMACS) {
        SkipToken(file_in);
      }
// convert out of 1-indexing
      int row = ReadInt(file_in, format) - 1;
      int column = ReadInt(file_in, format) - 1;
// read float or int depending on format
      FLOAT weight = ReadFloat(file_in, format);

// convert to resistances
      FLOAT resistance = FLOAT(1.0) / weight;

// add in undirected edge, even if the input is directed
      result.AddEdge(row, column, resistance);

// fprintf(stderr, "%d %d %lf\n", row, column, resistance); fflush(stderr);
    }

    CloseInputFile(file_in);
    return result;
  }

  void WriteGraph(string file_name, char format, const Graph &graph) {
    FILE* file_out = OpenAsWrite(file_name);

    int num_edges = 0;
    for (int u = 0; u < graph.n; ++u) {
      num_edges += graph.neighbor_list[u].size();
    }
    num_edges /= 2;

    if(format == DIMACS) {
      fprintf(file_out, "c DIMACS format generated from graphToGraph\n");
      fprintf(file_out, "c   missing source / sink vertices\n");
      fprintf(file_out, "p max ");
    }

    WriteInt(file_out, format, graph.n);
    WriteSpace(file_out, format);
    WriteInt(file_out, format, num_edges);
    WriteNewLine(file_out, format);

    for (int u = 0; u < graph.n; ++u) {
      for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
           it != graph.neighbor_list[u].end(); ++it) {
        if (format == ASCII || format == BINARY){
          if(u < it -> v) {
            WriteInt(file_out, format, u + 1);
            WriteSpace(file_out, format);
            WriteInt(file_out, format, it -> v + 1);
            WriteSpace(file_out, format);
//convert back to weights
            FLOAT weight = FLOAT(1) / it -> resistance;
            WriteFloat(file_out, format, weight);
            WriteNewLine(file_out, format);
          }
        } else if (format == DIMACS) {
          fprintf(file_out, "a %d %d %.16lf\n", u + 1, it -> v + 1,
              PrintFloat(FLOAT(1) / it -> resistance));
        }
      }
    }
    CloseOutputFile(file_out);
  }

 
  Matrix ReadMMMatrix(string file_name) {
    FILE* file_in = OpenAsRead(file_name);
    SkipHeader(file_in, ASCII);

    int n = ReadInt(file_in, 0);
    int m = ReadInt(file_in, 0);
    int nnz = ReadInt(file_in, 0);

    Matrix result = Matrix(n, m);
    for (int i = 0; i < nnz; ++i) {
      int row = ReadInt(file_in, 0) - 1;
      int column = ReadInt(file_in, 0) - 1;
      FLOAT value = ReadFloat(file_in, 0);
      result.addNonZero(row, column, value);
    }
    CloseInputFile(file_in);

    result.sortAndCombine();
    return result;
  }

  Vec ReadMMVec(string file_name) {
    FILE* file_in = fopen(file_name.c_str(), "r");
    SkipHeader(file_in, '%');

    int n = ReadInt(file_in, 0);
    int m = ReadInt(file_in, 0);
    assert(m == 1);

    Vec result(n);

    for (int i = 0; i < n; ++i) {
      result[i] = ReadFloat(file_in, 0);
    }
    CloseInputFile(file_in);
    return result;
  }

  void WriteMMVec(const Vec *vec, string file_name) {
    FILE* file_out = OpenAsWrite(file_name);

    fprintf(file_out, "%%%%MatrixMarket matrix array real general\n");
    fprintf(file_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");
    fprintf(file_out, "%d %d\n", vec -> n, 1);
    for (int i =0 ; i < vec -> n; ++i) {
      fprintf(file_out, "%.16lf\n", PrintFloat((*vec)[i]));
    }
    CloseOutputFile(file_out);
  }

  void WriteMMMatrix(string file_name, Matrix &matrix) {
  //warning: this `reorders' the matrix
    matrix.sortAndCombine();
    FILE* file_out = OpenAsWrite(file_name);

    fprintf(file_out, "%%%%MatrixMarket matrix coordinate real symmetric\n");
    fprintf(file_out, "%% Generated \?\?-\?\?-\?\?\?\?\n");

    fprintf(file_out, "%d %d %d\n", matrix.n,
      matrix.m, int(matrix.non_zero.size()));
    for (vector<MatrixElement>::iterator it = matrix.non_zero.begin();
        it != matrix.non_zero.end(); ++it) {
      fprintf(file_out, "%d %d ", it -> column + 1, it -> column + 1);
      fprintf(file_out, "%.16lf\n", PrintFloat(it -> value));
    }
    CloseOutputFile(file_out);
  }

  vector<int> ReadVectorInt(string file_name, char format, int n) {
    FILE* file_in = OpenAsRead(file_name);
    vector<int> result(n);
    for (int i = 0; i < n; ++i) {
      // Converting from 1-indexing
      result[i] = ReadInt(file_in, format) - 1;
    }
    CloseOutputFile(file_in);
    return result;
  }


  void WriteVectorInt(string file_name, char format, const vector<int> v) {
    FILE* file_out = OpenAsWrite(file_name);
    for (int i = 0; i < int(v.size()); ++i) {
      // Converting to 1-indexing
      WriteInt(file_out, format, v[i] + 1);
      WriteNewLine(file_out, format); 
    }
    CloseOutputFile(file_out);
  }

}      // namespace IO
#endif  // GENERATORS_IO_H_
