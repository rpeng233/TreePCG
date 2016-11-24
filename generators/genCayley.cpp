/*
 * NOTE: this is a stand-alone code that does't use the graph
 * and IO libraries
 *
 * generates Cayley graphs on n vertices with a number of generators.
 *
 * Input parameters:
 *   first: output file name 
 *   second: n
 *   then 2n pairs, each containing (skip amount, weight)
 *
 * Output 
 *   Written by Richard Peng on Oct 27, 16
 *
 *
*/

#include <cstdio>
#include <cstdlib>
#include "../common/io.h"
#include "../common/graph.h"

int n, k;

int main(int argt, char **args) {
  if (argt <= 4 || (argt - 4) % 2 == 1) {
    fprintf(stderr, "INCORRECT # OF ARGS,");
    fprintf(stderr, "should be FILENAME, format,");
    fprintf(stderr, "then n followed by 2k info about generators\n");
    exit(0);
  }

  string out_file = args[1];
  char out_format = IO::ParseFormat(args[2]);

  k = (argt - 3) / 2;
  sscanf(args[3], "%d", &n);
  Graph result(n);

  fprintf(stderr, "Making Cayley graph on %d vertices\n", n);
  fprintf(stderr, "with %d generators\n", k);

  for (int i = 0; i < k; ++i) {
    int skip;
    double resistance;

    sscanf(args[4 + i * 2], "%d", &skip);
    sscanf(args[5 + i * 2], "%lf", &resistance);

    fprintf(stderr, "edge type %d:\n", i);
    fprintf(stderr, "  skip = %d, resistance = %0.16lf\n", skip, resistance);
    for (int j = 0; j < n; ++j) {
      result.AddEdge(j, (j + skip) % n, resistance);
    }
  }
  IO::WriteGraph(out_file, out_format, result);
  return 0;
}
