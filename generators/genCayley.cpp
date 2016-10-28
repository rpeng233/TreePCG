/*
 * generates Cayley graphs on n vertices with a number of generators.
 *
 * Input parameters:
 *   first: output file name 
 *   second: n
 *   then 2n pairs, each containing (skip amount, weight)
 *
 * Output 
 * Written by Richard Peng on Oct 27, 16
 *
 *
*/

#include <cstdio>
#include <cstdlib>

int n, k;


int main(int argt, char **args) {
  if (argt <= 3 || (argt - 3) % 2 == 1) {
    fprintf(stderr, "INCORRECT # OF ARGS,");
    fprintf(stderr, "should be FILENAME,");
    fprintf(stderr, "then n followed by 2k info about generators\n");
    exit(0);
  }

  FILE *f_out = fopen(args[1], "w");
  int buffer[3];

  k = (argt - 3) / 2;
  sscanf(args[2], "%d", &n);


  buffer[0] = n;
  buffer[1] = k * n;
  fwrite(buffer, sizeof(int), 2, f_out);

  fprintf(stderr, "MAKING CAYLEY GRAPH ON %d VERTICES", n);
  fprintf(stderr, "WITH %d GENERATORS\n", k);

  for (int i = 0; i < k; ++i) {
    int skip, weight;

    sscanf(args[3 + i * 2], "%d", &skip);
    sscanf(args[4 + i * 2], "%d", &weight);

    fprintf(stderr, "EDGE TYPE %d: SKIP = %d, WEIGHT = %d\n", i, skip, weight);


    for (int j = 0; j < n; ++j) {
      buffer[0] = j + 1;
      buffer[1] = (j + skip) % n + 1;
      buffer[2] = weight;
      fwrite(buffer, sizeof(int), 3, f_out);
    }
  }
  fclose(f_out);
  return 0;
}
