/*
 * generates Cayley graphs on n vertices with a number of generators.
 *
 * Input parameters:
 *   first one: n
 *   then 2n pairs, each containing (skip amount, weight)
 *
 *
 * Written by Richard on Oct 26, 16
 *
 *
*/

#include <cstdio>
#include <cstdlib>

int n, k;


int main(int argt, char **args) {
  if (argt <= 1 || (argt - 2) % 2 == 1) {
    fprintf(stderr, "INCORRECT # OF ARGS,");
    fprintf(stderr, "should be n followed by 2k info about generators\n");
    exit(0);
  }

  k = (argt - 2) / 2;
  sscanf(args[1], "%d", &n);
  printf("%d %d\n", n, k * n);

  fprintf(stderr, "MAKING CAYLEY GRAPH ON %d VERTICES");
  fprintf(stderr, "WITH %d GENERATORS\n", n, k);

  for (int i = 0; i < k; ++i) {
    int skip, weight;

    sscanf(args[2 + i * 2], "%d", &skip);
    sscanf(args[3 + i * 2], "%d", &weight);

    fprintf(stderr, "EDGE TYPE %d: SKIP = %d, WEIGHT = %d\n", i, skip, weight);


    for (int j = 0; j < n; ++j) {
      printf("%d %d %d\n", j + 1, (j + skip) % n + 1, weight);
    }
  }
  return 0;
}
