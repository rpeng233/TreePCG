#include <cstdio>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  int n = atoi(argv[1]);
  double e = atof(argv[2]);
  double v = 1;
  printf("%%%% %d-by-%d diagonal matrxi with entries set to %lf^i\n", n, n, e);
  printf("%d %d %d\n", n, n, n);
  for(int i = 0; i < n; ++i) {
    v = v * e;
    printf("%d %d %lf\n", i + 1, i + 1, v);
  }
  return 0;
}
