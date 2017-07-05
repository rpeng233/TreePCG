#include <cstdio>

const int n = 100;
const double e = 1.2;

int main() {
  double v = 1;
  printf("%%%% %d-by-%d diagonal matrxi with entries set to %lf^i\n", n, n, e);
  printf("%d %d %d\n", n, n, n);
  for(int i = 0; i < n; ++i) {
    v = v * e;
    printf("%d %d %lf\n", i + 1, i + 1, v);
  }
  return 0;
}
