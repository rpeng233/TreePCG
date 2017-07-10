#include <cstdio>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  int n = atoi(argv[1]);
  printf("%%%% length %d all ones vector\n", n);
  printf("%d %d\n", n, 1);
  for(int i = 0; i < n; ++i) {
    printf("1\n");
  }
  return 0;
}
