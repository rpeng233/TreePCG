#include "graph.h"
#include "io.h"

int main(void) {
  EdgeListR es;
  
  ReadEdgeList(stdin, es);
  WriteMtx(stdout, es);

  return 0;
}
