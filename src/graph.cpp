#include <vector>
#include "graph.h"


void cayley(size_t n,
            const std::vector<size_t>& skips,
            const std::vector<FLOAT>& resistances,
            EdgeList& es) {
  es.Clear();
  es.n = n;
  assert(skips.size() == resistances.size());
  for (size_t i = 0; i < skips.size(); i++) {
    size_t skip = skips[i];
    FLOAT r = resistances[i];
    for (size_t j = 0; j < n; j++) {
      es.AddEdge(j, (j + skip) % n, r);
    }
  }
}

void line(size_t n,
          const std::vector<FLOAT>& resistances,
          EdgeList& es) {
  es.Clear();
  es.n = n;
  assert(n - 1 <= resistances.size());
  for (size_t i = 0; i < n - 1; i++) {
    es.AddEdge(i, i + 1, resistances[i]);
  }
}

void line(size_t n,
          EdgeList& es) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n - 1; i++) {
    es.AddEdge(i, i + 1, (FLOAT) 1);
  }
}

void cycle(size_t n,
           const std::vector<FLOAT>& resistances,
           EdgeList& es) {
  assert(n == resistances.size());
  line(n, resistances, es);
  es.AddEdge(0, n - 1, resistances[n - 1]);
}

void cycle(size_t n,
           EdgeList& es) {
  line(n, es);
  es.AddEdge(0, n - 1, (FLOAT) 1);
}

void torus(size_t n, size_t m, EdgeList& es) {
  es.Clear();
  es.n = n * m;
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < m; j++) {
      es.AddEdge(i * m + j, i * m + (j + 1) % m, (FLOAT) 1);
      es.AddEdge(i * m + j, ((i + 1) % n) * m + j, (FLOAT) 1);
    }
  }
}

void gnp(size_t n, double p, EdgeList& es) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if ((double) rand() / RAND_MAX < p) {
        es.AddEdge(i, j, (FLOAT) 1);
      }
    }
  }
}
