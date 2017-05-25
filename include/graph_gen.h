#ifndef INCLUDE_GRAPH_GEN_H_
#define INCLUDE_GRAPH_GEN_H_

struct One {
  FLOAT operator() () const {
    return FLOAT(1);
  }
};

template <typename EdgeT, typename WeightGen>
void line(size_t n, EdgeList<EdgeT>& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n - 1; i++) {
    es.AddEdge(EdgeT(i, i + 1, wgen()));
  }
}

template <typename EdgeT>
void line(size_t n, EdgeList<EdgeT>& es) {
  line(n, es, One());
}

template <typename EdgeT, typename WeightGen>
void cycle(size_t n, EdgeList<EdgeT>& es, WeightGen wgen=WeightGen()) {
  line(n, es, wgen);
  es.AddEdge(EdgeT(0, n - 1, wgen()));
}

template <typename EdgeT>
void cycle(size_t n, EdgeList<EdgeT>& es) {
  cycle(n, es, One());
}

template <typename EdgeT, typename WeightGen>
void grid2(size_t n, size_t m, EdgeList<EdgeT>& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n * m;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      if (j < m - 1) {
        es.AddEdge(EdgeT(i * m + j, i * m + j + 1, wgen()));
      }
      if (i < n - 1) {
        es.AddEdge(EdgeT(i * m + j, (i + 1) * m + j, wgen()));
      }
    }
  }
}

template <typename EdgeT>
void grid2(size_t n, size_t m, EdgeList<EdgeT>& es) {
  grid2(n, m, es, One());
}

template <typename EdgeT, typename WeightGen>
void grid3(size_t n1, size_t n2, size_t n3,
           EdgeList<EdgeT>& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n1 * n2 * n3;
  for (size_t i = 0; i < n1; i++) {
    for (size_t j = 0; j < n2; j++) {
      for (size_t k = 0; k < n3; k++) {
        if (i < n1 - 1) {
          es.AddEdge(EdgeT(
              i * n2 * n3 + j * n3 + k,
              (i + 1) * n2 * n3 + j * n3 + k,
              wgen()
          ));
        }
        if (j < n2 - 1) {
          es.AddEdge(EdgeT(
              i * n2 * n3 + j * n3 + k,
              i * n2 * n3 + (j + 1) * n3 + k,
              wgen()
          ));
        }
        if (k < n3 - 1) {
          es.AddEdge(EdgeT(
              i * n2 * n3 + j * n3 + k,
              i * n2 * n3 + j * n3 + k + 1,
              wgen()
          ));
        }
      }
    }
  }
}

template <typename EdgeT>
void grid3(size_t n1, size_t n2, size_t n3, EdgeList<EdgeT>& es) {
  grid3(n1, n2, n3, es, One());
}

template <typename EdgeT, typename WeightGen>
void gnp(size_t n, double p, EdgeList<EdgeT>& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if ((double) rand() / RAND_MAX < p) {
        es.AddEdge(EdgeT(i, j, wgen()));
      }
    }
  }
}

template <typename EdgeT>
void gnp(size_t n, double p, EdgeList<EdgeT>& es) {
  gnp(n, p, es, One());
}

inline
void cayley(size_t n,
            const std::vector<size_t>& skips,
            const std::vector<FLOAT>& resistances,
            EdgeList<EdgeR>& es) {
  es.Clear();
  es.n = n;
  assert(skips.size() == resistances.size());
  for (size_t i = 0; i < skips.size(); i++) {
    size_t skip = skips[i];
    FLOAT r = resistances[i];
    for (size_t j = 0; j < n; j++) {
      es.AddEdge(EdgeR(j, (j + skip) % n, r));
    }
  }
}

template <typename EdgeT>
void recursive_c_helper(EdgeList<EdgeT>& es,
                        size_t N, size_t i, size_t j, size_t n, size_t m) {
  if (n == 0 && m == 0) return;

  if (m == 0) {
    size_t ii = i + n / 2 + 1;
    recursive_c_helper(es, N, i, j, n / 2, m);
    recursive_c_helper(es, N, ii, j, n - n / 2 - 1, m);
    es.AddEdge(EdgeT((ii - 1) * N + j, ii * N + j, 1));
    return;
  }

  if (n == 0) {
    size_t jj = j + m / 2 + 1;
    recursive_c_helper(es, N, i, j, n, m / 2);
    recursive_c_helper(es, N, i, jj, n, m - m / 2 - 1);
    es.AddEdge(EdgeT(i * N + jj - 1, i * N + jj, 1));
    return;
  }

  size_t ii = i + n / 2 + 1;
  size_t jj = j + m / 2 + 1;
  recursive_c_helper(es, N, i, j, n / 2, m / 2);
  recursive_c_helper(es, N, ii, j, n - n / 2 - 1, m / 2);
  recursive_c_helper(es, N, i, jj, n / 2, m - m / 2 - 1);
  recursive_c_helper(es, N, ii, jj, n - n / 2 - 1, m - m / 2 - 1);
  es.AddEdge(EdgeT((ii - 1) * N + (jj - 1), (ii - 1) * N + jj, 1));
  es.AddEdge(EdgeT((ii - 1) * N + (jj - 1), ii * N + jj - 1, 1));
  es.AddEdge(EdgeT(ii * N + jj, ii * N + jj - 1, 1));
}

template <typename EdgeListType>
void recursive_c(size_t n, size_t m, EdgeListType& es) {
  es.Clear();
  es.n = n * m;
  recursive_c_helper(es, n, 0, 0, n - 1, m - 1);
}

#endif  // INCLUDE_GRAPH_GEN_H_
