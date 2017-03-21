/********************************************************************************
 *
 * Graph Definitions
 *
 * API Interface:
 *   struct Graph: definition of a graph (0-INDEXED).
 *     int .n: the size of graph
 *     vector<ArcR> .neighbor[i]:
 *         the neighbors of vertex i, 0 <= i <n
 *     .freeMemory(): destroys graph and frees memory
 *
 *   struct TreePlusEdges: definition of a connected graph with
 *        a spanning tree (0-INDEXED), rooted at .root.
 *     int .n: the size of graph
 *     vector<ArcR> .children[i]:
 *        the children of vertex i in spanning tree, 0 <= i < n
 *        the tree is DIRECTED, and 0 is always the root
 *
 *     vector<EdgeR> .off_tree_edges: the list off-tree edges
 *     .freeMemory(): destroys graph and frees memory
 *
 * NOTE:
 *   The `edge weights' in Graph and TreePlusEdges
 *   are always RESISTANCE, not WEIGHT!!!
 *
 *   Graph and GraphSP behave like objects in python or javascript. That is, if you do
 *     Graph A; Graph B=A;
 *   Then A and B will be actually pointing to the same object,
 *   i.e. modifying B will result in A being modified as well!
 *
 *   However, C++ does not have an automatic garbage collection system,
 *   so you have to run .freeMemory() to free memory of graphs that are no longer needed.
 *
 *   In Graph, every edge appears twice, in both direction.
 *   In TreePlusEdges, however, every edge appears in only one direction.
 *
 ********************************************************************************/

#ifndef INCLUDE_GRAPH_H_
#define INCLUDE_GRAPH_H_

#include <limits>
#include <map>
#include <queue>
#include <vector>
#include "common.h"

struct EdgeC;

struct EdgeR {
  size_t u, v;
  FLOAT resistance;

  EdgeR() {
    u = 0;
    v = 0;
    resistance = 0;
  }

  EdgeR(size_t u_, size_t v_, FLOAT resistance_) {
    u = u_;
    v = v_;
    resistance = resistance_;
  }

  EdgeR& operator=(const EdgeR& e) {
    u = e.u;
    v = e.v;
    resistance = e.resistance;
    return *this;
  }

  EdgeR& operator=(const EdgeC& e);

  FLOAT Resistance() const {
    return resistance;
  }

  FLOAT Conductance() const {
    return 1.0 / resistance;
  }

  void SetResistance(FLOAT r) {
    resistance = r;
  }

  void SetConductance(FLOAT c) {
    resistance = 1.0 / c;
  }

  bool operator <(const EdgeR &o) const {
    if (this->u != o.u) {
      return this->u < o.u;
    }

    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->resistance < o.resistance;
  }
};

struct EdgeC {
  size_t u, v;
  FLOAT conductance;

  EdgeC() {
    u = 0;
    v = 0;
    conductance = 0;
  }

  EdgeC(size_t u_, size_t v_, FLOAT conductance_) {
    u = u_;
    v = v_;
    conductance = conductance_;
  }

  EdgeC& operator=(const EdgeR& e);

  FLOAT Resistance() const {
    return 1.0 / conductance;
  }

  FLOAT Conductance() const {
    return conductance;
  }

  void SetResistance(FLOAT r) {
    conductance = 1.0 / r;
  }

  void SetConductance(FLOAT c) {
    conductance = c;
  }

  bool operator <(const EdgeC &o) const {
    if (this->u != o.u) {
      return this->u < o.u;
    }

    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->conductance < o.conductance;
  }
};

inline EdgeR& EdgeR::operator=(const EdgeC& e) {
  u = e.u;
  v = e.v;
  resistance = 1 / e.conductance;
  return *this;
}

inline EdgeC& EdgeC::operator=(const EdgeR& e) {
  u = e.u;
  v = e.v;
  conductance = 1 / e.resistance;
  return *this;
}

template <typename EdgeT>
struct EdgeList {
  size_t n;
  std::vector<EdgeT> edges;

  EdgeList() {
    n = 0;
  }

  template <typename E2>
  EdgeList(const EdgeList<E2>& es) {
    n = es.n;
    edges.resize(es.edges.size());

    for (size_t i = 0; i < edges.size(); i++) {
      edges[i] = es.edges[i];
    }
  }

  size_t Size() const {
    return edges.size();
  }

  // void AddEdge(size_t u, size_t v, FLOAT weight) {
  //   edges.push_back(EdgeT(u, v, weight));
  // }

  void AddEdge(const EdgeT& e) {
    edges.push_back(e);
  }

  void Clear() {
    n = 0;
    edges.clear();
  }

  void Reserve(size_t n) {
    edges.reserve(n);
  }

  void Resize(size_t n) {
    edges.resize(n);
  }

  void Free() {
    n = 0;
    edges = std::vector<EdgeT>();
  }

  EdgeT& operator[](const size_t i) {
    return edges[i];
  }

  const EdgeT& operator[](const size_t i) const {
    return edges[i];
  }
};

typedef EdgeList<EdgeC> EdgeListC;
typedef EdgeList<EdgeR> EdgeListR;

struct ArcR {
  size_t v;
  FLOAT resistance;

  ArcR() {
    v = 0;
    resistance = 0;
  }

  ArcR(int v_, FLOAT resistance_) {
    v = v_;
    resistance = resistance_;
  }

  ArcR& operator=(const EdgeR& e) {
    resistance = e.resistance;
    return *this;
  }

  FLOAT Resistance() const {
    return resistance;
  }

  FLOAT Conductance() const {
    return 1.0 / resistance;
  }

  void SetResistance(FLOAT r) {
    resistance = r;
  }

  void SetConductance(FLOAT c) {
    resistance = 1.0 / c;
  }

  bool operator <(const ArcR &o) const {
    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->resistance < o.resistance;
  }
};

struct ArcC {
  size_t v;
  FLOAT conductance;

  ArcC() {
    v = 0;
    conductance = 0;
  }

  ArcC(int v_, FLOAT conductance_) {
    v = v_;
    conductance = conductance_;
  }

  ArcC& operator=(const EdgeC& e) {
    conductance = e.conductance;
    return *this;
  }

  FLOAT Resistance() const {
    return 1.0 / conductance;
  }

  FLOAT Conductance() const {
    return conductance;
  }

  void SetResistance(FLOAT r) {
    conductance = 1.0 /r;
  }

  void SetConductance(FLOAT c) {
    conductance = 1.0 / c;
  }

  bool operator <(const ArcC &o) const {
    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->conductance < o.conductance;
  }
};

struct TreeVertexR {
  size_t parent;
  FLOAT parent_resistance;

  FLOAT ParentResistance() const {
    return parent_resistance;
  }

  FLOAT ParentConductance() const {
    return 1.0 / parent_resistance;
  }

  void SetParentR(size_t p, FLOAT r) {
    parent = p;
    parent_resistance = r;
  }

  void SetParentC(size_t p, FLOAT c) {
    parent = p;
    parent_resistance = 1.0 / c;
  }
};

template <typename VertexT>
struct Tree {
  size_t n;
  std::vector<VertexT> vertices;

  Tree() {
    n = 0;
  }

  Tree(size_t n_) {
    n = n_;
    vertices.resize(n);

    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  void Resize(size_t n_) {
    n = n_;
    vertices.clear();
    vertices.resize(n);

    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  void SetParentR(size_t v, size_t p, FLOAT r) {
    vertices[v].SetParentR(p, r);
  }

  void SetParentC(size_t v, size_t p, FLOAT c) {
    vertices[v].SetParentC(p, c);
  }

  VertexT& operator[](const size_t i) {
    return vertices[i];
  }

  const VertexT& operator[](const size_t i) const {
    return vertices[i];
  }
};

typedef Tree<TreeVertexR> TreeR;

template <typename ArcT>
struct AdjacencyArray {
  size_t n;
  std::vector<ArcT> arcs;
  std::vector<size_t> first_arc;

  AdjacencyArray() {
    n = 0;
  }

  template <typename EdgeT>
  AdjacencyArray(const EdgeList<EdgeT>& es) {
    BuildGraph(es);
  }

  template <typename EdgeT>
  void BuildGraph(const EdgeList<EdgeT>& es) {
    n = es.n;
    size_t m = es.edges.size();

    arcs.resize(m * 2);
    first_arc.resize(n + 1);
    std::vector<size_t> degrees(n);

    for (typename std::vector<EdgeT>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      degrees[it->u]++;
      degrees[it->v]++;
    }

    first_arc[0] = 0;
    for (size_t i = 0; i < n - 1; i++) {
      first_arc[i + 1] = degrees[i] + first_arc[i];
    }
    first_arc[n] = m * 2;

    size_t tmp_index;
    for (typename std::vector<EdgeT>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      degrees[it->u]--;
      tmp_index = first_arc[it->u] + degrees[it->u];
      arcs[tmp_index].v = it->v;
      arcs[tmp_index] = *it;

      degrees[it->v]--;
      tmp_index = first_arc[it->v] + degrees[it->v];
      arcs[tmp_index].v = it->u;
      arcs[tmp_index] = *it;
    }
  }

  void FreeMemory() {
    n = 0;
    arcs = std::vector<ArcT>();
    first_arc = std::vector<size_t>();
  }
};

struct AdjacencyMap {
  size_t n;
  std::vector<std::map<size_t, FLOAT> > neighbor_map;

  void AddEdgeR(size_t u, size_t v, const FLOAT& r) {
    neighbor_map[u][v] += 1 / r;
    neighbor_map[v][u] += 1 / r;
  }

  void AddEdgeR(const EdgeR& e) {
    neighbor_map[e.u][e.v] += 1 / e.resistance;
    neighbor_map[e.v][e.u] += 1 / e.resistance;
  }

  AdjacencyMap(const EdgeList<EdgeR>& es) {
    n = es.n;
    neighbor_map.resize(n);
    for (size_t i = 0; i < es.edges.size(); i++) {
      AddEdgeR(es[i]);
    }
  }

  template <typename VertexType>
  AdjacencyMap(const Tree<VertexType>& tree) {
    n = tree.n;
    neighbor_map.resize(n);
    for (size_t i = 0; i < n; i++) {
      size_t p = tree.vertices[i].parent;
      if (i == p) continue;
      const FLOAT& r = tree.vertices[i].ParentResistance();
      AddEdgeR(i, p, r);
    }
  }

};

template <typename TreeType, typename ArcT>
void DijkstraTree(const AdjacencyArray<ArcT>& graph, size_t root, TreeType& tree) {
  size_t n = graph.n;

  tree.Resize(n);

  std::vector<bool> finished(n, false);
  std::vector<double> dist(n, std::numeric_limits<double>::max());
  std::priority_queue<std::pair<double, size_t> > queue;

  dist[root] = 0;
  queue.push(std::make_pair(0, root));
  size_t count = 0;
  while (count < n) {
    size_t u = queue.top().second;
    queue.pop();
    if (finished[u]) {
      continue;
    }
    finished[u] = true;
    count++;
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const ArcR& arc = graph.arcs[i];
      double tmp = dist[u] + (double) arc.Resistance();
      if (tmp < dist[arc.v]) {
        dist[arc.v] = tmp;
        tree.SetParentR(arc.v, u, arc.Resistance());
        queue.push(std::make_pair(-tmp, arc.v));
      }
    }
  }
}

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
              i * n1 * n2 + j * n2 + k,
              (i + 1) * n1 * n2 + j * n2 + k,
              wgen()
          ));
        }
        if (j < n2 - 1) {
          es.AddEdge(EdgeT(
              i * n1 * n2 + j * n2 + k,
              i * n1 * n2 + (j + 1) * n2 + k,
              wgen()
          ));
        }
        if (k < n3 - 1) {
          es.AddEdge(EdgeT(
              i * n1 * n2 + j * n2 + k,
              i * n1 * n2 + j * n2 + k + 1,
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

// struct Graph {
//   int n;
//   std::vector<ArcR> *neighbor_list;
// 
//   Graph() {
//     n = 0;
//     neighbor_list = NULL;
//   }
// 
//   Graph(int n_) {
//     n = n_;
//     neighbor_list = new std::vector<ArcR>[n];
//   }
// 
//   Graph(const Graph &o) {
//     n = o.n;
//     neighbor_list = o.neighbor_list;
//   }
// 
//   Graph &operator =(const Graph & o) {
//     n = o.n;
//     neighbor_list = o.neighbor_list;
//     return (*this);
//   }
// 
//   void AddEdge(int u, int v, FLOAT resistance) {
//     neighbor_list[u].push_back(ArcR(v, resistance));
//     neighbor_list[v].push_back(ArcR(u, resistance));
//   }
// 
//   void AddEdge(EdgeR e) {
//     AddEdge(e.u, e.v, e.resistance);
//   }
// 
//   void SortAndCombine() const {
//     for (int u = 0; u < n; ++u) {
//       int new_degree = 0;
//       sort(neighbor_list[u].begin(), neighbor_list[u].end());
// 
//       int last_vertex = -1;
//       FLOAT last_weight = 0;
// 
//       for (std::vector<ArcR>::iterator it = neighbor_list[u].begin();
//           it != neighbor_list[u].end(); ++it) {
//         if (it -> v != last_vertex) {
//           if (last_vertex >= 0) {
//             neighbor_list[u][new_degree] =
//               ArcR(last_vertex, FLOAT(1) / last_weight);
//             new_degree++;
//           }
//           last_vertex = it -> v;
//         }
//         last_weight += FLOAT(1) / it -> resistance;
//       }
// 
//       if (last_vertex >= 0) {
//         neighbor_list[u][new_degree] =
//           ArcR(last_vertex, FLOAT(1) / last_weight);
//         new_degree++;
//       }
//       neighbor_list[u].resize(new_degree);
//     }
//   }
// 
// 
//   void FreeMemory() const {
//     delete neighbor_list;
//   }
// };

#endif  // INCLUDE_GRAPH_H_
