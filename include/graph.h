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

#include <map>
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

  EdgeR& operator=(const EdgeC& e);

  FLOAT Resistance() const {
    return resistance;
  }

  FLOAT Conductance() const {
    return 1 / resistance;
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

  FLOAT Resistance() {
    return 1 / conductance;
  }

  FLOAT Conductance() {
    return conductance;
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

  bool operator <(const ArcC &o) const {
    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->conductance < o.conductance;
  }
};

/*
struct Vertex {
  std::vector<ArcR> nghbrs;

  void AddArc(size_t v, FLOAT resistance) {
    nghbrs.push_back(ArcR(v, resistance));
  }

  void SortAndCombine() {
    if (nghbrs.empty()) return;
    sort(nghbrs.begin(), nghbrs.end());
    size_t new_degree = 0;
    size_t last_nghbr = nghbrs[0].v;
    FLOAT conductance = 0;

    for(std::vector<ArcR>::iterator ii = nghbrs.begin();
         ii != nghbrs.end(); ++ii) {
      if (ii -> v == last_nghbr) {
        conductance += 1 / ii -> resistance;
      } else {
        nghbrs[new_degree].v = last_nghbr;
        nghbrs[new_degree].resistance = 1 / conductance;
        new_degree++;
        last_nghbr = ii -> v;
        conductance = 1 / ii -> resistance;
      }
    }
    nghbrs[new_degree].v = last_nghbr;
    nghbrs[new_degree].resistance = 1 / conductance;
    new_degree++;
    nghbrs.resize(new_degree);
  }
};
*/

struct Graph {
  int n;
  std::vector<ArcR> *neighbor_list;

  Graph() {
    n = 0;
    neighbor_list = NULL;
  }

  Graph(int n_) {
    n = n_;
    neighbor_list = new std::vector<ArcR>[n];
  }

  Graph(const Graph &o) {
    n = o.n;
    neighbor_list = o.neighbor_list;
  }

  Graph &operator =(const Graph & o) {
    n = o.n;
    neighbor_list = o.neighbor_list;
    return (*this);
  }

  void AddEdge(int u, int v, FLOAT resistance) {
    neighbor_list[u].push_back(ArcR(v, resistance));
    neighbor_list[v].push_back(ArcR(u, resistance));
  }

  void AddEdge(EdgeR e) {
    AddEdge(e.u, e.v, e.resistance);
  }

  void SortAndCombine() const {
    for (int u = 0; u < n; ++u) {
      int new_degree = 0;
      sort(neighbor_list[u].begin(), neighbor_list[u].end());

      int last_vertex = -1;
      FLOAT last_weight = 0;

      for (std::vector<ArcR>::iterator it = neighbor_list[u].begin();
          it != neighbor_list[u].end(); ++it) {
        if (it -> v != last_vertex) {
          if (last_vertex >= 0) {
            neighbor_list[u][new_degree] =
              ArcR(last_vertex, FLOAT(1) / last_weight);
            new_degree++;
          }
          last_vertex = it -> v;
        }
        last_weight += FLOAT(1) / it -> resistance;
      }

      if (last_vertex >= 0) {
        neighbor_list[u][new_degree] =
          ArcR(last_vertex, FLOAT(1) / last_weight);
        new_degree++;
      }
      neighbor_list[u].resize(new_degree);
    }
  }


  void FreeMemory() const {
    delete neighbor_list;
  }
};

struct TreeVertexR {
  size_t parent;
  size_t degree; // TODO: right now the degree doesn't count the parent
  FLOAT parent_resistance;
  bool eliminated;

  TreeVertexR()
  : parent(0), degree(0), parent_resistance(0.0), eliminated(false)
  { }
};

struct TreeR {
  size_t n;
  std::vector<TreeVertexR> vertices;

  TreeR() {
    n = 0;
  }

  TreeR(size_t n_) {
    n = n_;
    vertices.resize(n);
    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  void SetParent(size_t v, size_t p, FLOAT r) {
    size_t old_p = vertices[v].parent;
    if (old_p != v) {
      vertices[old_p].degree--;
    }

    vertices[v].parent = p;
    vertices[v].parent_resistance = r;
    vertices[p].degree++;
  }
};

/*
struct TreePlusEdgesR {
  size_t n;
  std::vector<TreeVertexR> vertices;
  std::vector<EdgeR> off_tree_edges;

  TreePlusEdgesR() {
    n = 0;
  }

  TreePlusEdgesR(size_t n_) {
    n = n_;
    vertices.resize(n);
    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  TreePlusEdgesR(const TreePlusEdgesR &o) {
    n = o.n;
    vertices = o.vertices;
    off_tree_edges = o.off_tree_edges;
  }

  TreePlusEdgesR &operator =(const TreePlusEdgesR &o) {
    n = o.n;
    vertices = o.vertices;
    off_tree_edges = o.off_tree_edges;
    return (*this);
  }

  void SetParent(size_t v, size_t p, FLOAT r) {
    vertices[v].parent = p;
    vertices[v].parent_resistance = r;
    vertices[p].degree++;
  }

  void AddEdge(size_t u, size_t v, FLOAT r) {
    off_tree_edges.push_back(EdgeR(u, v, r));
    vertices[u].degree++;
    vertices[v].degree++;
  }
};
*/

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

  EdgeT& operator[](const size_t i) {
    return edges[i];
  }

  const EdgeT& operator[](const size_t i) const {
    return edges[i];
  }
};

/*
struct EdgeListC;

struct EdgeListR {
  size_t n;
  std::vector<EdgeR> edges;

  EdgeListR() {
    n = 0;
  }

  EdgeListR(const EdgeListC& es);

  size_t size() {
    return edges.size();
  }

  void AddEdge(size_t u, size_t v, FLOAT resistance) {
    edges.push_back(EdgeR(u, v, resistance));
  }

  void Clear() {
    n = 0;
    edges.clear();
  }

  EdgeR& operator[](const size_t i) {
    return edges[i];
  }

  const EdgeR& operator[](const size_t i) const {
    return edges[i];
  }
};

struct EdgeListC {
  size_t n;
  std::vector<EdgeC> edges;

  EdgeListC() {
    n = 0;
  }

  EdgeListC(const EdgeListR& es);

  size_t size() {
    return edges.size();
  }

  void AddEdge(size_t u, size_t v, FLOAT conductance) {
    edges.push_back(EdgeC(u, v, conductance));
  }

  void Clear() {
    n = 0;
    edges.clear();
  }

  EdgeC& operator[](const size_t i) {
    return edges[i];
  }

  const EdgeC& operator[](const size_t i) const {
    return edges[i];
  }
};

inline EdgeListR::EdgeListR(const EdgeListC& es) {
  n = es.n;
  edges.resize(es.edges.size());

  for (size_t i = 0; i < edges.size(); i++) {
    edges[i] = es.edges[i];
  }
}

inline EdgeListC::EdgeListC(const EdgeListR& es) {
  n = es.n;
  edges.resize(es.edges.size());

  for (size_t i = 0; i < edges.size(); i++) {
    edges[i] = es.edges[i];
  }
}
*/

template <typename ArcT>
struct AdjacencyArray {
  size_t n;
  std::vector<ArcT> arcs;
  std::vector<size_t> first_arc;

  AdjacencyArray() {
    n = 0;
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

  template <typename EdgeT>
  AdjacencyArray(const EdgeList<EdgeT>& es) {
    BuildGraph(es);
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

  /*
  AdjacencyMap(const TreePlusEdgesR& tree) {
    n = tree.n;
    neighbor_map.resize(n);
    for (size_t i = 0; i < tree.vertices.size(); i++) {
      size_t p = tree.vertices[i].parent;
      if (i == p) continue;
      AddEdgeR(i, p, tree.vertices[i].parent_resistance);
    }
    for (size_t i = 0; i < tree.off_tree_edges.size(); i++) {
      AddEdgeR(tree.off_tree_edges[i]);
    }
  }
  */

  AdjacencyMap(const TreeR& tree) {
    n = tree.n;
    neighbor_map.resize(n);
    for (size_t i = 0; i < n; i++) {
      size_t p = tree.vertices[i].parent;
      if (i == p) continue;
      const FLOAT& r = tree.vertices[i].parent_resistance;
      AddEdgeR(i, p, r);
    }
  }

};

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
    for(size_t j = 0; j < m; j++) {
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

#endif  // INCLUDE_GRAPH_H_
