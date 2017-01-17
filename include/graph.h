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

  EdgeR(size_t _u, size_t _v, FLOAT _resistance) {
    u = _u;
    v = _v;
    resistance = _resistance;
  }

  EdgeR(const EdgeR &o) {
    u = o.u;
    v = o.v;
    resistance = o.resistance;
  }

  EdgeR &operator = (const EdgeR &o) {
    u = o.u;
    v = o.v;
    resistance = o.resistance;
    return (*this);
  }

  EdgeR& operator=(const EdgeC& e);

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

  EdgeC(size_t _u, size_t _v, FLOAT _conductance) {
    u = _u;
    v = _v;
    conductance = _conductance;
  }

  EdgeC(const EdgeC &o) {
    u = o.u;
    v = o.v;
    conductance = o.conductance;
  }

  EdgeC &operator = (const EdgeC &o) {
    u = o.u;
    v = o.v;
    conductance = o.conductance;
    return (*this);
  }

  EdgeC& operator=(const EdgeR& e);

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

  ArcR(int _v, FLOAT _resistance) {
    v = _v;
    resistance = _resistance;
  }

  ArcR(const ArcR &o) {
    v = o.v;
    resistance = o.resistance;
  }

  ArcR &operator = (const ArcR &o) {
    v = o.v;
    resistance = o.resistance;
    return (*this);
  }

  bool operator <(const ArcR &o) const {
    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->resistance < o.resistance;
  }
};

inline bool CompareByResistance(const ArcR &a1, const ArcR &a2) {
  return a1.resistance < a2.resistance;
}

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

struct Graph {
  int n;
  std::vector<ArcR> *neighbor_list;

  Graph() {
    n = 0;
    neighbor_list = NULL;
  }

  Graph(int _n) {
    n = _n;
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

struct TreeVertex {
  size_t parent;
  size_t children_count;
  FLOAT parent_resistance;
  bool eliminated;

  TreeVertex()
  : parent(0), children_count(0), parent_resistance(0.0), eliminated(false)
  { }
};

struct Tree {
  size_t n;
  std::vector<TreeVertex> vertices;

  Tree() {
    n = 0;
  }

  Tree(size_t _n) {
    n = _n;
    vertices.resize(n);
    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  void SetParent(size_t v, size_t p, FLOAT r) {
    size_t old_p = vertices[v].parent;
    if (old_p != v) {
      vertices[old_p].children_count--;
    }

    vertices[v].parent = p;
    vertices[v].parent_resistance = r;
    vertices[p].children_count++;
  }
};

struct EdgeListC;

struct EdgeListR {
  size_t n;
  std::vector<EdgeR> edges;

  EdgeListR() {
    n = 0;
  }

  EdgeListR(const EdgeListC& es);

  void AddEdge(size_t u, size_t v, FLOAT resistance) {
    edges.push_back(EdgeR(u, v, resistance));
  }

  void Clear() {
    n = 0;
    edges.clear();
  }
};

struct EdgeListC {
  size_t n;
  std::vector<EdgeC> edges;

  EdgeListC() {
    n = 0;
  }

  EdgeListC(const EdgeListR& es);

  void AddEdge(size_t u, size_t v, FLOAT conductance) {
    edges.push_back(EdgeC(u, v, conductance));
  }

  void Clear() {
    n = 0;
    edges.clear();
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

struct Graph2 {
  size_t n;
  std::vector<ArcR> arcs;
  std::vector<size_t> first_arc;

  Graph2(const EdgeListR& es) {
    n = es.n;
    size_t m = es.edges.size();

    arcs.resize(m * 2);
    first_arc.resize(n + 1);
    std::vector<size_t> degrees(n);

    for (std::vector<EdgeR>::const_iterator it = es.edges.begin();
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
    for (std::vector<EdgeR>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      degrees[it->u]--;
      tmp_index = first_arc[it->u] + degrees[it->u];
      arcs[tmp_index].v = it->v;
      arcs[tmp_index].resistance = it->resistance;

      degrees[it->v]--;
      tmp_index = first_arc[it->v] + degrees[it->v];
      arcs[tmp_index].v = it->u;
      arcs[tmp_index].resistance = it->resistance;
    }
  }
};

struct Graph3 {
  size_t n;
  std::vector<std::map<size_t, FLOAT> > neighbor_map;

  Graph3(const EdgeListR& es) {
    n = es.n;
    neighbor_map.resize(n);
    for (std::vector<EdgeR>::const_iterator it = es.edges.begin();
         it != es.edges.end();
         ++it) {
      neighbor_map[it->u][it->v] = it->resistance;
      neighbor_map[it->v][it->u] = it->resistance;
    }
  }
};

class One {
public:
  FLOAT operator() () const {
    return FLOAT(1);
  }
};

template <class EdgeListType, class WeightGen>
void line(size_t n, EdgeListType& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n - 1; i++) {
    es.AddEdge(i, i + 1, wgen());
  }
}

template <class EdgeListType>
void line(size_t n, EdgeListType& es) {
  line(n, es, One());
}

template <class EdgeListType, class WeightGen>
void cycle(size_t n, EdgeListType& es, WeightGen wgen=WeightGen()) {
  line(n, es, wgen);
  es.AddEdge(0, n - 1, wgen());
}

template <class EdgeListType>
void cycle(size_t n, EdgeListType& es) {
  cycle(n, es, One());
}

template <class EdgeListType, class WeightGen>
void torus(size_t n, size_t m, EdgeListType& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n * m;
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < m; j++) {
      es.AddEdge(i * m + j, i * m + (j + 1) % m, wgen());
      es.AddEdge(i * m + j, ((i + 1) % n) * m + j, wgen());
    }
  }
}

template <class EdgeListType>
void torus(size_t n, size_t m, EdgeListType& es) {
  torus(n, m, es, One());
}

template <class EdgeListType, class WeightGen>
void gnp(size_t n, double p, EdgeListType& es, WeightGen wgen=WeightGen()) {
  es.Clear();
  es.n = n;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if ((double) rand() / RAND_MAX < p) {
        es.AddEdge(i, j, wgen());
      }
    }
  }
}

template <class EdgeListType>
void gnp(size_t n, double p, EdgeListType& es) {
  gnp(n, p, es, One());
}

inline
void cayley(size_t n,
            const std::vector<size_t>& skips,
            const std::vector<FLOAT>& resistances,
            EdgeListR& es) {
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

void gnp(size_t n, double p, EdgeListR& es);

/*
struct TreePlusEdges {
  size_t n;
  // size_t root;
  std::vector<TreeVertex> vertices;
  std::vector<EdgeR> off_tree_edges;

  TreePlusEdges() {
    n = 0;
    // root = 0;
  }

  TreePlusEdges(size_t _n) {
    n = _n;
    // root = 0;
    vertices.resize(n);
    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  TreePlusEdges(const TreePlusEdges &o) {
    n = o.n;
    // root = o.root;
    vertices = o.vertices;
    off_tree_edges = o.off_tree_edges;
  }

  TreePlusEdges &operator =(const TreePlusEdges &o) {
    n = o.n;
    // root = o.root;
    vertices = o.vertices;
    off_tree_edges = o.off_tree_edges;
    return (*this);
  }

  void setParent(size_t v, size_t p, FLOAT r) {
    vertices[v].parent = p;
    vertices[v].parent_resistance = r;
    vertices[p].ref_count++;
  }

  void addEdge(size_t u, size_t v, FLOAT r) {
    off_tree_edges.push_back(EdgeR(u, v, r));
  }
};
*/

#endif  // INCLUDE_GRAPH_H_
