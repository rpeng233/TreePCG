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

#endif  // INCLUDE_GRAPH_H_
