/********************************************************************************
 *
 * Graph Definitions
 *
 * API Interface:
 *   struct Graph: definition of a graph (0-INDEXED).
 *     int .n: the size of graph
 *     vector<Arc> .neighbor[i]:
 *         the neighbors of vertex i, 0 <= i <n
 *     .freeMemory(): destroys graph and frees memory
 *
 *   struct TreePlusEdges: definition of a connected graph with
 *        a spanning tree (0-INDEXED), rooted at .root.
 *     int .n: the size of graph
 *     vector<Arc> .children[i]:
 *        the children of vertex i in spanning tree, 0 <= i < n
 *        the tree is DIRECTED, and 0 is always the root
 *
 *     vector<Edge> .off_tree_edges: the list off-tree edges
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

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <vector>
#include "common.h"

struct Edge {
  size_t u, v;
  FLOAT resistance;

  Edge() {
    u = 0;
    v = 0;
    resistance = 0;
  }

  Edge(size_t _u, size_t _v, FLOAT _resistance) {
    u = _u;
    v = _v;
    resistance = _resistance;
  }

  Edge(const Edge &o) {
    u = o.u;
    v = o.v;
    resistance = o.resistance;
  }

  Edge &operator = (const Edge &o) {
    u = o.u;
    v = o.v;
    resistance = o.resistance;
    return (*this);
  }

  bool operator <(const Edge &o) const {
    if (this->u != o.u) {
      return this->u < o.u;
    }

    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->resistance < o.resistance;
  }

};

struct Arc {
  size_t v;
  FLOAT resistance;

  Arc() {
    v = 0;
    resistance = 0;
  }

  Arc(int _v, FLOAT _resistance) {
    v = _v;
    resistance = _resistance;
  }

  Arc(const Arc &o) {
    v = o.v;
    resistance = o.resistance;
  }

  Arc &operator = (const Arc &o) {
    v = o.v;
    resistance = o.resistance;
    return (*this);
  }

  bool operator <(const Arc &o) const {
    if (this->v != o.v) {
      return this->v < o.v;
    }

    return this->resistance < o.resistance;
  }
};

inline bool CompareByResistance(const Arc &a1, const Arc &a2) {
  return a1.resistance < a2.resistance;
}

struct Vertex {
  std::vector<Arc> nghbrs;

  void addArc(size_t v, FLOAT resistance) {
    nghbrs.emplace_back(v, resistance);
  }

  void sortAndCombine() {
    if (nghbrs.empty()) return;
    sort(nghbrs.begin(), nghbrs.end());
    size_t new_degree = 0;
    size_t last_nghbr = nghbrs[0].v;
    FLOAT conductance = 0;

    for (const auto& a : nghbrs) {
      if (a.v == last_nghbr) {
        conductance += 1 / a.resistance;
      } else {
        nghbrs[new_degree].v = last_nghbr;
        nghbrs[new_degree].resistance = 1 / conductance;
        new_degree++;
        last_nghbr = a.v;
        conductance = 1 / a.resistance;
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
  std::vector<Arc> *neighbor_list;

  Graph() {
    n = 0;
    neighbor_list = NULL;
  }

  Graph(int _n) {
    n = _n;
    neighbor_list = new std::vector<Arc>[n];
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
    neighbor_list[u].push_back(Arc(v, resistance));
    neighbor_list[v].push_back(Arc(u, resistance));
  }

  void AddEdge(Edge e) {
    AddEdge(e.u, e.v, e.resistance);
  }

  void sortAndCombine() const {
    for(int u = 0; u < n; ++u) {
      int new_degree = 0;
      sort(neighbor_list[u].begin(), neighbor_list[u].end());

      int last_vertex = -1;
      FLOAT last_weight = 0;

      for (std::vector<Arc>::iterator it = neighbor_list[u].begin();
          it != neighbor_list[u].end(); ++it) {
        if (it -> v != last_vertex) {
          if (last_vertex >= 0) {
            neighbor_list[u][new_degree] =
              Arc(last_vertex, FLOAT(1) / last_weight);
            new_degree++;
          }
          last_vertex = it -> v;
        }
        last_weight += FLOAT(1) / it -> resistance;
      }

      if (last_vertex >= 0) {
        neighbor_list[u][new_degree] =
          Arc(last_vertex, FLOAT(1) / last_weight);
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

  Tree(size_t _n)
  {
    n = _n;
    vertices.resize(n);
    for (size_t i = 0; i < n; i++) {
      vertices[i].parent = i;
    }
  }

  void setParent(size_t v, size_t p, FLOAT r) {
    vertices[v].parent = p;
    vertices[v].parent_resistance = r;
    vertices[p].children_count++;
  }
};

/*
struct TreePlusEdges {
  size_t n;
  // size_t root;
  std::vector<TreeVertex> vertices;
  std::vector<Edge> off_tree_edges;

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
    off_tree_edges.emplace_back(u, v, r);
  }
};
*/

#endif
