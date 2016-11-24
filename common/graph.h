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
  int u, v;
  FLOAT resistance;

  Edge() {
    u = 0;
    v = 0;
    resistance = 0;
  }

  Edge(int _u, int _v, FLOAT _resistance) {
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
};

struct Arc {
  int v;
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
};

bool CompareByResistance(const Arc &a1, const Arc &a2) {
  return a1.resistance < a2.resistance;
}


struct Graph {
  int n;
  vector<Arc> *neighbor_list;

  Graph() {
    n = 0;
    neighbor_list = NULL;
  }

  Graph(int _n) {
    n = _n;
    neighbor_list = new vector<Arc>[n];
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


  void FreeMemory() const {
    delete neighbor_list;
  }
};

struct TreePlusEdges {
  int n;
  int root;
  vector<Arc> *children;
  vector<Edge> off_tree_edge;

  TreePlusEdges() {
    n = 0;
    root = 0;
    children = NULL;
    off_tree_edge.clear();
  }

  TreePlusEdges(int _n) {
    n = _n;
    root = 0;
    children = new vector<Arc>[n];
    off_tree_edge.clear();
  }

  TreePlusEdges(const TreePlusEdges &o) {
    n = o.n;
    root = o.root;
    children = o.children;
    off_tree_edge = o.off_tree_edge;
  }

  TreePlusEdges &operator =(const TreePlusEdges &o) {
    n = o.n;
    root = o.root;
    children = o.children;
    off_tree_edge = o.off_tree_edge;
    return (*this);
  }

  void freeMemory() const {
    delete children;
  }
};

#endif
