/********************************************************************************
 * 
 * Graph Definitions
 * 
 * API Interface:
 *   struct Graph: definition of a graph (1-INDEXED).
 *     int .n: the size of graph
 *     vector< pair<int,FLOAT> > .neighbor[i]:
 *         the neighbors of vertex i, 1<=i<=n
 *     .freeMemory(): destroys graph and frees memory
 * 
 *   struct TreePlusEdges: definition of a connected graph with a spanning tree (1-INDEXED).
 *     int .n: the size of graph
 *     vector< pair<int,FLOAT> > .e[i]: the children of vertex i in spanning tree, 1<=i<=n
 *                                       the tree is DIRECTED, and 1 is always the root
 *     vector< tuple<int,int,FLOAT> > .o: the off-tree edges
 *     .freeMemory(): destroys graph and frees memory
 * 
 * NOTE:
 *   The resistance in Graph and GraphSP objects are always RESISTANCE, not WEIGHT!!!
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
 *   In GraphSP, however, every edge appears in only one direction.
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

  void FreeMemory() const {
    delete neighbor_list;
  }
};

struct TreePlusEdges {
  int n;
  vector<Arc> *children;
  vector<Edge> *off_tree_edge;

  TreePlusEdges() {
    n = 0;
    children = NULL;
    off_tree_edge = NULL;
  }

  TreePlusEdges(int _n) {
    n = _n;
    children = new vector<Arc>[n];
    off_tree_edge = new vector<Edge>;
  }

  TreePlusEdges(const TreePlusEdges &o) {
    n = o.n;
    children = o.children;
    off_tree_edge = o.off_tree_edge;
  }

  TreePlusEdges &operator =(const TreePlusEdges &o) {
    n = o.n;
    children = o.children;
    off_tree_edge = o.off_tree_edge;
    return (*this);
  }

  void freeMemory() const {
    delete children;
    delete off_tree_edge;
  }
};

#endif
