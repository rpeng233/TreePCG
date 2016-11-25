/********************************************************************************
 * 
 * Low Stretch Spanning Tree utility
 * 
 * Public API:
 *   vector<FLOAT> StretchCalculator::CalculatePathResistance
 *        (const TreePlusEdges &graph):
 *     Input: graph in tree plus off-tree edges format.
 *     Returns: 
 *       The tree-path resistance corresponding to off-tree edges of g, in the same order as g.o.
 *      Stretch of an off-tree edge is by definition (tree-path resistance)/(off-tree edge resistance).
 *   
 *   FLOAT CalculateTotalStretch(const GraphSP &g):
 *     Input: g, the graph.
 *     Returns: FLOAT with the total stretch of g.
 * 
 *   vector<int> DijkstraTree(const Graph &g):
 *     Input: g, the graph.
 *     Returns:  vector<int> containing all parents
 *
 *   TreePlusEdges GraphToTreePlusEdges(const Graph &g, vector<int> parents)
 *     Input: g, the graph
 *            parents, all parents of vertices in the graph
 *     Returns: TreePlusEdges with the tree marked as special 
 * 
 ********************************************************************************/

#ifndef __TREEFINDER_H__
#define __TREEFINDER_H__

#include <algorithm>
#include <cstdlib>
#include <stack>
#include <vector>
#include <queue>
#include <utility>
#include "../common/common.h"
#include "../common/graph.h"

// TOFIX: move some of these utils into common
//
//

struct UnionFind{
// union find, without rank, just path compression
// this is not on the runtime bottleneck
  int n;
  int *parent;

  UnionFind(int n) {
    parent = new int[n];
    for (int i = 0; i < n; ++i) {
      parent[i] = i;
    }
  }

  int FindSet(int x) {
  // written iteratively since we don't do union-by-rank
  // therefore might blow up stack
    vector<int> list;
    while(parent[parent[x]] != parent[x]) {
      list.push_back(x);
      x = parent[x];
    }
    for(vector<int>::iterator it = list.begin(); it != list.end(); ++it) {
      parent[*it] = parent[x];
    }
    return parent[x];
  }

  void Union(int u, int v) {
  // important that u's parent is set to v
  // offline LCA crucially relies on this
    u = FindSet(u);
    v = FindSet(v);
    parent[u] = parent[v];
  }

  ~UnionFind() {
    delete parent;
  } 
};

void DFS(int n, int root, vector<Arc> *arc, int *event, int *parent) {
//
// TOFIX: this assumes that things are connected
//
// This works on both the arc list of trees
// and the adjacency list of a graph, aka. it doesn't
// assume that arc is symmetric
//
// returns:
//   list of events, with 2 * u = add u onto stack
//                        2 * u + 1 = remove u from stack
//
//   parents of all nodes in the DFS trees
//
  char *visited = new char[n];
  int *i1 = new int[n];
  for (int i = 0; i < n; ++i) {
    visited[i] = 0;
    i1[i] = 0;
  }

  std::stack<int> dfs_stack;
  dfs_stack.push(root);
  parent[root] = -1;
  visited[root] = 1;
  int t = 1;
  event[0] = root * 2;

  while (!dfs_stack.empty()) {
    int u = dfs_stack.top();
    if (i1[u] == int(arc[u].size())) {
      event[t++] = u * 2 + 1;
      dfs_stack.pop();
    } else {
      int v = arc[u][i1[u]].v;
      i1[u]++;
      if (!visited[v]) {
        visited[v] = 1;
        parent[v] = u;
        event[t++] = v * 2;
        dfs_stack.push(v);
      }
    }
  }

  delete visited;
  delete i1;
}

TreePlusEdges GraphToTreePlusEdges(const Graph &graph,
                                   const vector<int> parent) {
  TreePlusEdges h(graph.n);
  int parents_found = 0;
  for (int u = 0; u < graph.n; ++u) {
    if (parent[u] == -1) {
      h.root = u;
    }
    for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
         it != graph.neighbor_list[u].end(); ++it) {
      if (u == parent[it -> v]) {
        parents_found++;
// printf("%d -> %d %lf\n", u, it -> v, it -> resistance);
        h.children[u].push_back(Arc(it -> v, it -> resistance));
      } else if (parent[u] != it -> v && u < it -> v) {
// printf("%d == %d %lf\n", u, it -> v, it -> resistance);
        h.off_tree_edge.push_back(Edge(u, it -> v, it -> resistance));
      }
    }
  }
  assert(parents_found == graph.n - 1);
  return h;
}


namespace TreeFinder {
  vector<double> GetStretch(const TreePlusEdges graph) {
    int *event = new int[graph.n * 2];
    int *parent = new int[graph.n];
    DFS(graph.n, graph.root, graph.children, event, parent);

    double *distance_to_root = new double[graph.n];
    int *discover_time = new int[graph.n];
    distance_to_root[graph.root] = 0;
    for (int i = 0; i < graph.n * 2; ++i) {
      if (event[i] % 2 == 0) {
        int u = event[i] / 2;
        discover_time[u] = i;
        for (vector<Arc>::iterator it = graph.children[u].begin();
            it != graph.children[u].end(); ++it) {
          if (parent[it -> v] == u) {
            distance_to_root[it -> v] =
              distance_to_root[u] + it -> resistance;
          }
        }
      }
    }
/*
*    for(int u = 0; u < graph.n; ++u) {
*      fprintf(stderr, "%d-->>%d %lf\n", u, parent[u], distance_to_root[u]);
*    }
*/
    vector<pair<int, int> > *query = new vector<pair<int, int> > [graph.n];
    for (int i = 0; i < int(graph.off_tree_edge.size()); ++i) {
      Edge edge = graph.off_tree_edge[i];
      int u = edge.u;
      int v = edge.v;
      if (discover_time[u] > discover_time[v]) {
        swap(u, v);
      }
      query[v].push_back(make_pair(u, i));
    }

    UnionFind uf(graph.n);
    vector<double> result(graph.off_tree_edge.size());
    for (int i = 0; i < graph.n * 2; ++i) {
// printf("(%d)%d\n", event[i] % 2, event[i] / 2); fflush(stdout);
      int u = event[i] / 2;
      if (event[i] % 2 == 0) {
        for (vector<pair<int, int> >::iterator it = query[u].begin();
            it != query[u].end(); ++it) {
          int v = it -> first;
// printf("%d %d %d\n", u, v, uf.FindSet(v));
          result[it -> second] = (distance_to_root[u] + distance_to_root[v]
              - distance_to_root[uf.FindSet(v)] * double(2))
                / graph.off_tree_edge[it -> second].resistance;
        }
      } else if (parent[u] != -1) {
        uf.Union(u, parent[u]);
      }
    }

    for(int i = 0; i < int(graph.off_tree_edge.size()); ++i) {
      int u = graph.off_tree_edge[i].u;
      int v = graph.off_tree_edge[i].v;
      double tree_resistance = distance_to_root[u] + distance_to_root[v];
      while (u != v) {
        if(distance_to_root[u] < distance_to_root[v]) {
          swap(u, v);
        }
        u = parent[u];
      }
      tree_resistance -= distance_to_root[u] * double(2);
      double tree_stretch = tree_resistance
                          / graph.off_tree_edge[i].resistance;
      assert(fabs(result[i] - tree_stretch) < 1e-9);
//      fprintf(stderr, "%lf %lf\n", result[i], tree_resistance
    }

    delete event;
    delete discover_time;
    delete distance_to_root;
    for (int i = 0; i < graph.n; ++i) {
      query[i].clear();
    }
    delete[] query;

    return result;
  }

  double TotalStretch(const TreePlusEdges graph) {
    vector<double> stretch = GetStretch(graph);

    double total = 0;
    for (int i = 0; i < int(graph.off_tree_edge.size()); ++i) {
// fprintf(stderr, "%lf\n", stretch[i]);
      total += stretch[i];
    }
    return total;
  }

  vector<int> DijkstraTree(const Graph &graph) {
// Dijkstra tree from a random vertex
//
//
    int root = rand() % graph.n;

    fprintf(stderr, "Running Dijkstra on %d Vertices\n", graph.n);
    fflush(stdout);

    double *dist = new double[graph.n];
    char *checked = new char[graph.n];
    vector<int> parent(graph.n);
    for (int i = 0; i < graph.n; ++i) {
      dist[i] = 1e100;  // TOFIX: should be global CONST
      parent[i] = -1;
      checked[i] = 0;
    }

// TOFIX: move Dijkstra into its own method (probably
// necessary for exponential start time clustering
    priority_queue<pair<double, int> > Q;
    dist[root] = 0;
    Q.push(make_pair(0, root));

    while (!Q.empty()) {
      int u = Q.top().second;
// fprintf(stderr, "%d %lf\n", u, dist[u]);
      Q.pop();
      if (!checked[u]) {
        checked[u] = 1;
        for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
             it != graph.neighbor_list[u].end(); ++it) {
          double temp = dist[u] + it -> resistance;
          if (temp < dist[it -> v]) {
            dist[it -> v] = temp;
            parent[it -> v] = u;
            Q.push(make_pair(-temp, it -> v));
          }
        }
      }
    }

    delete checked;
    delete dist;
    return parent;
  }

  bool CompareByResistance(const Edge &e1, const Edge &e2) {
    return e1.resistance < e2.resistance;
  }

  vector<int> MinimumSpanningTree(const Graph &graph) {
    vector<Edge> edge_list;
    for (int u = 0; u < graph.n; ++u) {
      for (vector<Arc>::iterator it = graph.neighbor_list[u].begin();
           it != graph.neighbor_list[u].end(); ++it) {
        if (u < it -> v) {
          edge_list.push_back(Edge(u, it -> v, it -> resistance));
        }
      }
    }
    sort(edge_list.begin(), edge_list.end(), CompareByResistance);

    Graph h(graph.n);
    UnionFind uf(graph.n);

    int tree_edges = 0;
    for (vector<Edge>::iterator it = edge_list.begin();
         it != edge_list.end(); ++it) {
// fprintf(stderr, "%d==%d %lf\n", it -> u, it -> v, it -> resistance);
      if (uf.FindSet(it -> u) != uf.FindSet(it -> v)) {
        tree_edges++;
        uf.Union(it -> u, it -> v);
        h.AddEdge(*it);
      }
    }
    assert(tree_edges == graph.n - 1);

    vector<int> parent(graph.n);
    int *event = new int[graph.n * 2];
    DFS(graph.n, 0, h.neighbor_list, event, &(parent[0]));
    return parent;
  }
}  // namespace TreeFinder

/*
namespace TreeFinder
{
  GraphSP findLowStretchTree(const Graph &g)
  {
    static FLOAT *vis=new FLOAT[MAXN];
    rep(i,1,g.n) vis[i]=0;
    Graph g2(g.n);
    rep(i,1,g.n)
    {
      rept(it,g.e[i]) vis[it->first]+=1.0/it->second;
      rept(it,g.e[i])
        if (vis[it->first]!=0)
        {
          g2.e[i].push_back(make_pair(it->first,1.0/vis[it->first]));
          vis[it->first]=0;
        }
    }
    GraphSP h1=DijkstraTreeFinder::solve(g2);
    GraphSP h2=MSTTreeFinder::solve(g2);
    g2.freeMemory();
    FLOAT st1=StretchCalculator::calculateTotalStretch(h1);
    FLOAT st2=StretchCalculator::calculateTotalStretch(h2);
    if (st1<st2)
    {
      h2.freeMemory();
      //printf("findLowStretchTree (Dijkstra): average stretch = %.16lf\n",(double)(st1/h1.o.size()));
      return h1;
    }
    else
    {
      h1.freeMemory();
      //printf("findLowStretchTree (MST): average stretch = %.16lf\n",(double)(st2/h2.o.size()));
      return h2;
    }
  }
  
  GraphSP findLowStretchTree(const GraphSP &g)
  {
    Graph g2=Graph(g.n);
    rep(i,1,g.n) 
      rept(it,g.e[i])
      {
        g2.e[i].push_back(*it);
        g2.e[it->first].push_back(make_pair(i,it->second));
      }
    rept(it,g.o)
    {
      int x=get<0>(*it), y=get<1>(*it); FLOAT z=get<2>(*it);
      g2.e[x].push_back(make_pair(y,z));
      g2.e[y].push_back(make_pair(x,z));
    }
    GraphSP h=findLowStretchTree(g2);
    g2.freeMemory();
    return h;
  }
}
*/

#endif
