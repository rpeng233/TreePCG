#include <limits>
#include <random>
#include <vector>
#include "binary_heap.h"
#include "graph.h"

using std::vector;

struct DijkstraVtx {
  size_t parent;
  size_t ancestor;
  size_t parent_edge_id;
  double distance;

  void Initialize(size_t i) {
    parent = ancestor = i;
    distance = std::numeric_limits<double>::infinity();
  }
};

struct DistanceGreater {
  const vector<DijkstraVtx>& vs;

  DistanceGreater(const vector<DijkstraVtx>& vs_) : vs(vs_) { }

  bool operator() (size_t u, size_t v) const {
    return vs[u].distance > vs[v].distance;
  }
};

struct Edge {
  size_t u;
  size_t v;
  size_t original_id;
  FLOAT resistance;

  Edge& operator=(const EdgeR& e) {
    u = e.u;
    v = e.v;
    resistance = e.resistance;
  }
};

struct Arc {
  size_t v;
  size_t original_id;
  FLOAT resistance;

  Arc& operator=(const Edge& e) {
    original_id = e.original_id;
    resistance = e.resistance;

    return *this;
  }
};

void Dijkstra(const Graph2<Arc>& graph,
              BinaryHeap<DistanceGreater>& queue,
              std::vector<DijkstraVtx>& vs) {
  while (queue.size != 0) {
    size_t u = queue.Top();
    queue.Pop();
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const Arc& a = graph.arcs[i];
      double new_distance = vs[u].distance + (double) a.resistance;
      if (vs[a.v].distance <= new_distance) {
        continue;
      }
      vs[a.v].distance = new_distance;
      vs[a.v].parent = u;
      vs[a.v].parent_edge_id = a.original_id;
      vs[a.v].ancestor = vs[u].ancestor;
      queue.BubbleUp(a.v);
    }
  }
}

void AKPW(const EdgeList<EdgeR>& es_, EdgeList<EdgeR>& tree) {
  EdgeList<Edge> es(es_);
  tree.n = es_.n;

  for (size_t i = 0; i < es.Size(); i++) {
    es[i].original_id = i;
  }

  size_t n = es.n;
  size_t m = es.Size();
  vector<DijkstraVtx> vs(n);
  DistanceGreater cmp(vs);
  BinaryHeap<DistanceGreater> queue(0, &cmp);
  Graph2<Arc> g;

  std::mt19937 rng(std::random_device{}());
  std::exponential_distribution<> exponential(0.5);
  // maps vertex ID in quotient graph to vertex ID in original graph
  vector<size_t> original_id(n);
  vector<size_t> tmp_map(n);
  // maps cluster ID -> {1, 2, 3, ...}
  vector<size_t> new_id(n);

  for (size_t i = 0; i < n; i++) {
    original_id[i] = i;
  }

  size_t count = 0;
  for (;;) {
    for (size_t i = 0; i < vs.size(); i++) {
      vs[i].Initialize(i);
      vs[i].distance = exponential(rng);
    }
    queue.Reset(n, &cmp);
    g.BuildGraph(es);
    Dijkstra(g, queue, vs);

    n = 0;
    for (size_t i = 0; i < vs.size(); i++) {
      assert(vs[vs[i].ancestor].ancestor == vs[i].ancestor);
      if (i == vs[i].ancestor) {
        // vertex i is a center, assign it a new ID for next round
        new_id[i] = n;
        tmp_map[n] = original_id[i];
        n++;
      }
      if (i != vs[i].parent) {
        tree.AddEdge(es_[vs[i].parent_edge_id]);
      }
    }

    m = 0;
    for (size_t i = 0; i < vs.size(); i++) {
      for (size_t j = g.first_arc[i]; j < g.first_arc[i + 1]; j++) {
        const Arc& a = g.arcs[j];
        if (a.v < i || vs[i].ancestor == vs[a.v].ancestor) {
          continue;
        }
        es[m].u = new_id[vs[i].ancestor];
        es[m].v = new_id[vs[a.v].ancestor];
        es[m].resistance = a.resistance;
        es[m].original_id = a.original_id;
        m++;
      }
    }
    // std::cout << count++ << ' ' << n << ' ' << m << ' ' << tree.edges.size() << '\n';
    if (n == 1 || m == 0) break;
    vs.resize(n);
    es.n = n;
    es.edges.resize(m);
    original_id.swap(tmp_map);
    // for (size_t i = 0; i < es.size(); i++) {
    //   es[i].u = new_id[es[i].u];
    //   es[i].v = new_id[es[i].v];
    // }
  }
}
