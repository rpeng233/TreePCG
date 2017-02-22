#include <cstddef>
#include <limits>
#include <random>
#include <vector>
#include "binary_heap.h"
#include "graph.h"
#include "pairing_heap.h"

using std::vector;

struct DijkstraVtx {
  PairingNode pairing_node;
  // size_t parent;
  size_t center;
  size_t parent_edge_id;
  double distance;

  void Initialize(size_t i) {
    // parent = center = i;
    center = i;
    distance = std::numeric_limits<double>::infinity();
  }
};

struct DistanceLess {
  bool operator () (const DijkstraVtx& u, const DijkstraVtx& v) const {
    return u.distance < v.distance;
  }
};

typedef
PairingHeap<DijkstraVtx, offsetof(DijkstraVtx, pairing_node), DistanceLess>
DijkstraQueue;

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

  Edge() : u(0), v(0), original_id(0) { }

  Edge(size_t u_, size_t v_, size_t id)
    : u(u_), v(v_), original_id(id)
  { }

  Edge& operator=(const EdgeR& e) {
    u = e.u;
    v = e.v;
  }
};

struct Arc {
  size_t v;
  size_t original_id;

  Arc& operator=(const Edge& e) {
    original_id = e.original_id;
    return *this;
  }
};

void Dijkstra(const AdjacencyArray<Arc>& graph,
              BinaryHeap<DistanceGreater>& queue,
              std::vector<DijkstraVtx>& vs) {
  while (queue.size != 0) {
    const size_t u = queue.Top();
    queue.Pop();
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const Arc& a = graph.arcs[i];
      double new_distance = vs[u].distance + 1;
      if (vs[a.v].distance <= new_distance) {
        continue;
      }
      vs[a.v].distance = new_distance;
      // vs[a.v].parent = u;
      vs[a.v].parent_edge_id = a.original_id;
      vs[a.v].center = vs[u].center;
      queue.BubbleUp(a.v);
    }
  }
}

void Dijkstra(const AdjacencyArray<Arc>& graph,
              vector<DijkstraVtx>& vs,
              DijkstraQueue& queue) {
  while (queue.Size() > 0) {
    const size_t u = &queue.Get_min() - &vs[0];
    // std::cout << u << '\n';
    queue.Delete_min();
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const Arc& a = graph.arcs[i];
      double new_distance = vs[u].distance + 1;
      // std::cout << a.v << ' ' << vs[a.v].distance << "   " << new_distance << '\n';
      if (vs[a.v].distance <= new_distance) {
        continue;
      }
      vs[a.v].distance = new_distance;
      // vs[a.v].parent = u;
      vs[a.v].parent_edge_id = a.original_id;
      vs[a.v].center = vs[u].center;
      queue.Decrease_key(vs[a.v]);
    }
  }
}

void AKPW_unweighted(const EdgeList<EdgeR>& es_, EdgeList<EdgeR>& tree) {
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
  AdjacencyArray<Arc> g;

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
      assert(vs[vs[i].center].center == vs[i].center);
      if (i == vs[i].center) {
        // vertex i is a center, assign it a new ID for next round
        new_id[i] = n;
        tmp_map[n] = original_id[i];
        n++;
      } else {
        tree.AddEdge(es_[vs[i].parent_edge_id]);
      }
    }

    m = 0;
    for (size_t i = 0; i < vs.size(); i++) {
      for (size_t j = g.first_arc[i]; j < g.first_arc[i + 1]; j++) {
        const Arc& a = g.arcs[j];
        if (a.v < i || vs[i].center == vs[a.v].center) {
          continue;
        }
        es[m].u = new_id[vs[i].center];
        es[m].v = new_id[vs[a.v].center];
        // es[m].resistance = a.resistance;
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

size_t inline Find(vector<size_t>& centers, size_t u) {
  if (centers[u] != u) {
    centers[u] = Find(centers, centers[u]);
  }
  return centers[u];
}

void AKPW(const EdgeList<EdgeR>& es, EdgeList<EdgeR>& tree) {
  size_t n = es.n;
  size_t m = es.Size();
  size_t count = 0;

  if (m == 0) return;

  tree.Clear();
  tree.n = n;
  tree.Reserve(n - 1);

  vector<DijkstraVtx> vs(n);
  vector<size_t> centers(n);
  vector<size_t> remaining(n);
  vector<size_t> tmp;
  tmp.reserve(n);
  EdgeList<Edge> es2;
  DistanceGreater cmp(vs);
  DijkstraQueue queue;
  AdjacencyArray<Arc> g;

  std::mt19937 rng(std::random_device{}());
  std::exponential_distribution<> exponential(0.5);

  es2.Reserve(m);
  es2.n = n;

  for (size_t i = 0; i < vs.size(); i++) {
    vs[i].Initialize(i);
    vs[i].distance = exponential(rng);
    queue.Insert(vs[i]);
    centers[i] = i;
    remaining[i] = i;
  }

  FLOAT bucket = es[0].resistance;
  FLOAT bucket_factor = 10;
  size_t idx = 0;
  m = 0;

  for (;;) {
    bucket *= bucket_factor;
    while (idx < es.Size() && es[idx].resistance < bucket) {
      size_t u = Find(centers, es[idx].u);
      size_t v = Find(centers, es[idx].v);
      es2.AddEdge(Edge(u, v, idx));
      idx++;
    }

    if (es2.Size() == 0) break;

    // for (size_t i = 0; i < es2.Size(); i++) {
    //   std::cout << es2[i].u << ' ' << es2[i].v << ' ' << es2[i].original_id << '\n';
    // }
    g.BuildGraph(es2);
    Dijkstra(g, vs, queue);
    queue.Clear();
    es2.Clear();
    es2.n = n;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      assert(vs[vs[u].center].center == vs[u].center);
      if (u == vs[u].center) {
        // std::cout << u << '\n';
        tmp.push_back(u);
        vs[u].Initialize(u);
        vs[u].distance = exponential(rng);
        queue.Insert(vs[u]);
      } else {
        // assert(u != vs[u].parent);
        tree.AddEdge(es[vs[u].parent_edge_id]);
        count++;
        // std::cout << "yeah\n";
        centers[u] = vs[u].center;
      }
    }

    if (tmp.size() <= 1) break;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      for (size_t j = g.first_arc[u]; j < g.first_arc[u + 1]; j++) {
        const Arc& a = g.arcs[j];
        if (a.v < u || vs[u].center == vs[a.v].center) {
          continue;
        }
        size_t v1 = Find(centers, u);
        size_t v2 = Find(centers, a.v);
        es2.AddEdge(Edge(v1, v2, a.original_id));
      }
    }


    remaining.swap(tmp);
    tmp.clear();
  }
  std::cout << tree.Size() << '\n';
}
