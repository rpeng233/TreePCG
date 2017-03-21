#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>
#include "binary_heap.h"
#include "graph.h"
#include "pairing_heap.h"

using std::vector;
using std::pair;

struct DijkstraVtx {
  PairingNode pairing_node;
  size_t center;
  size_t parent_edge_id;
  double distance;

  void Initialize(size_t i) {
    center = i;
    distance = std::numeric_limits<double>::infinity();
  }
};

struct BFSVtx {
  size_t center;
  size_t parent_edge_id;
  double distance;
  bool done;

  void Initialize(size_t i) {
    center = i;
    distance = std::numeric_limits<double>::infinity();
    done = false;
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

struct AKPWEdge {
  size_t u;
  size_t v;
  size_t original_id;

  AKPWEdge() : u(0), v(0), original_id(0) { }

  AKPWEdge(size_t u_, size_t v_, size_t id)
    : u(u_), v(v_), original_id(id)
  { }

  AKPWEdge& operator=(const EdgeR& e) {
    u = e.u;
    v = e.v;
    return *this;
  }
};

struct AKPWArc {
  size_t v;
  size_t original_id;

  AKPWArc& operator=(const AKPWEdge& e) {
    original_id = e.original_id;
    return *this;
  }
};

void Dijkstra(const AdjacencyArray<AKPWArc>& graph,
              vector<DijkstraVtx>& vs,
              DijkstraQueue& queue) {
  while (queue.Size() > 0) {
    const size_t u = &queue.Get_min() - &vs[0];
    // std::cout << u << '\n';
    queue.Delete_min();
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const AKPWArc& a = graph.arcs[i];
      double new_distance = vs[u].distance + 1;
      // std::cout << a.v << ' ' << vs[a.v].distance << "   " << new_distance << '\n';
      if (vs[a.v].distance <= new_distance) {
        continue;
      }
      vs[a.v].distance = new_distance;
      vs[a.v].parent_edge_id = a.original_id;
      vs[a.v].center = vs[u].center;
      queue.Decrease_key(vs[a.v]);
    }
  }
}

void BFS(const AdjacencyArray<AKPWArc>& graph,
         vector<BFSVtx>& vs,
         vector<size_t>& q1,
         vector<size_t>& q2) {
  size_t i = 0;
  size_t j = 0;

  for (;;) {
    size_t u;
    while (i < q1.size() && vs[q1[i]].done) {
      i++;
    }
    if (i < q1.size() && j < q2.size()) {
      if (vs[q1[i]].distance < vs[q2[j]].distance) {
        u = q1[i++];
      } else {
        u = q2[j++];
      }
    } else if (i < q1.size()) {
      u = q1[i++];
    } else if (j < q2.size()) {
      u = q2[j++];
    } else {
      break;
    }
    for (size_t i = graph.first_arc[u]; i < graph.first_arc[u + 1]; i++) {
      const AKPWArc& a = graph.arcs[i];
      double new_distance = vs[u].distance + 1;
      if (vs[a.v].distance <= new_distance) {
        continue;
      }
      vs[a.v].distance = new_distance;
      vs[a.v].parent_edge_id = a.original_id;
      vs[a.v].center = vs[u].center;
      vs[a.v].done = true;
      q2.push_back(a.v);
    }
  }
}

void AddExponentialShift(vector<BFSVtx>& vs, vector<size_t>& q, double lambda) {
  std::random_shuffle(q.begin(), q.end());

  std::mt19937 rng(std::random_device{}());
  std::uniform_real_distribution<double> unif(0, 1);
  size_t n = q.size();
  double sum = 0;
  for (size_t i = n; i > 0; i--) {
    sum += std::log(1 - unif(rng)) / (lambda * i);
    vs[q[i - 1]].distance = sum;
  }
}

static inline size_t Find(vector<size_t>& centers, size_t u) {
  if (centers[u] != u) {
    centers[u] = Find(centers, centers[u]);
  }
  return centers[u];
}

void AKPW(const EdgeList<EdgeR>& es, EdgeList<EdgeR>& tree) {
  const size_t n = es.n;
  const size_t m = es.Size();
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
  EdgeList<AKPWEdge> es2;
  DijkstraQueue queue;
  AdjacencyArray<AKPWArc> g;

  std::mt19937 rng(std::random_device{}());
  std::exponential_distribution<> exponential(0.01);

  es2.Reserve(m);
  es2.n = n;

  for (size_t i = 0; i < vs.size(); i++) {
    vs[i].Initialize(i);
    vs[i].distance = -exponential(rng);
    queue.Insert(vs[i]);
    centers[i] = i;
    remaining[i] = i;
  }

  FLOAT bucket = es[0].resistance;
  FLOAT bucket_factor = 2;
  size_t idx = 0;

  for (;;) {
    bucket *= bucket_factor;
    while (idx < es.Size() && es[idx].resistance < bucket) {
      size_t u = Find(centers, es[idx].u);
      size_t v = Find(centers, es[idx].v);
      es2.AddEdge(AKPWEdge(u, v, idx));
      idx++;
    }

    if (es2.Size() == 0 && idx == es.Size()) break;

    g.BuildGraph(es2);
    Dijkstra(g, vs, queue);
    queue.Clear();
    es2.Clear();
    es2.n = n;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      assert(vs[vs[u].center].center == vs[u].center);
      if (u == vs[u].center) {
        tmp.push_back(u);
        vs[u].Initialize(u);
        vs[u].distance = -exponential(rng);
        queue.Insert(vs[u]);
      } else {
        tree.AddEdge(es[vs[u].parent_edge_id]);
        count++;
        centers[u] = vs[u].center;
      }
    }

    if (tmp.size() <= 1 && idx == es.Size()) break;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      for (size_t j = g.first_arc[u]; j < g.first_arc[u + 1]; j++) {
        const AKPWArc& a = g.arcs[j];
        if (a.v < u || vs[u].center == vs[a.v].center) {
          continue;
        }
        size_t v1 = Find(centers, u);
        size_t v2 = Find(centers, a.v);
        es2.AddEdge(AKPWEdge(v1, v2, a.original_id));
      }
    }


    remaining.swap(tmp);
    tmp.clear();
    // std::cout << tree.Size() << '\n';
  }

  // std::cout << tree.Size() << '\n';
}

void AKPW2(const EdgeList<EdgeR>& es, EdgeList<EdgeR>& tree) {
  const size_t n = es.n;
  const size_t m = es.Size();
  size_t count = 0;

  if (m == 0) return;

  tree.Clear();
  tree.n = n;
  tree.Reserve(n - 1);

  vector<BFSVtx> vs(n);
  vector<size_t> centers(n);
  vector<size_t> remaining(n);
  vector<size_t> tmp;
  tmp.reserve(n);
  EdgeList<AKPWEdge> es2;
  vector<size_t> q1, q2;
  AdjacencyArray<AKPWArc> g;

  std::mt19937 rng(std::random_device{}());

  es2.Reserve(m);
  es2.n = n;

  q1.reserve(n);
  q2.reserve(n);
  for (size_t i = 0; i < vs.size(); i++) {
    vs[i].Initialize(i);
    q1.emplace_back(i);
    centers[i] = i;
    remaining[i] = i;
  }

  FLOAT bucket = es[0].resistance;
  FLOAT bucket_factor = 2;
  size_t idx = 0;

  for (;;) {
    bucket *= bucket_factor;
    while (idx < es.Size() && es[idx].resistance < bucket) {
      size_t u = Find(centers, es[idx].u);
      size_t v = Find(centers, es[idx].v);
      es2.AddEdge(AKPWEdge(u, v, idx));
      idx++;
    }

    if (es2.Size() == 0 && idx == es.Size()) break;

    g.BuildGraph(es2);
    AddExponentialShift(vs, q1, 0.05);
    BFS(g, vs, q1, q2);

    q1.clear();
    q2.clear();
    es2.Clear();
    es2.n = n;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      assert(vs[vs[u].center].center == vs[u].center);
      if (u == vs[u].center) {
        tmp.push_back(u);
        vs[u].Initialize(u);
        q1.emplace_back(u);
      } else {
        tree.AddEdge(es[vs[u].parent_edge_id]);
        count++;
        centers[u] = vs[u].center;
      }
    }

    if (tmp.size() <= 1 && idx == es.Size()) break;

    for (size_t i = 0; i < remaining.size(); i++) {
      const size_t u = remaining[i];
      for (size_t j = g.first_arc[u]; j < g.first_arc[u + 1]; j++) {
        const AKPWArc& a = g.arcs[j];
        if (a.v < u || vs[u].center == vs[a.v].center) {
          continue;
        }
        size_t v1 = Find(centers, u);
        size_t v2 = Find(centers, a.v);
        es2.AddEdge(AKPWEdge(v1, v2, a.original_id));
      }
    }

    remaining.swap(tmp);
    tmp.clear();
    // std::cout << tree.Size() << '\n';
  }

  // std::cout << tree.Size() << '\n';
}
