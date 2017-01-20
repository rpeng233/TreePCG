#include <iostream>

template <class Compare>
class BinaryHeap {
public:
  size_t size;
  std::vector<size_t> heap;
  std::vector<size_t> location;
  const Compare& less;

  BinaryHeap(size_t n, const Compare& less_)
    : size(n), heap(n), location(n), less(less_) {
    for (size_t i = 0; i < n; i++) {
      heap[i] = i;
    }
    std::make_heap(heap.begin(), heap.end(), less);
    for (size_t i = 0; i < n; i++) {
      location[heap[i]] = i;
    }
  }

  void BubbleDown(size_t key) {
    size_t l, r, max_child;
    size_t i = location[key];
    assert(i < size);
    do {
      l = i * 2 + 1;
      r = l + 1;
      if (l >= size) break;
      if (r < size && less(heap[l], heap[r])) {
        max_child = r;
      } else {
        max_child = l;
      }
      if (less(heap[i], heap[max_child])) {
        std::swap(location[heap[i]], location[heap[max_child]]);
        std::swap(heap[i], heap[max_child]);
        i = max_child;
      } else {
        break;
      }
    } while (true);
    // for (size_t ii = 0; ii < heap.size(); ii++) {
    //   assert(location[heap[ii]] == ii);
    //   assert(heap[location[ii]] == ii);
    // }
    // assert(std::is_heap(heap.begin(), heap.begin() + size, less));
  }

  void BubbleUp(size_t key) {
    size_t i = location[key];
    assert(i < size);
    size_t p = (i - 1) / 2;
    while (i != 0 && less(heap[p], heap[i])) {
      std::swap(location[heap[p]], location[heap[i]]);
      std::swap(heap[p], heap[i]);
      i = p;
      p = (i - 1) / 2;
    }
    // for (size_t ii = 0; ii < heap.size(); ii++) {
    //   assert(location[heap[ii]] == ii);
    //   assert(heap[location[ii]] == ii);
    // }
    // assert(std::is_heap(heap.begin(), heap.begin() + size, less));
  }

  size_t Top() {
    return heap[0];
  }

  void Pop() {
    std::swap(location[heap[0]], location[heap[size - 1]]);
    std::swap(heap[0], heap[size - 1]);
    size--;
    BubbleDown(heap[0]);
  }
};

