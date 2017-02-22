/* Pairing heap for POD */

#ifndef PAIRING_HEAP_H
#define PAIRING_HEAP_H

#include <cassert>
#include <iostream>
#include <vector>

class PairingNode {

public:
  PairingNode *prev;
  PairingNode *next;
  PairingNode *first_child;

  PairingNode() {
    prev = next = first_child = nullptr;
  }

  void reset() {
    prev = next = first_child = nullptr;
  }
};


template <typename T, std::size_t Offset, typename Comp>
class PairingHeap {

public:
  PairingHeap();
  std::size_t Size();
  void Insert(T& n);
  void Decrease_key(T& n);
  T& Get_min();
  void Delete_min();
  void Clear();

private:
  std::size_t hsize;
  PairingNode *root;

  PairingNode* link(PairingNode *x, PairingNode *y);
};


template <typename T, std::size_t Offset, typename Comp>
PairingHeap<T, Offset, Comp>::PairingHeap() {
  hsize = 0;
  root = nullptr;
}


template <typename T, std::size_t Offset, typename Comp>
std::size_t PairingHeap<T, Offset, Comp>::Size() {
  return hsize;
}


template <typename T, std::size_t Offset, typename Comp>
void PairingHeap<T, Offset, Comp>::Insert(T& t) {
  PairingNode *n = (PairingNode *) ((char *) &t + Offset);
  if (n == nullptr) return;
  assert(n->prev == nullptr);
  assert(n->next == nullptr);
  assert(n->first_child == nullptr);
  root = link(root, n);
  hsize++;
}


template <typename T, std::size_t Offset, typename Comp>
void PairingHeap<T, Offset, Comp>::Decrease_key(T& t) {
  PairingNode *n = (PairingNode *) ((char *) &t + Offset);
  if (n == nullptr || n == root) return;
  assert(n->prev != nullptr);

  if (n->next != nullptr) {
    n->next->prev = n->prev;
  }

  if (n->prev->first_child == n) {
    n->prev->first_child = n->next;
  } else {
    n->prev->next = n->next;
  }

  n->prev = n->next = nullptr;
  root = link(root, n);
}


template <typename T, std::size_t Offset, typename Comp>
T& PairingHeap<T, Offset, Comp>::Get_min() {
  return *((T *) ((char *) root - Offset));
}


template <typename T, std::size_t Offset, typename Comp>
void PairingHeap<T, Offset, Comp>::Delete_min() {
  if (hsize == 0) return;

  PairingNode *n = root->first_child;
  std::vector<PairingNode *> tmp_list;

  while (n != nullptr && n->next != nullptr) {
    PairingNode *n1 = n;
    PairingNode *n2 = n->next;
    n = n2->next;
    n1->prev = n1->next = nullptr;
    n2->prev = n2->next = nullptr;
    tmp_list.push_back(link(n1, n2));
  }

  if (n != nullptr) {
    n->prev = nullptr;
  }

  for (PairingNode *m : tmp_list) {
    assert(m->prev == nullptr);
    n = link(n, m);
  }

  root->prev = root->next = root->first_child = nullptr;
  root = n;
  hsize--;
}

template <typename T, std::size_t Offset, typename Comp>
void PairingHeap<T, Offset, Comp>::Clear() {
  hsize = 0;
  root = nullptr;
}

template <typename T, std::size_t Offset, typename Comp>
PairingNode *PairingHeap<T, Offset, Comp>::link(PairingNode *x, PairingNode *y) {
  if (x == nullptr) return y;
  if (y == nullptr) return x;

  assert(x->prev == nullptr);
  assert(x->next == nullptr);
  assert(y->prev == nullptr);
  assert(y->next == nullptr);

  T *tx = (T *)((char *) x - Offset);
  T *ty = (T *)((char *) y - Offset);

  Comp cmp;

  if (cmp(*tx, *ty)) {
    y->prev = x;
    y->next = x->first_child;
    if (x->first_child != nullptr) {
      x->first_child->prev = y;
    }
    x->first_child = y;
    assert(x->prev == nullptr);
    assert(x->next == nullptr);
    return x;
  } else {
    x->prev = y;
    x->next = y->first_child;
    if (y->first_child != nullptr) {
      y->first_child->prev = x;
    }
    y->first_child = x;
    assert(y->prev == nullptr);
    assert(y->next == nullptr);
    return y;
  }
}

#endif
