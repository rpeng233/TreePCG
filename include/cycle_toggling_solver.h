/* A Cycle Toggling Solver
 *
 * Given as input is a rooted tree plus the off-tree edges. Tree edges are
 * oriented from the child to the parent, and off-tree edges are oriented from
 * the endpoint with the smaller ID to the endpoint with bigger ID.
 *
 * In order to add a one unit of cycle flow formed by an off-tree edge, we
 * simply add one unit of flow from one endpoint of the edge all the way up to
 * the root, then add one negative unit of flow from the other point to the
 * root.
 *
 * The tree is decomposed into a number of chains using the heavy-light
 * decomposition. We let each internal vertex to own the edge to its parent.
 * Due to our edge orientation, each chain points upward to the head of the
 * chain. Thus a flow on the tree can be decomposed into a union of flows on
 * the chains.
 *
 * The flow on a chain is represented by a BST, where each BST vertex
 * corresponds to a vertex in the chain. The BST is sorted according the
 * position in the chain from the head. Each BST vertex maintains flow
 * information from itself to its left most descendent (including that
 * descendent's edge pointing up the chain). This gives a decomposition of any
 * flow on the chain. See below for an example. The BST of different chains are
 * then linked together according to the original tree, forming a virtual tree.
 *
 * The BST:                            4
 *                                    / \
 *                                   /   \
 *                                  /     \
 *                                 /       \
 *                                /         \
 *                               /           \
 *                              /             \
 *                             2               6
 *                            / \             / \
 *                           /   \           /   \
 *                          /     \         /     \
 *                         1       3       5       7
 *
 *
 * The chain:           <--1<--2<--3<--4<--5<--6<--7
 *
 * Flow mained at each BST vertex:
 *
 *                      <--1
 *                      <------2
 *                              <--3
 *                      <--------------4
 *                                      <--5
 *                                      <------6
 *                                              <--7
 *
 * Since all we ever do is to push flow up to the head of the chain, the number
 * of updates we need is porportional to the height of the tree. For instance,
 * if we want to push flow from vertex 6 up to the, we only * need to update
 * the flows starting from 6 and 4.
 * 
 * In addition to the decomposed chain flow, each BST vertex also maintain the sum
 * of resistance-flow product on the range its responsible fore. Unlike the
 * decomposed chain flow which are kept separate (in a way, maintained lazily),
 * this value needs to be kept up to date, since we want to query it by only
 * walking up the tree.
 *
 * This is done by propogating the delta in the sum up the tree whenever we
 * change the flows. Whenever we walk up from a left child, we need to update
 * the sum at current vertex by whatever happened in its left subtree. In the
 * above example, if we push flow from vertex 3, we not only need to update the
 * flow values at 2 and 3, but we also need to update the resistance*flow sum
 * at vertex 4.
 *
 * To query the resistance*flow sum, we walk up the tree collecting the sum at
 * each vertex. However since flows are lazily kept seperate, we need to also
 * add in the contribution from the upper level whenever we came up from a left
 * child. Once again in the example above, if we query from vertex 3, we
 * combine the sum values at veretx 2 and 3, but vertex 4 might be sending a
 * non-zero flow.
 *
 * To combine the decomposed flow into a single flow, we perform a pre-order
 * traversal on the BST, pushing down the flow to the left subtree.
 *
 */

#ifndef INCLUDE_CYCLE_TOGGLING_SOLVER_H__
#define INCLUDE_CYCLE_TOGGLING_SOLVER_H__

#include "common.h"
#include "graph.h"

class CycleTogglingSolver {
public:
  CycleTogglingSolver(const TreeR& tree_, const EdgeListR& off_tree_es_);

private:
  struct Node {
    size_t parent;
    double resistance;
    double flow;
    double resistance_to_root;
  };

  static const char LEFT = 0;
  static const char RIGHT = 1;
  static const char VIRTUAL = 2;

  struct VirtualNode {
    size_t parent;
    size_t left_child;
    size_t right_child;
    double rf_sum;                // the sum of resistance*flow to the left most descendent
    double flow_to_left;          // flow to the left most descendent
    double resistance_to_left;    // resistance to the left most descendent
    char type;                    // One of LEFT, RIGHT and VIRTUAL
  };

  struct OffTreeEdge {
    size_t u;
    size_t v;
    size_t lca;
    double flow;
    double resistance;
  }

  std::vector<Node> tree;
  std::vector<VirtualNode> hld;
  std::vector<OffTreeEdge> es;
};

#endif  // INCLUDE_CYCLE_TOGGLING_SOLVER_H__
