#ifndef INCLUDE_AKPW_H
#define INCLUDE_AKPW_H

#include "graph.h"
void AKPW(const EdgeList<EdgeR>& es, EdgeList<EdgeR>& tree);
void AKPW_unweighted(const EdgeList<EdgeR>& es_, EdgeList<EdgeR>& tree);

#endif // INCLUDE_AKPW_H
