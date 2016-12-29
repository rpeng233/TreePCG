#ifndef PCHOL_H
#define PCHOL_H

#include <vector>
#include <unordered_map>
#include "common.h"
#include "graph.h"


class PartialElim {
public:
	size_t v;
	size_t nghbr1;
	size_t nghbr2;
	FLOAT resistance1;
	FLOAT resistance2;
  FLOAT rhs;
	char deg;

	PartialElim(size_t v_, size_t p, FLOAT r_)
	{
		v = v_;
		nghbr1 = p;
		resistance1 = r_;
		rhs = 0;
		deg = 1;
	}

	PartialElim(size_t v_, size_t n1, FLOAT r1, size_t n2, FLOAT r2)
	{
		v = v_;
		nghbr1 = n1;
		nghbr2 = n2;
		resistance1 = r1;
		resistance2 = r2;
		rhs = 0;
		deg = 2;
	}
};


class PartialCholesky {
public:
  std::vector<PartialElim> elims;
  std::unordered_map<size_t, size_t> relabeling;
};

std::pair<TreePlusEdges, PartialCholesky>
partial_cholesky(const TreePlusEdges& tree);

std::vector<FLOAT> eliminate_rhs(const std::vector<FLOAT>& rhs_,
                                 PartialCholesky& pchol);

void back_substitution(const PartialCholesky& pchol,
                       const std::vector<FLOAT>& partial_x,
                       std::vector<FLOAT>& x);

#endif
