import scipy as sp
import numpy as np
import math as math

from decimal import *
from scipy import io
from math import sqrt

def pcg(mat, b, pre, lhs, maxits, verbose):
	n = mat.shape[0]

	x = np.array([Decimal(0) for i in range(n)])

	r = np.copy(b)
	z = pre(r)
	p = np.copy(z)

	rho = np.dot(r, z)

	print "Running pcg with ", getcontext().prec, " precision"

	itcnt = 0
	while itcnt < maxits:
		itcnt += 1

		q = mulMatVec(mat,p)

		al = rho / np.dot(p,q)

		x = x + al * p
		r -= al * q

		errMN = (np.linalg.norm((lhs - x) * mulMatVec(mat, lhs - x)) / np.linalg.norm(lhs * mulMatVec(mat, lhs))).sqrt()
		err2 = np.linalg.norm(mulMatVec(mat, x) - b) / np.linalg.norm(b)

		if verbose:
			print "iteration ", itcnt
			print '%.100f' % errMN
			print '%.100f' % err2
			print ""

			# print al

		z = pre(r)

		oldrho = rho
		rho = np.dot(z,r)

		beta = rho/oldrho
		p = z + beta * p

	return x

def treeSolver(tree, diag):
	n = tree.shape[0]

	perm = bfsOrd(tree)
	invperm = np.copy(perm)

	perm = np.argsort(perm)

	permTree = permMat(tree, perm)
	permLapTree = lap(permTree) + sp.sparse.csc_matrix(np.diag(permArray(diag,perm)))

	father = [0 for i in range(n)]

	for u in range(1,n):
		nbrs = sp.sparse.find(permTree[u,:])[1]

		for v in nbrs:
			if v < u:
				father[u] = v
				break

	def f(b):
		geld = np.array([Decimal(permLapTree[i,i]) for i in range(n)])
		aux = np.array([Decimal(el) for el in permArray(b, perm)])
		if np.linalg.norm(diag) == 0:
			aux = aux - np.mean(aux)

		for u in reversed(range(n)):
			if u == 0:
				break

			nbrs = sp.sparse.find(permTree[u,:])[1]
			weights = sp.sparse.find(permTree[u,:])[2]

			for i in range(len(nbrs)):
				v = nbrs[i]
				w = Decimal(weights[i])

				if v < u:
					fact = w / geld[u]

					geld[v] -= w * fact
					aux[v] += aux[u] * fact

		res = np.array([Decimal(1) for i in range(n)])
		if np.linalg.norm(diag) != 0:
			res[0] = res[0] / geld[0]

		for i in range(1,n):
			res[i] = (aux[i] + Decimal(permTree[father[i],i]) * res[father[i]]) / geld[i]

		if np.linalg.norm(diag) == 0:
			res = res - np.mean(res)
		res = permArray(res,invperm)

		return res

	return f

# return the bfs order traversal of a tree
def bfsOrd(tree):
	n = tree.shape[0]

	Q = [0]

	vis = [0 for i in range(n)]
	vis[0] = 1

	left = -1
	right = 0
	while left < right:
		left = left + 1

		u = Q[left]

		neigh = sp.sparse.find(tree[u,:])[1]

		for v in neigh:
			if vis[v] == 0:
				vis[v] = 1
				right = right + 1
				Q.append(v)

	return np.array(Q)

def permArray(v, perm):
	res = np.copy(v)
	for i in range(len(v)):
		res[perm[i]] = v[i]
	return res

# read a matrix in matrix market format and turn it into a CSC matrix with a certain precision
def readMatrix(fileName):
	A = sp.io.mmread(fileName)
	A = A.tocsc()

	return A

def readArray(fileName):
	x = sp.io.mmread(fileName)

	newx = []
	for elem in x:
		newx.append(Decimal(elem[0]))
	newx = np.array(newx);

	return newx - np.mean(newx)

def lap(M):
	n = M.shape[0]

	row,col,v = sp.sparse.find(M)

	diag = [0 for i in range(n)]
	for i in range(len(v)):
		diag[row[i]] += v[i]
		diag[col[i]] += v[i]
	diag = np.array(diag)

	row = np.append(row, [i for i in range(n)])
	col = np.append(col, [i for i in range(n)])
	v = np.append(-v, diag / 2)

	return sp.sparse.coo_matrix((v, (row, col))).tocsc()

def permMat(M,perm):
	row,col,v = sp.sparse.find(M)
	return sp.sparse.coo_matrix((v, (perm[row], perm[col]))).tocsc()


def mulMatVec(M, x):
	res = np.array([Decimal(0) for i in range(M.shape[0])])

	row,col,v = sp.sparse.find(M)

	for ind in range(len(v)):
		res[row[ind]] += Decimal(v[ind]) * x[col[ind]]

	return res


# example code

getcontext().prec = 1000

path = "/Users/serbanstan/git/TreePCG/graphs/pathDisjoint_1000_exp30/"

A = readMatrix(path + "graph.mtx");
tree = readMatrix(path + "tree.mtx");
truex = readArray(path + "x.vec");
b = mulMatVec(lap(A), truex)
f = treeSolver(tree, np.array([0 for i in range(tree.shape[1])]))

ss = np.sum(truex)

# myx = f(b)
# print np.linalg.norm(mulMatVec(lap(tree), myx) - b) / np.linalg.norm(b)

# print sum(truex)

for x in mulMatVec(lap(A), mulMatVec(lap(A), truex)):
	print '%.100f' % x


# pcg(lap(A), b, f, truex, 30, True)