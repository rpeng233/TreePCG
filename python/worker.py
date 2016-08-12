import scipy as sp
import numpy as np
import math as math
import time
import csv
import sys

from decimal import *
from scipy import io
from math import sqrt

def pcg(mat, b, pre, lhs, maxits, verbose):
	print "Running pcg with ", getcontext().prec, " precision"

	n = len(b)

	x = np.array([Decimal('0') for i in range(n)])

	r = np.copy(b)
	z = pre(r)
	p = np.copy(z)

	rho = np.dot(r, z)

	allErrMN = np.array([])
	allErr2 = np.array([])

	t = time.time()

	itcnt = 0
	while itcnt < maxits:
		itcnt += 1

		q = mulIJVVec(mat,p)

		al = rho / np.dot(p,q)

		x = x + al * p
		r -= al * q

		errMN = (np.dot((lhs - x), mulIJVVec(mat, lhs - x)) / np.dot(lhs, mulIJVVec(mat, lhs))).sqrt()
		err2 = np.linalg.norm(mulIJVVec(mat, x) - b) / np.linalg.norm(b)

		allErrMN = np.append(allErrMN, float(errMN))
		allErr2 = np.append(allErr2, float(err2))

		z = pre(r)

		oldrho = rho
		rho = np.dot(z,r)

		beta = rho/oldrho
		p = z + beta * p

		if verbose:
			print "iteration ", itcnt
			print '%.100f' % errMN
			print '%.100f' % err2
			print 'total runtime so far ', time.time() - t
			print ""

		if itcnt % 100 == 0:
			print 'saving results after ', itcnt, ' iterations'

			with open('_log_python_tree_' + str(treeInd) + '_' + str(getcontext().prec) + '_' + str(len(allErr2)) + 'iters.csv', 'w') as csvfile:
			    fieldnames = ['A-norm', '2-norm']
			    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

			    writer.writeheader()
			    for i in range(len(allErr2)):
				    writer.writerow({'A-norm': allErrMN[i], '2-norm': allErr2[i]})



	return x,allErrMN,allErr2

def treeSolver(tree, diag):
	n = tree.shape[0]

	perm = bfsOrd(tree)
	invperm = np.copy(perm)
	perm = np.argsort(perm)

	permTree = permMat(tree, perm)
	father = [0 for i in range(n)]

	for u in range(1,n):
		nbrs = sp.sparse.find(permTree[u,:])[1]

		for v in nbrs:
			if v < u:
				father[u] = v
				break

	# compute the stuff on the diagonal
	auxgeld = [Decimal(el) for el in permArray(diag,perm)]
	row,col,weight = sp.sparse.find(permTree)
	for i in range(len(weight)):
		auxgeld[row[i]] += Decimal(weight[i])

	def f(b):
		geld = np.copy(auxgeld)

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
					fact = Decimal(w / geld[u])

					geld[v] -= w * fact
					aux[v] += aux[u] * fact

		res = np.array([Decimal('1') for i in range(n)])
		if np.linalg.norm(diag) != 0:
			res[0] = res[0] / geld[0]

		for i in range(1,n):
			res[i] = (aux[i] + Decimal(permTree[father[i],i]) * res[father[i]]) / geld[i]

		# print "solve accuracy ", '%.100f' % (np.linalg.norm(mulIJVVec(lap(permTree), res) - permArray(b,perm)) / np.linalg.norm(b))

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

	newx = np.array([Decimal(el[0]) for el in x])

	return newx - np.mean(newx)

def lap(M):
	n = M.shape[0]

	u,v,w = sp.sparse.find(M)

	u = np.array(u)
	v = np.array(v)
	w = np.array([-Decimal(el) for el in w])

	diag = np.array([Decimal('0') for i in range(n)])
	for i in range(len(u)):
		diag[u[i]] -= w[i]

	u = np.append(u, np.array([i for i in range(n)]))
	v = np.append(v, np.array([i for i in range(n)]))
	w = np.append(w, diag)

	return (u,v,w)

def permMat(M,perm):
	row,col,v = sp.sparse.find(M)
	return sp.sparse.coo_matrix((v, (perm[row], perm[col]))).tocsc()

def mulIJVVec(IJV, x):
	res = np.array([Decimal('0') for i in range(len(x))])

	nnz = len(IJV[0])
	for ind in range(nnz):
		res[IJV[0][ind]] += Decimal(IJV[2][ind]) * x[IJV[1][ind]]

	return res


# example code



# for curPrecision in precision:

precision = int(sys.argv[1])
treeInd = int(sys.argv[2])

print 'solving for tree ', treeInd, ' with precision ', precision

getcontext().prec = precision

path = "/Users/serbanstan/git/TreePCG/graphs/rand_1000_u1000/"

A = readMatrix(path + "graph.mtx");
tree = readMatrix(path + "tree" + str(treeInd) + ".mtx");
truex = readArray(path + "x.vec");

b = mulIJVVec(lap(A), truex)
b = b - np.mean(b) # with the newest changes this line should be redundant

f = treeSolver(tree, np.array([0 for i in range(tree.shape[1])]))

# myx = f(b)
# print "per solve precision ", '%.20f' % (np.linalg.norm(mulIJVVec(lap(tree), myx) - b) / np.linalg.norm(b))

myx,myErrMN,myErr2 = pcg(lap(A), b, f, truex, 100, True)

with open('_log_python_tree_' + str(treeInd) + '_' + str(getcontext().prec) + '_' + str(len(myErr2)) + 'iters.csv', 'w') as csvfile:
    fieldnames = ['A-norm', '2-norm']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for i in range(len(myErr2)):
	    writer.writerow({'A-norm': myErrMN[i], '2-norm': myErr2[i]})



