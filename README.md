This code implements two solvers based on low stretch spanning trees:
* directly precondition by the low stretch tree,
* sample a number of edges based on stretch, use that solve as a preconditioner

Both of these methods only have one layer of iterative methods. So preconditioned conjugate gradient works, and is the method of choice.


The subdirectories are:
* graphs: folders containing the test data
     each folder contains a graph.bin/asc/mtx file with the graph.
     tree.?.mtx: spanning trees containing the indices of edges that form a low stretch spanning tree.

* common: `standard' structs used, and routines for reading/writing from them.

     common.h: common header files

     graph.h: Graphs stored as adjacency lists, and
		TreePlusEdges with tree + list of off tree edges

     matrix.h: matrix / vectors
           compatible with matrix market formats. 

     io.h: methods for reading / writing these. Supported formats include:
           ASC: ASCII
           BIN: binary, 4 bit ints and etc
           DIMACS: DIMACS format with comments/flags to indicate edges

     formats are picked via a flag in the function call.


* generators: various generators for graphs including:
        Cayley graphs
        RMF???

     Tools for converting between the various graph / matrix formats


* lowStretchTrees
       treeFinder.h: methods for finding low stretch spanning trees.
       stretchStatistics.cpp: takes a graph and a tree, outputs stats about the distribution of stretches.
       buildTree.cpp: builds tree from a graph, outputs to a tree file.


* solvers:
       solver.h: define the standard API of a solver.
       PCG.h: preconditioned conjugate gradient code.