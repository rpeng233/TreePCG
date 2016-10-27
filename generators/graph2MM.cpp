/*
 * Converts .txt files in graph sparse format to matrix market format
 * 
 * takes parameters:
 *     input file name
 *     output file name
 *
 * Copied from graphSp2MM written by Haoran,
 *     editted by Richard on Oct 26, 16
 */

#include <iostream>
#include <cstdio>
#include <string>
#include "io.h"
#include "graph.h"
#include "matrix.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Please specify graph file and output file" << endl;
        return -1;
    }

    string graphFile = argv[1];
    string outFile = argv[2];
    GraphSP g = IO::readGraph(graphFile);
    Mat A = IO::constructMatrixFromGraph(g);
    IO::saveMM(A, outFile);
    return 0;
}
