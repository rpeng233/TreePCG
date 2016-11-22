/*
 * Converts .txt files in graph sparse format to matrix market format
 * 
 * takes parameters:
 *     input file name
 *     output file name
 *
 * Copied from graphSp2MM written by Haoran,
 *     editted by Richard on Nov 22, 16
 */

#include <iostream>
#include <cstdio>
#include <string>
#include "../common/io.h"
#include "../common/graph.h"
#include "../common/matrix.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Please specify graph file and output file" << endl;
        return -1;
    }

    string graph_file = argv[1];
    string out_file = argv[2];

    Graph g = IO::ReadGraph(graph_file, 0);
    Matrix laplacian_g = IO::GraphToMatrix(g);
    IO::SaveMMMatrix(laplacian_g, out_file);
    return 0;
}
