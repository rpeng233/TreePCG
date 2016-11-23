/*
 * Converts .txt files in graph sparse format to matrix market format
 * 
 * takes parameters:
 *     input file name
 *     output file name
 *
 * Copied from graphSp2MM written by Haoran,
 *     editted by Richard 
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
        cerr << "FORMAT: graph file and output file" << endl;
        return -1;
    }

    fprintf(stderr, "CONVERTING BINARY FILE:\n");
    fprintf(stderr, "%s\n", argv[1]);
    fprintf(stderr, "TO MTX FILE:\n");
    fprintf(stderr, "%s\n", argv[2]);

    string graphFile = argv[1];
    string outFile = argv[2];

    Graph graph = IO::ReadGraph(graphFile, 1);
    fprintf(stderr, "FINISHED READING FILE\n");
    Matrix matrix = IO::GraphToMatrix(graph);
    fprintf(stderr, "CONSTRUCTED MATRIX\n");
    IO::SaveMMMatrix(matrix, outFile);
    fprintf(stderr, "DONE SAVING TO FILE\n");
    return 0;
}
