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

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "FORMAT: graph file, format, output file" << endl;
        return -1;
    }

    string graph_file = argv[1];
    char format = IO::ParseFormat(argv[2]);
    string out_file = argv[3];

    fprintf(stderr, "CONVERTING %s in FORMAT %s:\n",
      graph_file.c_str(), argv[2]);
    fprintf(stderr, "TO MTX FILE: %s\n", out_file.c_str());

    Graph graph = IO::ReadGraph(graph_file, format);
    fprintf(stderr, "FINISHED READING FILE\n");
    Matrix matrix = IO::GraphToMatrix(graph);
    fprintf(stderr, "CONSTRUCTED MATRIX\n");
    IO::WriteMMMatrix(out_file, matrix);
    fprintf(stderr, "DONE SAVING TO FILE\n");
    return 0;
}
