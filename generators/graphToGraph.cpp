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

int main(int argc, char *argv[]) {
    if (argc < 5) {
        cerr << "FORMAT: input, format, output, format" << endl;
        return -1;
    }

    fprintf(stderr, "CONVERTING GRAPH %s in FORMAT %s\n", argv[1], argv[2]);
    fprintf(stderr, "TO GRAPH: %s in FORAMT %s\n", argv[3], argv[4]);

    string in_file = argv[1];
    char in_format = IO::ParseFormat(argv[2]);
    string out_file = argv[3];
    char out_format = IO::ParseFormat(argv[4]);

    Graph graph = IO::ReadGraph(in_file, in_format);
    fprintf(stderr, "FINISHED READING FILE\n");
    IO::WriteGraph(out_file, out_format, graph);
    fprintf(stderr, "DONE SAVING TO FILE\n");
    return 0;
}
