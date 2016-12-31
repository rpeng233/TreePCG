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
#include <vector>
#include "../common/io.h"
#include "../common/graph.h"
#include "../lowStretchTrees/treeFinder.h"

int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr, "FORMAT: tree type, graph file, format\n");
        fprintf(stderr, "        output tree file, format\n");
        return -1;
    }

    string tree_type = argv[1];
    string graph_file = argv[2];
    char graph_format = IO::ParseFormat(argv[3]);
    string out_file = argv[4];
    char out_format = IO::ParseFormat(argv[5]);

    Graph graph = IO::ReadGraph(graph_file, graph_format);
    graph.SortAndCombine();
    fprintf(stderr, "Read in graph\n");

// IO::WriteGraph("__stdout", IO::ASCII, graph);

    vector<int> tree;
    if (tree_type[0] == 'D') {
      tree = TreeFinder::DijkstraTree(graph);
    } else if (tree_type[0] == 'M') {
      tree = TreeFinder::MinimumSpanningTree(graph);
    }


    IO::WriteVectorInt(out_file, out_format, tree);
    fprintf(stderr, "Saved tree to file %s\n", argv[4]);

    TreePlusEdges g1 = GraphToTreePlusEdges(graph, tree);
    fprintf(stderr, "Total Stretch = %lf\n",
      TreeFinder::TotalStretch(g1));
    return 0;
}
