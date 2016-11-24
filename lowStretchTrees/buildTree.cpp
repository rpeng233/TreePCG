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

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "FORMAT: tree type, graph file, output tree file" << endl;
        return -1;
    }

    string tree_type = argv[1];
    string graph_file = argv[2];
    string out_file = argv[3];

    Graph graph = IO::ReadGraph(graph_file, 1);
    fprintf(stderr, "FINISHED READING FILE\n");

    vector<int> tree;
    if (tree_type[0] == 'D') {
      tree = TreeFinder::DijkstraTree(graph);
    } else if (tree_type[0] == 'M') {
      tree = TreeFinder::MinimumSpanningTree(graph);
    }


    IO::WriteVector(out_file, 0, tree);
    fprintf(stderr, "DONE SAVING TO FILE\n");

    TreePlusEdges g1 = GraphToTreePlusEdges(graph, tree);
    fprintf(stderr, "Total Stretch = %lf\n",
      TreeFinder::TotalStretch(g1));
    return 0;
}
