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
    if (argc < 5) {
        fprintf(stderr, "FORMAT: graph file, format\n");
        fprintf(stderr, "        tree file, format\n");
        return -1;
    }

    string graph_file = argv[1];
    char graph_format = IO::ParseFormat(argv[2]);
    string tree_file = argv[3];
    char tree_format = IO::ParseFormat(argv[4]);

    Graph graph = IO::ReadGraph(graph_file, graph_format);
    graph.SortAndCombine();
    fprintf(stderr, "FINISHED READING GRAPH\n");

// IO::WriteGraph("__stdout", IO::ASCII, graph);

    vector<int> tree = IO::ReadVectorInt(tree_file, tree_format, graph.n);
//for (int i = 0; i < graph.n; ++i) fprintf(stderr, "%d\n", tree[i]);
    fprintf(stderr, "FINISHED READING TREE\n");

    TreePlusEdges g1 = GraphToTreePlusEdges(graph, tree);
    fprintf(stderr, "Total Stretch = %lf\n",
      TreeFinder::TotalStretch(g1));
    vector<FLOAT> stretch = TreeFinder::GetStretch(g1);

    sort(stretch.begin(), stretch.end());
    double min_stretch = stretch[0];
    fprintf(stderr, "Minimum Stretch = %0.6lf\n", PrintFloat(min_stretch));
    double max_stretch = stretch[stretch.size() - 1];
    fprintf(stderr, "Maximum Stretch = %0.6lf\n", PrintFloat(max_stretch));
    int l = 0, r = 0;
    for(double x = min_stretch; x <= max_stretch; x *= double(2)) {
      while(l < (int)stretch.size() && stretch[l] < x) {
        l++;
      }
      while(r < (int)stretch.size() && stretch[r] < x * double(2)) {
        r++;
      }
      printf("# of edges with stretch between [%0.6lf, %0.6lf]: %d\n",
        x, x * double(2), r - l);
    }
    return 0;
}
