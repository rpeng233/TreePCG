rm *.bin *.txt
./../../generators/genCayley graph.bin B 11 2 4 5 8
./../../generators/graphToGraph graph.bin B graph.asc.txt A
./../../generators/graphToGraph graph.bin B graph.dimacs.txt D
./../../generators/graphToMatrix graph.bin B laplacian.MM.txt
./../../lowStretchTrees/buildTree D graph.bin B tree_D.bin B
./../../lowStretchTrees/buildTree D graph.bin B tree_D.asc.txt A
./../../lowStretchTrees/stretchStatistics graph.bin B tree_D.asc.txt A
./../../lowStretchTrees/buildTree M graph.bin B tree_M.bin B
./../../lowStretchTrees/buildTree M graph.bin B tree_M.asc.txt A
./../../lowStretchTrees/stretchStatistics graph.bin B tree_M.asc.txt A
