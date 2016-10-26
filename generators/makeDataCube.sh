g++ genCayley.cpp -o genCayley -O2
g++ graph2MM.cpp -o graph2MM -O2 -std=c++11

./genCayley 125000 1 1 50 1 2500 1 > temp
./graph2MM temp ../graphs/cayley_cube_125000_uniform/graph.mtx

./genCayley 1000000 1 1 100 100 10000 100000 > temp
./graph2MM temp ../graphs/cayley_cube_1000000_geometric/graph.mtx
