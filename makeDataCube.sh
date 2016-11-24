g++ genCayley.cpp -o genCayley -O2
g++ graphBin2MM.cpp -o graphBin2MM -O2 -std=c++11 -Wl,--stack,1000000000

./genCayley ../graphs/cayley_cube_125000_uniform/graph.bin 125000 1 1 50 1 2500 1
./graphBin2MM ../graphs/cayley_cube_125000_uniform/graph.bin ../graphs/cayley_cube_125000_uniform/graph.mtx

./genCayley ../graphs/cayley_cube_1000000_geometric/graph.bin 1000000 1 1 100 100 100000 100000
./graphBin2MM ../graphs/cayley_cube_1000000_geometric/graph.bin ../graphs/cayley_cube_1000000_geometric/graph.mtx

./genCayley ../graphs/cayley_cube_8000000_uniform/graph.bin 8000000 1 1 200 1 400000 1
./graphBin2MM ../graphs/cayley_cube_8000000_uniform/graph.bin ../graphs/cayley_cube_8000000_uniform/graph.mtx
