make
./genCayley temp.bin B 11 2 7 4 15
./graphToGraph temp.bin B temp.asc.txt A
./graphToGraph temp.bin B temp.dimacs.txt D
./graphToMatrix temp.bin B temp.MM.txt
