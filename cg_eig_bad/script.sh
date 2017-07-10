#!/bin/bash
for d in ../diagonals/*/
do
    echo $d
    ./cg_solver $d > $d/log$1.txt
done
