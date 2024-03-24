#!/usr/bin/env bash

python3 FOAM_to_readable.py
g++ -O2 -c -std=c++17 Mesh.cpp -o Mesh.o
g++ -O2 -c -std=c++17 SparseAddress.cpp -o SparseAddress.o
g++ -O2 -c -std=c++17 SparseMat.cpp -o SparseMat.o
g++ -O2 -c -std=c++17 SLS.cpp -o SLS.o
g++ -O2 -c -std=c++17 GSSmooth.cpp -o GSSmooth.o
g++ -O2 -c -std=c++17 GSSolver.cpp -o GSSolver.o
g++ -O2 -c -std=c++17 SIMPLE_Solve.cpp -o SIMPLE_Solve.o
g++ SparseAddress.o GSSmooth.o GSSolver.o Mesh.o SparseMat.o SLS.o SIMPLE_Solve.o -o SIMPLE_Solve_execute
./SIMPLE_Solve_execute