CXX = icc
CXXFLAGS += -I./include -I/home/shen/suiteSparse/include
CXXFLAGS += -L/home/shen/suiteSparse/lib
CXXFLAGS += -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis
CXXFLAGS += -Wall -pedantic
CXXFLAGS += -O3
# CXXFLAGS += -g
# CXXFLAGS += -ffast-math


all : bin/test # bin/gen_cayley bin/graph_to_matrix bin/graph_to_graph

bin/test : include/*.h
bin/test : obj/tree_solver.o obj/cholesky.o obj/stretch.o
bin/test : obj/akpw.o obj/aug_tree_precon.o obj/sparse_cholesky.o
bin/test : obj/incomplete_cholesky.o
bin/test : src/test.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< obj/*.o -o $@

obj/tree_solver.o : include/graph.h
obj/tree_solver.o : src/tree_solver.cpp include/tree_solver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/cholesky.o : include/binary_heap.h
obj/cholesky.o : src/cholesky.cpp include/cholesky.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/stretch.o : include/graph.h include/disjoint_set.h
obj/stretch.o : src/stretch.cpp include/stretch.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/akpw.o : include/graph.h include/binary_heap.h include/pairing_heap.h
obj/akpw.o : src/akpw.cpp include/akpw.h
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/aug_tree_precon.o : include/cholesky.h include/graph.h
obj/aug_tree_precon.o : src/aug_tree_precon.cpp include/aug_tree_precon.h
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/sparse_cholesky.o : include/cholesky.h include/graph.h
obj/sparse_cholesky.o : src/sparse_cholesky.cpp include/sparse_cholesky.h
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/incomplete_cholesky.o : include/cholesky.h include/graph.h
obj/incomplete_cholesky.o : src/incomplete_cholesky.cpp include/incomplete_cholesky.h
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

bin/gen_cayley : generators/gen_cayley.cpp include/*.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_graph: generators/graph_to_graph.cpp include/*.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_matrix: generators/graph_to_matrix.cpp include/*h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	rm -rf bin/*
	rm -rf obj/*
