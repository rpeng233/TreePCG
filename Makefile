CXX = g++
CXXFLAGS += -I./include
# CXXFLAGS += -Wall -pedantic
CXXFLAGS += -O3
# CXXFLAGS += -g
# CXXFLAGS += --std=c++11
# CXXFLAGS += -ffast-math


all : bin/test bin/gen_cayley bin/graph_to_matrix bin/graph_to_graph

bin/test : include/*.h
bin/test : obj/tree_solver.o obj/min_degree_solver.o obj/low_stretch_tree.o
bin/test : src/test.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< obj/*.o -o $@

obj/tree_solver.o : include/graph.h
obj/tree_solver.o : src/tree_solver.cpp include/tree_solver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/min_degree_solver.o : include/min_degree_solver.h include/binary_heap.h
obj/min_degree_solver.o : src/min_degree_solver.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/low_stretch_tree.o : include/graph.h include/disjoint_set.h
obj/low_stretch_tree.o : src/low_stretch_tree.cpp include/low_stretch_tree.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

bin/gen_cayley : generators/gen_cayley.cpp include/*.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_graph: generators/graph_to_graph.cpp include/*.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_matrix: generators/graph_to_matrix.cpp include/*h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	rm -rf bin/*
	rm -rf obj/*
