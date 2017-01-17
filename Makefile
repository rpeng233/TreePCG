CXX = g++
CXXFLAGS += -I./include
# CXXFLAGS += -Wall -pedantic
CXXFLAGS += -O3
# CXXFLAGS += -g
# CXXFLAGS += --std=c++11 -ffast-math


all : bin/test bin/gen_cayley bin/graph_to_matrix bin/graph_to_graph

bin/test : include/*.h
bin/test : src/test.cpp obj/tree_solver.o obj/min_degree_solver.o
	$(CXX) $(CXXFLAGS) --std=c++11 $< obj/*.o -o $@

obj/tree_solver.o : src/tree_solver.cpp include/tree_solver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/min_degree_solver.o : include/min_degree_solver.h include/binary_heap.h
obj/min_degree_solver.o : src/min_degree_solver.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

bin/gen_cayley : generators/gen_cayley.cpp include/*.h
	$(CXX) $(CXXFLAGS) $(CFLAGS) $< -o $@

bin/graph_to_graph: generators/graph_to_graph.cpp include/*.h
	$(CXX) $(CXXFLAGS) $(CFLAGS) $< -o $@

bin/graph_to_matrix: generators/graph_to_matrix.cpp include/*h
	$(CXX) $(CXXFLAGS) $(CFLAGS) $< -o $@

clean :
	rm -f bin/*
	rm -f obj/*
