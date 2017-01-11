CXX = g++
CXXFLAGS += -I./include
CXXFLAGS += -Wall -pedantic
# CXXFLAGS += -O3
# CXXFLAGS += -g
CXXFLAGS += 


all : bin/test bin/gen_cayley bin/graph_to_matrix bin/graph_to_graph

bin/test : include/*.h
bin/test : src/test.cpp obj/tree_solver.o obj/graph.o
	$(CXX) $(CXXFLAGS) $< obj/*.o -o $@

obj/tree_solver.o : src/tree_solver.cpp include/tree_solver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/graph.o : src/graph.cpp include/graph.h
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
