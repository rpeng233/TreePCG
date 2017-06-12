CXX = g++
CXXFLAGS += -I./include
CXXFLAGS += -Wall -pedantic
CXXFLAGS += -O2
# CXXFLAGS += -g
# CXXFLAGS += -ffast-math
LDFLAGS += -lcholmod -lamd -lcolamd -lcamd -lccolamd -lmetis

all : bin/test bin/es2mtx

bin/test : include/*.h
bin/test : obj/tree_solver.o obj/cholesky.o obj/stretch.o
bin/test : obj/akpw.o obj/aug_tree_precon.o obj/sparse_cholesky.o
bin/test : obj/cholmod_solver.o obj/incomplete_cholesky.o
bin/test : obj/flow_gradient_solver.o obj/cycle_toggling_solver.o
bin/test : obj/partial_cholesky.o obj/aug_tree_chain.o obj/precon_chain_solver.o
bin/test : obj/mmio.o
bin/test : src/test.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< obj/*.o -o $@ $(LDFLAGS)

bin/es2mtx : include/*.h
bin/es2mtx : obj/mmio.o
bin/es2mtx : src/es2mtx.cpp
	$(CXX) $(CXXFLAGS) $< obj/mmio.o -o $@

obj/tree_solver.o : include/*.h
obj/tree_solver.o : src/tree_solver.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/cholesky.o : include/*.h
obj/cholesky.o : src/cholesky.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/stretch.o : include/*.h
obj/stretch.o : src/stretch.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/akpw.o : include/*.h
obj/akpw.o : src/akpw.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/aug_tree_precon.o : include/*.h
obj/aug_tree_precon.o : src/aug_tree_precon.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/cholmod_solver.o : include/*.h
obj/cholmod_solver.o : src/cholmod_solver.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/sparse_cholesky.o : include/*.h
obj/sparse_cholesky.o : src/sparse_cholesky.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/incomplete_cholesky.o : include/*.h
obj/incomplete_cholesky.o : src/incomplete_cholesky.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/partial_cholesky.o : include/*.h
obj/partial_cholesky.o : src/partial_cholesky.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/aug_tree_chain.o : include/*.h
obj/aug_tree_chain.o : src/aug_tree_chain.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/precon_chain_solver.o : include/*.h
obj/precon_chain_solver.o : src/precon_chain_solver.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/flow_gradient_solver.o : include/*.h
obj/flow_gradient_solver.o : src/flow_gradient_solver.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

obj/cycle_toggling_solver.o : include/*.h
obj/cycle_toggling_solver.o : src/cycle_toggling_solver.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 $< -c -o $@

obj/mmio.o : src/mmio.c include/mmio.h
	$(CXX) $(CXXFLAGS) -w $< -c -o $@

bin/gen_cayley : generators/gen_cayley.cpp include/graph.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_graph: generators/graph_to_graph.cpp include/graph.h
	$(CXX) $(CXXFLAGS) $< -o $@

bin/graph_to_matrix: generators/graph_to_matrix.cpp include/graph.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	rm -rf bin/*
	rm -rf obj/*
