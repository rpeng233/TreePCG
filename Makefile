CXX = g++ --std=c++11
CXXFLAGS += -I./include
CXXFLAGS += -Wall -pedantic
# CXXFLAGS += -O3
# CXXFLAGS += -g
CXXFLAGS += 


all : bin/test

bin/test : include/*.h
bin/test : src/test.cpp obj/tree_solver.o
	$(CXX) $(CXXFLAGS) $< obj/*.o -o $@

obj/tree_solver.o : src/tree_solver.cpp include/tree_solver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean :
	rm -f bin/*
	rm -f obj/*
