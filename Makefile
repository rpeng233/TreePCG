CXX = g++ --std=c++11
CXXFLAGS += -I./common
CXXFLAGS += -Wall
# CXXFLAGS += -O3
# CXXFLAGS += -g


all : bin/test

bin/test : common/graph.h common/common.h common/matrix.h
bin/test : src/test.cpp obj/TreeSolver.o
	$(CXX) $(CXXFLAGS) $^ -o $@

obj/TreeSolver.o : src/TreeSolver.cpp common/TreeSolver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean :
	rm -f bin/*
	rm -f obj/*
