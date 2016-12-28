CXX = g++ --std=c++11
CXXFLAGS += -I./common
CXXFLAGS += -I./solvers
CXXFLAGS += -Wall
# CXXFLAGS += -O3
# CXXFLAGS += -g


all : bin/test

bin/test : common/matrix.h solvers/PCGSolver.h
bin/test : solvers/test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	rm -f bin/*
