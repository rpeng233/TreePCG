CXX = g++ --std=c++11
CXXFLAGS += -I./common
CXXFLAGS += -Wall -pedantic
# CXXFLAGS += -O3
# CXXFLAGS += -g
CXXFLAGS += -lmpfr


all : bin/test

bin/test : common/*.h
bin/test : src/test.cpp obj/TreeSolver.o
	$(CXX) $(CXXFLAGS) $< obj/*.o -o $@

obj/TreeSolver.o : src/TreeSolver.cpp common/TreeSolver.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean :
	rm -f bin/*
	rm -f obj/*
