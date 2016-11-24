CCOMP = g++
CFLAGS = -Wall

all: buildTree
buildTree: buildTree.cpp treeFinder.h ../common/io.h ../common/graph.h
	$(CCOMP) $(CFLAGS) buildTree.cpp -o buildTree

clean:
	rm -f buildTree
