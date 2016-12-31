CCOMP = g++
CFLAGS = -Wall

all: buildTree stretchStatistics

buildTree: buildTree.cpp treeFinder.h ../common/io.h ../common/graph.h
	$(CCOMP) $(CFLAGS) buildTree.cpp -o buildTree

stretchStatistics: stretchStatistics.cpp treeFinder.h ../common/io.h ../common/graph.h
	$(CCOMP) $(CFLAGS) stretchStatistics.cpp -o stretchStatistics

clean:
	rm -f buildTree stretchStatistics
