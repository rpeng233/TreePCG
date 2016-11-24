CCOMP = g++
CFLAGS = -Wall

all: genCayley graphToMatrix graphToGraph
genCayley: genCayley.cpp ../common/io.h ../common/graph.h
	$(CCOMP) $(CFLAGS) genCayley.cpp -o genCayley
graphToGraph: graphToGraph.cpp ../common/io.h ../common/graph.h
	$(CCOMP) $(CFLAGS) graphToGraph.cpp -o graphToGraph
graphToMatrix: graphToMatrix.cpp ../common/io.h ../common/graph.h ../common/matrix.h
	$(CCOMP) $(CFLAGS) graphToMatrix.cpp -o graphToMatrix

clean:
	rm -f genCayley graphToMatrix graphToGraph
