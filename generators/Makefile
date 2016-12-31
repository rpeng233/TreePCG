CCOMP = g++
CFLAGS = -Wall

all: gen_cayley graph_to_matrix graph_to_graph
gen_cayley: gen_cayley.cpp ../include/io.h ../include/graph.h
	$(CCOMP) $(CFLAGS) gen_cayley.cpp -o gen_cayley
graph_to_graph: graph_to_graph.cpp ../include/io.h ../include/graph.h
	$(CCOMP) $(CFLAGS) graph_to_graph.cpp -o graph_to_graph
graph_to_matrix: graph_to_matrix.cpp ../include/io.h ../include/graph.h ../include/matrix.h
	$(CCOMP) $(CFLAGS) graph_to_matrix.cpp -o graph_to_matrix

clean:
	rm -f gen_cayley graph_to_matrix graph_to_graph
