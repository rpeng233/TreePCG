CCOMP = g++
CFLAGS = -Wall

all: genCayley graphAsc2MM graphBin2MM
genCayley: genCayley.cpp
	$(CCOMP) $(CFLAGS) genCayley.cpp -o genCayley
graphAsc2MM: graphAsc2MM.cpp
	$(CCOMP) $(CFLAGS) graphAsc2MM.cpp -o graphAsc2MM
graphBin2MM: graphBin2MM.cpp
	$(CCOMP) $(CFLAGS) graphBin2MM.cpp -o graphBin2MM

clean:
	rm -f genCayley graph2MM graphAsc2MM graphBin2MM
