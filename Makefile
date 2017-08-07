#Makefile for findAlmost

#change these to the locations of PEBBL and nauty respectively
PEBBL_DIR = /opt/acro-pebbl
NAUTY_DIR = /opt/nauty25r9


MPICPP       = mpic++ -std=c++0x -O3 
NAUTYGRAPH   = NautyGraph.cpp
PEBBL   = -I$(PEBBL_DIR)/include -L$(PEBBL_DIR)/lib -lpebbl -lutilib
NAUTY   = -I$(NAUTY_DIR) $(NAUTY_DIR)/nauty.a

findAlmost : findAlmost.cpp NautyGraph.cpp
	$(MPICPP) findAlmost.cpp $(NAUTYGRAPH) $(NAUTY) $(PEBBL) -o findAlmost

clean :
	rm -f findAlmost
