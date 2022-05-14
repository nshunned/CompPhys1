ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) 
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3

all: Problem_1 Problem_2

Problem_1:  Problem_1.cpp RKn.o
	$(GXX) -o Problem_1 Problem_1.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS) RKn.o

Problem_2:  Problem_2.cpp RKn.o
	$(GXX) -o Problem_2 Problem_2.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS) RKn.o

RKn.o: RKn.cpp RKn.hpp
	$(GXX) -c -o RKn.o RKn.cpp $(ROOTCFLAGS) 

clean:
	rm -f Problem_1 Problem_2 RKn.o *~
