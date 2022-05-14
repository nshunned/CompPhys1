CPPFLAGS=-O3 -Wall
#CPPFLAGS=-g -Wall

all: compareinteg

compareinteg: compareinteg.cpp integrators.o
	g++ $(CPPFLAGS) -o compareinteg compareinteg.cpp integrators.o

integrators.o: integrators.hpp integrators.cpp
	g++ $(CPPFLAGS) -c integrators.cpp

clean: 
	rm -f *.o compareinteg *.png *.ps *.pdf *.dat *~ *.d *.pcm

