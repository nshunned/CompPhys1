# here we access the root configuration, include files, and libraries
ROOTCFLAGS=$(shell root-config --cflags)
ROOTINC=$(shell root-config --incdir)
ROOTLIBDIR=$(shell root-config --libdir)
ROOTLIBS=$(shell root-config --libs) -lMinuit
ROOTLDFLAGS=$(shell root-config --ldflags)

ROOTC=$(ROOTCFLAGS) 
#-I$(ROOTINC)
ROOTLINK=-L$(ROOTLIBDIR) $(ROOTLIBS) $(ROOTLDFLAGS)




a11: myrand drand48 badrandom drand48_test badrandom_test entropy generators

myrand_plot: myrand.cpp
	g++ -O -Wall -omyrand myrand.cpp

drand48_plot: drand48.cpp
	g++ -O -Wall -odrand48 drand48.cpp

badrandom_plot: badrandom.cpp
	g++ -O -Wall -obadrandom badrandom.cpp

drand48_test: drand48_test.cpp
	g++ -O -Wall -odrand48_test drand48_test.cpp

badrandom_test: badrandom_test.cpp
	g++ -O -Wall -obadrandom_test badrandom_test.cpp

entropy: entropy.cpp
	g++ -O -Wall -oentropy entropy.cpp

generators: generators.cpp
	g++ -O -Wall $(ROOTC) -o generators generators.cpp $(ROOTLINK)


clean:
	rm -f drand48 badrandom drand48_test entropy 
	rm -f badrandom_test generators *.dat *.bin *.gpg *.so *.d *~
