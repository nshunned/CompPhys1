all: onespin twospin ising1d ising2d ising2d_vs_T ising1d_vs_T
onespin: onespin.cpp
	g++ -Wall -O -o onespin onespin.cpp
ising1d: ising1d.cpp
	g++ -Wall -O -o ising1d ising1d.cpp
ising2d: ising2d.cpp
	g++ -Wall -O -o ising2d ising2d.cpp
ising2d_vs_T: ising2d_vs_T.cpp
	g++ -Wall -O -o ising2d_vs_T ising2d_vs_T.cpp
ising1d_vs_T: ising1d_vs_T.cpp
	g++ -Wall -O -o ising1d_vs_T ising1d_vs_T.cpp
twospin: twospin.cpp
	g++ -Wall -o twospin twospin.cpp
clean:
	rm -f onespin twospin ising1d ising2d 
	rm -f ising2d_vs_T ising1d_vs_T
	rm -f a.out *.dat *~
