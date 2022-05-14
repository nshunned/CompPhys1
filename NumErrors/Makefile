# We build two versions of the shared library
# Nothing special is required to load the C++ libraries into python using ROOT
# If we use ctypes and load the functions by hand in python, they must be
# linked using the extern "C" setting.  Otherwise python will not be able to
# recognize the names of the functions (c.f. name mangling)


all: mexp mexpso mexpsoc

mexp: mexp.cpp mexp.h
	g++ -O -Wall -omexp mexp.cpp    # build executable
mexpso: mexp.cpp mexp.h
	g++ -O -Wall -fPIC -shared -omexp.so mexp.cpp # build shared library
mexpsoc: mexp.cpp mexp.h  # build "so" w/ C style linkage
	g++ -O -Wall -fPIC -DEXTERNC -shared -omexp.soc mexp.cpp


clean:
	rm -f mexp mexp.so mexp.soc *~
