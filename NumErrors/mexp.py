# Here we will use the functions in our mexp.cpp file inside of python!
# **version w/o ROOT libraries**

from __future__ import print_function  # ensure same behaviour on python2(3)

from ctypes import cdll  # pythonic interfaces for C(++) libraries
from ctypes import c_double
from math import exp
import sys
import argparse

parser = argparse.ArgumentParser(description='mexp program checker')
parser.add_argument('-p', '--doplot', default=None, action="store_true", help='Make a plot')
parser.add_argument('-x', '--xtest', type=float, default=0.1, help='Value of x for exp(-x)')
parser.add_argument('-m', '--xmin', type=float, default=0.001, help='minimum x for scan')
parser.add_argument('-M', '--xmax', type=float, default=100, help='maximum x for scan')
args = parser.parse_args()


# stuff for making plots
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 18})  # make fonts readable!

# load our library and define the input/output interfaces
mylib=cdll.LoadLibrary('./mexp.soc')
# we have to define input/ouput data types for EVERY function we want to use
# in a particular library
mylib.mexp.argtypes=[c_double]  # input data type[s]
mylib.mexp.restype=c_double     # return value type

x=args.xtest
val=mylib.mexp(x)

print("exp(",x,")=",val,"rel_error=",abs(val-exp(-x))/exp(-x))

# plot results for a range of x values - note plot will truncate if results are NAN 
if not args.doplot: sys.exit()

x=args.xmin
delta=1.01
xary=[]
yary=[]
while x<args.xmax:  # loop over x-range and test out mexp fcn
    val=mylib.mexp(x)
    err=abs(val-exp(-x))/exp(-x)
    if np.isfinite(err):  # truncate for large errors
        xary.append(x)
        yary.append(err)
    x=x*delta
    
# convert lists is np arrays
xary=np.array(xary)
yary=np.array(yary)


plt.plot(xary,yary)
plt.xlabel('x')
plt.yscale('log')
plt.ylabel('Rel. error')
plt.title('Relative error in exp(-x) calculation')
print("Close plot window to exit")
plt.show()





