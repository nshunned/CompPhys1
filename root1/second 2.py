# second.py
# Here we specify a named function for ROOT.
# The function is effectively added to the ROOT framework
# Usage in ROOT: .L second.C
#                randomHist() or randomHist(int)
#Note: you must reload the function after making changes!

import ROOT as r

def randomHist(entries=1000):
    my_hist = r.TH1F("my_hist", "My First Histogram", 100, -5, 5)
    my_hist.FillRandom("gaus", entries)
    return my_hist