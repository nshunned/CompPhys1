# first.py
# pythonic version

import ROOT as r
import sys

my_hist = r.TH1F("my_hist", "My First Histogram", 100, -5, 5);
my_hist.FillRandom("gaus", 10000);
my_hist.Draw();


#alternatively, we can make our TCanvas and subdivide it if we want
tc=r.TCanvas("tc","My Canvas")
tc.Divide(2,2)   # columns, rows  2x3 TPads in this case
tc.cd(1)         # TPads are numbered gogin left to right, top to bottom
my_hist.Draw()
tc.cd(2)
my_hist.Draw("E")           # draw with error bars
tc.cd(3)
my_hist.Draw("E")           # draw with errorbars
my_hist.Draw("L,same")      # overlay with lines through the points
tc.cd(4)
my_hist.Draw("HBAR")        # draw as horixontal bar chart
tc.Update()


print("Hit return to exit")
sys.stdin.readline()

