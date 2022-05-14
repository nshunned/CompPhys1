import ROOT as r
import sys

fa1 = r.TF1("fa1","sin(x)/x*exp(x/4)",1e-6,10);  # a very basic function
# now make a function with a parameter
fa1p = r.TF1("fa1p","[1]*sin(x)/x*exp(x/[0])",1e-6,10);  # add 2 parameters
fa1p.SetParameter(0,5);
fa1p.SetParameter(1,0.8);
fa1p.SetLineColor(r.kBlue);
fa1.Draw();
fa1p.Draw("same");

print("Hit return to exit")
sys.stdin.readline()

