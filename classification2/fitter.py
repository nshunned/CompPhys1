import ROOT as r
import sys

# cuts
tfSim=r.TFile("simulation_sel.root")
tfDat=r.TFile("dataSmSig_sel.root")
hBkg=tfSim.Get("hb_var3")
hSig=tfSim.Get("hs_var3")
hDat=tfDat.Get("hd_var3")

# loglikelihood
# tFile=r.TFile("logLike.root")
# hSig=tFile.Get("sll_2")
# hBkg=tFile.Get("bll_2")
# hDat=tFile.Get("dll_2")

hDat.Print()

# first plot the data and S/B distributions, all normalized to same # of events
tc1=r.TCanvas("tc1")
hd=r.TH1F(hDat) 
hs=r.TH1F(hSig) ; hs.SetLineColor(r.kRed)
hb=r.TH1F(hBkg) ; hb.SetLineColor(r.kBlue)
hs.Scale(1.0/hs.Integral() * hDat.Integral())
hb.Scale(1.0/hb.Integral() * hDat.Integral())
hst=r.THStack("hs","")
hst.SetTitle("S,B dists normalized to data")
hst.Add(hd)
hst.Add(hb)
hst.Add(hs)
hst.Draw("nostack")
#hd.Draw()
#hb.Draw("same,hist")
#hs.Draw("same,hist")



r.gROOT.ProcessLine(".L fitter.C")
r.fitter(hBkg,hSig,hDat)




print("Hit return to exit")
sys.stdin.readline()
