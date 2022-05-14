#include <TLatex.h>

TH1F* hdata;
TH1F* syst1;
TH1F* syst2;
TF1* fparam;

void setCanvas(TCanvas* canvas) {
  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
}

void problem1() {
    // read from file
    TFile* file = new TFile("mydata.root");
    hdata = new TH1F(*(TH1F*)file->Get("Data"));
    syst1 = new TH1F(*(TH1F*)file->Get("syst1"));
    syst2 = new TH1F(*(TH1F*)file->Get("syst2"));

    // fit original data with linear func & retrieve fit params
    hdata->Fit("pol1","Q0");
    TF1 *fitFunc = hdata->GetFunction("pol1");
    double par0Stat = fitFunc->GetParameter(0);
    double par1Stat = fitFunc->GetParameter(1);
    double par0StatErr = fitFunc->GetParError(0);
    double par1StatErr = fitFunc->GetParError(1);

    // syst1 up
    TH1F* data1Up = new TH1F(*hdata);
    data1Up->Add(syst1);
    data1Up->Fit("pol1","Q0");
    fitFunc = data1Up->GetFunction("pol1");
    double par0Syst1Up = fitFunc->GetParameter(0);
    double par1Syst1Up = fitFunc->GetParameter(1);
    
    // syst2up
    TH1F* data2Up = new TH1F(*hdata);
    data2Up->Add(syst2);
    data2Up->Fit("pol1","Q0");
    fitFunc = data2Up->GetFunction("pol1");
    double par0Syst2Up = fitFunc->GetParameter(0);
    double par1Syst2Up = fitFunc->GetParameter(1);

    // up errors
    double par0SystUpErr = TMath::Sqrt((par0Stat-par0Syst1Up)*(par0Stat-par0Syst1Up) + (par0Stat-par0Syst2Up)*(par0Stat-par0Syst2Up));
    double par1SystUpErr = TMath::Sqrt((par1Stat-par1Syst1Up)*(par1Stat-par1Syst1Up) + (par1Stat-par1Syst2Up)*(par1Stat-par1Syst2Up));

    // syst1 down
    TH1F* data1Down = new TH1F(*hdata);
    data1Down->Add(syst1,-1);
    data1Down->Fit("pol1","Q0");
    fitFunc = data1Down->GetFunction("pol1");
    double par0Syst1Down = fitFunc->GetParameter(0);
    double par1Syst1Down = fitFunc->GetParameter(1);

    // syst2 down
    TH1F* data2Down = new TH1F(*hdata);
    data2Down->Add(syst2,-1);
    data2Down->Fit("pol1","Q0");
    fitFunc = data2Down->GetFunction("pol1");
    double par0Syst2Down = fitFunc->GetParameter(0);
    double par1Syst2Down = fitFunc->GetParameter(1);

    // down errors
    double par0SystDownErr = TMath::Sqrt((par0Stat-par0Syst1Down)*(par0Stat-par0Syst1Down) + (par0Stat-par0Syst2Down)*(par0Stat-par0Syst2Down));
    double par1SystDownErr = TMath::Sqrt((par1Stat-par1Syst1Down)*(par1Stat-par1Syst1Down) + (par1Stat-par1Syst2Down)*(par1Stat-par1Syst2Down));

    cout << "fit parameters" << endl;
    cout << "p0 (y-intercept) = " << par0Stat << "; errors: stat = " << par0StatErr << ", syst up = " << par0SystUpErr << ", syst down = " << par0SystDownErr << endl;
    cout << "p1 (slope) = " << par1Stat << "; errors: stat = " << par1StatErr << ", syst up = " << par1SystUpErr << ", syst down = " << par1SystDownErr << endl;

    file->Close();
}

double pdf(double* xPtr, double par[]) {
    double x = *xPtr;
    double p0 = par[0];
    double p1 = par[1];
    double s1 = par[2];
    double s2 = par[3];
    double delta1 = syst1->GetBinContent(syst1->GetBin(x));
    double delta2 = syst2->GetBinContent(syst2->GetBin(x));
    return p0 + p1*x + s1*delta1 + s2+delta2;
}

double chi2(TH1F* h, TF1* f) {
    double chi2=0;
    for (int i=1; i<h->GetNbinsX(); i++) {
        int y = (int)(h->GetBinContent(i));
        double ey = h->GetBinError(i);
        if (ey == 0)
            continue;
        double yFit = f->Eval(h->GetBinCenter(i));
        chi2 += (y-yFit)/ey * (y-yFit)/ey;
    }
    double s1 = f->GetParameter(2);
    double s2 = f->GetParameter(3);
    return chi2 + s1*s1 + s2*s2;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag) {
    for (int i=0; i<npar; i++) {
        fparam->SetParameter(i,par[i]);
    }
    f = chi2(hdata,fparam);
}

void problem2() {
    // setup canvas
    TCanvas* canvas = new TCanvas();
    setCanvas(canvas);

    // read from file
    TFile* file = new TFile("mydata.root");
    hdata = new TH1F(*(TH1F*)file->Get("Data"));
    syst1 = new TH1F(*(TH1F*)file->Get("syst1"));
    syst2 = new TH1F(*(TH1F*)file->Get("syst2"));
    double xmin = hdata->GetBinCenter(hdata->FindFirstBinAbove()-1);
    double xmax = hdata->GetBinCenter(hdata->FindLastBinAbove()+1);

    const int npar = 4;
    TMinuit minuit(npar);
    minuit.SetFCN(fcn);

    // function to fit the data
    TF1* myfunc = new TF1("myfunc", pdf, xmin, xmax, npar);
    fparam = myfunc;

    // set up arrays for minimizer 1
    TString parName[npar];
    double par[npar];
    double stepSize[npar];
    double minVal[npar];
    double maxVal[npar];

    parName[0] = "y-intercept";
    parName[1] = "slope";
    parName[2] = "s1";
    parName[3] = "s2";
    par[0] = 66;
    par[1] = 26;
    par[2] = 1;
    par[3] = -1;
    for (int i=0; i<npar; i++) {
        stepSize[i] = TMath::Abs(par[i]*0.01);
        minVal[i] = 0;
        maxVal[i] = 0;
    }

    // initialize parameters
    for (int i=0; i<npar; i++) {
        minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
    }

    // minimize
    minuit.Migrad();
    double outpar[npar], err[npar];
    for (int i=0; i<npar; i++) {
        minuit.GetParameter(i,outpar[i],err[i]);
    }

    // plot the result
    hdata->SetTitle(";x;f(p,s)");
    hdata->Draw("e");
    myfunc->SetParameters(outpar);
    myfunc->Draw("same");

    myfunc->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
    myfunc->SetLineColor(1);             //  black (default)
    myfunc->SetLineWidth(1);
}

void makeFits() {
    problem1();
    problem2();
}