// fit_user.C
// example of fitting with a user defined function

void setCanvas() {
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetTitleSize(0.075,"t");
  gPad->SetLeftMargin(0.125);
  gPad->SetBottomMargin(0.125);
}

// A function producing two peaks on top of an exponentialy falling background
// Depends on several parameters
// Note: this need not be a 1D function
// Generic interface for fcn of n input-values and m parameters
// Functions with this interface may be used to construct a "TFunction"
Double_t myfunction(Double_t *xin, Double_t *par) {
     Double_t x=xin[0];

     Double_t bkgScale=par[0];

     Double_t alpha=par[1];
     Double_t beta=par[2];
     Double_t background = pow(x/beta,-1.0*alpha);

     Double_t A1=par[3];
     Double_t mu1=par[4];
     Double_t sig1=par[5];
     Double_t peak1=A1*TMath::Exp(-0.5*(x-mu1)*(x-mu1)/sig1/sig1);

     Double_t A2=par[6];
     Double_t mu2=par[7];
     Double_t sig2=par[8];
     Double_t peak2=A2*exp(-0.5*(x-mu2)*(x-mu2)/sig2/sig2);

     return bkgScale*background+peak1+peak2;
}

Double_t ning(Double_t *xin, Double_t *par) {
     Double_t x=xin[0];
     return par[0]*x*x*x*x + par[1]*x*x*x + par[2]*x*x + par[3]*x + par[4]*TMath::Sin(par[5]*x+par[6]) + par[7];
}


void fit_user(int entries=100000) {
     TFile *f = new TFile("data1.root","recreate");

     ////////////////////////////////////
     //            PART I
     ////////////////////////////////////
     
     TCanvas *canvas1 = new TCanvas("canvas1", "Problem 1", 700, 500);
     setCanvas();

     // create myfunction & set parameters
     int nParam1=9;
     TF1 *f1 = new TF1("f1",myfunction,300,1000,nParam1);
     f1->SetParameters(1e9,4.7,40,5000,500,2,1200,800,25);

     // fill a histogram with random data using f1 as a pdf
     TH1F *ranHist1 = new TH1F("ranHist1", "Random Histogram",500,300,1000);
     ranHist1->FillRandom("f1",entries);
     ranHist1->Draw("e");

     // all fits begin with initial guesses at the best parameter values
     f1->SetParameters(0.25e6,4,80,150,500,2,30,800,20);
     f1->Draw("same");
     
     // perform the fit here
     ranHist1->Fit("f1","EQ");

     // retrieve fitted f1
     TF1 *f1Fitted = ranHist1->GetFunction("f1");
     f1Fitted->Print();


     // retrieve central parameter values and errors
     cout << "fit parameters" << endl;
     for (int i=0; i<nParam1; i++) {
          cout << i << "\t" << f1Fitted->GetParameter(i) << " +- " << f1Fitted->GetParError(i) << endl;
     }

     // get chi2, pvalue
     Double_t chi2_1 = f1Fitted->GetChisquare();  
     cout << "chi2 / ndf: " << chi2_1 << " / " << f1Fitted->GetNDF() << " reduced chi2: " << chi2_1/f1Fitted->GetNDF() << endl;
     cout << " P-value: " << TMath::Prob(chi2_1,f1Fitted->GetNDF()) << endl;

     // print canvas
     canvas1->Print("fit1.pdf","pdf");

     ////////////////////////////////////
     //            PART II
     ////////////////////////////////////

     TCanvas *canvas2 = new TCanvas("canvas2", "Problem 2", 700, 500);
     setCanvas();

     // read & draw histogram
     TFile *file = new TFile("datadist.root");
     TH1F *ranHist2 = new TH1F(*(TH1F*)file->Get("h"));
     ranHist2->Draw("e");

     // create function & try to guess fit
     int nParam2=8;
     TF1 *f2 = new TF1("f2",ning,0,12,nParam2);
     f2->SetParameters(-0.0205613, 0.350682, -2.62072, 17.8655,
                       4.8857, -1.84536, 4.46804,
                       7.5547);
     f2->Draw("same");
     
     // perform the fit here
     ranHist2->Fit("f2","EQ");

     // retrieve fitted f2
     TF1 *f2Fitted = ranHist2->GetFunction("f2");
     f2Fitted->Print();

     // retrieve central parameter values and errors
     cout << "fit parameters" << endl;
     for (int i=0; i<nParam2; i++) {
          cout << i << "\t" << f2Fitted->GetParameter(i) << " +- " << f2Fitted->GetParError(i) << endl;
     }

     // get chi2, pvalue
     Double_t chi2_2 = f2Fitted->GetChisquare(); 
     cout << "chi2 / ndf: " << chi2_2 << " / " << f2Fitted->GetNDF() << " reduced chi2: " << chi2_1/f2Fitted->GetNDF() << endl;
     cout << " P-value: " << TMath::Prob(chi2_2,f2Fitted->GetNDF()) << endl;

     // print canvas
     canvas2->Print("fit2.pdf","pdf");

     f->Write();
}

