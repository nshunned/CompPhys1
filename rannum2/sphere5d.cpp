#include "Math/QuasiRandom.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"

#include <iostream>

using namespace TMath;
using namespace std;

// dimension
const int DIM = 5;
// volume of unit n-ball and the container n-cube
const double VOL_BALL = Power(Pi(), DIM/2.0)/Gamma(DIM/2.0+1);
const double VOL_CUBE = Power(2,5);

double distance(double x1, double x2, double x3, double x4, double x5) {
    return Sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5);
}

double fixedGrid(int total, bool rand) {
    // points for a single "edge" of the 5D grid
    int n = (int)Power(total, 1.0/5); // number of grid points
    double h = 2.0/n;                 // separation between points on a single "edge"
    double x[n];                      // coordinate value for each point on a single "edge"
    for (int i=0; i<n; i++) x[i] = i*h + h/2 - 1;
    
    // number of points inside the 5-ball
    int count=0;

    // lower-discrepency method
    // coordinate of each point is randomly shifted within
    //     twice the separation of the fixed grid
    if (rand) {
        TRandom3 generator(0);
        for (int i1=0; i1<n; i1++)
        for (int i2=0; i2<n; i2++)
        for (int i3=0; i3<n; i3++)
        for (int i4=0; i4<n; i4++)
        for (int i5=0; i5<n; i5++) {
            double r1 = (generator.Rndm()-0.5)*2*h;
            double r2 = (generator.Rndm()-0.5)*2*h;
            double r3 = (generator.Rndm()-0.5)*2*h;
            double r4 = (generator.Rndm()-0.5)*2*h;
            double r5 = (generator.Rndm()-0.5)*2*h;
            double d = distance(x[i1]+r1, x[i2]+r2, x[i3]+r3, x[i4]+r4, x[i5]+r5);
            if (d <= 1)
                count++;
        }
    }
    // uniform fixed grid
    else {
        for (int i1=0; i1<n; i1++)
        for (int i2=0; i2<n; i2++)
        for (int i3=0; i3<n; i3++)
        for (int i4=0; i4<n; i4++)
        for (int i5=0; i5<n; i5++) {
            double d = distance(x[i1], x[i2], x[i3], x[i4], x[i5]);
            if (d <= 1)
                count++;
        }
    }

    return VOL_CUBE*((double)count/total);
}

double random3(int total) {
    // number of points inside the 5-ball
    int count=0;

    // randomly put a point in the 5-cube
    TRandom3 generator(0);
    for(int i=0; i<total; i++) {
        double x[5];
        for (int i=0; i<5; i++) x[i] = (generator.Rndm()-0.5)*2;

        double d = distance(x[0], x[1], x[2], x[3], x[4]);
        if (d <= 1) count++;

        //cout << d << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << endl;
    }
    return VOL_CUBE*((double)count/total);
}

double sobol(int total) {
    using namespace ROOT::Math;

    // randomly put a point in the 5-cube and count the number of points inside the 5-ball
    int count=0;
    QuasiRandomSobol generator(5);
    for(int i=0; i<total; i++) {
        double x[5];
        generator.Next(x);
        double d = distance(x[0], x[1], x[2], x[3], x[4]);
        //cout << d << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << endl;
        if (d <= 1) count++;
    }
    // cout << "count: " << count << endl;
    return VOL_CUBE*((double)count/total);
}

void sphere5d() {
    // double N = Power(30,5);
    // cout << fixedGrid(N, 0) << endl;
    // cout << fixedGrid(N, 1) << endl;
    // cout << random3(N) << endl;
    // cout << sobol(N) << endl;

    // list of total points
    double totals[] = {1e1, 1e2, 1e3, 1e4, 1e5, 1.5e5, 1e6, 1.5e6, 1e7, 1.5e7, 1e8};
    int nTotals = sizeof(totals)/sizeof(totals[0]);

    // list of errors
    double eFixed[nTotals];
    double eLow[nTotals];
    double eRand3[nTotals];
    double eSobol[nTotals];

    // calculate error for each method
    for (int i=0; i<nTotals; i++) {
        eFixed[i] = Abs(fixedGrid(totals[i],0) - VOL_BALL);
        eLow[i] = Abs(fixedGrid(totals[i],1) - VOL_BALL);
        eRand3[i] = Abs(random3(totals[i]) - VOL_BALL);
        eSobol[i] = Abs(sobol(totals[i]) - VOL_BALL);
    }

    // set up canvas
    UInt_t dw = gClient->GetDisplayWidth();
    UInt_t dh = gClient->GetDisplayHeight();
    TCanvas *canvas = new TCanvas("canvas","Error vs. N",dw/2,dw/3*2);
    gPad->SetLogx();
    gPad->SetLogy();

    // TGraphs for each method
    TGraph *fixed = new TGraph(nTotals, totals, eFixed);
    TGraph *low = new TGraph(nTotals, totals, eLow);
    TGraph *rand3 = new TGraph(nTotals, totals, eRand3);
    TGraph *sobol = new TGraph(nTotals, totals, eSobol);

    // make them look nice
    fixed->SetNameTitle("Fixed Grid", ";Number of Points Used;Error");
    low->SetNameTitle("Low Discrepency", ";Number of Points Used;Error");
    rand3->SetNameTitle("TRandom3", ";Number of Points Used;Error");
    sobol->SetNameTitle("Sobol", ";Number of Points Used;Error");
    // line color
    fixed->SetLineColor(2);
    low->SetLineColor(3);
    rand3->SetLineColor(4);
    sobol->SetLineColor(5);
    // line width
    fixed->SetLineWidth(2);
    low->SetLineWidth(2);
    rand3->SetLineWidth(2);
    sobol->SetLineWidth(2);
    // marker style
    fixed->SetMarkerStyle(5);
    low->SetMarkerStyle(5);
    rand3->SetMarkerStyle(5);
    sobol->SetMarkerStyle(5);
    // axes
    fixed->GetXaxis()->SetTitleOffset(1.4);
    fixed->GetYaxis()->SetTitleOffset(1.4);
    fixed->GetYaxis()->SetRangeUser(1e-7,100);

    // draw them
    fixed->Draw("");
    low->Draw("SAMELP");
    rand3->Draw("SAMELP");
    sobol->Draw("SAMELP");

    // legend
    TLegend* legend = new TLegend(0.5, 0.32, 0.16, 0.15, "");
    legend->SetTextSize(0.035);
    legend->AddEntry(fixed, "Fixed Grid");
    legend->AddEntry(low, "Low Discrepency");
    legend->AddEntry(rand3, "TRandom3");
    legend->AddEntry(sobol, "Sobol");
    legend->Draw();

    canvas->Print("sphere5d.pdf","pdf");
}