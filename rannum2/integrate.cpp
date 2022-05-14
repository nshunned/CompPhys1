
#include "TRandom3.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <cstdio>
#include <iostream>

using namespace TMath;
using namespace std;

double const X_MIN=1;
double const X_MAX=4;
double const X_RANGE = X_MAX-X_MIN;
double const A=0.285369684759;

// Monte-Carlo integration of the following function
double func(double x) {
    return Exp(-Abs(x*Log(2*x*x)-15/x));
}

double darting(int total) {
    // let's get darting over a 3x1 rectangle
    TRandom3 generator(0);
    int count=0;
    for (int i=0; i<total; i++) {
        double x = generator.Rndm()*X_RANGE + X_MIN;
        double y = generator.Rndm();
        if(y <= func(x)) count++;
    }
    return X_RANGE*count/total;
}

double rectangles(int total) {
    // freopen("output.txt","w",stdout);
    
    // optimized container and its volume
    double points[12][2] = {
        {1.00, 0.05}, {2.0, 0.05}, {2.10, 0.10}, {2.2, 0.2}, {2.25, 0.3}, {2.39, 0.7},
        {2.45, 1.00}, {2.5, 1.00}, {2.55, 0.75}, {2.7, 0.6}, {2.80, 0.2}, {4.00, 0.1}
    };
    double vol=0;
    for(int i=0; i<11; i++) {
        vol += (points[i+1][0]-points[i][0]) * points[i+1][1];
        // cout << vol << endl;
    }
    
    // throw and count darts
    TRandom3 generator(0);
    int count=0;
    for (int in=0; in<total;) {
        double x = generator.Rndm()*X_RANGE + X_MIN;
        double y = generator.Rndm();
        for (int i=0; i<11; i++) {
            if (x >= points[i][0] && x <= points[i+1][0] && y <= points[i+1][1]) {
                in++;
                if(y <= func(x)) count++;
                break;
            }
        }
        // cout << x << " " << y << " " << func(x) << endl;
    }
    // fclose (stdout);
    return vol*count/total;
}

void integrate() {
    // list of total points
    double totals[] = {1e1, 1e2, 1e3, 1e4, 1e5, 1.5e5, 1e6, 1.5e6, 1e7, 1.5e7, 1e8};
    int nTotals = sizeof(totals)/sizeof(totals[0]);

    // list of errors
    double eDarting[nTotals];
    double eRectangles[nTotals];

    // calculate error for each method
    for (int i=0; i<nTotals; i++) {
        eDarting[i] = Abs(darting(totals[i]) - A)/A;
        eRectangles[i] = Abs(rectangles(totals[i]) - A)/A;
    }

    // set up canvas
    UInt_t dw = gClient->GetDisplayWidth();
    UInt_t dh = gClient->GetDisplayHeight();
    TCanvas *canvas = new TCanvas("canvas","Relative Errors vs. N",dw/2,dw/3*2);
    gPad->SetLogx();
    gPad->SetLogy();

    // TGraphs for each method
    TGraph *gDarting = new TGraph(nTotals, totals, eDarting);
    TGraph *gRectangles = new TGraph(nTotals, totals, eRectangles);

    // make them look nice
    gDarting->SetNameTitle("Darting", ";Number of Points Used;Relative Error");
    gRectangles->SetNameTitle("Rectangles", ";Number of Points Used;Relative Error");
    // line color
    gDarting->SetLineColor(2);
    gRectangles->SetLineColor(3);
    // line width
    gDarting->SetLineWidth(2);
    gRectangles->SetLineWidth(2);
    // marker style
    gDarting->SetMarkerStyle(5);
    gRectangles->SetMarkerStyle(5);
    // axes
    gDarting->GetXaxis()->SetTitleOffset(1.4);
    gDarting->GetYaxis()->SetTitleOffset(1.4);
    gDarting->GetYaxis()->SetRangeUser(1e-7,100);

    // draw them
    gDarting->Draw("");
    gRectangles->Draw("SAMELP");
    
    // legend
    TLegend* legend = new TLegend(0.5, 0.32, 0.16, 0.15, "");
    legend->SetTextSize(0.035);
    legend->AddEntry(gDarting, "Darting");
    legend->AddEntry(gRectangles, "Rectangles");
    legend->Draw();

    canvas->Print("random.pdf","pdf");
} 