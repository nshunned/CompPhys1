#include <iostream>

using namespace std;

double scan(TH1F *s, TH1F *b){
    int nbins=s->GetNbinsX();
    double sig=0;
    int bestcut=0;
    for (int i=1; i<nbins; i++){
        double event_s=s->Integral(i,nbins);
        double event_b=b->Integral(i,nbins);
        double sig_i=event_s/sqrt(event_s+event_b);
        if (sig_i>sig) {
            sig=sig_i;
            bestcut=i;
        }
    }
    cout << "Max significance cut for high cut on ";
    cout << s->GetName() << " = " << sig << " at bin " << bestcut << " : " << s->GetBinLowEdge(bestcut+1) << endl;
    return sig;
}

void scanner() {
    // 1. retrieve histograms from file
    TFile *file=new TFile("logLike.root");
    TH1F *hlls[3], *hllr[3], *hll2[3];
    TString lls[3]={"sll_s","bll_s","dll_s"};
    TString llr[3]={"sll_r","bll_r","dll_r"};
    TString ll2[3]={"sll_2","bll_2","dll_2"};
    for (int i=0; i<3; i++) {
        hlls[i]=(TH1F*)file->Get(lls[i]);
        hllr[i]=(TH1F*)file->Get(llr[i]);
        hll2[i]=(TH1F*)file->Get(ll2[i]);
    }

    // 2. scan and calculate expected significance
    scan(hlls[0],hlls[1]);
    scan(hllr[0],hllr[1]);
    scan(hll2[0],hll2[1]);
}