void plotMC(){
  gStyle->SetOptStat(0);
  TFile *in=new TFile("simulation.root");
  TTree *s=(TTree*)in->Get("TreeS");
  TTree *b=(TTree*)in->Get("TreeB");
  TFile *in2=new TFile("dataSmSig.root");
  TTree *d=(TTree*)in2->Get("TreeD");
  TFile *in3=new TFile("dataBigSig.root");
  TTree *d2=(TTree*)in3->Get("TreeD");

  b->SetLineColor(kBlue);
  b->SetMarkerColor(kBlue);
  s->SetLineColor(kRed);
  s->SetMarkerColor(kRed);
  d->SetLineColor(kBlack);
  d2->SetLineColor(kBlack);

  //  d->SetLineWidth(2);

  //  b->SetFillColor(2);
  TFile *out=new TFile("sbDists.root","recreate");
  TCanvas *s1 = new TCanvas("s1","s1",1200,300);
  s1->Divide(5,1);
  s1->cd(1); s->Draw("var1>>var1s");
  s1->cd(2); s->Draw("var2>>var2s"); 
  s1->cd(3); s->Draw("var3>>var3s",""); 
  s1->cd(4); s->Draw("var4>>var4s"); 
  s1->cd(5); s->Draw("var5>>var5s",""); 
  s1->SaveAs("s1.png");

  TCanvas *b1 = new TCanvas("b1","b1",1200,300);
  b1->Divide(5,1);
  b1->cd(1); b->Draw("var1>>var1b"); 
  b1->cd(2); b->Draw("var2>>var2b"); 
  b1->cd(3); b->Draw("var3>>var3b","");
  b1->cd(4); b->Draw("var4>>var4b"); 
  b1->cd(5); b->Draw("var5>>var5b",""); 
  b1->SaveAs("b1.png");
  out->Write();

  TCanvas *sb1 = new TCanvas("sb1","sb1",1200,300);
  sb1->Divide(5,1);
  sb1->cd(1); s->Draw("var1"); b->Draw("var1","","same"); d->Draw("var1","","e same");
  //sb1->cd(2); s->Draw("var2"); b->Draw("var2","","same");  d->Draw("var2","","e same");
  sb1->cd(2); b->Draw("var2"); s->Draw("var2","","same");  d->Draw("var2","","e same");
  sb1->cd(3); s->Draw("var3",""); b->Draw("var3","","same"); d->Draw("var3","","e same");
  sb1->cd(4); s->Draw("var4"); b->Draw("var4","","same"); d->Draw("var4","","e same"); 
  sb1->cd(5); s->Draw("var5",""); b->Draw("var5","","same");  d->Draw("var5","","e same");
  sb1->SaveAs("sb1.png");

  TCanvas *sb12=new TCanvas("sb12","sb12",1200,300);
  sb12->Divide(5,1);
  sb12->cd(1); s->Draw("var1"); b->Draw("var1","","same"); d2->Draw("var1","2","e same");
  sb12->cd(2); s->Draw("var2"); b->Draw("var2","","same");  d2->Draw("var2","2","e same");
  sb12->cd(3); s->Draw("var3",""); b->Draw("var3","","same"); d2->Draw("var3","","e same");
  sb12->cd(4); s->Draw("var4"); b->Draw("var4","","same"); d2->Draw("var4","2","e same"); 
  sb12->cd(5); s->Draw("var5",""); b->Draw("var5","","same");  d2->Draw("var5","","e same");
  
  sb12->SaveAs("sb12.png");



  TCanvas *s2=new TCanvas("s2","s2",1200,600);
  s2->Divide(5,2);  // 1,2 1,3 1,4 1,5 2,3 2,4 2,5 3,4 3,5 4,5
  s2->cd(1); s->Draw("var1:var2");
  s2->cd(2); s->Draw("var1:var3","var3<50");
  s2->cd(3); s->Draw("var1:var4");
  s2->cd(4); s->Draw("var1:var5","var5<80");
  s2->cd(5); s->Draw("var2:var3","var3<50");
  s2->cd(6); s->Draw("var2:var4");
  s2->cd(7); s->Draw("var2:var5","var5<80");
  s2->cd(8); s->Draw("var3:var4","var3<50");
  s2->cd(9); s->Draw("var5:var3","var5<80&&var3<50");
  s2->cd(10); s->Draw("var4:var5","var5<80");
  s2->SaveAs("s2.png");

  TCanvas *b2=new TCanvas("b2","b2",1200,600);
  b2->Divide(5,2);  // 1,2 1,3 1,4 1,5 2,3 2,4 2,5 3,4 3,5 4,5
  b2->cd(1); b->Draw("var1:var2");
  b2->cd(2); b->Draw("var1:var3","");
  b2->cd(3); b->Draw("var1:var4");
  b2->cd(4); b->Draw("var1:var5","");
  b2->cd(5); b->Draw("var2:var3","");
  b2->cd(6); b->Draw("var2:var4");
  b2->cd(7); b->Draw("var2:var5","var5<80");
  b2->cd(8); b->Draw("var3:var4","var3<50");
  b2->cd(9); b->Draw("var5:var3","var5<80&&var3<50");
  b2->cd(10); b->Draw("var4:var5","var5<80");
  b2->SaveAs("b2.png");


  TCanvas *sb2=new TCanvas("sb2","sb2",1200,600);
  TF1 *c1a=new TF1("c1a","3*x-0.2",-1,2);  c1a->SetLineColor(kGreen);
  TF1 *c1b=new TF1("c1b","3*x-2.9",-1,2);  c1b->SetLineColor(kGreen);
  TEllipse *c6 = new TEllipse(3,0.5,2,1.0); c6->SetLineColor(kGreen);  c6->SetFillStyle(0);  c6->SetLineWidth(3);
  TF1 *c9a=new TF1("c9a","7.0",0,55);  c9a->SetLineColor(kGreen);
  TF1 *c9b=new TF1("c9b","x+9",0,55);  c9b->SetLineColor(kGreen);

  sb2->Divide(5,2);  // 1,2 1,3 1,4 1,5 2,3 2,4 2,5 3,4 3,5 4,5
  sb2->cd(1);
  s->SetMarkerColor(kBlue);
  s->SetLineColor(kBlue);
  b->SetMarkerColor(kRed);
  s->Draw("var1:var2>>sv12");  b->Draw("var1:var2>>bv12");
  TH1F *sv12 = (TH1F*)gDirectory->Get("sv12");// bv12->SetMarkerColor(kRed);
  TH1F *bv12 = (TH1F*)gDirectory->Get("bv12"); 
  bv12->Draw(); sv12->Draw("box same");
  sb2->cd(2);
  s->Draw("var1:var3>>sv13","var3<50"); b->Draw("var1:var3>>bv13");
  TH1F *sv13 = (TH1F*)gDirectory->Get("sv13");
  TH1F *bv13 = (TH1F*)gDirectory->Get("bv13"); 
  bv13->Draw(); sv13->Draw("box same");
  sb2->cd(3);
  s->Draw("var1:var4>>sv14");   b->Draw("var1:var4>>bv14");
  TH1F *sv14 = (TH1F*)gDirectory->Get("sv14");
  TH1F *bv14 = (TH1F*)gDirectory->Get("bv14");
  bv14->Draw(); sv14->Draw("box same");
  sb2->cd(4);
  s->Draw("var1:var5>>sv15","var5<80"); b->Draw("var1:var5>>bv15");
  TH1F *sv15 = (TH1F*)gDirectory->Get("sv15");
  TH1F *bv15 = (TH1F*)gDirectory->Get("bv15");
  bv15->Draw(); sv15->Draw("box same");  
  sb2->cd(5);
  s->Draw("var2:var3>>sv23","var3<50"); b->Draw("var2:var3>>bv23");
  TH1F *sv23 = (TH1F*)gDirectory->Get("sv23");
  TH1F *bv23 = (TH1F*)gDirectory->Get("bv23");
  bv23->Draw(); sv23->Draw("box same"); 
  sb2->cd(6);
  s->Draw("var2:var4>>sv24"); b->Draw("var2:var4>>bv24");
  TH1F *sv24 = (TH1F*)gDirectory->Get("sv24");
  TH1F *bv24 = (TH1F*)gDirectory->Get("bv24");
  bv24->Draw(); sv24->Draw("box same"); 
  sb2->cd(7);
  s->Draw("var2:var5>>sv25","var5<80"); b->Draw("var2:var5>>bv25");
  TH1F *sv25 = (TH1F*)gDirectory->Get("sv25");
  TH1F *bv25 = (TH1F*)gDirectory->Get("bv25");
  bv25->Draw(); sv25->Draw("box same"); 
  sb2->cd(8);
  s->Draw("var3:var4>>sv34","var3<50"); b->Draw("var3:var4>>bv34");
  TH1F *sv34 = (TH1F*)gDirectory->Get("sv34");
  TH1F *bv34 = (TH1F*)gDirectory->Get("bv34");
  bv34->Draw(); sv34->Draw("box same");
  sb2->cd(9);
  s->Draw("var5:var3>>sv35","var5<80&&var3<50"); b->Draw("var5:var3>>bv35");
  TH1F *sv35 = (TH1F*)gDirectory->Get("sv35");
  TH1F *bv35 = (TH1F*)gDirectory->Get("bv35");
  bv35->Draw(); sv35->Draw("box same");
  sb2->cd(10);
  s->Draw("var4:var5>>sv45","var5<80");b->Draw("var4:var5>>bv45");
  TH1F *sv45 = (TH1F*)gDirectory->Get("sv45");
  TH1F *bv45 = (TH1F*)gDirectory->Get("bv45");
  bv45->Draw(); sv45->Draw("box same");

  sb2->SaveAs("sb2.png");

}

