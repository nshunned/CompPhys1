// generate data with systematics
// "True" (parent) distribution is a linear function
// First sample the parent pdf
// Then add a systematic that affects the shape
// Finally add a systematic that shifts the overall normaliztion

double fcn(Double_t *x, Double_t *par){
  //return par[0] + par[1]*TMath::Log(x[0]) + par[2]*TMath::Log(x[0])*TMath::Log(x[0]);
  return par[0] + par[1]*TMath::Sin(par[2]*x[0]);
}


void makeData(){
  gStyle->SetOptStat(0);  // turn off histogram stat box
  TFile * tf=new TFile("mydata.root","recreate");

  TF1 model("model",fcn,1,500,3);
  model.SetParameter(0,1.23);
  model.SetParameter(1,10);
  model.SetParameter(2,1);
  TH1F *data=new TH1F("data","data;x;y",100,0,500);
  data->Sumw2();   // calculate proper uncertianties w/ weighted data
  data->FillRandom("model",2000);
  

  data->Draw();
  data->Write();
  
  //  data->Fit(&model);
  return;
}

