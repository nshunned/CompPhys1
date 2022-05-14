#include <TMath.h>

using TMath::Exp;

// double gaus
double gaus2(double* xPtr, double par[]){        
  double x = *xPtr;
  double A=par[0], m1=par[1], s1=par[2];
  double B=par[3], m2=par[4], s2=par[5];
  double e1=(x-m1)/s1;
  double e2=(x-m2)/s2;
  return A*Exp(-0.5*e1*e1) + B*Exp(-0.5*e2*e2);
}

// 4-parameters (one can be eliminated)
// This is what you end up with if you just shift the
// peak of the function in Wikipedia and add a normalization
// parameter
double gumbel(double* xPtr, double par[]){        
  double x0 = *xPtr;
  double C=par[0], mu=par[1]; 
  double a=par[2], b=par[3]; 
  double x=x0-mu;
  double ex=-(b*Exp(-a*x)+a*x);
  return C*Exp(ex);
}

// 3-parameter version, a more sensible formulation
double gumbel2(double* xPtr, double par[]){        
  double x0 = *xPtr;
  double C=par[0]*1000, mu=par[1]; 
  double a=par[2]; 
  double x=x0-mu;
  double ex=-(Exp(-a*x)+a*x);
  return C*a*Exp(ex);
}


double calcCHI2(TH1F* h, TF1* f, int *dof=0){
  double chi2=0;
  if (dof!=0) *dof=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    double x=h->GetBinCenter(i);
    double n=(h->GetBinContent(i));
    double s=(h->GetBinError(i));
    double mu=f->Eval(x);
    if (s>0){
      chi2+=(n-mu)*(n-mu)/s/s;
      if (dof) (*dof)++;
    }
  }
  return chi2;  
}
