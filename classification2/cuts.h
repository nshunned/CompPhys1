#ifndef CUTS_H_
#define CUTS_H_

// use this function to define events passing cuts
// a simple example of one cut is given
bool pass(Float_t *v){
  Float_t var1=v[0]; //, var2=v[1], var3=v[2], var4=v[3], var5=v[4];
  Float_t var2=v[1];
  Float_t var3=v[2];
  Float_t var4=v[3];
  Float_t var5=v[4];
  bool pass= (var1>66 && var1<150) && var3>0.2 && (var4>80 && var4<170) && var5>-5;
  if (pass) return true;
  return false;
}

#endif
