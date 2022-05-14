//************************************************************************          
// fib.cpp 
//************************************************************************

#include <vector>
#include <iostream>
#include <stdlib.h> 

using namespace std;

void add(vector<int>* term, vector<int>* term_1, vector<int>* term_2);
void print_vec(vector<int>* term);

int main(int argc, char* argv[]) {
  int n_digits;
  if (argc == 2)
    n_digits = atoi(argv[1]);
  else
    cout << "Please enter the correct number of parameters!" << endl;
  
  vector<int>* term_2 = new vector<int>;
  vector<int>* term_1 = new vector<int>;
  vector<int>* term = new vector<int>;
  vector<int>* temp;

  term_1->push_back(1);
  term->push_back(1);

  int n_th = 2;
  while((int)term->size() != n_digits) {
    n_th++;
    
    temp = term_2;
    term_2 = term_1;
    term_1 = term;
    
    term = temp;
    term->clear();
    add(term, term_1, term_2);
  }

  cout << "F" << n_th << " = ";
  print_vec(term);
  
  return 1;
}

void add(vector<int>* term, vector<int>* term_1, vector<int>* term_2) {
  int s = 0;
  int r = 0;

  int d_1 = (int)term_1->size();
  int d_2 = (int)term_2->size();

  for(int i = 0; i < d_1 - d_2; i++)
    term_2->push_back(0);

  for(int i = 0; i < d_1; i++) {
    s = term_1->at(i) + term_2->at(i) + r;
    r = s/10;
    s = s%10;
    term->push_back(s);
  }
  if(d_1 == d_2 && r)
    term->push_back(1);
}

void print_vec(vector<int>* term) {
  for(int i = term->size()-1; i >= 0; i--)
    cout << term->at(i);
  cout << endl;
}
