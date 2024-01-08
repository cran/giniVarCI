#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Linearization(NumericVector y, double G, int n, double meany){ // based on the biased gini
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double sumz  = 0.0; // Sum z
  double sumz2 = 0.0; // Sum z2
  double Zi;
  double Fi; // Sum value
  double Fyi; // Sum value y*delta(yj<yi)
  bool Delta1;
  bool Delta2;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Fi  = 0.0;
    Fyi = 0.0;
    for(int j = 0; j < n; j++){
      Delta1 = y[j]<=y[i];
      Delta2 = y[i]<=y[j];
      if(Delta1){
        Fi += 1.0;
      }
      if(Delta2){
        Fyi += y[j];
      }
    }  // for j
    Zi = (2.0*y[i]*Fi/n - (G + 1.0) * (y[i] + meany) + 2.0*Fyi/n) / meany;
    sumz += Zi;
    sumz2 += pow(Zi,2);
  } // for i
  return ( sumz2 / n - pow(sumz/n, 2)) / (n - 1.0); // Obtain and return the variance
}
