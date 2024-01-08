#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double NormalityCpp(NumericVector y, double G, int n, double meany){ // based on the biased gini
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double sumu1 = 0.0; // Sum u1
  double sumu2 = 0.0; // Sum u12
  double Fi; // Sum value
  double Fyi; // Sum value y*delta(yj<yi)
  bool Delta1;
  bool Delta2;
  double u1i;
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
    u1i = 2.0*(y[i]*Fi + Fyi)/n - (G + 1.0)*y[i];
    sumu1 += u1i;
    sumu2 += pow(u1i,2);
  } // for i
  return (sumu2/n - pow(sumu1/n, 2)) /( pow(meany, 2) * (n - 1.0)) ; // Obtain and return the variance
}
