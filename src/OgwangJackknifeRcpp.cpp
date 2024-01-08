#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double OgwangJackknife(NumericVector SORTy, double G, int n, double sumy, double sumiy){
// To be applied on the biased expression (see Langel and Tille, 2013).
// ------------------------------------
// Declaration of vectors and variables
// ------------------------------------
  double sumyk  = 0.0;   // Sum y from 1 to k
  double sumGk  = 0.0;   // Sum GkG
  double sumGk2 = 0.0;   // Sum Gk^2
  double Gk;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    sumyk += SORTy[i];
    Gk = G + (2.0 / (sumy - SORTy[i])) * (SORTy[i] * sumiy / (sumy * n) + sumiy / (n * (n - 1.0)) - (sumy - sumyk + (i + 1.0)*SORTy[i])/ (n - 1.0)) - 1.0/(n * (n - 1.0));
    sumGk += Gk;
    sumGk2 += pow(Gk,2);
    } // for i
  return (n - 1.0) * (sumGk2/n - pow(sumGk/n, 2) ); // Obtain and return  Var(G)
}
