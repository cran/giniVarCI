#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double LinearizationA(NumericVector y, double G, int n, double meany){ // y must be ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector YbarHat(n);
  double sumz  = 0.0; // Sum z
  double sumz2 = 0.0; // Sum z2
  double Zi;
  double CumSumY = 0.0;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    CumSumY += y[i];
    YbarHat[i] = CumSumY/(i + 1.0);
  }
  for(int i = 0; i < n; i++){
    Zi = (2.0*(i + 1.0)*(y[i] - YbarHat[i])/n  + meany - y[i] - G * (meany + y[i]) ) / meany;
    sumz += Zi;
    sumz2 += pow(Zi,2);
  } // for i
  return ( sumz2 / n - pow(sumz/n,2)) / (n - 1.0); // Obtain and return the variance
}
