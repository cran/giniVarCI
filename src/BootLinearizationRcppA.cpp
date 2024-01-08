#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tlinearizationA(NumericVector y, int n, int B, double G) { // Based on the biased Gini. Bootstrapping the linearization method.
  NumericVector BootG(B);
  double Gb;
  for(int i = 0; i < B; i++) {
    NumericVector SubSample = y[ floor(runif(n, 0, n)) ];
    std::sort(SubSample.begin(), SubSample.end());
    // Calculating the estimator
    double sum0 = 0.0;
    double sum1 = 0.0;
    for(int j = 0; j < n; j++){
      sum0 += SubSample[j];
      sum1 += (j + 1.0) * SubSample[j]/n;
    }
    double meany = sum0/n;
    Gb = (2.0*sum1/meany - n - 1.0 ) / n;
// Loop for calculating variance for Subsample
double sumz  = 0.0; // Sum z
double sumz2 = 0.0; // Sum z2
double Zj = 0.0;
double CumSumY = 0.0;
NumericVector YbarHat(n);
for(int j = 0; j < n; j++){
  CumSumY += y[j];
  YbarHat[j] = CumSumY/(j + 1.0);
}
for(int j = 0; j < n; j++){
  Zj = (2.0*(j + 1.0)*(y[j] - YbarHat[j])/n  + meany - y[j] - G * (meany + y[j]) ) / meany;
  sumz += Zj;
  sumz2 += pow(Zj,2);
} // for i
double meanz = sumz/n;
double VarGb = ( sumz2 / n - pow(meanz,2)) / (n - 1.0);
BootG[i] = (Gb - G)/sqrt(VarGb);
  } // for - i
return BootG;
}
