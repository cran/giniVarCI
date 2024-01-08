#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tjackknife(NumericVector y, int n, int B, double G) { // Based on the biased Gini. Bootstrapping the linearization method.
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector BootG(B);
  NumericVector SubSample(n);
  double Gb;
  double sum0;
  double sum1;
  double sumyk;   // Sum y from 1 to k
  double sumGk;   // Sum GkG
  double sumGk2;  // Sum Gk^2
  double Gk;
  double VarGb;   // Obtain Var(Gb)
  // ------------------------------------
// Loop for bootrapping
  for(int b = 0; b < B; b++) {
    SubSample = y[ floor(runif(n, 0, n)) ];
    std::sort(SubSample.begin(), SubSample.end());
// Loop for estimating G for each subsample
    sum0 = 0.0;
    sum1 = 0.0;
    for(int i = 0; i < n; i++){
      sum0 += SubSample[i];
      sum1 += (i + 1.0) * SubSample[i]/n;
    }
    Gb = (2.0 * n * sum1 / sum0 - n - 1.0 ) / n;
// Loop for calculating variance for Subsample
    sumyk  = 0.0;   // Sum y from 1 to k
    sumGk  = 0.0;   // Sum GkG
    sumGk2 = 0.0;   // Sum Gk^2
    for(int i = 0; i < n; i++){
      sumyk += SubSample[i];
      Gk = Gb + (2.0 / (sum0 - SubSample[i])) * (SubSample[i] * sum1 / sum0 + (sum1*n) / (n * (n - 1.0)) - (sum0 - sumyk + (i + 1.0)*SubSample[i])/ (n - 1.0)) - 1.0/(n * (n - 1.0));
      sumGk += Gk;
      sumGk2 += pow(Gk,2);
    } // for i
   VarGb = (n - 1.0) * (sumGk2/n - pow(sumGk/n, 2) ); // Obtain Var(Gb)
   BootG[b] = (Gb - G)/sqrt(VarGb);
  } // for - b
  return BootG;
}
