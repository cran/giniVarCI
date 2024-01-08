#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pBootstrap(NumericVector y, int n, int B, double alpha) {
  // Computes the Percentile Bootstrap and the variance
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector BootG(B); //Bootstrap estimates
  NumericVector SubSample(n);
  NumericVector output(3); //Lower limit; Upper limit; variance
  double sumGb  = 0.0;
  double sumGb2 = 0.0;
  double sum0;
  double sum1;
  // ------------------------------------
  for(int b = 0; b < B; b++) {
    SubSample = y[ floor(runif(n, 0, n)) ];
    std::sort(SubSample.begin(), SubSample.end());
    sum0 = 0.0;
    sum1 = 0.0;
    for(int i = 0; i < n; i++){
      sum0 += SubSample[i];
      sum1 += (i + 1.0) * SubSample[i]/ n;
    }
    BootG[b] = (2.0 * n * sum1 / sum0 - n - 1.0)  / n;
    sumGb  += BootG[b];
    sumGb2 += pow(BootG[b],2);
  }
  std::sort(BootG.begin(), BootG.end());
  output[0] = BootG[ floor(B * alpha/2.0)];
  output[1] = BootG[ floor(B - B*alpha/2.0)];
  output[2] = (sumGb2 - B * pow(sumGb/B, 2) )/(B - 1.0);
  return output;
}
