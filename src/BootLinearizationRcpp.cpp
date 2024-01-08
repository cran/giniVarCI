#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tlinearization(NumericVector y, int n, int B, double G) { // Based on the biased Gini. Bootstrapping the linearization method.
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector BootG(B);
  NumericVector SubSample(n);
  double Gb;
  double sum0;
  double sum1;
  double sumz; // Sum z
  double sumz2; // Sum z2
  double Fi;
  double Fyi; // Sum value y*delta(yj<yi)
  bool Delta1;
  bool Delta2;
  double Zi;
  double VarGb;
  // ------------------------------------
  for(int b = 0; b < B; b++) {
    SubSample = y[ floor(runif(n, 0, n)) ];
    std::sort(SubSample.begin(), SubSample.end());
    sum0 = 0.0;
    sum1 = 0.0;
    for(int i = 0; i < n; i++){
      sum0 += SubSample[i];
      sum1 += (i + 1.0) * SubSample[i]/n;
    }
    Gb = (2.0 * n * sum1 / sum0 - n - 1.0 ) / n;
    // Loop for calculating variance for Subsample
    sumz  = 0.0;
    sumz2 = 0.0;
    for(int i = 0; i < n; i++){
      Fi  = 0.0;
      Fyi = 0.0;
      for(int j = 0; j < n; j++){
        Delta1 = SubSample[j]<=SubSample[i];
        Delta2 = SubSample[i]<=SubSample[j];
        if(Delta1){
          Fi += 1.0;
        }
        if(Delta2){
          Fyi += SubSample[j];
        }
      }  // for j
      Zi = (2.0*SubSample[i]*Fi/n - (Gb + 1.0)*(SubSample[i] + sum0/n) + 2.0*Fyi/n) * n / sum0;
      sumz += Zi;
      sumz2 += pow(Zi, 2);
    } // for i
    VarGb = ( sumz2/n - pow(sumz/n, 2) ) / (n - 1.0);
    BootG[b] = (Gb - G)/sqrt(VarGb);
  } // for - b
return BootG;
}
