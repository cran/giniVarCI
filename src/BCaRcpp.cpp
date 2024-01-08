#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector BCaCpp(NumericVector y, double sumy, int n, int B, double G, double alpha) {
  // Computes the BCa confidence limits
  // Vector y must be ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector BootG(B); //Bootstrap estimates for the Gini index
  NumericVector SubSample(n);
  double sumGb  = 0.0;
  double sumGb2 = 0.0;
  double sumz0  = 0.0;
  double sum0;
  double sum1;
  double meany;
  double z0;
  double sumu1;
  double sumu2;
  double sumu3;
  NumericVector x(n);
  int ni = n - 1; // Sample size for the Jackknife procedure
  double sumjy; // Sum value for the product i*x(i)/n
  double sumyi;
  double Gi;
  double sumdif2;
  double acc;
  double zalpha1;
  double tt1;
  double zalpha2;
  double tt2;
  NumericVector output(3); //Lower limit; Upper limit; variance
  // ------------------------------------
  for(int i = 0; i < B; i++) {
    SubSample = y[ floor(runif(n, 0, n)) ];
// ------------------------------------------------------
// Estimation of the Gini index for each bootstrap sample
// Computation of the bias correction factor z0__________
// ------------------------------------------------------
    std::sort(SubSample.begin(), SubSample.end());
    sum0 = 0.0;
    sum1 = 0.0;
    for(int j = 0; j < n; j++){
      sum0 += SubSample[j];
      sum1 += (j + 1.0) * SubSample[j]/ n;
    }
    meany = sum0/n;
    BootG[i] = (2.0*sum1/meany - n - 1.0)  / n;
    if(BootG[i] < G){
      sumz0 = sumz0 + 1.0;
    }
    sumGb  += BootG[i];          // This is for the variance
    sumGb2 += pow(BootG[i],2);   // This is for the variance
  } // end for i=0; i<B
  z0 = R::qnorm(sumz0/B, 0.0, 1.0, TRUE, FALSE);
// ------------------------------------------
// Computation of the acceleration factor acc
// ------------------------------------------
  sumu1 = 0.0;
  sumu2 = 0.0;
  sumu3 = 0.0;
 for(int i = 0; i < n; i++) {
   x = clone(y);
   x.erase(i);
   sumjy = 0.0; // Sum value for the product i*x(i)/n
   for(int j = 0; j < ni; j++){
     sumjy += (j + 1.0) * x[j] / ni;
   }
   sumyi = sumy - y[i];
   Gi = 2.0*sumjy/sumyi - 1.0 - 1.0/ni;
   sumu1 += Gi;
   sumu2 += pow(Gi,2);
   sumu3 += pow(Gi,3);
 }
sumdif2 = sumu2 - pow(sumu1,2)/n;
acc = (-2.0*pow(sumu1,3)/pow(n,2) + 3*sumu1*sumu2/n - sumu3)/(6*pow(sumdif2,1.5));
zalpha1 = R::qnorm(alpha/2, 0.0, 1.0, TRUE, FALSE);
tt1 = R::pnorm(z0 + (z0 + zalpha1) / (1.0 - acc * (z0 + zalpha1)), 0.0, 1.0, TRUE, FALSE);
zalpha2 = R::qnorm(1-alpha/2, 0.0, 1.0, TRUE, FALSE);
tt2 = R::pnorm(z0 + (z0 + zalpha2) / (1.0 - acc * (z0 + zalpha2)), 0.0, 1.0, TRUE, FALSE);
// -------------------------------------
// Quantiles for the bootstrap estimates
// -------------------------------------
 std::sort(BootG.begin(), BootG.end());
 output[0] = BootG[ floor(B * tt1)];
 output[1] = BootG[ floor(B * tt2)];
 output[2] = (sumGb2 - B * pow(sumGb/B, 2) )/(B - 1.0);
 return output;
}
