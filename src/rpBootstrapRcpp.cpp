#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector rpBootstrap(NumericVector y, NumericVector w, int n, int B, double alpha, int method) {
  // Computes the Percentile Bootstrap and the variance using the Rescaled Bootstrap
  // y and w must be ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector BootG(B); //Bootstrap estimates
  NumericVector Labels(n); // Also used for Fwi in Method == 5
  NumericVector output(3); //Lower limit; Upper limit; variance
  double sumGb  = 0.0;
  double sumGb2 = 0.0;
  double Nhatb; // Sum for Nhat
  double yhatb; // Sum for mean
  double Sum1;
  double Sum2;
  int K;
  int n1 = n;
  if (method == 5)
  {
    n1 = n + 1;
  }
  // ------------------------------------
  for(int b = 0; b < B; b++) {
    std::vector<double> wstar(n1, 0.0);
    Labels = floor(runif(n, 0, n));
    // Calculate w^star for each bootstrap sample
    for(int i = 0; i < n; i++){
      K = Labels[i];
      if (method == 5)
      {
        K = Labels[i] + 1;
      }
      wstar[K] += w[Labels[i]]*n/(n - 1.0);
    }
    if (method <= 3)
    {
      Nhatb = 0.0;
      yhatb = 0.0;
      Sum1  = 0.0;
      Sum2  = 0.0;
      for(int i = 0; i < n; i++){
        Nhatb += wstar[i];
        yhatb += wstar[i]*y[i];
        Sum1 += pow(wstar[i],2)*y[i];
        Sum2 += wstar[i]*y[i]*Nhatb; // This is accumulative Nhatb
      }
      BootG[b] = ( 2.0 * Sum2 - Sum1 ) / ( Nhatb *  yhatb ) - 1.0;
    }
    if (method == 4)
    {
      Nhatb = 0.0;
      yhatb = 0.0;
      for(int i = 0; i < n; i++){
        Nhatb += wstar[i];
        yhatb += wstar[i]*y[i];
      } // for i
      Sum2 = 0.0;
      for(int i = 0; i < n; i++){
        Sum1 = 0.0;
        for(int j = 0; j < n; j++){
          if (i!=j){
            Sum1 += wstar[j]*std::min(y[i],y[j]);
          }
        } // for j
        Sum2 += wstar[i]*Sum1/(Nhatb-wstar[i]);
      } // for i
      BootG[b] = 1.0 - Sum2/yhatb;
    }
    if (method == 5)
    {
      Nhatb = 0.0;
      yhatb = 0.0;
      Sum1 = 0.0; // Sum for w * Fwi
      Sum2 = 0.0; // cumulative w0
      for(int i = 0; i < n; i++){
        K = i + 1;
        Nhatb += wstar[K];
        yhatb += wstar[K]*y[i];
        Sum2 += wstar[i];
        Labels[i] = Sum2 + wstar[K]/2.0;
        Sum1 += wstar[K]*Labels[i];
      }
      Sum2 = 0.0;
      for(int i = 0; i < n; i++){
        K = i + 1;
        Sum2 += wstar[K]*(y[i]-yhatb/Nhatb)*(Labels[i]/Nhatb - Sum1/pow(Nhatb,2) );
      }
      BootG[b] = 2.0*Sum2/yhatb;
    }
    sumGb  += BootG[b];
    sumGb2 += pow(BootG[b],2);
  }
  std::sort(BootG.begin(), BootG.end());
  output[0] = BootG[ floor(B * alpha/2.0)];
  output[1] = BootG[ floor(B - B*alpha/2.0)];
  output[2] = (sumGb2 - B * pow(sumGb/B, 2) )/(B - 1.0);
  return output;
}

