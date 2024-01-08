#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex4Rcpp(NumericVector y, int n, bool bc){ //y is ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector CumSum(n + 1);
  CumSum[0] = 0.0;
  NumericVector pi(n + 1);
  pi[0] = 0.0;
  NumericVector qi(n + 1);
  double Sumy = 0.0;
  double Sum;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    CumSum[i + 1] = Sumy;
    pi[i + 1] = (i + 1.0)/n;
  }
  qi = CumSum/Sumy;
  Sum = 0.0;
  for(int i = 0; i < n; i++){
     Sum += (qi[i + 1] + qi[i])*(pi[i + 1] - pi[i]);
    }
if (bc == FALSE)
{
  Output  = 1.0 - Sum;
}
else
{
  Output  = n*(1.0 - Sum)/(n - 1.0);
}
return Output;
}
