#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex3Rcpp(NumericVector y, int n, bool bc){ // y is ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector CumSum(n);
  NumericVector qi(n);
  double Sumy = 0.0;
  double SumNum;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    CumSum[i] = Sumy;
  }
  qi = CumSum/Sumy;
  SumNum = 0.0;
  for(int i = 0; i < (n-1); i++){
     SumNum += qi[i];
    }
if (bc == FALSE)
{
  Output  = (n - 1.0)/n - 2.0*SumNum/n;
}
else
{
  Output  = 1.0 - 2.0*SumNum/(n - 1.0);
}
return Output;
}
