#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex2Rcpp(NumericVector y, int n, bool bc){ // y is ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector CumSum(n);
  NumericVector pi(n);
  NumericVector qi(n);
  double Sumy = 0.0;
  double SumNum;
  double SumDen;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    CumSum[i] = Sumy;
    pi[i] = (i + 1.0)/n;
  }
  qi = CumSum/Sumy;
  SumNum = 0.0;
  SumDen = 0.0;
  for(int i = 0; i < (n-1); i++){
     SumNum += pi[i]-qi[i];
     SumDen += pi[i];
    }
if (bc == FALSE)
{
  Output  = (n - 1.0)*SumNum/(n*SumDen);
}
else
{
  Output  = SumNum/SumDen;
}
return Output;
}
