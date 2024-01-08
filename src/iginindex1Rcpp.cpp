#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex1Rcpp(NumericVector y, int n, bool bc){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double Sumy = 0.0;
  double Sum = 0.0;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    for(int j = 0; j < n; j++){
      Sum += std::abs(y[i]-y[j]);
    }  // for j
  } // for i
if (bc == FALSE)
{
  Output  = Sum/(2.0*Sumy*n);
}
else
{
  Output  = Sum/(2.0*Sumy*(n - 1.0));
}
return Output;
}
