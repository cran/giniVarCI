#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex8Rcpp(NumericVector y, int n, bool bc){
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
      Sum += std::min(y[i],y[j]);
    }  // for j
  } // for i
if (bc == FALSE)
{
  Output  = 1.0 - Sum/(Sumy*n);
}
else
{
  Output  = n/(n - 1.0) - Sum/(Sumy*(n - 1.0));
}
return Output;
}
