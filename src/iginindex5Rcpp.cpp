#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex5Rcpp(NumericVector y, int n, bool bc){ // y is ordered
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double Sumy = 0.0;
  double Sumiy = 0.0;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    Sumiy += (i + 1.0) * y[i]/n;
    }
if (bc == FALSE)
{
  Output  = 2.0*Sumiy/Sumy - (n + 1.0) / n;
}
else
{
  Output  = 2.0*n*Sumiy/(Sumy*(n-1.0)) - (n + 1.0) / (n - 1.0);
}
return Output;
}
