#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex10Rcpp(NumericVector y, int n, bool bc, NumericMatrix Matrix, int NumCol){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double Sumy = 0.0;
  double Sum = 0.0;
  double Output;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += y[i];
  }
 Sum = 0.0;
 for(int i = 0; i < NumCol; i++){
    Sum += std::abs(Matrix(0,i) - Matrix(1,i));
  }
if (bc == FALSE)
{
  Output  = (n - 1.0)*Sum/(2.0*Sumy*NumCol);
}
else
{
  Output  = n * Sum/(2.0*Sumy*NumCol);
}

return Output;
}
