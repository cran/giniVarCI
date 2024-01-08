#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex9Rcpp(NumericVector y, int n, bool bc){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Fn(n);
  double Sumy = 0.0;
  double Sum = 0.0;
  double SumF;
  double Output;
  // ------------------------------------
  for(int j = 0; j < n; j++){
      SumF = 0.0;
      for(int i = 0; i < n; i++){
       if (y[i]<y[j])   SumF += 1;
       if (y[i]==y[j])  SumF += 0.5;
      }
      Fn[j] = SumF/n;
    Sumy += y[j];
    Sum += y[j]*Fn[j];
  }
if (bc == FALSE)
{
  Output  = 2.0*Sum/Sumy - 1.0;
}
else
{
  Output  = 2.0 * n * Sum/(Sumy*(n - 1.0)) - n/(n - 1.0);
}
return Output;
}
