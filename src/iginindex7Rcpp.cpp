#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
double iginindex7Rcpp(NumericVector y, int n, bool bc){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Fn(n);
  double SumF;
  double Sumy;
  double Sum;
  double Output;
  // ------------------------------------
    for(int j = 0; j < n; j++){
      SumF = 0.0;
      for(int i = 0; i < n; i++){
       if (y[i]<y[j])   SumF += 1;
       if (y[i]==y[j])  SumF += 0.5;
      }
      Fn[j] = SumF/n;
  }
  Sumy = 0.0;
  Sum = 0.0;
  for(int i = 0; i < n; i++){
    Sumy += y[i];
    for(int j = 0; j < n; j++){
      Sum += std::abs(y[i] - y[j]) * std::abs(Fn[i] - Fn[j]);
    }  // for j
  } // for i
if (bc == FALSE)
{
  Output  = Sum/(Sumy*n);
}
else
{
  Output  = Sum/(Sumy*(n - 1.0));
}
return Output;
}
