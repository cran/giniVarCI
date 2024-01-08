#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List fginindex1Rcpp(NumericVector y, NumericVector w, int n){ //bool bc Ignore bc
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double Sumy = 0.0;
  double Sum = 0.0;
  double Nhat = 0.0;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Sumy += w[i]*y[i];
    Nhat += w[i];
    for(int j = 0; j < n; j++){
      Sum += w[i]*w[j]*std::abs(y[i]-y[j]);
    }  // for j
  } // for i
//  double Output;
//  if (bc == FALSE)
//  {
//    Output  = Sum/(2.0*Sumy*Nhat);
//  }
//  else
//  {
//    Output  = Nhat*Sum/(2.0*Sumy*(pow(Nhat,2)-Sumw2) );
//  }
List res = List::create(Named("Ghat") = Sum/(2.0*Sumy*Nhat),
                        Named("Nhat") = Nhat,
                        Named("MeanW") = Sumy/Nhat);

return res;  //Output;
}
