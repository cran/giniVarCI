#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
List fginindex5Rcpp(NumericVector y, NumericVector w0, int n){ // y, and w0 must be ordered; bool bc
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Fwi(n); // Vector for Fwi
  double Nhat = 0.0; // Sum for Nhat
  double yhat = 0.0; // Sum for mean
  double wFwi = 0.0; // Sum for w * Fwi
  double Cumw0 = 0.0; // cumulative w0
  double SumG = 0.0;
  int k;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    k = i + 1;
    Nhat += w0[k];
    yhat += w0[k]*y[i];
    Cumw0 += w0[i];
    Fwi[i] = Cumw0 + w0[k]/2.0; // Here, Cumw0 is the accumulative of w.
    wFwi += w0[k]*Fwi[i];
  }
  for(int i = 0; i < n; i++){
    k = i+1;
    SumG += w0[k]*(y[i]-yhat/Nhat)*(Fwi[i]/Nhat - wFwi/pow(Nhat,2) );
      }
//  double Output;
//  if (bc == FALSE)
//  {
//    Output  = 2.0 * SumG * (pow(Nhat,2) - Sumw2) / (pow(Nhat,2) * yhat)  ;
//  }
//  else
//  {
//    Output  = 2.0*SumG/yhat;
//  }
//return 2.0*SumG/yhat; //Output;
List res = List::create(Named("Ghat") = 2.0*SumG/yhat,
                          Named("Nhat") = Nhat,
                          Named("MeanW") = yhat/Nhat);
return res;  //Output;
}
