#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
List fginindex2Rcpp(NumericVector y, NumericVector w, int n){ //bool bc
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double Nhat = 0.0; // Sum for Nhat
  double yhat = 0.0; // Sum for mean
  double wwy = 0.0; // Sum for w2*y
  double wyCS = 0.0; // Sum for sum w*y*CS
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Nhat += w[i];
    yhat += w[i]*y[i];
    wwy += pow(w[i],2)*y[i];
    wyCS += w[i]*y[i]*Nhat; // This is accumulative Nhat
  }
//  double Output;
//  if (bc == FALSE)
//  {
//    Output  = ( 2.0 * wyCS - wwy ) / ( Nhat *  yhat ) - 1.0;
//  }
//  else
//  {
//    Output  = Nhat * (( 2.0 * wyCS - wwy ) / yhat  - Nhat)/(pow(Nhat,2) - Sumw2);
//  }
//return ( 2.0 * wyCS - wwy ) / ( Nhat *  yhat ) - 1.0; //Output;
List res = List::create(Named("Ghat") = ( 2.0 * wyCS - wwy ) / ( Nhat *  yhat ) - 1.0,
                          Named("Nhat") = Nhat,
                          Named("MeanW") = yhat/Nhat);
return res;  //Output;
  }
