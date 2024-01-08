#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
List fginindex4Rcpp(NumericVector y, NumericVector w, int n){ //bool bc
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double sumj;
  double Nhat = 0.0; // Sum for Nhat
  double yhat = 0.0; // Sum for mean
  double sumi;
  // ------------------------------------
  for(int i = 0; i < n; i++){
    Nhat += w[i];
    yhat += w[i]*y[i];
  } // for i
  sumi = 0.0;
  for(int i = 0; i < n; i++){
    sumj = 0.0;
    for(int j = 0; j < n; j++){
      if (i!=j){
        sumj += w[j]*std::min(y[i],y[j]);
      }
    } // for j
  sumi += w[i]*sumj/(Nhat-w[i]);
  } // for i
//  double Output;
//  if (bc == FALSE)
//  {
//    Output = (pow(Nhat,2) - Sumw2)*(1.0 - sumi/yhat)/pow(Nhat,2);
//  }
//  else
//  {
//    Output = 1.0 - sumi/yhat;
//  }
//return 1.0 - sumi/yhat; //Output;
List res = List::create(Named("Ghat") = 1.0 - sumi/yhat,
                          Named("Nhat") = Nhat,
                          Named("MeanW") = yhat/Nhat);
return res;  //Output;
}
