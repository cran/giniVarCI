#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List giniMeansRcpp(NumericVector SORTx){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  int n = SORTx.size(); // Size of vector
  double sum0 = 0.0; // Sum value for the mean
  double sum1 = 0.0; // Sum value for the product i*x(i)/n
  // ------------------------------------
  for(int i = 0; i < n; i++){
    sum0 += SORTx[i];
    sum1 += (i + 1.0) * SORTx[i] / n;
  }
  List res = List::create(Named("Mean") = sum0/n,
                          Named("Meaniy") = sum1,
                          Named("Sample.size") = n);
  return res;  // Obtain and return a list with sums and the sample size
}
