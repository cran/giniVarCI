#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector Order(NumericVector x, NumericVector y) {
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  return x[idx];
}
