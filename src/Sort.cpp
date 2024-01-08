#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector Sort(NumericVector x) {
std::sort(x.begin(), x.end());
    return x;
}
