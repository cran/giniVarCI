#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double wLinearization(NumericVector y, NumericVector w, NumericMatrix Delta, int n, double Ghat, double Nhat, double MeanW, String varformula, NumericVector PiU, int N){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Z(n);
  double Fi  = 0.0; // Sum value
  double Fyi = 0.0; // Sum value y*delta(yj<yi)
  bool Delta1;
  bool Delta2;
  double output = 0.0;
// We obtain the vector Z
  for(int i = 0; i < n; i++){
    Fi  = 0.0; // Sum value
    Fyi = 0.0; // Sum value y*delta(yj<yi)
    for(int j = 0; j < n; j++){ // We calculate z[i], sums are on label "j"
      Delta1 = y[j]<=y[i];
      Delta2 = y[i]<=y[j];
      if(Delta1){
        Fi += w[j];
      }
      if(Delta2){
        Fyi += w[j]*y[j];
      }
    }  // for j
    Z[i] = (2.0*y[i]*Fi/Nhat - (Ghat + 1.0) * (y[i] + MeanW) + 2.0*Fyi/Nhat) / (Nhat * MeanW);
  } // for i
// We calculate the Variance
if(varformula == "SYG")
{
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      output += - Delta( i , j )*pow(Z[i]*w[i] - Z[j]*w[j],2)/2.0;
    }
  }
}
if(varformula == "HT")
{
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      output += Delta( i , j )*w[i]*w[j]*Z[i]*Z[j];
    }
  }
}
if(varformula == "HR")
{
  double SumPiU = 0.0;
  for(int i = 0; i < N; i++){
    SumPiU += PiU[i];
  }
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      output += (1-0 - 1.0/w[i] - 1.0/w[j] + SumPiU/n)*pow(Z[i]*w[i] - Z[j]*w[j],2);
    }
  }
  output = output / (n - 1.0);
}
return output; // Return the variance
}
