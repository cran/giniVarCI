#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double wLinearizationA(NumericVector reOrder, NumericVector w, NumericVector ySORT, NumericVector wSORT, NumericMatrix Delta, int n, double Ghat, double Nhat, double MeanW, String varformula, NumericVector PiU, int N){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Z(n);
  double Nhati  = 0.0; // Cum Sum w ordered
  double Sumwy; // Sum wSORT*ySORT
  double output;
  // We obtain the vector Z
for(int i = 0; i < n; i++){
    Nhati  += wSORT[i];
    Sumwy = 0.0;
    for(int j = 0; j < i; j++){ // We calculate z[i], sums are on label "j"
      Sumwy += wSORT[j]*ySORT[j];
    }  // for j
    Z[i] = (2.0*Nhati*(ySORT[i] - Sumwy/Nhati) + Nhat*(MeanW - ySORT[i] - Ghat*(MeanW + ySORT[i]) ) )/  (pow(Nhat,2) * MeanW);
  } // for i
// We calculate the Variance
Z = Z[reOrder]; // original order according y, for using Delta correctly
output = 0.0;
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
