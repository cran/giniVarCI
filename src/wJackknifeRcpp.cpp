#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double wJackknife(NumericVector y, NumericVector w, NumericMatrix Delta, int n, double Ghat, double Nhat, String varformula, NumericVector PiU, int N, NumericVector ORDER, NumericVector reOrder, int method){
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector Z(n);
  NumericVector ySORT = y[ORDER];
  NumericVector wSORT = w[ORDER];
  NumericVector wSORT0(n);
  double Gk;
  double Nhatk; // Sum for Nhat
  double yhatk; // Sum for mean
  double output;
  // ------------------------------------
  if (method <=3)
  {
    double wwyk; // Sum for w2*y
    double wyCSk; // Sum for sum w*y*CS
    for(int k = 0; k < n; k++){ // Remove element k for Jackknife
      Nhatk = 0.0;
      yhatk = 0.0;
      wwyk  = 0.0;
      wyCSk = 0.0;
      for(int i = 0; i < n; i++){
        if (k != i)
        {
          Nhatk += wSORT[i];
          yhatk += wSORT[i]*ySORT[i];
          wwyk += pow(wSORT[i],2)*ySORT[i];
          wyCSk += wSORT[i]*ySORT[i]*Nhatk; // This is accumulative Nhat
        }
      }
      Gk = (2.0 * wyCSk - wwyk)/(Nhatk * yhatk ) - 1.0;
      Z[k] = (1.0 - wSORT[k]/Nhat) * (Ghat - Gk)/wSORT[k];
    } // for i
    Z = Z[reOrder]; // original order according y, for using Delta correctly
  } // method = 1 or 2 or 3
  if (method == 4)
  {
    double sumj;
    double sumi;
    for(int k = 0; k < n; k++){ // Remove element k for Jackknife
      Nhatk = 0.0; // Sum for Nhat
      yhatk = 0.0; // Sum for mean
      for(int i = 0; i < n; i++){
        if (k != i)
        {
          Nhatk += w[i];
          yhatk += w[i]*y[i];
        }
      } // for i
      sumi = 0.0;
      for(int i = 0; i < n; i++){
        if (k!=i)
        {
          sumj = 0.0;
          for(int j = 0; j < n; j++)
          {
            if (k!=j)
            {
              if (i!=j){
                sumj += w[j]*std::min(y[i],y[j]);
              }
            }
          } // for j
          sumi += w[i]*sumj/(Nhatk - w[i]);
        }
      } // for i
      Gk = 1.0 - sumi/yhatk;
      Z[k] = (1.0 - w[k]/Nhat) * (Ghat - Gk)/w[k];
    }
  } // end method = 4
  if (method == 5)
  {
    int nk = n - 1; // Sample size for the Jackknife procedure
    NumericVector Fwi(nk); // Vector for Fwi
    NumericVector yJ(n); // Vector y for Jackknife
    NumericVector wJ(n); // Vector w for Jackknife
    double wFwik; // Sum for w * Fwi
    double Cumw0k; // cumulative w0
    int i1;
    double MeanF;
    double SumG = 0.0;
    for(int k = 0; k < n; k++)
    { // Remove element k for Jackknife
      yJ = clone(ySORT);
      yJ.erase(k);
      wJ = clone(wSORT);
      wJ.erase(k);
      wJ.insert(wJ.begin(),0.0);
      Nhatk = 0.0;
      yhatk = 0.0;
      wFwik = 0.0;
      Cumw0k = 0.0;
      for(int i = 0; i < nk; i++){
        i1 = i + 1;
        Nhatk += wJ[i1];
        yhatk += wJ[i1]*yJ[i];
        Cumw0k += wJ[i];
        Fwi[i] = Cumw0k + wJ[i1]/2.0; // Here, Cumw0 is the accumulative of w.
        wFwik += wJ[i1]*Fwi[i];
      }
      MeanF = wFwik/pow(Nhatk,2);
      SumG = 0.0;
      for(int i = 0; i < nk; i++){
        i1 = i+1;
        SumG += wJ[i1]*(yJ[i]-yhatk/Nhatk)*(Fwi[i]/Nhatk - MeanF);
      }
      Gk = 2.0*SumG/yhatk;
      Z[k] = (1.0 - wSORT[k]/Nhat) * (Ghat - Gk)/wSORT[k];
    }
    Z = Z[reOrder]; // original order according y, for using Delta correctly
  }
// --------------------
// Variance computation
// --------------------
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
        output += (1.0 - 1.0/w[i] - 1.0/w[j] + SumPiU/n)*pow(Z[i]*w[i] - Z[j]*w[j],2);
      }
    }
    output = output / (n - 1.0);
  }
  return output; // Return the variance
}

