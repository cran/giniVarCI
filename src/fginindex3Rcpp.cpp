#include <Rcpp.h>
using namespace Rcpp;
// Place the export tag right above function declaration.
// [[Rcpp::export]]
List fginindex3Rcpp(NumericVector y, NumericVector w, int n){ // bool bc
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  double SumF;
  double Sumy = 0.0;
  double Sum = 0.0;
  int Delta1 = 0;
  int Delta2 = 0;
  double Nhat = 0.0;
  // ------------------------------------
  for(int j = 0; j < n; j++){
     SumF = 0.0;
    for(int i = 0; i < n; i++){
      Delta1 = 0;
      Delta2 = 0;
      if (y[i]<y[j]){
        Delta1 = 1;
      }
      if (y[i]==y[j]){
        Delta2 = 1;
      }
      SumF += w[i]*(Delta1 + 0.5*Delta2);
    }
    Nhat += w[j];
    Sumy += w[j]*y[j];
    Sum  += w[j]*y[j]*SumF;
  }
//  double Output;
//  if (bc == FALSE)
//  {
//    Output  = 2.0*Sum/(Sumy*Nhat) - 1.0;
//  }
//  else
//  {
//    Output  = pow(Nhat,2)*(2.0*Sum/(Sumy*Nhat) - 1.0)/(pow(Nhat,2) - Sumw2);
//  }

//return 2.0*Sum/(Sumy*Nhat) - 1.0; //Output;
List res = List::create(Named("Ghat") = 2.0*Sum/(Sumy*Nhat) - 1.0,
                        Named("Nhat") = Nhat,
                        Named("MeanW") = Sumy/Nhat);

return res;  //Output;
}
