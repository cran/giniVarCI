#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector emplikCpp(NumericVector y, double sumy, double G, int n, double criticalvalue, double precisionEL, int maxiterEL){ // based on the biased gini
  NumericVector output(3); //Contains: Lowerlimit; Upperlimit; Variance
// -------------------------------------------------------------------
// Computation of the distribution function for each y_i
// Computation of u1 and u2 for Sigma22 and Sigma32 (Qin et al., 2010)
// -------------------------------------------------------------------
// ------------------------------------
// Declaration of vectors and variables
// ------------------------------------
  NumericVector DFi(n);
  NumericVector ForLOG(n);
  double sumu1; // Sum u1
  double sumu12; // Sum u12
  double sumu2; // Sum u2
  double sumu22; // Sum u22
  double pert   = 1e-8;
  double rtol   = 1e-6;
  double atol   = 1e-8;
  double ctol   = 1e-8;
  double start  = 0.0;
  double Fi; // Sum delta(yi<yj)
  double Fyi; // Sum y*delta(yj<yi)
  double u1i;
  double u2i;
  double Sigma32;
  double k;
  double L;
  double U;
  double Dif;
  double x;
  double refx;
  double ewt;
  double reffx;
  double Z;
  double delt;
  double fx;
  double jacob;
  double relchange;
  double Lambda;
  double SumLog;

  // ------------------------------------
  sumu1  = 0.0;
  sumu12 = 0.0;
  sumu2  = 0.0;
  sumu22 = 0.0;
for(int i = 0; i < n; i++){
    Fi  = 0.0;
    Fyi = 0.0;
    for(int j = 0; j < n; j++){
      if(y[j]<=y[i]){
        Fi += 1.0;
      }
      if(y[i]<=y[j]){
        Fyi += y[j];
      }
    }  // End for j
    DFi[i] = Fi/n;
    u1i = 2.0*(y[i]*Fi + Fyi)/n - (G + 1.0)*y[i];
    u2i = 2.0*y[i]*Fi/n  - (G + 1.0)*y[i];
    sumu1  += u1i;
    sumu12 += pow(u1i,2);
    sumu2  += u2i;
    sumu22 += pow(u2i,2);
  } // for i
  Sigma32 =  (sumu12 - n*pow(sumu1/n, 2))/(n-1);
  // double Sigma22 =  (sumu22 - n*pow(sumu2/n, 2))/(n-1);
  output[2] = n*Sigma32/pow(sumy, 2);
  k = (sumu22 - n*pow(sumu2/n, 2))/(Sigma32*(n-1)); // Sigma22/Sigma32;
// -------------------------------
// Finding an initial bound: [L,U]
// -------------------------------
  L = G - 0.05;
  U = G + 0.05;
  if(L<0){
    L = 0;
  }
  if (U>1){
    U = 1;
  }
  if (L!=0)
  {
    Dif = -1.0;
    while(Dif <= 0.0)
    {
      L = L - 0.01;
//   Note: Theta = L
// Search root (lambda)
  x = start;
  for(int i = 0; i < maxiterEL; i++){
  refx = x;
  ewt  = rtol*std::abs(x) + atol;
  reffx = 0.0;
  for(int j = 0; j < n; j++){
    Z = (2*DFi[j] - 1.0)*y[j] - L*y[j];
    reffx += Z/(1.0 + x*Z);
  }
  reffx = reffx/n;
if (std::abs(reffx/ewt)<1.0){
    break;
  }
  delt = std::max(std::abs(x) * pert, pert);
  x    = x + delt;
  fx = 0.0;
  for(int j = 0; j < n; j++){
    Z = (2.0*DFi[j] - 1.0)*y[j] - L*y[j];
    fx += Z/(1.0 + x*Z);
  }
  fx = fx/n;
  jacob   = (fx - reffx) / delt;
      x  = refx;
  relchange =  -1.0 * reffx/jacob;
    if (std::abs(relchange) < ctol){
      break;
    }
      x  = x + relchange;
}
Lambda = x;
ForLOG = 1.0 + Lambda*( (2.0*DFi - 1.0)*y - L*y );
NumericVector::iterator it = std::min_element(ForLOG.begin(), ForLOG.end());
      if (*it <= 0.0){
        Dif = 1.0;
      }
      else
      {
        SumLog = 0.0;
        for(int i = 0; i < n; i++){
          SumLog += log(ForLOG[i]);
        }
        // double RTheta = (-1.0)*SumLog;
        Dif  = 2.0 * SumLog - pow(k,-1) * criticalvalue;
      }
      if (L<=0.0){
        Dif = 1.0;
      }
    }                // end while (Dif <=0)
  }                  //end if (L!=0)
  if (U!=1)
  {
    Dif = -1.0;
    while(Dif <= 0.0)
    {
      U = U + 0.01;
      //  Theta == U
      // Search root (lambda)
      x = start;
      for(int i = 0; i < maxiterEL; i++){
        refx = x;
        ewt  = rtol*std::abs(x) + atol;
        reffx = 0.0;
        for(int j = 0; j < n; j++){
          Z = (2*DFi[j] - 1.0)*y[j] - U*y[j];
          reffx += Z/(1.0 + x*Z);
        }
        reffx = reffx/n;
        if (std::abs(reffx/ewt)<1.0){
          break;
        }
        delt = std::max(std::abs(x) * pert, pert);
        x    = x + delt;
        fx = 0.0;
        for(int j = 0; j < n; j++){
          Z = (2.0*DFi[j] - 1.0)*y[j] - U*y[j];
          fx += Z/(1.0 + x*Z);
        }
        fx = fx/n;
        jacob   = (fx - reffx) / delt;
        x  = refx;
        relchange =  -1.0 * reffx/jacob;
        if (std::abs(relchange) < ctol){
          break;
        }
        x  = x + relchange;
      }      // end search root
      Lambda = x;
      ForLOG = 1.0 + Lambda*( (2.0*DFi - 1.0)*y - U*y );
      NumericVector::iterator it = std::min_element(ForLOG.begin(), ForLOG.end());
      if (*it <= 0.0){
        Dif = 1.0;
      }
      else
      {
        SumLog = 0.0;
        for(int i = 0; i < n; i++){
          SumLog += log(ForLOG[i]);
        }
        //double RTheta = (-1.0)*SumLog;
        Dif  = 2.0 * SumLog - pow(k,-1) * criticalvalue;
      }
      if (U>=1.0){
        Dif = 1.0;
      }
    }              // while (Dif <=0)
  }                // end if (U!=1)
// --------------------------------------------------------
// Start the search for the final Confidence interval [L,U]
// --------------------------------------------------------
  L = L - precisionEL;
  Dif =  1.0;
  while(Dif > 0.0)
  {
    L = L + precisionEL;
    // Search root (lambda)
    x = start;
    for(int i = 0; i < maxiterEL; i++){
      refx = x;
      ewt  = rtol*std::abs(x) + atol;
      reffx = 0.0;
      for(int j = 0; j < n; j++){
        Z = (2*DFi[j] - 1.0)*y[j] - L*y[j];
        reffx += Z/(1.0 + x*Z);
      }
      reffx = reffx/n;
      if (std::abs(reffx/ewt)<1.0){
        break;
      }
      delt = std::max(std::abs(x) * pert, pert);
      x    = x + delt;
      fx = 0.0;
      for(int j = 0; j < n; j++){
        Z = (2.0*DFi[j] - 1.0)*y[j] - L*y[j];
        fx += Z/(1.0 + x*Z);
      }
      fx = fx/n;
      jacob   = (fx - reffx) / delt;
      x  = refx;
      relchange =  -1.0 * reffx/jacob;
      if (std::abs(relchange) < ctol){
        break;
      }
      x  = x + relchange;
    }
    Lambda = x;
    ForLOG = 1.0 + Lambda*( (2.0*DFi - 1.0)*y - L*y );
    NumericVector::iterator it = std::min_element(ForLOG.begin(), ForLOG.end());
    if (*it <= 0.0){
      Dif = 1.0;
    }
    else
    {
      SumLog = 0.0;
      for(int i = 0; i < n; i++){
        SumLog += log(ForLOG[i]);
      }
      //double RTheta = (-1.0)*SumLog;
      Dif  = 2.0 * SumLog - pow(k,-1) * criticalvalue;
    }
    if (L>=U){
      Dif = -1.0;
    }
  } // end while (Dif > 0)
  if (L<U){
    U = U + precisionEL;
    Dif =  1.0;
    while(Dif > 0.0)
    {
      U = U - precisionEL;
        // Search root (lambda)
      x = start;
      for(int i = 0; i < maxiterEL; i++){
        refx = x;
        ewt  = rtol*std::abs(x) + atol;
        reffx = 0.0;
        for(int j = 0; j < n; j++){
          Z = (2*DFi[j] - 1.0)*y[j] - U*y[j];
          reffx += Z/(1.0 + x*Z);
        }
        reffx = reffx/n;
        if (std::abs(reffx/ewt)<1.0){
          break;
        }
        delt = std::max(std::abs(x) * pert, pert);
        x    = x + delt;
        fx = 0.0;
        for(int j = 0; j < n; j++){
          Z = (2.0*DFi[j] - 1.0)*y[j] - U*y[j];
          fx += Z/(1.0 + x*Z);
        }
        fx = fx/n;
        jacob   = (fx - reffx) / delt;
        x  = refx;
        relchange =  -1.0 * reffx/jacob;
        if (std::abs(relchange) < ctol){
          break;
        }
        x  = x + relchange;
      }
      Lambda = x;
      ForLOG = 1.0 + Lambda*( (2.0*DFi - 1.0)*y - U*y );
      NumericVector::iterator it = std::min_element(ForLOG.begin(), ForLOG.end());
    if (*it <= 0.0){
        Dif = 1.0;
      }
      else
      {
        SumLog = 0.0;
        for(int i = 0; i < n; i++){
          SumLog += log(ForLOG[i]);
        }
        //double RTheta = (-1.0)*SumLog;
        Dif  = 2.0 * SumLog - pow(k,-1) * criticalvalue;
      }
      if (U<=L){
        Dif = -1.0;
      }
    } // while (Dif > 0)
  }
  else {
    L = G;
    U = G;
  }
  output[0] = L;
  output[1] = U;
  return output; // Obtain and return the Confidence interval
}
