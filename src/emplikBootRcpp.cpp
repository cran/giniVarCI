#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector emplikBootCpp(NumericVector y, double G, int n, int B, double precisionEL, int maxiterEL, double alpha){
  // CI using EL and critical value using Bootstrap. It assumes the biased Gini index.
  // Bootstrap samples are used to give an estimator of the variance.
  NumericVector output(3); // Lowerlimit; UpperLimit; variance
  // ------------------------------------
  // Declaration of vectors and variables
  // ------------------------------------
  NumericVector DFi(n);
  NumericVector BootR(B); //It contains R(theta) values
  NumericVector SubSample(n);
  NumericVector DFb(n);
  NumericVector ForLOG(n);
  double pert   = 1e-8;
  double rtol   = 1e-6;
  double atol   = 1e-8;
  double ctol   = 1e-8;
  double start  = 0.0;
  double Fi;
  double sumGb; // sum of Gb   for variance computation
  double sumGb2; // sum of Gb^2 for variance computation
  double sum0;
  double sum1;
  double Gb;
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
  int B2; // The numbers of replicates after removing NAs
  double criticalvalue;
  double L;
  double U;
  double Dif;
  // ------------------------------------
// ------------------------------------------------------------------
// Computation of DFi (distribution function) for the original sample
// ------------------------------------------------------------------
for(int i = 0; i < n; i++){
    Fi  = 0.0;
    for(int j = 0; j < n; j++){
      if(y[j]<=y[i]){
        Fi += 1.0;
      }
    }   // end for j
    DFi[i] = Fi/n;
  }     // end for i
// -----------------------------------------------------------
// Critical value using bootstrap. Only in this case theta = G
// -----------------------------------------------------------
// ----------------------
// Loop for bootstrapping
// ----------------------
sumGb  = 0.0;
sumGb2 = 0.0;
for(int b = 0; b < B; b++) {
  SubSample = y[ floor(runif(n, 0, n)) ];
  std::sort(SubSample.begin(), SubSample.end());
  // --------------------------------------------------------
  // Loop for estimating DF and Gini index for each subsample
  // --------------------------------------------------------
  sum0 = 0.0;
  sum1 = 0.0;
  for(int i = 0; i < n; i++){
    Fi  = 0.0; // Sum value for calculating Fi (Only the SUM)
    sum0 += SubSample[i];
    sum1 += (i + 1.0) * SubSample[i]/ n;
    for(int j = 0; j < n; j++){
    if(SubSample[j]<=SubSample[i]){
        Fi += 1.0;
      }
    }  // end for j
    DFb[i] = Fi/n;
  }    // end for i
  Gb = (2.0 * n * sum1 / sum0 - n - 1.0)  / n;
  sumGb  += Gb;
  sumGb2 += pow(Gb,2);
  // ------------------------------------------------------------------
  // Loop for searching Lambda. We use theta = G (see Qin et al., 2010)
  // ------------------------------------------------------------------
  x = start;
  for(int i = 0; i < maxiterEL; i++){
    refx = x;
    ewt  = rtol*std::abs(x) + atol;
    reffx = 0.0;
    for(int j = 0; j < n; j++){
      Z = (2*DFb[j] - 1.0)*SubSample[j] - G*SubSample[j];
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
      Z = (2.0*DFb[j] - 1.0)*SubSample[j] - G*SubSample[j];
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
  ForLOG = 1.0 + Lambda*( (2.0*DFb - 1.0)*SubSample - G*SubSample );
  NumericVector::iterator it = std::min_element(ForLOG.begin(), ForLOG.end());
  if (*it > 0.0){
    SumLog = 0.0;
    for(int i = 0; i < n; i++){
      SumLog += log(ForLOG[i]);
    }
    BootR[b]= SumLog; // we remove -1. Then, other -1 was required.
  }
  else {
    BootR[b] = NA_REAL;
  }
}     // end  for b
// -------------------------------------------
// Variance estimation using Bootstrap samples
// -------------------------------------------
output[2] = (sumGb2 - B * pow(sumGb/B, 2) )/(B - 1.0);
// --------------------------------------
// Critical value using bootstrap samples
// --------------------------------------
NumericVector BootR2 = BootR[!is_na(BootR)]; //NAs are removed
B2 = BootR2.size(); // The numbers of replicates after removing NAs
std::sort(BootR2.begin(), BootR2.end());
criticalvalue = BootR2[ floor(B2*(1-alpha))];
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
      // double RTheta = SumLog; // -1 is omitted. The next line also had -1.
      Dif  = SumLog - criticalvalue;
    }
    if (L<=0.0){
      Dif = 1.0;
    }
  }               // end while (dif <=0)
}                 // end if (L!=0)
if (U!=1)
{
  Dif = -1.0;
  while(Dif <= 0.0)
  {
    U = U + 0.01;
    //    Theta == U
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
    } // end search root
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
      //double RTheta = SumLog; //-1 is ommited. The next line also had -1.
      Dif = SumLog - criticalvalue;
    }
    if (U>=1.0){
      Dif = 1.0;
    }
  }      // end while (Dif <= 0.0)
}        // end if (U!=1)
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
    //double RTheta = SumLog; //-1 is ommited. The next line also had -1.
    Dif = SumLog - criticalvalue;
  }
  if (L>=U){
    Dif = -1.0;
  }
}                          // end while (Dif > 0)
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
      // double RTheta = SumLog; //-1 is ommited. The next line also had -1.
      Dif  = SumLog - criticalvalue;
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
// ------
// Output
// ------
  output[0] = L;
  output[1] = U;
return output; // Obtain and return the Lower limit, Upper limit, and variance
}
