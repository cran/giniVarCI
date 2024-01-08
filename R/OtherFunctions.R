
GiniDagum <- function(shape2.p, gini, shape1.a)
{
  gamma(shape2.p)*gamma(2*shape2.p + 1/shape1.a)/(gamma(2*shape2.p)*gamma(shape2.p + 1/shape1.a)) - 1 - gini
}


GiniGamma <- function(k,gini)
{
  gamma((2*k + 1)/2)/(k*gamma(k)*sqrt(pi)) - gini
}


GiniChisq <- function(k,gini)
{
  2*gamma((1 + k)/2)/(k*gamma(k/2)*sqrt(pi)) - gini
}


SampleDagum <- function(shape2.p, gini)
{
  SOL.root <- stats::uniroot(GiniDagum, c(0.01, 10000), gini = gini, shape2.p = shape2.p)
  SOL.root$root
}


SampleGamma <- function(gini)
{
  SOL.root <- stats::uniroot(GiniGamma, c(0.000001,171), gini=gini)
  SOL.root$root
}


SampleChisq <- function(gini)
{
  SOL.root <- stats::uniroot(GiniChisq, c(1e-20, 341), gini=gini)
  SOL.root$root
}




iginindex1 <- function(y, Sample.size, bias.correction)
{
  ABSi <- sapply(y, function(Yi, Y) sum(abs(Y-Yi)), Y=y)
  output <- ifelse (bias.correction,  sum(ABSi)/(2*sum(y)*(Sample.size - 1)),
                    sum(ABSi)/(2*sum(y)*Sample.size) )
  return(output)
}


iginindex2 <- function(y, Sample.size, bias.correction)
{
  CumSum <- cumsum(y) # y must be ordered
  qi <- CumSum/CumSum[Sample.size]
  pi <- (1:Sample.size)/Sample.size
  Num <- sum(pi[1:(Sample.size - 1L)]-qi[1:(Sample.size - 1L)])
  Denom <- sum(pi[1:(Sample.size - 1L)])
  output <- ifelse (bias.correction,  Num/Denom,
                    (Sample.size - 1)*Num/(Denom*Sample.size) )
  return(output)
}



iginindex3 <- function(y, Sample.size, bias.correction)
{
  CumSum <- cumsum(y) # y must be ordered
  qi <- CumSum/CumSum[Sample.size]
  output <- ifelse (bias.correction,  1 - (2/(Sample.size - 1))*sum(qi[1:(Sample.size - 1L)]),
                    1 - 1/Sample.size - (2/Sample.size)*sum(qi[1:(Sample.size - 1L)]) )
  return(output)
}





iginindex4 <- function(y, Sample.size,  bias.correction)
{
  pi <- c(0, (1:Sample.size)/Sample.size)
  CumSum <- cumsum(y) # y must be ordered
  qi <- c(0, CumSum/CumSum[Sample.size])
  output <- ifelse (bias.correction,  Sample.size*(1 - sum((qi[2L:(Sample.size + 1L)] + qi[1L:Sample.size])*(pi[2L:(Sample.size + 1L)] - pi[1L:Sample.size]) ))/(Sample.size - 1),
                    1 - sum( (qi[2L:(Sample.size + 1L)] + qi[1L:Sample.size])*(pi[2L:(Sample.size + 1L)] - pi[1L:Sample.size]) ) )
  return(output)
}





iginindex5 <- function(y, Sample.size, bias.correction)
{
 #  y must be ordered
  output <- ifelse(bias.correction,
                   (2/((Sample.size-1)*mean(y)))*sum(((1:Sample.size)/Sample.size)*y) - (Sample.size+1)/(Sample.size-1), (2/sum(y))*sum(((1:Sample.size)/Sample.size)*y) - 1 - 1/Sample.size)
  return(output)
}



iginindex6 <- function(y, Sample.size, bias.correction)
{
  # y must be ordered
  output <- ifelse(bias.correction,
                   2 * stats::cov( (1:Sample.size), y ) / sum(y) ,
                   2 * (Sample.size - 1 ) * stats::cov( (1:Sample.size), y )/( Sample.size * sum(y) ) )
  return(output)
}


iginindex7 <- function(y, Sample.size, bias.correction)
{
  Fucnt.MidpointDF <- function(Yi,Y)  mean((Y<Yi) + 0.5*(Y==Yi))
  MidpointDFi <- sapply(y, Fucnt.MidpointDF, Y=y)
  Funct.AbsDFi <- function(i, Y) sum(abs(Y-Y[i])*abs(MidpointDFi-MidpointDFi[i]) )
  AbsDFi <- sapply(1:Sample.size, Funct.AbsDFi, Y=y)
  output <- ifelse(bias.correction,
                   mean(AbsDFi)/(mean(y)*(Sample.size-1)),
                   mean(AbsDFi)/sum(y) )
  return(output)
}



iginindex8 <- function(y, Sample.size, bias.correction)
{
  if (bias.correction == TRUE)
  {
    Funct.Min <- function(i, Y)  sum(pmin(Y[-i],Y[i]))
    zi <- sapply(1:Sample.size, Funct.Min, Y=y)
    output <- 1- sum(zi)/((Sample.size-1)*sum(y))
  }
  else if (bias.correction == FALSE)
  {
    Funct.Min <- function(i, Y)  sum(pmin(Y,Y[i]))
    zi <- sapply(1:Sample.size, Funct.Min, Y=y)
    output <- 1 - sum(zi)/(Sample.size*sum(y))
  }
  return(output)
}


iginindex9 <- function(y, Sample.size, bias.correction)
{
  Fucnt.MidpointDF <- function(Yi,Y)  mean((Y<Yi) + 0.5*(Y==Yi))
  MidpointDFi <- sapply(y, Fucnt.MidpointDF, Y=y)
  output <- ifelse(bias.correction,
                   2*(Sample.size-1)*sum(y*MidpointDFi)/mean(y) - Sample.size/(Sample.size-1),
                   2*mean(y*MidpointDFi)/mean(y) - 1 )
  return(output)
}


iginindex10 <- function(y, Sample.size, bias.correction, Matrix, NumCol)
{
  SUM <- sum(abs(Matrix[1L,]-Matrix[2L,]))
  output <- ifelse(bias.correction, (1/(2*mean(y)))*SUM/NumCol,
                   ((Sample.size - 1)/(2*sum(y)))*SUM/NumCol )
  return(output)
}



iinterval <- function(y, interval = c("zjackknife", "tjackknife", "zalinearization", "zblinearization", "talinearization", "tblinearization", "pbootstrap", "BCa", "ELchisq", "ELboot"), alpha, GiniBiased, Sample.size, sumy, sumiy, B, BC, bias.correction, precisionEL, maxiterEL)
{
  interval <- match.arg(interval)
  if (interval == "ELchisq")  criticalvalue <- stats::qchisq(alpha, df=1, lower.tail=FALSE)
  ResInt <- switch(interval,
                   zjackknife     = OgwangJackknife(y, GiniBiased, Sample.size, sumy, sumiy),
                   zalinearization = LinearizationA(y, GiniBiased, Sample.size, sumy/Sample.size),
                   zblinearization = Linearization(y, GiniBiased, Sample.size, sumy/Sample.size),
                   tjackknife     = tjackknife(y, Sample.size, B, GiniBiased),       # File BootJackknifeRcpp
                   talinearization = tlinearizationA(y, Sample.size, B, GiniBiased), # File BootLinearizationRcpp
                   tblinearization = tlinearization(y, Sample.size, B, GiniBiased),  # File BootLinearizationRcppA
                   pbootstrap     = pBootstrap(y, Sample.size, B, alpha),
                   BCa            = BCaCpp(y, sumy, Sample.size, B, GiniBiased, alpha),
                   ELchisq        = emplikCpp(y, sumy, GiniBiased, Sample.size, criticalvalue, precisionEL, maxiterEL),
                   ELboot        = emplikBootCpp(y, GiniBiased, Sample.size, B, precisionEL, maxiterEL, alpha)
  )
  if( (interval == "zjackknife") || (interval == "zalinearization") || (interval == "zblinearization") )
  {
    Z <- stats::qnorm(1 - alpha/2)
    M <- matrix(c(GiniBiased-Z*sqrt(ResInt), GiniBiased+Z*sqrt(ResInt) ), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    outputBiased   <- list(Gini = GiniBiased,    Interval = M,    Variance = ResInt        )
    outputUnbiased <- list(Gini = GiniBiased/BC, Interval = M/BC, Variance = ResInt/(BC^2L) )
  }
  if(interval == "tjackknife")
  {
    VarG <- OgwangJackknife(y, GiniBiased, Sample.size, sumy, sumiy)
    t2   <- stats::quantile(ResInt, alpha/2,   names = FALSE)
    t1   <- stats::quantile(ResInt, 1 - alpha/2, names = FALSE)
    M    <- matrix(c(GiniBiased - t1*sqrt(VarG), GiniBiased - t2*sqrt(VarG)), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    outputBiased   <- list(Gini = GiniBiased,    Interval = M,    Variance = VarG        )
    outputUnbiased <- list(Gini = GiniBiased/BC, Interval = M/BC, Variance = VarG/(BC^2L) )
  }
  if(interval == "talinearization")
  {
    VarG <- LinearizationA(y, GiniBiased, Sample.size, sumy/Sample.size)
    t2   <- stats::quantile(ResInt, alpha/2,   names = FALSE)
    t1   <- stats::quantile(ResInt, 1 - alpha/2, names = FALSE)
    M    <- matrix(c(GiniBiased - t1*sqrt(VarG), GiniBiased - t2*sqrt(VarG)), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    outputBiased   <- list(Gini = GiniBiased,    Interval = M,    Variance = VarG        )
    outputUnbiased <- list(Gini = GiniBiased/BC, Interval = M/BC, Variance = VarG/(BC^2) )
  }
  if(interval == "tblinearization")
  {
    VarG <- Linearization(y, GiniBiased, Sample.size, sumy/Sample.size)
    t2   <- stats::quantile(ResInt, alpha/2,   names = FALSE)
    t1   <- stats::quantile(ResInt, 1-alpha/2, names = FALSE)
    M    <- matrix(c(GiniBiased - t1*sqrt(VarG), GiniBiased - t2*sqrt(VarG)), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    outputBiased   <- list(Gini = GiniBiased,    Interval = M,    Variance = VarG        )
    outputUnbiased <- list(Gini = GiniBiased/BC, Interval = M/BC, Variance = VarG/(BC^2) )
  }
  if( (interval == "pbootstrap") || (interval == "BCa") || (interval == "ELchisq") || (interval == "ELboot") )
  {
    M    <- matrix(c(ResInt[1L], ResInt[2L]), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    outputBiased   <- list(Gini = GiniBiased,    Interval = M,    Variance = ResInt[3L] )
    outputUnbiased <- list(Gini = GiniBiased/BC, Interval = M/BC, Variance = ResInt[3L]/(BC^2L) )
  }
  ifelse (bias.correction,  output <- outputUnbiased, output <- outputBiased)
  # Trunc the limits at [0,1]
  if (output$Interval[1L,1L]<0) output$Interval[1L,1L] <- 0
  if (output$Interval[1L,2L]>1) output$Interval[1L,2L] <- 1
  return(output)
}




## Using abs over differences
fginindex1 <- function(y, w) #, bias.correction)
{
  N.Hat <- sum(w)
  Funct.ABSi <- function(Yi, Y, W) sum(W*abs(Y-Yi))
  ABSi <- sapply(y, Funct.ABSi, Y=y, W=w)
  Mean.y <- sum(w*y)/N.Hat
#  output <- ifelse (bias.correction,  sum(w*ABSi)/(2*(N.Hat^2-sum(w^2))*Mean.y),
#                    sum(w*ABSi)/(2*N.Hat^2*Mean.y) )
#  return(output)
  return(sum(w*ABSi)/(2*N.Hat^2*Mean.y))
}



# Using the Eurostat definition
fginindex2 <- function(y, w) #, bias.correction)
{
  CSw <- cumsum(w)
  N.Hat <- sum(w)
  Num <-  2*sum(w*y*CSw) - sum(w^2*y)
  Den <- N.Hat*sum(w*y)
#  output <- ifelse (bias.correction,  N.Hat^2*(Num/Den - 1)/(N.Hat^2 - sum(w^2)),
#                    Num/Den - 1 )
#return(output)
return(Num/Den - 1)
}


# Using the smooth distribution function
fginindex3 <- function(y, w) #, bias.correction)
{
  N.Hat <- sum(w)
  Funct.FDi <- function(Yi, Y, W, N.Hat) sum(W*((Y<Yi) + 0.5*(Y==Yi)) )/N.Hat
  FDi <- sapply(y, Funct.FDi, Y=y, W=w, N.Hat=N.Hat)
  Mean.y <- sum(w*y)/N.Hat
  #  output <- ifelse (bias.correction,  N.Hat*(2*sum(w*y*FDi)/Mean.y - N.Hat)/(N.Hat^2-sum(w^2)) ,
  #                    2*sum(w*y*FDi)/(N.Hat*Mean.y) - 1 )
  #  return(output)
  return(2*sum(w*y*FDi)/(N.Hat*Mean.y) - 1)
}



# using Berger and Gedek-Balay (2020)
fginindex4 <- function(y, w) #, bias.correction)
{
  N.Hat <- sum(w)
  Sample.size <- length(y)
  Funct.Min <- function(i, Y, W, N.Hat)  sum(W[-i]*pmin(Y[-i],Y[i]))/(N.Hat - W[i])
  zi <- sapply(1:Sample.size, Funct.Min, Y=y, W=w, N.Hat=N.Hat)
#  output <- ifelse (bias.correction,  1- sum(w*zi)/sum(w*y), (N.Hat^2 - sum(w^2))*(1- sum(w*zi)/sum(w*y) )/N.Hat^2
#                  )
#  return(output)
  return(1- sum(w*zi)/sum(w*y))
}


fginindex5 <- function(y, w) # , bias.correction)
{
  N.Hat <- sum(w)
  F.hat <- ( w / 2 + c(0, utils::head(cumsum(w), -1)) ) / N.Hat
  Mean.F <- sum(w*F.hat)/N.Hat
  Mean.y <- sum(w*y)/N.Hat
#  output <- ifelse (bias.correction,  2*sum( w * (y - Mean.y) * (F.hat - Mean.F) )/(N.Hat*Mean.y),
#                    2*sum( w * (y - Mean.y) * (F.hat - Mean.F) ) * (N.Hat^2 - sum(w^2))/(N.Hat^3*Mean.y) )
#  return(output)
  return(2*sum( w * (y - Mean.y) * (F.hat - Mean.F) )/(N.Hat*Mean.y))
}



finterval <- function(y, w, interval = c("zjackknife", "zalinearization", "zblinearization", "pbootstrap"), Pi, Pij, PiU, alpha, Sample.size, Ghatmethod, Ghat, Nhat, MeanW, varformula, ORDER, method, B)
{
# ghatmethod: Gini index estimation using the method selected by the user. Used for Jackknife
# Ghat: Gini index estimation with expression a or b. Used for Linearization.
  interval <- match.arg(interval)
  if (is.null(Pi)) # w is used for point estimation. Pi is required and checked for variance estimation.
  {
    Pi <- 1/w
    if ( (min(Pi) < 0L) || (max(Pi) > 1L) )   stop("values of the object 'Pi' must be in the interval [0,1]")
  }
  #        if ( (interval == "zjackknife") || (interval == "tjackknife") || (interval == "zlinearization") || (interval == "tlinearization")  )
  if ( (interval == "zjackknife") || (interval == "zalinearization") || (interval == "zblinearization")  )
  {
    MPi <- matrix(Pi, nrow = Sample.size, ncol = 1L)
    if (!is.null(Pij))
    {
      dimPij <- dim(Pij)
      if (dimPij[1L] != dimPij[2L])  stop("object 'Pij' is not a square matrix")
      if (dimPij[1L] != Sample.size) stop("object 'Pij' is not a (n by n) square matrix, where n is the lenght of the object 'y'")
    }
    else
    {
      Pij <- MPi%*%t(MPi)*( 1 - (1 - MPi)%*%t(1 - MPi)/ sum(1 - MPi) )
      diag(Pij) <- Pi
    }
    Delta <- (Pij - MPi%*%t(MPi))/Pij
  } # End intervals that require Pij
  if (missing(PiU) && (varformula == "HR")) stop("object 'PiU' must be specified when 'varformula = HR'")
  if (!missing(PiU)) {
    PiU <- PiU[!is.na(PiU)]
    if ( (min(PiU) < 0L) || (max(PiU) > 1L) )   stop("inclusion probabilities must be in the interval [0,1]")
    N <- length(PiU)
    if ( N <= Sample.size) stop("length of 'PiU' must be larger than length of objects 'Pi' or 'w'")
  }
   else
  { # Any value because it is required in Rcpp
     PiU <- 0
     N <- 0
   }
  ResInt <- switch(interval,
                   zjackknife      = wJackknife(y, w, Delta, Sample.size, Ghatmethod, Nhat, varformula, PiU, N, (ORDER - 1L), (order(ORDER) - 1L), method),
                   zalinearization = wLinearizationA((order(ORDER)-1), w, y[ORDER], w[ORDER], Delta, Sample.size, Ghat, Nhat, MeanW, varformula, PiU, N),
                   zblinearization = wLinearization(y, w, Delta, Sample.size, Ghat, Nhat, MeanW, varformula, PiU, N),
                   pbootstrap     = rpBootstrap(y[ORDER], w[ORDER], Sample.size, B, alpha, method)
                   )
if(interval == "zjackknife")
  {
    if(ResInt>=0)
    {
      Z <- stats::qnorm(1 - alpha/2)
      M <- matrix(c(Ghatmethod - Z*sqrt(ResInt), Ghatmethod + Z*sqrt(ResInt) ), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    }
    else
    {
      M <- matrix(c(NA, NA), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
      warning("the variance estimate is negative. Use other 'varformula' or an alternative variance estimation method")
    }
    output   <- list(Gini = Ghat,    Interval = M,    Variance = ResInt        )
  }
if( (interval == "zalinearization") || (interval == "zblinearization") )
  {
  if(ResInt>=0)
  {
    Z <- stats::qnorm(1 - alpha/2)
    M <- matrix(c(Ghat - Z*sqrt(ResInt), Ghat + Z*sqrt(ResInt) ), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
  }
  else
  {
    M <- matrix(c(NA, NA), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    warning("the variance estimate is negative. Use other 'varformula' or an alternative variance estimation method")
  }
  output   <- list(Gini = Ghat,    Interval = M,    Variance = ResInt        )
  }
  if(interval == "pbootstrap")
  {
    M    <- matrix(c(ResInt[1L], ResInt[2L]), nrow = 1L, dimnames = list(c(), c("lower","upper")) )
    output   <- list(Gini = Ghatmethod,    Interval = M,    Variance = ResInt[3L] )
  }
return(output)
}

