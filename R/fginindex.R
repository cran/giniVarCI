fginindex <- function(y, w, method = 2L, Pi = NULL, na.rm = TRUE, useRcpp = TRUE )
{
  method <- as.integer(method)
  if ( (method < 1L) || (method > 5L) )   stop("object 'method' must be an integer between 1 and 5")
  if (missing(y))  stop("object 'y' is not specified")
  if (missing(w) && is.null(Pi))  stop("objects 'w' or 'Pi' must be specified. Use igini() for samples values without survey weights/inclusion probabilities")
  if (!is.numeric(y))  stop("object 'y' must be numeric")
  if (!missing(w) && !is.numeric(w))  stop("object 'w' must be numeric")
  if (!is.null(Pi) && !is.numeric(Pi))  stop("object 'Pi' must be numeric")
  if ( na.rm == TRUE )
  {
    y <- y[!is.na(y)]
    if (min(y) >= 0L)
    {
      Sample.size <- length(y)
      if (Sample.size > 1L)
      {
        if (!missing(w))
        {
          w <- w[!is.na(w)]
          LengthW <- length(w)
          if (LengthW == 1L)
          {
            LengthW <- Sample.size
            w <- rep(w, Sample.size)
            warning("length of object 'w' is 1 and has been replicated n times (same survey weights)")
          }
          if (Sample.size != LengthW) stop("object 'y' and 'w' have different length")
        }
        if (!is.null(Pi))
        {
          Pi <- Pi[!is.na(Pi)]
          if ( (min(Pi) < 0L) || (max(Pi) > 1L) )   stop("inclusion probabilities must be in the interval [0,1]")
          LengthPi <- length(Pi)
          if (LengthPi == 1L)
          {
            LengthPi <- Sample.size
            Pi <- rep(Pi, Sample.size)
            warning("length of object 'Pi' is 1 and has been replicated n times (same inclusion probabilities)")
          }
        }
        if (!missing(w) && !is.null(Pi) )
        {
          if(LengthW != LengthPi)              stop("specify 'w' or 'Pi' but not both")
          if ( sum(abs(w-1/Pi))  < 1e-15)   warning("specify 'w' or 'Pi' but not both")
          else                                 stop("specify 'w' or 'Pi' but not both")
        }
##############################################
#### Computing the estimator of the Gini index
##############################################
        if (missing(w)) w <- 1/Pi
        if (method == 1L)
        {
          output <- ifelse(useRcpp, fginindex1Rcpp(y, w, Sample.size)$Ghat, fginindex1(y, w) )
        }
        if (method == 2L)
        {
          ORDER <- order(y)
          output <- ifelse(useRcpp, fginindex2Rcpp(y[ORDER], w[ORDER], Sample.size)$Ghat, fginindex2(y[ORDER], w[ORDER]))
        }
        if (method == 3L)
        {
          output <- ifelse(useRcpp, fginindex3Rcpp(y, w, Sample.size)$Ghat, fginindex3(y, w) )
        }
        if (method == 4L)
        {
          output <- ifelse(useRcpp, fginindex4Rcpp(y, w, Sample.size)$Ghat, fginindex4(y, w) )
        }
        if (method == 5L)
        {
          ORDER <- order(y)
          output <- ifelse(useRcpp, fginindex5Rcpp(y[ORDER], c(0,w[ORDER]), Sample.size)$Ghat, fginindex5(y[ORDER], w[ORDER]) )
        }
      }
      else stop ("lenght of object 'y' must be larger than 1")
    }
    else
    {
      output <- NA
      warning("sample values must be non-negative")
    }
  }
  else output <- NA
  return(output)
}
