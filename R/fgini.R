fgini <- function(y, w, method = 2L, interval = NULL, Pi = NULL, Pij = NULL, PiU, alpha = 0.05, B = 1000L, na.rm = TRUE, varformula = "SYG", large.sample = FALSE) #, precisionEL=1e-4, maxiterEL=100)
{
  method <- as.integer(method)
  if ( (method < 1L) || (method > 5L) )   stop("object 'method' must be an integer between 1 and 5")
  if (missing(y))  stop("object 'y' is not specified")
  if (missing(w) && is.null(Pi))  stop("objects 'w' or 'Pi' must be specified. Use igini() for samples values without survey weights/inclusion probabilities")
  if (!is.numeric(y))  stop("object 'y' must be a numeric vector")
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
        if (Sample.size != LengthW) stop("objects 'y' and 'w' have different length")
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
          ResG <- fginindex1Rcpp(y, w, Sample.size)
          output <- ResG$Ghat
        }
        if (method == 2L)
        {
          if (large.sample == FALSE)
          {
            ORDER <- Order(1:Sample.size,y)
          }
          else
          {
            ORDER <- order(y)
          }
          ResG <- fginindex2Rcpp(y[ORDER], w[ORDER], Sample.size)
          output <- ResG$Ghat
        }
        if (method == 3L)
        {
          ResG <- fginindex3Rcpp(y, w, Sample.size)
          output <- ResG$Ghat
        }
        if (method == 4L)
        {
          ResG <- fginindex4Rcpp(y, w, Sample.size)
          output <- ResG$Ghat
        }
        if (method == 5L)
        {
          if (large.sample == FALSE)
          {
            ORDER <- Order(1:Sample.size,y)
          }
          else
          {
            ORDER <- order(y)
          }
          ResG <- fginindex5Rcpp(y[ORDER], c(0,w[ORDER]), Sample.size)
          output <- ResG$Ghat
        }
##############################################
##### Computing the Confidence Interval ######
##############################################
        if (!is.null(interval))
          {
        if ( (method == 1L) || (method == 3L) || (method == 4L) )
          {
          if (large.sample == FALSE)
          {
            ORDER <- Order(1:Sample.size,y)
          }
          else
          {
            ORDER <- order(y)
          }
          }
        if ( (method == 4L ) || (method == 5L))
        {
          ResG <- fginindex2Rcpp(y[ORDER], w[ORDER], Sample.size)
        }
          outputCI <- finterval(y, w, interval = interval, Pi, Pij, PiU, alpha, Sample.size, output, ResG$Ghat, ResG$Nhat, ResG$MeanW, varformula, ORDER, method, B)
          # Gini estimation in output is the one selected by the user
          outputCI$Gini <- output
          output <- outputCI
        } # is.null(interval)
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
