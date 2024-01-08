iginindex <- function(y, method = 5L, bias.correction = TRUE, cum.sums = NULL, na.rm = TRUE, useRcpp = TRUE )
  # It is not required to give the cumulative sums based on the non-decreasing order of the variable y.
{
  method <- as.integer(method)
  if ( (method < 1L) || (method > 10L) )   stop("object 'method' must be an integer between 1 and 10")
  if (missing(y) && is.null(cum.sums))     stop("objects 'y' or 'cum.sums' must be specified")
  if (!missing(y) && !is.numeric(y))       stop("object 'y' must be numeric")
  if (!is.null(cum.sums) && !is.numeric(cum.sums))    stop("object 'cum.sums' must be numeric")
  if ( na.rm == TRUE )
    {
  if (!missing(y))
      {
        y <- y[!is.na(y)]
        y.sort    <- Sort(y)
      }
  if (!is.null(cum.sums))
      {
        cum.sums <- cum.sums[!is.na(cum.sums)]
        y.sort.CS <- Sort(c(cum.sums[1L], diff(cum.sums)))
      }
  if (!missing(y) && !is.null(cum.sums) )
      {
        if (length(y)!=length(cum.sums))                 stop("specify 'y' or 'cum.sums' but not both.")
        if ( sum(abs(y.sort - y.sort.CS))  < 1e-15)   warning("specify 'y' or 'cum.sums' but not both.")
        else                                             stop("specify 'y' or 'cum.sums' but not both.")
      }
 if (!is.null(cum.sums)) y.sort <- y.sort.CS
  if (min(y.sort) >= 0L)
      {
        Sample.size <- length(y.sort)
        # NOTE: Cumulative sums are calculated because can be unordered
        if (Sample.size > 1L)
        {
          if (method == 1L)
          {
            output <- ifelse(useRcpp, iginindex1Rcpp(y.sort, Sample.size, bias.correction), iginindex1(y.sort, Sample.size, bias.correction) )
          } # End method 1
          if (method == 2L)
          {
            output <- ifelse(useRcpp, iginindex2Rcpp(y.sort, Sample.size, bias.correction), iginindex2(y.sort, Sample.size, bias.correction) )
          }
          if (method == 3L)
          {
            output <- ifelse(useRcpp, iginindex3Rcpp(y.sort, Sample.size, bias.correction), iginindex3(y.sort, Sample.size, bias.correction) )
          }
          if (method == 4L)
          {
            output <- ifelse(useRcpp, iginindex4Rcpp(y.sort, Sample.size, bias.correction), iginindex4(y.sort, Sample.size, bias.correction) )
          }
          if (method == 5L)
          {
            output <- ifelse(useRcpp, iginindex5Rcpp(y.sort, Sample.size, bias.correction), iginindex5(y.sort, Sample.size, bias.correction) )
          }
          if (method == 6L)
          {
            # The output is expression 5 if we calculate the covariance by hand.
            output <- ifelse(useRcpp, iginindex5Rcpp(y.sort, Sample.size, bias.correction), iginindex6(y.sort, Sample.size, bias.correction) )
          }
          if (method == 7L)
          {
            output <- ifelse(useRcpp, iginindex7Rcpp(y.sort, Sample.size, bias.correction), iginindex7(y.sort, Sample.size, bias.correction) )
          }
          if (method == 8L)
          {
            output <- ifelse(useRcpp, iginindex8Rcpp(y.sort, Sample.size, bias.correction), iginindex8(y.sort, Sample.size, bias.correction) )
          }
          if (method == 9L)
          {
            output <- ifelse(useRcpp, iginindex9Rcpp(y.sort, Sample.size, bias.correction), iginindex9(y.sort, Sample.size, bias.correction) )
          }
          if (method == 10L)
          {
            Funct.apply <- function(x, y) y[x]
            Matrix	<- apply(utils::combn(Sample.size,2L), 2L, Funct.apply, y = y.sort)
            NumCol <- as.integer(choose(Sample.size, 2L))
            output <- ifelse(useRcpp, iginindex10Rcpp(y.sort, Sample.size, bias.correction, Matrix, NumCol), iginindex10(y.sort, Sample.size, bias.correction, Matrix, NumCol) )
          }
        }
        else stop ("lenght of object 'y' must be larger than 1.")
      }
      else
      {
        output <- NA
        warning("sample values must be non-negative")
      }
    }
    else output <- NA # if (na.rm == TRUE)
  return(output)
}
