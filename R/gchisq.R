gchisq <- function(df)
{
  if (is.numeric(df) == TRUE)
  {
    output <- vector(length = length(df))
    pos.NA <- df<=0
    output[pos.NA] <- NA
    output[!pos.NA] <- 2*gamma((1 + df[!pos.NA])/2)/(df[!pos.NA]*gamma(df[!pos.NA]/2)*sqrt(pi))
    sum.NA <- sum(pos.NA)
    if (sum.NA == 1L)
    {
      warning("an element of the object 'df' is not positive")
    }
    else if (sum.NA >= 1L)
    {
      warning(paste(sum.NA, " elements of the object 'df' are not positive"))
    }
  }
  else
  {
    output <- NA
    warning("object 'df' must be numeric")
  }
  return(output)
}
