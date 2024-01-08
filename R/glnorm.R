glnorm <- function(sdlog)
{
  if (is.numeric(sdlog) == TRUE)
  {
    output <- 2 * stats::pnorm(sdlog/sqrt(2)) - 1
    pos.NA <- sdlog <= 0
    output[pos.NA] <- NA
    sum.NA <- sum(pos.NA)
    if (sum.NA == 1L)
    {
      warning("an element of the object 'sdlog' is not positive")
    }
    else if (sum.NA >= 1L)
    {
      warning(paste(sum.NA, " elements of the object 'sdlog' are not positive"))
    }
  }
  else
  {
    output <- NA
    warning("object 'sdlog' must be numeric")
  }
return(output)
}
