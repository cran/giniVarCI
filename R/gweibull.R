gweibull <- function(shape)
{
  if (is.numeric(shape) == TRUE)
  {
    output <- vector(length = length(shape))
    pos.NA <- shape <= 0
    output[pos.NA] <- NA
    output[!pos.NA] <-  1 - 2^(-1/shape[!pos.NA])
    sum.NA <- sum(pos.NA)
    if (sum.NA == 1L)
    {
      warning("an element of the object 'shape' is not positive")
    }
    else if (sum.NA >= 1L)
    {
      warning(paste(sum.NA, " elements of the object 'shape' are not positive"))
    }
  }
  else
  {
    output <- NA
    warning("object 'shape' must be numeric")
  }
  return(output)
}
