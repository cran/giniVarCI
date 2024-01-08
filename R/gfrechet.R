gfrechet <- function(shape)
{
  if (is.numeric(shape) == TRUE)
  {
    output <- vector(length = length(shape))
    pos.NA <- shape < 1
    output[pos.NA] <- NA
    output[!pos.NA] <-  2^(1/shape[!pos.NA]) - 1
    sum.NA <- sum(pos.NA)
    if (sum.NA == 1L)
    {
      warning("an element of the object 'shape' is smaller than 1")
    }
    else if (sum.NA >= 1L)
    {
      warning(paste(sum.NA, " elements of the object 'shape' are smaller than 1"))
    }
  }
  else
  {
    output <- NA
    warning("object 'shape' must be numeric")
  }
  return(output)
}
