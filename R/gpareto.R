gpareto <- function(shape)
{
  if (is.numeric(shape) == TRUE)
  {
    output <- 1/(2*shape - 1)
    pos.NA <- shape <= 0
    output[pos.NA] <- NA
    output[(shape > 0) & (shape < 1)] <- 1
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
