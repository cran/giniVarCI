gfisk <- function(shape1.a)
{
  if (is.numeric(shape1.a) == TRUE)
  {
    output <- 1/shape1.a
    pos.NA <- shape1.a <= 0
    output[pos.NA] <- NA
    output[(shape1.a > 0) & (shape1.a < 1)] <- 1
    sum.NA <- sum(pos.NA)
    if (sum.NA == 1L)
    {
      warning("an element of the object 'shape1.a' is not positive")
    }
    else if (sum.NA >= 1L)
    {
      warning(paste(sum.NA, " elements of the object 'shape1.a' are not positive"))
    }
  }
  else
  {
    output <- NA
    warning("object 'shape1.a' must be numeric")
  }
return(output)
}
