gdagum <- function(shape1.a, shape2.p)
{
  if ( (is.numeric(shape1.a) == FALSE) || (is.numeric(shape2.p) == FALSE) )
  {
    output <- NA
    warning("objects 'shape1.a' and 'shape2.p' must be numeric")
  }
  else if ( (length(shape1.a) != 1L) || (length(shape2.p) != 1L) )
  {
    output <- NA
    warning("objects 'shape1.a' and 'shape2.p' must have length 1")
  }
  else if ( (shape1.a <= 0L) || (shape2.p <= 0L) )
  {
    output <- NA
    warning("objects 'shape1.a' and 'shape2.p' must be positive")
  }
  else
  {
    output <- gamma(shape2.p)*gamma(2*shape2.p+1/shape1.a)/(gamma(2*shape2.p)*gamma(shape2.p+1/shape1.a)) - 1
    if( (output < 0) || (output > 1) )
    {
      output <- NA
      warning("the selected parameters 'shape1.a' and 'shape2.p' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
