gbeta <- function(shape1, shape2)
{
  if ( (is.numeric(shape1) == FALSE) || (is.numeric(shape2) == FALSE) )
  {
    output <- NA
    warning("objects 'shape1' and 'shape2' must be numeric")
  }
  else if ( (length(shape1) != 1L) || (length(shape2) != 1L) )
  {
    output <- NA
    warning("objects 'shape1' and 'shape2' must have length 1")
  }
  else if ( (shape1 <= 0) || (shape2 <= 0) )
  {
    output <- NA
    warning("objects 'shape1' and 'shape2' must be positive")
  }
  else
  {
    output <- (2/shape1)*beta(shape1+shape2,shape1+shape2)/(beta(shape1,shape1)*beta(shape2,shape2))
    if( (output < 0) || (output > 1) )
    {
      output <- NA
      warning("the selected parameters 'shape1' and 'shape2' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
