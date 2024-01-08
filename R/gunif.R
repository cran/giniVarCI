gunif <- function(min = 0, max = 1)
{
  if ( (is.numeric(min)==FALSE) || (is.numeric(max)==FALSE)  )
  {
    output <- NA
    warning("objects 'min' and 'max' must be numeric")
  }
  else if ( (length(min) != 1L) || (length(max) != 1L) )
  {
    output <- NA
    warning("objects 'min' and 'max' must have length 1")
  }
  else if ( (min < 0) || (max <= 0) )
  {
    output <- NA
    warning("objects 'min' and 'max' must be non-negative")
  }
  else if (min>=max)
  {
    output <- NA
    warning("object 'min' must be smaller than the object 'max'")
  }
  else
  {
    output <- (max - min)/(3*(min + max))
    if( (output < 0) || (output > 1) )
    {
      output <- NA
      warning("the selected lower and upper limits 'min' and 'max' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
