gparetoIV <- function(location = 0, scale = 1, inequality = 1, shape = 1)
{
  if ( (is.numeric(location) == FALSE) || (is.numeric(scale) == FALSE) || (is.numeric(inequality) == FALSE) || (is.numeric(shape) == FALSE) )
  {
    output <- NA
    warning("objects 'location', 'scale', 'inequality' and 'shape' must be numeric")
  }
  else if ( (length(location) != 1L) || (length(scale) != 1L) || (length(inequality) != 1L) || (length(shape) != 1L) )
  {
    output <- NA
    warning("objects 'location', 'scale', 'inequality' and 'shape' must have length 1")
  }
  else if ( (location < 0L) || (scale <= 0L) || (inequality <= 0L) || (shape <= 0L) )
  {
    output <- NA
    warning("objects 'scale', 'inequality' and 'shape' must be positive, and object 'location' non-negative")
  }
#  else if (inequality > shape)
#  {
#    output <- NA
#    warning("object 'shape' must be smaller than the object 'inequality'")
#  }
  else
  {
    InnerFunc <- function(y)
    {
      y * VGAM::dparetoIV(y, location, scale, inequality, shape)
    }
    InnerIntegral <- Vectorize(function(y) { stats::integrate(InnerFunc, location, VGAM::qparetoIV(y, location, scale, inequality, shape))$value})
    Integral2 <- stats::integrate(InnerIntegral, 0, 1)
    integral.MEAN <- stats::integrate(InnerFunc, lower = location, upper = Inf, subdivisions=2000)
    mean.ParetoIV <- integral.MEAN$value
    output <- 2*(0.5 - (1/mean.ParetoIV)*Integral2$value)
    if( (output < 0) || (output >1) )
    {
      output <- NA
      warning("the selected parameters 'location', 'scale', 'ineqaulity' and 'shape' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
