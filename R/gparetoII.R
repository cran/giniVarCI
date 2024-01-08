gparetoII <- function(location = 0, scale = 1, shape = 1)
{
  if ( (is.numeric(location) == FALSE) || (is.numeric(scale) == FALSE) || (is.numeric(shape) == FALSE) )
  {
    output <- NA
    warning("objects 'location', 'scale' and 'shape' must be numeric")
  }
  else if ( (length(location) != 1L) || (length(scale) != 1L) || (length(shape) != 1L) )
  {
    output <- NA
    warning("objects 'location', 'scale' and 'shape' must have length 1")
  }
  else if ( (location < 0L) || (scale <= 0L) || (shape <= 0L) )
  {
    output <- NA
    warning("objects 'scale' and 'shape' must be positive, and object 'location' non-negative")
  }
  else
  {
    InnerFunc <- function(y)
    {
      y * VGAM::dparetoII(y, location, scale, shape)
    }
    InnerIntegral <- Vectorize(function(y) { stats::integrate(InnerFunc, location, VGAM::qparetoII(y, location, scale, shape))$value})
    Integral2 <- stats::integrate(InnerIntegral, 0, 1)
    integral.MEAN <- stats::integrate(InnerFunc, lower = location, upper = Inf, subdivisions=2000)
    mean.ParetoII <- integral.MEAN$value
    output <- 2*(0.5 - (1/mean.ParetoII)*Integral2$value)
    if(output<0||output>1)
    {
      output <- NA
      warning("the selected parameters 'location', 'scale' and 'shape' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
