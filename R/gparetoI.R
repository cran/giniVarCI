gparetoI <- function(scale = 1, shape = 1)
{
  if ( (is.numeric(scale) == FALSE) || (is.numeric(shape) == FALSE) )
  {
    output <- NA
    warning("objects 'scale' and 'shape' must be numeric")
  }
  else if ( (length(scale) != 1L) || (length(shape) != 1L) )
  {
    output <- NA
    warning("objects 'scale' and 'shape' must have length 1")
  }
  else if ( (scale <= 0) || (shape <= 0) )
  {
    output <- NA
    warning("objects 'scale' and 'shape' must be positive")
  }
  else
  {
    location <- scale
    InnerFunc <- function(y)
    {
      y * VGAM::dparetoI(y, scale, shape)
    }
    InnerIntegral <- Vectorize(function(y) { stats::integrate(InnerFunc, location, VGAM::qparetoI(y, scale, shape))$value})
    Integral2 <- stats::integrate(InnerIntegral, 0, 1)
    integral.MEAN <- stats::integrate(InnerFunc, lower = location, upper = Inf, subdivisions=2000)
    mean.ParetoI <- integral.MEAN$value
    output <- 2*(0.5-(1/mean.ParetoI)*Integral2$value)
    if(output<0||output>1)
    {
      output <- NA
      warning("the selected parameters 'scale' and 'shape' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
