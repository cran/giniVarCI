ggompertz <- function(scale = 1, shape)
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
    InnerFunc <- function(y)
    {
      y * VGAM::dgompertz(y, scale, shape)
    }
    InnerIntegral <- Vectorize(function(y) { stats::integrate(InnerFunc, 0, VGAM::qgompertz(y, scale, shape))$value})
    Integral2 <- stats::integrate(InnerIntegral, 0, 1)

    integral.MEAN <- stats::integrate(InnerFunc, lower = 0, upper = Inf)
    mean.Gom <- integral.MEAN$value
    output <- 2*(0.5 - (1/mean.Gom)*Integral2$value)
    if((output<0)||(output>1))
    {
      output <- NA
      warning("the selected parameters 'scale' and 'shape' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
