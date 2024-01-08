gf <- function(df1, df2)
{
  if ( (is.numeric(df1) == FALSE) || (is.numeric(df2) == FALSE) )
  {
    output <- NA
    warning("objects 'df1' and 'df2' must be numeric")
  }
  else if ( (length(df1) != 1L) || (length(df2) != 1L) )
  {
    output <- NA
    warning("objects 'df1' and 'df2' must have length 1")
  }
  else if ( (df1 <= 0L) || (df2 <= 0L) )
  {
    output <- NA
    warning("objects 'df1' and 'df2' must be positive")
  }
  else if (df2<=2)
  {
    output <- NA
    warning("the expectation of a F distribution F(df1,df2) exists and is finite only when df2 > 2")
  }
  else
  {
    InnerFunc <- function(y)
    {
      y * stats::df(y,df1, df2)
    }
    InnerIntegral <- Vectorize(function(y) { stats::integrate(InnerFunc, 0, stats::qf(y, df1, df2))$value})
    Integral2 <- stats::integrate(InnerIntegral, 0, 1)
    mean.f <- df2/(df2-2)
    output <- 2*(0.5-(1/mean.f)*Integral2$value)
    if( (output < 0) || (output > 1) )
    {
      output <- NA
      warning("the selected degrees of freedom 'df1' and 'df2' give a Gini index out of the interval [0,1]")
    }
  }
  return(output)
}
