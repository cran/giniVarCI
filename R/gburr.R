gburr <- function(scale = 1, shape.g = 1, shape.s = 1)
{
  return(gparetoIV(location = 0, scale = scale, inequality = 1/shape.g, shape = shape.s))
}
