gparetoIII <- function(inequality = 1)
{
if (is.numeric(inequality) == TRUE)
{
  output <- inequality
  pos.NA0 <- output < 0
  pos.NA1 <- output > 1
  output[pos.NA0] <- NA
  output[pos.NA1] <- NA
  sum.NA0 <- sum(pos.NA0)
  if (sum.NA0 == 1L)
  {
    warning("an element of the object 'inequality' is smaller than 0")
  }
  else if (sum.NA0 >= 1L)
  {
    warning(paste(sum.NA0, " elements of the object 'inequality' are smaller than 0"))
  }
  sum.NA1 <- sum(pos.NA1)
  if (sum.NA1 == 1L)
  {
    warning("an element of the object 'inequality' is larger than 1")
  }
  else if (sum.NA1 >= 1L)
  {
    warning(paste(sum.NA1, " elements of the object 'inequality' are larger than 1"))
  }
}
else
{
  output <- NA
  warning("object 'inequality' must be numeric")
}
return(output)
}
