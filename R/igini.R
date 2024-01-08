igini <- function(y, bias.correction = TRUE, interval = NULL, B = 1000L, alpha = 0.05, cum.sums = NULL, na.rm = TRUE, precisionEL=1e-4, maxiterEL=100L, large.sample = FALSE)
{
  if (missing(y) && is.null(cum.sums) )  stop("object 'y' or 'cum.sums' must be specified")
  if (!missing(y) && !is.numeric(y))       stop("object 'y' must be numeric")
  if (!is.null(cum.sums) && !is.numeric(cum.sums))    stop("object 'cum.sums' must be numeric")
  if (na.rm == TRUE)
  {
    if (!missing(y))
    {
      y <- y[!is.na(y)]
      if (large.sample == FALSE)
      {
        y.sort <- Sort(y)
      }
      else
      {
        y.sort <- sort(y)
      }
    }
    if (!is.null(cum.sums))
    {
      cum.sums <- cum.sums[!is.na(cum.sums)]
      if (large.sample == FALSE)
      {
        y.sort.CS <- Sort(c(cum.sums[1L], diff(cum.sums)))
      }
      else
      {
        y.sort.CS <- sort(c(cum.sums[1L], diff(cum.sums)))
      }
    }
    if (!missing(y) && !is.null(cum.sums) )
    {
      if (length(y)!=length(cum.sums))                 stop("specify 'y' or 'cum.sums' but not both")
      if ( sum(abs(y.sort - y.sort.CS))  < 1e-15)   warning("specify 'y' or 'cum.sums' but not both")
      else                                             stop("specify 'y' or 'cum.sums' but not both")
    }
    if (!is.null(cum.sums)) y.sort <- y.sort.CS

   if (y.sort[1L] >= 0)
    {
      RES <- giniMeansRcpp(y.sort)
      # Gini = Unbiased version
      Gini <- 2 * RES$Meaniy / (RES$Mean * (RES$Sample.size - 1)) - (RES$Sample.size + 1 ) / (RES$Sample.size - 1 )
      BC <- (RES$Sample.size - 1) / RES$Sample.size  # Transformation to obtain the biased version
      GiniBiased <- BC * Gini
      output <- ifelse (bias.correction,  Gini, GiniBiased)
      # Gini index + Confidence interval
      if (!is.null(interval)) {
        output <- iinterval(y.sort, interval = interval, alpha, GiniBiased, RES$Sample.size, RES$Mean*RES$Sample.size, RES$Meaniy*RES$Sample.size, B, BC, bias.correction, precisionEL, maxiterEL)
      } # is.null(interval)
    } # (y.sort[1] >= 0L)
    else
    {
      output <- NA
      warning("sample values must be non-negative")
    }
  } # if (na.rm)
  else output <- NA
  return(output)
}
