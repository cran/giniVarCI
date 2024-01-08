gsample <- function(n, gini, distribution=c("pareto", "dagum", "lognormal", "fisk", "weibull",  "gamma", "chisq", "frechet"), scale = 1, meanlog = 0, shape2.p = 1, location = 0)
{
  if (!is.numeric(gini))   stop("object 'gini' must be numeric")
  if (!is.numeric(n))      stop("object 'n' must be numeric")
  if (length(n) > 1L)      stop("object 'n' must have length 1")
  distribution <- match.arg(distribution)
  COL <- length(gini)
  output <- matrix(nrow = n, ncol = COL, dimnames = list(c(), gini))
  for (i in 1L:COL)
  {
    if( (gini[i] < 0) || (gini[i] > 1)) stop("elements of the object 'gini' must be in the interval [0,1]")
    if((gini[i] == 0) && (distribution == "dagum"))   stop("the Dagum distribution requires a gini index larger than 0")
    if((gini[i] == 0) && (distribution == "weibull")) stop("the Weibull distribution requires a gini index larger than 0")
    if((gini[i] == 1) && (distribution == "weibull")) stop("the Weibull distribution requires a gini index smaller than 1")
    if((gini[i] == 1) && (distribution == "gamma"))   stop("the Gamma distribution requires a gini index smaller than 1")
    if((gini[i] == 0) && (distribution == "fisk"))    stop("the Fisk distribution requires a gini index larger than 0")
    if((gini[i] == 0) && (distribution == "frechet")) stop("the Frechet distribution requires a gini index larger than 0")
    output[,i] <- switch(distribution,
                         pareto    = VGAM::rpareto(n, scale, shape = (1 + gini[i])/(2*gini[i])),
                         lognormal = stats::rlnorm(n, meanlog = meanlog, sdlog = stats::qnorm( (gini[i] + 1)/2)*sqrt(2)),
                         dagum     = VGAM::rdagum(n, shape1.a = SampleDagum(shape2.p, gini[i]), shape2.p = shape2.p),
                         weibull   = stats::rweibull(n, shape= 1/log2(1/(1 - gini[i])), scale = scale),
                         gamma     = stats::rgamma(n, shape = SampleGamma(gini[i]), scale = scale),
                         chisq     = stats::rchisq(n, df = SampleChisq(gini[i])),
                         fisk      = VGAM::rfisk(n, scale, shape = 1/gini[i]),
                         frechet   = VGAM::rfrechet(n, location, scale, shape = 1/log2(gini[i] + 1)  )
    )
  }
  if (COL == 1L) output <- as.vector(output)
  return(output)
}


