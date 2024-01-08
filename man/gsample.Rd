\name{gsample}
\alias{gsample}
\title{Samples from a set of continuous probability distributions with user-defined Gini indices}
\usage{
gsample(
  n,
  gini,
  distribution = c("pareto", "dagum", "lognormal", "fisk", "weibull", "gamma",
  "chisq", "frechet"),
  scale = 1,
  meanlog = 0,
  shape2.p = 1,
  location = 0
)
}
\arguments{
\item{n}{An integer specifying the sample(s) size.}

\item{gini}{A numeric vector of values between 0 and 1, indicating the Gini indices for the continuous distribution from which samples are generated.}

\item{distribution}{A character string specifying the continuous probability distribution to be used to generate the sample. Possible values are \code{"pareto"}, \code{"dagum"}, \code{"lognormal"}, \code{"fisk"}, \code{"weibull"}, \code{"gamma"}, \code{"chisq"} and \code{"frechet"} for the Pareto, Dagum, logNormal, Fisk (Log-logistic), Weibull, Gamma, Chi-Squared and Frechet distributions, respectively.}

\item{scale}{The scale parameter for the Pareto, Dagum, Fisk, Weibull, Gamma and Frechet distributions. The default value is \code{scale = 1}.}

\item{meanlog}{The mean for the logNormal distribution on the log scale. The default value is \code{meanlog = 0}.}

\item{shape2.p}{The scale parameter \code{p} for the Dagum distribution. The default value is \code{shape2.p = 1}.}

\item{location}{The location parameter for the Frechet distribution. The default value is \code{location = 0}.}

}
\value{
A numeric vector (or matrix of order \eqn{n} \eqn{\times} size(\eqn{gini})) with the samples by columns extracted from the continuous probability distribution stated in \code{distribution} and the Gini indices corresponding to the vector \code{gini}.
}
\description{
Draws samples from a continuous probability distribution with Gini indices set by the user.
}
\details{
For each continuous probability distribution, parameters involved in the theoretical formulation of the Gini index (\eqn{G}) are selected such that \eqn{G} takes the values set in the argument  \code{gini}. Additional parameters required in the distribution can be set by the user, and default values are provided. \code{scale} is the scale parameter for the Pareto, Dagum, Fisk, Weibull, Gamma and Frechet distributions, \code{meanlog} is the mean for the Lognormal distribution on the log scale, \code{shape2.p} is the scale parameter \code{p}  for the Dagum distribution, and \code{location} is the location parameter for the Frechet distribution. Additional information for the continuous probability distributions used by this function can be seen in Kleiber and Kotz (2003), Johnson et al. (1995) and Yee (2022).
}
\note{
Underestimation problems may appear for large heavy-tailed distributions (Pareto, Dagum, Lognormal, Fisk and Frechet) and large values of \code{gini}. A larger sample size may solve/minimize this problem.
}

\examples{
# Sample from the Pareto distribution and parameter selected such that the Gini index is 0.3.
gsample(n = 10, gini = 0.3, "pareto")

# Samples from the Pareto distribution and gini indices 0.2 and 0.5.
gsample(n = 10, gini = c(0.2,0.5), "par", scale = 2)

# Samples from the Lognormal distribution and gini indices 0.2 and 0.5.
gsample(n = 10, gini = c(0.2,0.5), "lognormal", meanlog = 5)

# Samples from the Dagum distribution and gini indices 0.2 and 0.5.
gsample(n = 10, gini = c(0.2,0.5), "dagum")

# Samples from the Fisk (Log-logistic) distribution and gini indices 0.3 and 0.6.
gsample(n = 10, gini = c(0.3,0.6), "fisk")

# Sample from the Weibull distribution and parameter selected such that the Gini index is 0.2.
gsample(n = 10, gini = 0.2, "weibull")

# Sample from the Gamma distribution and parameter selected such that the Gini index is 0.3.
gsample(n = 10, gini = 0.2, "gamma")

# Samples from the Chi-Squared distribution and gini indices 0.3 and 0.6..
gsample(n = 10, gini = c(0.3,0.6), "chi")

# Samples from the Frechet distribution and gini indices 0.3 and 0.6.
gsample(n = 10, gini = c(0.3,0.6), "fre")
}

\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.
}
\seealso{
\code{\link{gpareto}}, \code{\link{gdagum}}, \code{\link{glnorm}}, \code{\link{gfisk}}, \code{\link{gweibull}}, \code{\link{ggamma}}, \code{\link{gchisq}}, \code{\link{gfrechet}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}

