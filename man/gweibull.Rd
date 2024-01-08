\name{gweibull}
\alias{gweibull}
\title{Gini index for the Weibull distribution with user-defined shape parameters}
\usage{
gweibull(shape)
}
\arguments{
\item{shape}{A vector of positive real numbers specifying shape parameters \eqn{a} of the Weibull distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculate the Gini indices for the Weibull distribution with \code{shape} parameters \eqn{a}.
}
\details{The Weibull distribution with scale parameter \eqn{\sigma}, \code{shape} parameter \eqn{a}, and denoted as \eqn{Weibull(\sigma, a)}, where \eqn{\sigma>0} and \eqn{a>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y) = \displaystyle \frac{a}{\sigma}\left(\frac{y}{\sigma}\right)^{a-1}e^{-(y/\sigma)^{a}},} and a cumulative distribution function given by
\deqn{F(y) = \displaystyle 1 - e^{-(y/\sigma)^{a}},}
where \eqn{y \geq 0}.

The Gini index can be computed as
\deqn{G = 1-2^{-1/a}.}
}

\note{The Gini index of the Weibull distribution does not depend on its scale parameter.}

\examples{
# Gini index for the Weibull distribution with 'shape = 1'.
gweibull(shape = 1)

# Gini indices for the Weibull distribution and different shape parameters.
gweibull(shape = 1:10)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{gbeta}}, \code{\link{ggamma}}, \code{\link{gchisq}}, \code{\link{gunif}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
