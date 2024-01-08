\name{ggamma}
\alias{ggamma}
\title{Gini index for the Gamma distribution with user-defined shape parameter}
\usage{
ggamma(shape)
}
\arguments{
\item{shape}{A vector of positive real numbers specifying the shape parameters \eqn{\alpha} of the Gamma distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini indices for the Gamma distribution with \code{shape} parameters \eqn{\alpha}.
}
\details{The Gamma distribution with \code{shape} parameter  \eqn{\alpha}, scale parameter \eqn{\sigma} and denoted as \eqn{Gamma(\alpha, \sigma)}, where \eqn{\alpha>0} and \eqn{\sigma>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995)
\deqn{f(y) = \displaystyle \frac{1}{\sigma^{\alpha}\Gamma(\alpha)}y^{\alpha-1}e^{-y/\sigma},}
and a cumulative distribution function given by
\deqn{F(y) = \frac{\gamma\left(\alpha, \frac{y}{\sigma}\right)}{\Gamma(\alpha)},}
where \eqn{y \geq 0}, the gamma function is defined by
\deqn{\Gamma(\alpha) = \int_{0}^{\infty}t^{\alpha-1}e^{-t}dt,}
and the lower incomplete gamma function is given by
\deqn{\gamma(\alpha,y) = \int_{0}^{y}t^{\alpha-1}e^{-t}dt.}

The Gini index can be computed as
\deqn{G = \displaystyle \frac{\Gamma\left(\frac{2\alpha+1}{2}\right)}{\alpha\Gamma(\alpha)\sqrt{\pi}}.} The Gamma distribution is related to the Chi-squared distribution: \eqn{Gamma(n/2, 2) = \chi_{n}^2}.
}

\note{The Gini index of the Gamma distribution does not depend on its scale parameter.}

\examples{
# Gini index for the Gamma distribution with 'shape = 1'.
ggamma(shape = 1)

# Gini indices for the Gamma distribution and different shape parameters.
ggamma(shape = 1:10)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{gchisq}}, \code{\link{gf}}, \code{\link{gbeta}}, \code{\link{gweibull}}, \code{\link{glnorm}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
