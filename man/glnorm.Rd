\name{glnorm}
\alias{glnorm}
\title{Gini index for the Log Normal distribution with user-defined standard deviations}
\usage{
glnorm(sdlog)
}
\arguments{
\item{sdlog}{A vector of positive real numbers specifying standard deviations \eqn{\sigma} of the Log Normal distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a standard deviation is non-numeric or non-positive.
}
\description{
Calculates the Gini indices for the Log Normal distribution with standard deviations \eqn{\sigma} (\code{sdlog}).
}
\details{The Log Normal distribution with mean \eqn{\mu}, standard deviation \eqn{\sigma} on the log scale (argument \code{sdlog}) and denoted as \eqn{logNormal(\mu, \sigma)}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995)
\deqn{f(y)=\displaystyle \frac{1}{\sqrt{2\pi}\sigma y}\exp\left[- \frac{(\ln(x) - \mu)^2}{2\sigma^2} \right],}
and a cumulative distribution function given by
\deqn{F(y)=\displaystyle \Phi\left(\frac{\ln(x) - \mu}{\sigma}\right),}
where \eqn{y > 0} and
\deqn{\Phi(y) = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^{y} e^{-t^{2}/2}dt}
is the cumulative distribution function of a standard Normal distribution.

The Gini index can be computed as
\deqn{G = 2\Phi\left( \displaystyle \frac{\sigma}{\sqrt{2}}\right) - 1.}
}

\note{The Gini index of the logNormal distribution does not depend on the mean parameter.}

\examples{
# Gini index for the Log Normal distribution with standard deviation 'sdlog = 2'.
glnorm(sdlog = 2)

# Gini indices for the Log Normal distribution with different standard deviations.
glnorm(sdlog = c(0.2, 0.5, 1:3))
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{ggamma}}, \code{\link{gpareto}}, \code{\link{gchisq}}, \code{\link{gweibull}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
