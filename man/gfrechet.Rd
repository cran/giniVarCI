\name{gfrechet}
\alias{gfrechet}
\title{Gini index for the Frechet distribution with user-defined shape parameters}
\usage{
gfrechet(shape)
}
\arguments{
\item{shape}{A vector of positive real numbers higher or equal than 1 specifying shape parameters \eqn{s} of the Frechet distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a shape parameter is non-numeric or smaller than 1.
}
\description{
Calculates the Gini indices for the Frechet distribution with \code{shape} parameters \eqn{s}.
}
\details{The Frechet distribution with location parameter \eqn{a}, scale parameter \eqn{b},  \code{shape} parameter \eqn{s} and denoted as \eqn{Frechet(a,b,s)}, where \eqn{a>0}, \eqn{b>0} and \eqn{s>0}, has a
probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995)
\deqn{f(y) = \displaystyle \frac{sb}{(y-a)^{2}} \left(\frac{b}{y-a}\right)^{s-1} \exp\left[- \displaystyle \left(\frac{b}{y-a}\right)^{s} \right],}
and a cumulative distribution function given by
\deqn{F(y)= \displaystyle \exp\left[- \displaystyle \left(\frac{b}{y-a}\right)^{s} \right],}
where \eqn{y > a}.

The Gini index, for \eqn{s \geq 1}, can be computed as
\deqn{G = 2^{1/s} -1.}
}

\note{The Gini index of the Frechet distribution does not depend on its location and scale parameters and only is defined when its shape parameter is at least 1.}

\examples{
# Gini index for the Frechet distribution with a shape parameter 's = 1'.
gfrechet(shape = 1)

# Gini indices for the Frechet distribution and different shape parameters.
gfrechet(shape = 1:10)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gfisk}}, \code{\link{gpareto}}, \code{\link{ggompertz}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
