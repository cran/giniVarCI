\name{gchisq}
\alias{gchisq}
\title{Gini index for the Chi-Squared distribution with user-defined degrees of freedom}
\usage{
gchisq(df)
}
\arguments{
\item{df}{A vector of positive real numbers specifying degrees of freedom of the Chi-Squared distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when degrees of freedom are non-numeric or non-positive.
}
\description{
Calculates Gini indices for the Chi-Squared distribution with degrees of freedom \eqn{n} (\code{df}).
}
\details{The Chi-Squared distribution with degrees of freedom \eqn{n} (argument \code{df}) and denoted as \eqn{\chi_{n}^2}, where \eqn{n>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995)
\deqn{f(y)= \displaystyle \frac{1}{2^{n/2}\Gamma\left(\frac{n}{2}\right)}y^{n/2-1}e^{-y/2},}
and a cumulative distribution function given by
\deqn{F(y) = \frac{\gamma\left(\frac{n}{2}, \frac{y}{2}\right)}{\Gamma(\alpha)},}
where \eqn{y \geq 0}, the gamma function is defined by
\deqn{\Gamma(\alpha) = \int_{0}^{\infty}t^{\alpha-1}e^{-t}dt,}
and the lower incomplete gamma function is given by
\deqn{\gamma(\alpha,y) = \int_{0}^{y}t^{\alpha-1}e^{-t}dt.}

The Gini index can be computed as
\deqn{G=\displaystyle \frac{2\Gamma\left( \frac{1+n}{2}\right)}{n\Gamma\left(\frac{n}{2}\right)\sqrt{\pi}}.} The Chi-Squared distribution is related to the Gamma distribution: \eqn{\chi_{n}^2 = Gamma(n/2, 2)}.
}

\examples{
# Gini index for the Chi-Squared distribution with degrees of freedom equal to 2.
gchisq(df = 2)

# Gini indices for the Chi-Squared distribution and different degrees of freedom.
gchisq(df = 5:10)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{ggamma}}, \code{\link{gf}}, \code{\link{gbeta}}, \code{\link{glnorm}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
