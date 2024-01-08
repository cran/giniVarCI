\name{gdagum}
\alias{gdagum}
\title{Gini index for the Dagum distribution with user-defined shape parameters}
\usage{
gdagum(shape1.a, shape2.p)
}
\arguments{
\item{shape1.a}{A positive real number specifying the shape1 parameter \eqn{a} of the Dagum distribution.}
\item{shape2.p}{A positive real number specifying the shape parameter \eqn{p} of the Dagum distribution.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini index for the Dagum distribution with shape parameters \eqn{a} (\code{shape1.a}) and \eqn{p} (\code{shape2.p}).
}
\details{The Dagum distribution with scale parameter \eqn{b}, shape parameters \eqn{a} (argument \code{shape1.a}) and \eqn{p} (argument \code{shape2.p}) and denoted as \eqn{Dagum(b,a,p)} , where \eqn{b>0}, \eqn{a>0} and \eqn{p>0},
has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Rodriguez, 1977; Yee, 2022)
\deqn{f(y) = \displaystyle \frac{ap}{y}\frac{\left(\frac{y}{b}\right)^{ap}}{ \left[\left(\frac{y}{b} \right)^{a} + 1 \right]^{p+1} },}
and a cumulative distribution function given by
\deqn{F(y)= \left[1 + \displaystyle \left( \frac{y}{b}\right)^{-a} \right]^{-p},}
where \eqn{y > 0}.

The Gini index can be computed as
\deqn{G = \displaystyle \frac{\Gamma(p)\Gamma(2p+1/a)}{\Gamma(2p)\Gamma(p+1/a)}-1,}
where the gamma function is defined as
\deqn{\Gamma(\alpha) = \int_{0}^{\infty}t^{\alpha-1}e^{-t}dt.} The Dagum distribution is also known the Burr III, inverse Burr, beta-K, or 3-parameter kappa distribution. The Dagum distribution is related to the Fisk (Log Logistic) distribution: \eqn{Dagum(b,a,1) = Fisk(b,a)}. The Dagum distribution is also related to the inverse Lomax distribution and the inverse paralogistic distribution (see Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022).
}

\note{The Gini index of the Dagum distribution does not depend on its scale parameter.}

\examples{
# Gini index for the Dagum distribution with shape parameters 'a = 2' and 'p = 20'.
gdagum(shape1.a = 2, shape2.p = 20)

}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gburr}}, \code{\link{gpareto}}, \code{\link{gfisk}}, \code{\link{ggompertz}}, \code{\link{gfrechet}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
