\name{gbeta}
\alias{gbeta}
\title{Gini index for the Beta distribution with user-defined shape parameters}
\usage{
gbeta(shape1, shape2)
}
\arguments{
\item{shape1}{A positive real number specifying the shape1 parameter \eqn{a} of the Beta distribution.}
\item{shape2}{A positive real number specifying the shape2 parameter \eqn{b} of the Beta distribution.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini index for the Beta distribution with shape parameters \eqn{a} (\code{shape1}) and  \eqn{b} (\code{shape2}).
}
\details{The Beta distribution with shape parameters \eqn{a} (argument \code{shape1}) and \eqn{b} (argument \code{shape2}) and denoted as \eqn{Beta(a,b)}, where \eqn{a>0} and \eqn{b>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y) = \displaystyle \frac{1}{B(a,b)}y^{a-1}(1-y)^{b-1},}
and a cumulative distribution function given by
\deqn{F(y)= \displaystyle \frac{B(y;a,b)}{B(a,b)} }
where \eqn{0 \leq y \leq 1},
\deqn{B(a,b) = \displaystyle \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}}
is the beta function,
\deqn{\Gamma(\alpha) = \int_{0}^{\infty}t^{\alpha-1}e^{-t}dt} is the gamma function, and
\deqn{B(y;a,b) = \displaystyle  \int_{0}^{y}t^{a-1}(1-t)^{b-1}dt} is the incomplete beta function.

The Gini index can be computed as
\deqn{G = \displaystyle \frac{2}{a}\frac{B(a+b,a+b)}{B(a,a)B(b,b)}.}
}

\examples{
# Gini index for the Beta distribution with shape parameters 'a = 2' and 'b = 1'.
gbeta(shape1 = 2, shape2 = 1)

# Gini index for the Beta distribution with shape parameters 'a = 1' and 'b = 2'.
gbeta(shape1 = 1, shape2 = 2)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gf}}, \code{\link{gunif}}, \code{\link{gweibull}}, \code{\link{ggamma}}, \code{\link{gchisq}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
