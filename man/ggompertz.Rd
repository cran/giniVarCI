\name{ggompertz}
\alias{ggompertz}
\title{Gini index for the Gompertz distribution with user-defined scale and shape parameters}
\usage{
ggompertz(
 scale = 1,
 shape
)
}
\arguments{
\item{scale}{A positive real number specifying the scale parameter \eqn{\beta} of the Gompertz distribution. The default value is \code{scale = 1}.}
\item{shape}{A positive real number specifying the shape parameter \eqn{\alpha} of the Gompertz distribution.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when a parameter is non-numeric or non-positive.
}
\description{
Calculate the Gini index for the Gompertz distribution with \code{scale} parameter \eqn{\beta} and \code{shape} parameter \eqn{\alpha}.
}
\details{The Gompertz distribution with \code{scale} parameter \eqn{\beta}, \code{shape} parameter \eqn{\alpha} and denoted as \eqn{Gompertz(\beta, \alpha)}, where \eqn{\beta>0} and \eqn{\alpha>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Rodriguez, 1977; Yee, 2022)
\deqn{f(y)= \alpha e^{\beta y}  \exp\left[ - \displaystyle \frac{\alpha}{\beta}\left(e^{\beta y} - 1 \right) \right],} and a cumulative distribution function given by
\deqn{F(y)= 1 -\exp\left[ - \displaystyle \frac{\alpha}{\beta}\left(e^{\beta y} - 1 \right) \right],} where \eqn{y \geq 0}.

The Gini index can be computed as
\deqn{G = 2\left(0.5 - \displaystyle \frac{1}{E[y]}\int_{0}^{1}\int_{0}^{Q(y)}yf(y)dy\right),}
where \eqn{Q(y)} is the quantile function of the Gompertz distribution, and \eqn{E[y]} is the expectation of the distribution. If \code{scale} is not specified it assumes the default value of 1.
}

\examples{
# Gini index for the Gompertz distribution with 'scale = 1' and 'shape = 3'.
ggompertz(scale = 1, shape = 3)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{ggamma}}, \code{\link{gbeta}}, \code{\link{gchisq}}, \code{\link{gpareto}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
