\name{gburr}
\alias{gburr}
\title{Gini index for the Burr Type XII (Singh-Maddala) distribution with user-defined scale and shape parameters}
\usage{
gburr(
 scale = 1,
 shape.g = 1,
 shape.s = 1
)
}
\arguments{
\item{scale}{A positive real number specifying the scale parameter \eqn{b} of the Burr Type XII (Singh-Maddala) distribution. The default value is \code{scale = 1}.}
\item{shape.g}{A positive real number specifying the shape parameter \eqn{g} of the Burr Type XII (Singh-Maddala) distribution. The default value is \code{shape.g = 1}.}
\item{shape.s}{A positive real number specifying the shape parameter \eqn{s} of the Burr Type XII (Singh-Maddala) distribution. The default value is \code{shape.s = 1}.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when any of the parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini index for the Burr Type XII (Singh-Maddala) distribution with  \code{scale} parameter \eqn{b} and shape parameters \eqn{g} (\code{shape.g}) and \eqn{s} (\code{shape.s}).
}
\details{The Burr Type XII (Singh-Maddala) distribution with \code{scale} parameter \eqn{b}, shape parameters \eqn{g} (argument \code{shape.g}) and \eqn{s} (argument \code{shape.s}) and denoted as \eqn{BurrXII(b,g,s)}, where \eqn{b>0}, \eqn{g>0} and \eqn{s>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Rodriguez, 1977; Yee, 2022)
\deqn{f(y) = \displaystyle \frac{gs}{b}\left(\frac{y}{b}\right)^{g-1}\left[1 + \left(\frac{y}{b}\right)^{g}\right]^{-(s+1)},}
and a cumulative distribution function given by
\deqn{F(y)=1-\left[1 + \displaystyle \left( \frac{y}{b}\right)^{g} \right]^{-s},}
where \eqn{y>0}.

The Gini index can be computed as
\deqn{G = 2\left(0.5 - \displaystyle \frac{1}{E[y]}\int_{0}^{1}\int_{0}^{Q(y)}yf(y)dy\right),}
where
\eqn{Q(y)} is the quantile function of the Burr Type XII (Singh-Maddala) distribution, and \eqn{E[y]} is the expectation of the distribution. The Burr Type XII (Singh-Maddala) distribution is related to the Pareto (IV) distribution: \eqn{BurrXII(b,g,s) = ParetoIV(0,b,1/g,s)}.}

\examples{
# Gini index for the Burr Type XII distribution with 'scale = 1', 'shape.g = 2', 'shape.s = 1'.
gburr(scale = 1, shape.g = 2, shape.s = 1)

# Gini index for the Burr Type XII distribution with 'scale = 1', 'shape.g = 5', 'shape.s = 3'.
gburr(scale = 1, shape.g = 5, shape.s = 3)

}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Rodriguez, R. N. (1977). A guide to the Burr type XII distributions. Biometrika, 64(1), 129-134.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gparetoIV}}, \code{\link{gpareto}}, \code{\link{gparetoI}}, \code{\link{gparetoII}}, \code{\link{gparetoIII}}, \code{\link{gfisk}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
