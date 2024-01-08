\name{gparetoII}
\alias{gparetoII}
\title{Gini index for the Pareto (II) distribution with user-defined location, scale and shape parameters}
\usage{
gparetoII(
 location = 0,
 scale = 1,
 shape = 1
)
}
\arguments{
\item{location}{A positive real number specifying the location parameter \eqn{a} of the Pareto (II) distribution. The default value is \code{location = 0}.}
\item{scale}{A positive real number specifying the scale parameter \eqn{b} of the Pareto (II) distribution. The default value is \code{scale = 1}.}
\item{shape}{A positive real number specifying the shape parameter \eqn{s} of the Pareto (II) distribution. The default value is \code{shape = 1}.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when a parameter is non-numeric or positive, except the location parameter that can be equal to 0.
}
\description{
Calculates the Gini index for the Pareto (II) distribution with \code{location} parameter \eqn{a}, \code{scale} parameter \eqn{b} and \code{shape} parameter \eqn{s}.
}
\details{The Pareto (II) distribution with \code{location} parameter \eqn{a}, \code{scale} parameter \eqn{b}, \code{shape} parameter \eqn{s} and denoted as \eqn{ParetoII(a,b,s)}, where \eqn{a \geq 0}, \eqn{b>0} and \eqn{s>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y)= \displaystyle \frac{s}{b} \left[1 + \left( \frac{y-a}{b}\right)\right]^{-(s+1)},} and a cumulative distribution function given by
\deqn{F(y)=1-\left(1 + \displaystyle  \frac{y-a}{b} \right)^{-s},}
where \eqn{y>a}.

The Gini index can be computed as
\deqn{G = 2\left(0.5 - \displaystyle \frac{1}{E[y]}\int_{0}^{1}\int_{0}^{Q(y)}yf(y)dy\right),}
where \eqn{Q(y)} is the quantile function of the Pareto (II) distribution, and \eqn{E[y]} is the expectation of the distribution. If \code{location} is not specified it assumes the default value of 0, and \code{scale} and \code{shape} assume the default value of 1. The Pareto (II) distribution is related to the Pareto (IV) distribution: \eqn{ParetoII(a,b,s) = ParetoIV(a,b,1,s)}.
}

\examples{
# Gini index for the Pareto (II) distribution with parameters 'a = 1', 'b = 1' and 's = 3'.
gparetoII(location = 1, scale = 1, shape = 3)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gpareto}}, \code{\link{gparetoI}}, \code{\link{gparetoIII}}, \code{\link{gparetoIV}}, \code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gfisk}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
