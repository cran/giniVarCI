\name{gparetoI}
\alias{gparetoI}
\title{Gini index for the Pareto (I) distribution with user-defined scale and shape parameters}
\usage{
gparetoI(
 scale = 1,
 shape = 1
)
}
\arguments{
\item{scale}{A positive real number specifying the scale parameter \eqn{b} of the Pareto (I) distribution. The default value is \code{scale = 1}.}
\item{shape}{A positive real number specifying the shape parameter \eqn{s} of the Pareto (I) distribution. The default value is \code{shape = 1}.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when a parameter is non-numeric or non-positive.
}
\description{
Calculate the Gini index for the Pareto (I) distribution with \code{scale} parameter \eqn{b} and \code{shape} parameter \eqn{s}.
}
\details{The Pareto (I) distribution with \code{scale} parameter \eqn{b}, \code{shape} parameter \code{s} and denoted as \code{ParetoI(b,s)}, where \eqn{b>0} and \eqn{s>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y)= \displaystyle \frac{s}{b} \left(\frac{y}{b}\right)^{-(s+1)},} and a
cumulative distribution function given by
\deqn{F(y)=1 - \displaystyle \left(\frac{y}{b}\right)^{-s},}
where \eqn{y>b}.

The Gini index can be computed as
\deqn{G = 2\left(0.5 - \displaystyle \frac{1}{E[y]}\int_{0}^{1}\int_{0}^{Q(y)}yf(y)dy\right),}
where \eqn{Q(y)} is the quantile function of the Pareto (I) distribution, and \eqn{E[y]} is the expectation of the distribution. If \code{scale} or \code{shape} are not specified they assume the default value of 1. The Pareto (I) distribution is related to the Pareto (IV) distribution: \eqn{ParetoI(b,s) = ParetoIV(b,b,1,s)}
}

\examples{
# Gini index for the Pareto (I) distribution with scale 'b = 1' and shape 's = 3'.
gparetoI(scale = 1, shape = 3)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gpareto}}, \code{\link{gparetoII}}, \code{\link{gparetoIII}}, \code{\link{gparetoIV}}, \code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gfisk}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
