\name{gparetoIII}
\alias{gparetoIII}

\title{Gini index for the Pareto (III) distribution with user-defined inequality parameters}
\usage{
gparetoIII(
 inequality = 1
)
}
\arguments{
\item{inequality}{A vector of positive numbers in the \eqn{[0,1]} interval specifying inequality parameters \eqn{g} of the Pareto (III) distribution. The default value is \code{inequality = 1}.}
}

\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a inequality parameter is non-numeric or it is out of the interval \eqn{[0,1]}.
}
\description{
Calculate the Gini index for the Pareto (III) distribution with \code{inequality} parameters \eqn{g}.
}
\details{The Pareto (III) distribution with location parameter \eqn{a}, scale parameter \eqn{b}, \code{inequality} parameter \code{g} and denoted as \eqn{ParetoIII(a,b,g)}, where \eqn{a>0}, \eqn{b>0}, and \eqn{g \in [0,1]}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y)= \displaystyle \frac{1}{bg} \left( \frac{y-a}{b}\right)^{1/g-1} \left[1 + \left( \frac{y-a}{b}\right)^{1/g} \right]^{-2},} and a cumulative distribution function given by
\deqn{F(y)=1-\left[1 + \displaystyle \left( \frac{y-a}{b}\right)^{1/g} \right]^{-1},}
where \eqn{y>a}.

The Gini index is \eqn{G = g.}

If \code{inequality} is not specified it assumes the default value of 1. The Pareto (III) distribution is related to the Pareto (IV) distribution: \eqn{ParetoIII(a,b,g) = ParetoIV(a,b,g,1)}.}

\note{The Gini index of the Pareto (III) distribution does not depend on its location and scale parameters.}

\examples{
# Gini index for the Pareto (III) distribution with inequality parameter 'g = 0.3'.
gparetoIII(inequality = 0.3)

# Gini indices for the Pareto (III) distribution with different inequality parameters.
gparetoIII(inequality = seq(0.1, 0.9, by=0.1))

}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gpareto}}, \code{\link{gparetoI}}, \code{\link{gparetoII}}, \code{\link{gparetoIV}}, \code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gfisk}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
