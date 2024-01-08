\name{gfisk}
\alias{gfisk}
\title{Gini index for the Fisk (Log Logistic) distribution with user-defined shape parameters}
\usage{
gfisk(shape1.a)
}
\arguments{
\item{shape1.a}{A vector of positive real numbers specifying shape parameters \eqn{a} of the Fisk (Log Logistic) distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini indices for the Fisk (Log Logistic) distribution with shape parameters \eqn{a} (\code{shape1.a}).
}
\details{The Fisk (Log Logistic) distribution with scale parameter \eqn{b}, shape parameter \eqn{a} (argument \code{shape1.a}) and denoted as \eqn{Fisk(b,a)}, where \eqn{b>0} and \eqn{a>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y) = \displaystyle \frac{a}{y}\frac{\left(\frac{y}{b}\right)^{a}}{ \left[\left(\frac{y}{b} \right)^{a} + 1 \right]^{2} },} and a cumulative distribution function given by
\deqn{F(y)=1-\left[1 + \displaystyle \left( \frac{y}{b}\right)^{a} \right]^{-1},}
where \eqn{y \geq 0}.

The Gini index can be computed as
\deqn{G = \left\{
   \begin{array}{cl}
    1 , & 0< a <1; \\
     \displaystyle \frac{1}{a}, & a \geq 1.
   \end{array}
 \right.
} The Fisk (Log Logistic) distribution is related to the Dagum distribution: \eqn{Fisk(b,a) = Dagum(b,a,1)}.
}

\note{The Gini index of the Fisk (Log Logistic) distribution does not depend on its scale parameter.}

\examples{
# Gini index for the Fisk distribution with a shape parameter 'a = 2'.
gfisk(shape1.a = 2)

# Gini indices for the Fisk distribution and different shape parameters.
gfisk(shape1.a = 1:10)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.
}
\seealso{
\code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gpareto}}, \code{\link{ggompertz}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
