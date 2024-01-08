\name{gpareto}
\alias{gpareto}
\title{Gini index for the Pareto distribution with user-defined shape parameters}
\usage{
gpareto(shape)
}
\arguments{
\item{shape}{A vector of positive real numbers specifying shape parameters \eqn{\alpha} of the Pareto distribution.}
}
\value{
A numeric vector with the Gini indices. A \code{NA} is returned when a shape parameter is non-numeric or non-positive.
}
\description{
Calculates the Gini indices for the Pareto distribution with \code{shape} parameters \eqn{\alpha}.
}
\details{The Pareto distribution with scale parameter \eqn{k}, \code{shape} parameter \eqn{\alpha} and denoted as \eqn{Pareto(k, \alpha)}, where \eqn{k>0} and \eqn{\alpha>0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y)=\displaystyle \frac{\alpha k^{\alpha}}{y^{\alpha +1}}, }
and a cumulative distribution function given by
\deqn{F(y) = \displaystyle 1 - \left(\frac{k}{y}\right)^{\alpha},}
where \eqn{y \geq k}.

The Gini index can be computed as
\deqn{G = \left\{
   \begin{array}{cl}
    1 , & 0<\alpha <1; \\
     \displaystyle \frac{1}{2\alpha-1}, & \alpha \geq 1.
   \end{array}
 \right.
}
}

\note{The Gini index of the Pareto distribution does not depend on the shape parameter.}

\examples{
# Gini index for the Pareto distribution with 'shape = 2'.
gpareto(shape = 2)

# Gini indices for the Pareto distribution and different shape parameters.
gpareto(shape = 1:5)
}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.
}
\seealso{
\code{\link{gparetoI}}, \code{\link{gparetoII}}, \code{\link{gparetoIII}}, \code{\link{gparetoIV}}, \code{\link{gdagum}}, \code{\link{gburr}}, \code{\link{gfisk}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
