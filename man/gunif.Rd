\name{gunif}
\alias{gunif}
\title{Gini index for the Uniform distribution with user-defined lower and upper limits}
\usage{
gunif(
 min = 0,
 max = 1
)
}
\arguments{
\item{min}{A non-negative real number specifying the lower limit of the Uniform distribution. The default value is \code{min = 0}.}
\item{max}{A positive real number higher than \code{min} specifying the upper limit of the Uniform distribution. The default value is \code{max = 1}.}
}
\value{
A numeric value with the Gini index. A \code{NA} value is returned when a limit is non-numeric or non-negative, or \eqn{\min \geq \max}.
}
\description{
Calculates the Gini index for the Uniform distribution with lower limit \code{min} and upper limit \code{max}.
}
\details{The Uniform distribution with lower and upper limits \eqn{min}  and \eqn{max}, and denoted as \eqn{U(min,max)}, where \eqn{\min \geq 0}, \eqn{\max >0}, \eqn{\min < \max} and both must be finite, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995; Yee, 2022)
\deqn{f(y)= \displaystyle \frac{1}{\max - \min},} where \eqn{y \in [\min, \max]}. The cumulative distribution function is given by
\deqn{F(y) = \left\{
   \begin{array}{cl}
    0 , &  y < \min; \\
     \displaystyle \frac{y-\min}{\max - \min}, & y \in [\min, \max]; \\
    1 , &  y > \max.
   \end{array}
 \right.
}
The Gini index can be computed as
\deqn{G = \displaystyle \frac{\max - \min}{3(\min + \max)}.}
If \code{min} or \code{max} are not specified they assume the default values of 0 and 1, respectively.
}

\examples{
# Gini index for the Uniform distribution with lower limit 0 and upper limit 1.
gunif()

# Gini index for the Uniform distribution with lower limit 10 and upper limit 190.
gunif(min = 10, max = 190)


}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

Yee, T. W. (2022). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-7, https://CRAN.R-project.org/package=VGAM.

}
\seealso{
\code{\link{gbeta}}, \code{\link{ggamma}}, \code{\link{gchisq}}, \code{\link{gf}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
