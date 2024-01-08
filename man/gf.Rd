\name{gf}
\alias{gf}
\title{Gini index for the F distribution with user-defined degrees of freedom}
\usage{
gf(df1, df2)
}
\arguments{
\item{df1}{A positive real number specifying the degrees of freedom \eqn{\nu_1} of the F distribution.}
\item{df2}{A positive real number higher or equal than two specifying the degrees of freedom \eqn{\nu_2} of the F distribution.}
}
\value{
A numeric value with the Gini index. A \code{NA} is returned when degrees of freedom are non-numeric or \eqn{df1 \leq 0} or \eqn{df2 < 2} .
}
\description{
Calculates the Gini index for the F distribution with degrees of freedom \eqn{\nu_1} (\code{df1}) and \eqn{\nu_2} (\code{df2}).
}
\details{The F distribution with \eqn{\nu_1} (argument \code{df1}) and \eqn{\nu_2} (argument \code{df2}) degrees of freedom and denoted as \eqn{F_{\nu_1,\nu_2}}, where \eqn{\nu_1>0} and \eqn{\nu_2 > 0}, has a probability density function given by (Kleiber and Kotz, 2003; Johnson et al., 1995)
\deqn{f(y) = \displaystyle \frac{\Gamma\left(\frac{\nu_{1}}{2} + \frac{\nu_{2}}{2}\right)}{\Gamma\left(\frac{\nu_{1}}{2}\right)\Gamma\left(\frac{\nu_{2}}{2}\right)}\left( \frac{\nu_{1}}{\nu_{2}}\right)^{\nu_{1}/2}y^{\nu_{1}/2-1}\left(1 + \frac{\nu_{1}y}{\nu_{2}}\right)^{-(\nu_{1}+\nu_{2})/2},}
and a cumulative distribution function given by
\deqn{F(y)= \displaystyle I_{\nu_{1}y/(\nu_{1}y + \nu_{2})}\left( \frac{\nu_{1}}{2}, \frac{\nu_{2}}{2} \right),}
where \eqn{y \geq 0},
\deqn{\Gamma(\alpha) = \int_{0}^{\infty}t^{\alpha-1}e^{-t}dt}
is the gamma function,
\deqn{I_{y}(a,b)=\displaystyle \frac{B(y;a,b)}{B(a,b)}}
is the regularized incomplete beta function,
\deqn{B(a,b) = \displaystyle \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}}
is the beta function,
and
\deqn{B(y;a,b) = \displaystyle  \int_{0}^{y}t^{a-1}(1-t)^{b-1}dt} is the incomplete beta function.

The Gini index, for \eqn{\nu_2 \geq 2}, can be computed as
\deqn{G = 2\left(0.5 - \displaystyle \frac{\nu_{2} - 2}{ \nu_{2}}\int_{0}^{1}\int_{0}^{Q(y)}yf(y)dy\right),}
where \eqn{Q(y)} is the quantile function of the F distribution.
}

\examples{
# Gini index for the F distribution with 'df1 = 10' and 'df2 = 20' degrees of freedom.
gf(df1 = 10, df2 = 20)

}
\references{
Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) Continuous Univariate Distributions, volume 1, chapter 14. Wiley, New York.

}
\seealso{
\code{\link{gchisq}}, \code{\link{ggamma}}, \code{\link{ggompertz}}, \code{\link{glnorm}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
