\name{igini}
\alias{igini}
\title{Gini index, variances and confidence intervals in infinite populations}
\usage{
igini(
  y,
  bias.correction = TRUE,
  interval = NULL,
  B = 1000L,
  alpha = 0.05,
  cum.sums = NULL,
  na.rm = TRUE,
  precisionEL = 1e-04,
  maxiterEL = 100L,
  large.sample = FALSE
)
}
\arguments{
\item{y}{A vector with the non-negative real numbers to be used for estimating the Gini index. This argument can be missing if argument \code{cum.sums} is provided.}

\item{bias.correction}{A `TRUE/FALSE` logical value indicating whether the bias correction should be applied to the estimation of the Gini index. The default value is \code{bias.correction = TRUE}.}

\item{interval}{A character string specifying the type of variance estimation and confidence interval to be used, or \code{NULL} (the default value) to omit the computation of both variance and confidence interval. Possible values are \code{"zjackknife"}, \code{"tjackknife"}, \code{"zalinearization"}, \code{"zblinearization"}, \code{"talinearization"}, \code{"tblinearization"}, \code{"pbootstrap"}, \code{"BCa"}, \code{"ELchisq"} and \code{"ELboot"}. The default value is \code{interval = NULL}.}

\item{B}{A single integer specifying the number of bootstrap replicates. The default value is \code{B = 1000L}.}

\item{alpha}{A single numeric value between 0 and 1. If \code{interval} is not \code{NULL}, the confidence level to be used for computing  the confidence interval for the Gini is \code{1-alpha}. Some authors call \code{alpha} the significance level. The default value is \code{alpha = 0.05}.}

\item{cum.sums}{A vector with the non-negative real numbers specifying the cumulative sums of the variable used to estimate the Gini index. This argument can be \code{NULL} if argument \code{y} is provided. The default value is \code{cum.sums = NULL}.}

\item{na.rm}{A `TRUE/FALSE` logical value indicating whether \code{NA}'s should be removed before the computation proceeds. The default value is \code{na.rm = TRUE}.}

\item{precisionEL}{A single numeric value specifying the precision for the confidence interval based on the empirical likelihood method. The default value is \code{precisionEL = 1e-4}, i.e., limits of the confidence interval have a total of 4 decimal places.}

\item{maxiterEL}{A single integer specifying the maximal number of iterations allowed for the convergene of the empirical likelihood method. The default value is \code{maxiterEL = 100L}.}

\item{large.sample}{A `TRUE/FALSE` logical value indicating whether the sample is large to apply a faster algorithm to sort the sample values. The default value is \code{large.sample = FALSE}.}
}
\value{
When \code{interval = NULL}, a single numeric value between 0 and 1, containing the estimation of the Gini index based on the vector \code{y} or the vector \code{cum.sums}.
When \code{interval} is not \code{NULL}, a list of 3 components: a single numeric value with the estimation of the Gini index; a single numeric value with the variance estimation of the Gini index; and a numeric matrix with 1 row and 2 columns containing the lower and upper limits of the confidence intervals for the Gini index.
}
\description{
Estimation of the Gini index and computation of variances and confidence interval for infinite populations.
}
\details{
For a sample \eqn{S}, with size \eqn{n}, derived from an infinite population, the Gini index is estimated by

\deqn{\widehat{G}      = \displaystyle \frac{2}{\overline{y}n^{2}}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n} }
when \code{bias.correction = FALSE}, and by

\deqn{\widehat{G}^{bc} = \displaystyle \frac{2}{\overline{y}n(n-1)}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n-1} }
when \code{bias.correction = TRUE}. For more details, see Muñoz et al. (2023). The table below sumarises the various types of variances and confidence intervals that computes this function.
Methods based on the jackknife technique use the fast algorithm suggested by Ogwang (2000). The linearization technique for variance estimation (Deville, 1999) has been applied to the following estimators of the Gini index (Berger, 2008; Langel and Tille, 2013):
\deqn{\widehat{G}^{a} = \displaystyle \frac{1}{2\overline{y}n^{2}}\sum_{i \in S}\sum_{j\in S} |y_i-y_j|}
and
\deqn{\widehat{G}^{b}      = \displaystyle \frac{2}{\overline{y}n}\sum_{i \in S}y_{i}\widehat{F}_{n}(y_{i}) - 1,}
where
\deqn{\widehat{F}_{n}(y_i)=\frac{1}{n}\sum_{j \in S}\delta(y_j \leq y_i).}

\code{zalinearization} and \code{zblinearization} linearizate, respectively, the estimators \eqn{\widehat{G}^{a}} and \eqn{\widehat{G}^{b}}. The percentile bootstrap (see Qin et al., 2010) is computed using \code{pbootstrap}. \code{Bca} is the bias corrected bootstrap confidence interval (Efron and Tibshirani, 1993). \code{ELchisq} and \code{ELboot} are the confidence intervals based on the empirical likelihood method.
The \bold{vignette} \code{vignette("GiniVarInterval")} contains a detailed description of the various methods for variance estimation and confidence intervals for the Gini index.

\tabular{llll}{
\bold{Interval}          \tab \bold{Variance} \tab \bold{Critical values}  \tab \bold{References}\cr
_______________     \tab ____________  \tab __________________ \tab __________________________\cr
\code{zjackknife}        \tab Jackknife     \tab Normal                \tab Berger (2008)\cr
\code{tjackknife}        \tab Jackknife     \tab Studentized bootstrap \tab Biewen (2002); Berger (2008)\cr
\code{zalinearization}    \tab Linearization \tab Normal                \tab Langel and Tille (2013)\cr
\code{zblinearization}    \tab Linearization \tab Normal                \tab Berger (2008)\cr
\code{talinearization}    \tab Linearization \tab Studentized bootstrap \tab Langel and Tille (2013)\cr
\code{tblinearization}    \tab Linearization \tab Studentized bootstrap \tab Biewen (2002); Berger (2008)\cr
\code{pBootstrap}        \tab Bootstrap     \tab Percentile bootstrap  \tab Qin et al. (2010)\cr
\code{BCa}               \tab Bootstrap  \tab BCa bootstrap           \tab Davison and Hinkley (1997)\cr
\code{ELchisq}           \tab Linearization \tab Chi-Squared           \tab Qin et al. (2010)\cr
\code{ELboot}            \tab Bootstrap     \tab Percentile bootstrap  \tab Qin et al. (2010)\cr
}
}
\examples{
# Sample, with size 50, from a Lognormal distribution. The true Gini index is 0.5.
set.seed(123)
y <- gsample(n = 50, gini = 0.5, distribution = "lognormal")

# Bias corrected estimation of the Gini index.
igini(y)

# Estimation of the Gini index and confidence interval based on jackknife and studentized bootstrap.
igini(y, interval = "tjackknife")


}
\references{
Berger, Y. G. (2008). A note on the asymptotic equivalence of jackknife and linearization variance estimation for the Gini Coefficient. Journal of Official Statistics, 24(4), 541-555.

Biewen, M. (2002). Bootstrap inference for inequality, mobility and poverty measurement. Journal of Econometrics, 108(2), 317-342.

Davison, A. C., and Hinkley, D. V. (1997). Bootstrap Methods and Their Application (Cambridge Series in Statistical and Probabilistic Mathematics, No 1)–Cambridge University Press.

Deville, J.C. (1999). Variance Estimation for Complex Statistics and Estimators: Linearization and Residual Techniques. Survey Methodology, 25, 193--203.

Efron, B. and Tibshirani, R. (1993). An Introduction to the Bootstrap. Chapman and Hall, New York, London.

Langel, M., and Tille, Y. (2013). Variance estimation of the Gini index: revisiting a result several times published. Journal of the Royal Statistical Society: Series A (Statistics in Society), 176(2), 521-540.

Muñoz, J. F., Moya-Fernández, P. J., and Álvarez-Verdejo, E. (2023). Exploring and Correcting the Bias in the Estimation of the Gini Measure of Inequality. Sociological Methods & Research. https://doi.org/10.1177/00491241231176847

Ogwang, T. (2000). A convenient method of computing the Gini index and its standard error. Oxford Bulletin of Economics and Statistics, 62(1), 123-123.

Qin, Y., Rao, J. N. K., and Wu, C. (2010). Empirical likelihood confidence intervals for the Gini measure of income inequality. Economic Modelling, 27(6), 1429-1435.
}
\seealso{
\code{\link{icompareCI}}, \code{\link{iginindex}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
