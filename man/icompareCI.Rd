\name{icompareCI}
\alias{icompareCI}
\title{Comparisons of variance estimators and confidence intervals for the Gini index in infinite populations}
\usage{
icompareCI(
 y,
 B = 1000L,
 alpha = 0.05,
 plotCI = TRUE,
 digitsgini = 2L,
 digitsvar = 4L,
 cum.sums = NULL,
 na.rm = TRUE,
 precisionEL = 1e-4,
 maxiterEL = 100L,
 line.types = c(1L, 2L),
 colors = c("red", "green"),
 save.plot = FALSE
)
}
\arguments{
\item{y}{A vector with the non-negative real numbers to be used for estimating the Gini index. This argument can be missing if argument \code{cum.sums} is provided.}

\item{B}{A single integer specifying the number of bootstrap replicates. The default value is \code{B = 1000L}.}

\item{alpha}{A single numeric value between 0 and 1 specifying the confidence level 1-\code{alpha} to be used for computing the confidence interval for the Gini. Some authors call \code{alpha} the significance level. The default value is \code{alpha = 0.05}.}

\item{plotCI}{A `TRUE/FALSE` logical value indicating whether confidence intervals are compared using a plot. The default value is \code{plotCI = TRUE}.}

\item{digitsgini}{A single integer specifying the number of decimals used in the estimation of the Gini index and confidence intervals. The default value is \code{digitsgini = 2L}.}

\item{digitsvar}{A single integer specifying the number of decimals used in the variance estimation of the Gini index. The default value is \code{digitsvar = 4L}.}

\item{cum.sums}{A numeric vector of non-negative real numbers specifying the cumulative sums of the variable used to estimate the Gini index. This argument can be \code{NULL} if argument \code{y} is provided. The default value is \code{cum.sums = NULL}.}

\item{na.rm}{A `TRUE/FALSE` logical value indicating whether the \code{NA} should be removed before the computation proceeds. The default value is \code{na.rm = TRUE}.}

\item{precisionEL}{A single numeric value specifying the precision for the confidence interval based on the empirical likelihood method. The default value is \code{precisionEL = 1e-4}, i.e., limits of the confidence interval have a total of 4 decimal places.}

\item{maxiterEL}{A single integer specifying the maximum number of iterations allowed for the convergence in the empirical likelihood method. The default value is \code{maxiterEL = 100L}.}

\item{line.types}{A numeric vector with length equal 2 specifying the line types. See the function \code{plot} for the different line types. The default value is \code{lty = c(1L,2L)}.}

\item{colors}{A numeric vector with length equal 2 specifying the colors for lines of the plot. The default value is \code{colors = c("red", "green")}.}

\item{save.plot}{A `TRUE/FALSE` logical value indicating whether the ggplot object of the  plot comparing the confidence intervals should be saved in the output. The default value is \code{save.plot = FALSE}.}

}
\value{
If \code{save.plot = FALSE}, a data frame with columns:
\enumerate{
\item \code{interval}. The method used to construct the confidence interval.
\item \code{bc}. A `TRUE/FALSE` logical value indicating whether the  bias correction is applied.
\item \code{gini}. The estimation of the Gini index.
\item \code{lowerlimit}. The lower limit of the confidence interval.
\item \code{upperlimit}. The upper limit of the confidence interval.
\item \code{var.gini}. The variance estimation for the estimator of the Gini index.
}
If \code{save.plot = TRUE}, a list with two components: (i) `base.CI` a data frame of six columns as just described and (ii) `plot` a (ggplot) description of the plot, which is a list with components that contain the plot itself, the data, information about the scales, panels, etc. As a side-effect, a plot that compares the various methods for constructing confidence intervals for the Gini index is displayed. **ggplot2** is needed to be installed for this option to work.

If \code{plotCI = TRUE}, as a side-effect, a plot that compares the various methods for constructing confidence intervals for the Gini index is displayed. **ggplot2** is needed to be installed for this option to work.
}
\description{
Compares variance estimates and confidence intervals for the Gini index in infinite populations.
}
\details{
For a sample \eqn{S}, with size \eqn{n}, derived from an infinite population, the Gini index is estimated by two different versions (see Muñoz et al., 2023 for more details):

\deqn{\widehat{G}      = \displaystyle \frac{2}{\overline{y}n^{2}}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n}; }

\deqn{\widehat{G}^{bc} = \displaystyle \frac{2}{\overline{y}n(n-1)}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n-1}, }
where the label \eqn{bc} indicates that the bias correction is applied. The table below sumarises the various types of variances and confidence intervals that computes this function.
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

# Estimation of the Gini index and confidence intervals using different methods.
icompareCI(y)
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
\code{\link{igini}}, \code{\link{iginindex}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
