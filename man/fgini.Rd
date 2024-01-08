\name{fgini}
\alias{fgini}
\title{Gini index, variances and confidence intervals in finite populations}
\usage{
fgini(
  y,
  w,
  method = 2L,
  interval = NULL,
  Pi = NULL,
  Pij = NULL,
  PiU,
  alpha = 0.05,
  B = 1000L,
  na.rm = TRUE,
  varformula = "SYG",
  large.sample = FALSE
)

}
\arguments{
\item{y}{A vector with the non-negative real numbers to be used for estimating the Gini index.}

\item{w}{A numeric vector with the survey weights to be used for estimating the Gini index, the variance and the confidence interval. This argument can be missing if argument \code{Pi} is provided. }

\item{method}{An integer between 1 and 5 selecting one of the 5 methods detailed below for estimating the Gini index in finite populations. The default method is \code{method = 2L}.}

\item{interval}{A character string specifying the type of variance estimation and confidence interval to be used. Possible values are \code{"zjackknife"}, \code{"zalinearization"}, \code{"zblinearization"} and \code{"pbootstrap"}. \code{interval = NULL} omits the computation of both variance and confidence interval. The default value is \code{interval = NULL}.}

\item{Pi}{A numeric vector with the (sample) first inclusion probabilites to be used for estimating the Gini index, the variance and the confidence interval. This argument can be \code{NULL} if argument \code{w} is provided. The default value is \code{Pi = NULL}.}

\item{Pij}{A numeric square matrix with the (sample) second (joint) inclusion probabilites to be used for the variance estimation and the confidence interval. The Hajek approximation is used when \code{Pij = NULL}. This argument is used when  \code{interval={"zjackknife", "zalinearization", "zblinearization"}}. The default value is \code{Pij = NULL}.}

\item{PiU}{A numeric vector with the (population) first inclusion probabilites. This argument is only required when the Hartley-Rao expression for the variance estimation is selected (\code{varformula = "HR"}).}

\item{alpha}{A single numeric value between 0 and 1. If \code{interval} is not \code{NULL}, the confidence level to be used for computing  the confidence interval for the Gini is \code{1-alpha}. Some authors call \code{alpha} the significance level. The default value is \code{alpha = 0.05}.}

\item{B}{A single integer specifying the number of bootstrap replicates. This argument is required when \code{interval = "pbootsptrap"}. The default value is \code{B = 1000L}.}

\item{na.rm}{A `TRUE/FALSE` logical value indicating whether \code{NA}'s should be removed before the computation proceeds. The default value is \code{na.rm = TRUE}.}

\item{varformula}{A character string specifying the type of formula to be used for the variance estimator when \code{interval = {"zjackknife", "zalinearization", "zblinearization"}}. Possible values are \code{"HT"} (Hortvitz-Thompson), \code{"SYG"} (Sen-Yates-Grundy) and \code{"HR"} (Hartley-Rao). The default value is \code{varformula = "SYG"}.}

\item{large.sample}{A `TRUE/FALSE` logical value indicating indicating whether the sample is large to apply a faster algorithm to sort the sample values in the computation of the Gini index. The default value is \code{large.sample = FALSE}.}

}
\value{
When \code{interval = NULL}, the function returns a single numeric value between 0 and 1 informing about the estimation of the Gini index. When \code{interval} is not \code{NULL}, the function returns a list with 3 components: a single numeric value with the estimation of the Gini index; a single numeric value with the variance estimation of the Gini index; and a vector of length two containing the lower and upper limits of the confidence interval for the Gini index.
}
\description{
Estimates the Gini index and computes variances and confidence intervals in finite populations.
}
\details{
For a sample \eqn{S}, with size \eqn{n} and inclusion probabilities \eqn{\pi_i=P(i\in S)} (argument \code{Pi}), derived from a finite population \eqn{U}, with size \eqn{N}, different formulations of the Gini index have been proposed in the literature. his function estimates the Gini index, variances and confidence intervals using various formulations. The different methods for estimating the Gini index are (see also Muñoz et al., 2023):

\\
\code{Gini Index formulae.}

\code{method = 1} (Langel and Tillé, 2013)

\deqn{ \widehat{G}_{w1}= \displaystyle \frac{1}{2\widehat{N}^{2}\overline{y}_{w}}\sum_{i \in S}\sum_{j \in S}w_{i}w_{j}|y_{i}-y_{j}|,}
where \eqn{\widehat{N}=\sum_{i \in S}w_i}, \eqn{\overline{y}_{w}=\widehat{N}^{-1}\sum_{i \in S}w_{i}y_{i}}, and \eqn{w_i} are the survey weights. For example, the survey weights can be \eqn{w_i=\pi_{i}^{-1}}. \code{w} or \code{Pi} must be provided, but not both. It is required that \eqn{w_i = \pi_i^{-1}}, for \eqn{i \in S}, when both \code{w} and \code{Pi} are provided.



\code{method = 2} (Alfons and Templ, 2012; Langel and Tillé, 2013)

\deqn{ \widehat{G}_{w2} =\displaystyle \frac{2\sum_{i \in S}w_{(i)}^{+}\widehat{N}_{(i)}y_{(i)} - \sum_{i \in S}w_{i}^{2}y_{i}  }{\widehat{N}^{2}\overline{y}_{w}}-1,}
where \eqn{y_{(i)}} are the values \eqn{y_i} sorted in increasing order, \eqn{w_{(i)}^{+}} are the values \eqn{w_i} sorted according to the increasing order of the values \eqn{y_i}, and \eqn{\widehat{N}_{(i)}=\sum_{j=1}^{i}w_{(j)}^{+}}. Langel and Tillé (2013) show that \eqn{\widehat{G}_{w1} = \widehat{G}_{w2}}.


\code{method = 3} (Berger, 2008)

\deqn{ \widehat{G}_{w3} = \displaystyle \frac{2}{\widehat{N}\overline{y}_{w}}\sum_{i \in S}w_{i}y_{i}\widehat{F}_{w}^{\ast}(y_{i})-1, }
where
\deqn{\widehat{F}_{w}^{\ast}(t) = \displaystyle \frac{1}{\widehat{N}}\sum_{i \in S}w_{i}[\delta(y_i < t) + 0.5\delta(y_i = t)] }
is the smooth (mid-point) distribution function, and \eqn{\delta(\cdot)} is the indicator variable that takes the value 1 when its argument is true, and the value 0 otherwise. It can be seen that \eqn{\widehat{G}_{w2} = \widehat{G}_{w3}}.


\code{method = 4} (Berger and Gedik-Balay, 2020)

\deqn{\widehat{G}_{w4} = 1 - \displaystyle \frac{\overline{v}_{w}}{\overline{y}_{w}},}
where \eqn{\overline{v}_{w}=\widehat{N}^{-1}\sum_{i \in S}w_{i}v_{i}} and
\deqn{v_{i} = \displaystyle \frac{1}{\widehat{N} - w_{i}}\sum_{ \substack{j \in S\\ j\neq i}}\min(y_{i},y_{j}).}


\code{method = 5} (Lerman and Yitzhaki, 1989)

\deqn{\widehat{G}_{w5} = \displaystyle \frac{2}{\widehat{N}\overline{y}_{w}} \sum_{i \in S} w_{(i)}^{+}[y_{(i)} - \overline{y}_{w}]\left[ \widehat{F}_{w}^{LY}(y_{(i)}) - \overline{F}_{w}^{LY} \right], }
where
\deqn{\widehat{F}_{w}^{LY}(y_{(i)}) = \displaystyle \frac{1}{\widehat{N}}\left(\widehat{N}_{(i-1)} + \frac{w_{(i)}^{+}}{2} \right) }
and \eqn{\overline{F}_{w}^{LY}=\widehat{N}^{-1}\sum_{i \in S}w_{(i)}^{+}\widehat{F}_{w}^{LY}(y_{(i)})}.

\\
\code{Variances and confidence intervals.}

For a given estimator \eqn{\widehat{G}_{w}} and variable \eqn{z}, the Horvitz-Thompson type variance estimator (Hortvitz and Thompson, 1952)
\deqn{\widehat{V}_{HT}(\widehat{G}_{w}) = \displaystyle \sum_{i\in S}\sum_{j\in S}\breve{\Delta}_{ij}w_{i}w_{j}z_{i}z_{j}}
is computed when \code{varformula = "HT"}, where
\deqn{\breve{\Delta}_{ij}=\displaystyle \frac{\pi_{ij}-\pi_{i}\pi_{j}}{\pi_{ij}}}
and \eqn{\pi_{ij}} is the second (joint) inclusion probability of the individuals \eqn{i} and \eqn{j}, i.e., \eqn{\pi_{ij}=P\{(i,j)\in S)\} } (argument \code{Pij}).


The Sen-Yates-Grundy type variance estimator (Sen, 1953; Yates and Grundy, 1953)
\deqn{\widehat{V}_{SYG}(\widehat{G}_{w}) = - \displaystyle \frac{1}{2}\sum_{i\in S}\sum_{j\in S}\breve{\Delta}_{ij}(w_{i}z_i-w_{j}z_{j})^{2}} is computed when \code{varformula = "SYG"}, and the Hartley-Rao type variance estimator (Hartley and Rao, 1962) \deqn{\widehat{V}_{HR}(\widehat{G}_{w}) = \displaystyle \frac{1}{n-1}\sum_{i\in S}\sum_{\substack{j \in S\\ j < i}}\left(1-\pi_i-\pi_j + \frac{1}{n}\sum_{k\in U}\pi_{k}^{2} \right)(w_{i}z_i-w_{j}z_{j})^{2}} is computed when \code{varformula = "HR"}. Note that the The Horvitz-Thompson variance estimator can give negative values. We observe that both Horvitz-Thompson and Sen-Yates-Grundy variance estimators depend on second (joint) inclusion probabilities (argument \code{Pij}). The Hajek (1964) approximation \deqn{\pi_{ij}\cong \pi_{i}\pi_{j}\left[1- \displaystyle \frac{(1-\pi_{i})(1-\pi_{j})}{\sum_{i \in S}(1-\pi_{i})} \right]} is used when the second (joint) inclusion probabilities are not available (\code{Pij = NULL}). Note that the Hajek approximation is suggested for large-entropy sampling designs, large samples, and large populations (see Tille 2006; Berger and Tille, 2009; Haziza et al., 2008; Berger, 2011). For instance, this approximation is not recomended for highly-stratified samples (Berger, 2005). The Hartley-Rao variance estimator requires the first inclusion probabilities at the population level (argument \code{PiU}). \code{zjakknife} computes the confidence interval based on the jackknife technique with critical values based on the Normal approximation. \code{zalinearization} and \code{zblinearization} compute the confidence intervals based on the linearization technique applied to the estimators
\deqn{\widehat{G}_{w}^{a} = \widehat{G}_{w1}}
and
\deqn{\widehat{G}_{w}^{b} = \displaystyle \frac{2}{\widehat{N}\overline{y}_{w}}\sum_{i \in S}w_{i}y_{i}\widehat{F}_{w}(y_{i})-1,}
respectively, where
\deqn{\widehat{F}_{w}(t)=\frac{1}{\widehat{N}}\sum_{i \in S}w_i\delta(y_i \leq t).} Critical values are also based on the Normal approximation. \code{pbootstrap} computes the variance using the rescaled bootstrap, and the confidence interval is constructed using the percentile method. The \bold{vignette} \code{vignette("GiniVarInterval")} contains a detailed description of the various methods for variance estimation and confidence intervals for the Gini index.

The following table summarises the various types of variances and confidence intervals that the function  \code{fgini} computes. The argument \code{varformula} only applies for the jackknife and linearization techniques (see Berger, 2008; Langel and Tillé, 2013).

\tabular{llll}{
\bold{Interval}          \tab \bold{Variance} \tab \bold{Critical values}  \tab \bold{References}\cr
_______________     \tab ______________  \tab _________________ \tab _________________________\cr
\code{zjackknife}        \tab Jackknife     \tab Normal                \tab Berger (2008)\cr
\code{zalinearization}    \tab Linearization \tab Normal                \tab Langel and Tille (2013)\cr
\code{zblinearization}    \tab Linearization \tab Normal                \tab Berger (2008)\cr
\code{pBootstrap}        \tab Rescaled bootstrap     \tab Percentile bootstrap  \tab Berger and Gedik-Balay (2020) \cr
}
}
\examples{
# Income and weights (region 'Burgenland') from the 2006 Austrian EU-SILC (Package 'laeken').
data(eusilc, package="laeken")
y <- eusilc$eqIncome[eusilc$db040 == "Burgenland"]
w <- eusilc$rb050[eusilc$db040 == "Burgenland"]

# Estimation of the Gini index using 'method = 2' .
fgini(y, w)


y <- c(30428.83, 14976.54, 18094.09, 29476.79, 20381.93, 6876.17,
       10360.96, 8239.82, 29476.79, 32230.71)
w <- c(357.86, 480.99, 480.99, 476.01, 498.58, 498.58, 476, 498.58, 476.01, 476.01)

# Gini index estimation and confidence interval using:
 ## a: The method 2 for point estimation.
 ## b: The method 'zjackknife' for variance estimation.
 ## c: The Sen-Yates-Grundy type variance estimator.
 ## d: The Hajek approximation for the joint inclusion probabilities.
fgini(y, w, interval = "zjackknife")

# Gini index estimation and confidence interval using:
 ## a: The method 2 for point estimation.
 ## b: The method 'zalinearization' for variance estimation.
 ## c: The Sen-Yates-Grundy type variance estimator.
 ## d: The Hajek approximation for the joint inclusion probabilities.
fgini(y, w, interval = "zalinearization")

# Gini index estimation and confidence interval using:
 ## a: The method 3 for point estimation.
 ## b: The method 'zblinearization' for variance estimation.
 ## c: The Sen-Yates-Grundy type variance estimator.
 ## d: The Hajek approximation for the joint inclusion probabilities.
fgini(y, w, method = 3L, interval = "zblinearization")

# Gini index estimation and confidence interval using:
 ## a: The method 2 for point estimation.
 ## b: The method 'pbootstrap' for variance estimation.
 ## c: The percentile bootstrap method for the confidence interval.
fgini(y, w, interval = "pbootstrap")
}

\references{
Alfons, A., and Templ, M. (2012). Estimation of social exclusion indicators from complex surveys: The R package laeken. KU Leuven, Faculty of Business and Economics Working Paper.

Berger, Y. G. (2005). Variance estimation with highly stratified sampling designs with unequal probabilities. Australian & New Zealand Journal of Statistics, 47, 365–373.

Berger, Y. G. (2008). A note on the asymptotic equivalence of jackknife and linearization variance estimation for the Gini Coefficient. Journal of Official Statistics, 24(4), 541-555.

Berger, Y. G. (2011). Asymptotic consistency under large entropy sampling designs with unequal probabilities. Pakistan Journal of Statistics, 27, 407–426.

Berger, Y. G. and Tillé, Y. (2009). Sampling with unequal probabilities. In Sample Surveys: Design, Methods and Applications (eds. D. Pfeffermann and C. R. Rao), 39–54. Elsevier, Amsterdam

Berger, Y., and Gedik-Balay, I. (2020). Confidence intervals of Gini
coefficient under unequal probability sampling. Journal of Official Statistics, 36(2), 237-249.

Hajek, J. (1964). Asymptotic theory of rejective sampling with varying probabilities from a finite population. The Annals of Mathematical Statistics, 35, 4, 1491–1523.

Hartley, H. O., and Rao, J. N. K. (1962). Sampling with unequal probabilities and without replacement. The Annals of Mathematical Statistics, 350-374.

Haziza, D., Mecatti, F. and Rao, J. N. K. (2008). Evaluation of some approximate variance estimators under the Rao-Sampford unequal probability sampling design. Metron, LXVI, 91–108.

Horvitz, D. G. and Thompson, D. J. (1952). A generalization of sampling without replacement from a finite universe. Journal of the American Statistical Association, 47, 663–685.

Langel, M., and Tille, Y. (2013). Variance estimation of the Gini index: revisiting a result several times published. Journal of the Royal Statistical Society: Series A (Statistics in Society), 176(2), 521-540.

Lerman, R. I., and Yitzhaki, S. (1989). Improving the accuracy of estimates of Gini coefficients. Journal of econometrics, 42(1), 43-47.

Muñoz, J. F., Moya-Fernández, P. J., and Álvarez-Verdejo, E. (2023). Exploring and Correcting the Bias in the Estimation of the Gini Measure of Inequality. Sociological Methods & Research. https://doi.org/10.1177/00491241231176847

Sen, A. R. (1953). On the estimate of the variance in sampling with varying probabilities. Journal of the Indian Society of Agricultural Statistics, 5, 119–127.

Tillé, Y. (2006). Sampling Algorithms. Springer, New York.

Yates, F., and Grundy, P. M. (1953). Selection without replacement from within strata with probability proportional to size. Journal of the Royal Statistical Society B, 15, 253–261.

}

\seealso{
\code{\link{fginindex}}, \code{\link{fcompareCI}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
