\name{fginindex}
\alias{fginindex}
\title{Gini index for finite populations and different estimation methods.}
\usage{
fginindex(
 y,
 w,
 method = 2L,
 Pi = NULL,
 na.rm = TRUE,
 useRcpp = TRUE
)
}
\arguments{
\item{y}{A vector with the non-negative real numbers to be used for estimating the Gini index.}

\item{w}{A numeric vector with the survey weights to be used for estimating the Gini index. This argument can be missing if argument \code{Pi} is provided. }

\item{method}{An integer between 1 and 5 selecting one of the 5 methods detailed below for estimating the Gini index in finite populations. The default method is \code{method = 2L}.}

\item{Pi}{A numeric vector with the (sample) first inclusion probabilites to be used for estimating the Gini index. This argument can be \code{NULL} if argument \code{w} is provided. The default value is \code{Pi = NULL}.}

\item{na.rm}{A `TRUE/FALSE` logical value indicating whether \code{NA}'s should be removed before the computation proceeds. The default value is \code{na.rm = TRUE}.}

\item{useRcpp}{A `TRUE/FALSE` logical value indicating whether \code{Rcpp} (\code{useRcpp = TRUE}), or \code{R} (\code{useRcpp = FALSE}), is used for computation. The default value is \code{UseRcpp = TRUE}.}
}
\value{
A single numeric value between 0 and 1. The estimation of the Gini index.
}
\description{
Estimates the Gini index in finite populations, using different methods.
}
\details{
For a sample \eqn{S}, with size \eqn{n} and inclusion probabilities \eqn{\pi_i=P(i\in S)} (argument \code{Pi}), derived from a finite population \eqn{U}, with size \eqn{N}, different formulations of the Gini index have been proposed in the literature. This function estimates the Gini index using various formulations, and both \code{R} and \verb{C++} codes are implemented. This can be useful for research purposes, and speed comparisons can be made. The different methods for estimating the Gini index are (see also Muñoz et al., 2023):



\code{method = 1} (Langel and Tillé, 2013)

\deqn{ \widehat{G}_{w1}= \displaystyle \frac{1}{2\widehat{N}^{2}\overline{y}_{w}}\sum_{i \in S}\sum_{j \in S}w_{i}w_{j}|y_{i}-y_{j}|,}
where \eqn{\widehat{N}=\sum_{i \in S}w_i}, \eqn{\overline{y}_{w}=\widehat{N}^{-1}\sum_{i \in S}w_{i}y_{i}}, and \eqn{w_i} are the survey weights. For example, the survey weights can be \eqn{w_i=\pi_{i}^{-1}}. \code{w} or \code{Pi} must be provided, but not both. It is required that \eqn{w_i = \pi_i^{-1}}, for \eqn{i \in S}, when both \code{w} and \code{Pi} are provided.


\code{method = 2} (Alfons and Templ, 2012; Langel and Tillé, 2013)

\deqn{ \widehat{G}_{w2} =\displaystyle \frac{2\sum_{i \in S}w_{(i)}^{*}\widehat{N}_{(i)}y_{(i)} - \sum_{i \in S}w_{i}^{2}y_{i}  }{\widehat{N}^{2}\overline{y}_{w}}-1,}
where \eqn{y_{(i)}} are the values \eqn{y_i} sorted in increasing order, \eqn{w_{(i)}^{*}} are the values \eqn{w_i} sorted according to the increasing order of the values \eqn{y_i}, and \eqn{\widehat{N}_{(i)}=\sum_{j=1}^{i}w_{(j)}^{*}}. Langel and Tillé (2013) show that \eqn{\widehat{G}_{w1} = \widehat{G}_{w2}}.


\code{method = 3} (Berger, 2008)

\deqn{ \widehat{G}_{w3} = \displaystyle \frac{2}{\widehat{N}\overline{y}_{w}}\sum_{i \in S}w_{i}y_{i}\widehat{F}_{w}^{\ast}(y_{i})-1, }
where
\deqn{\widehat{F}_{w}^{\ast}(t) = \displaystyle \frac{1}{\widehat{N}}\sum_{i \in S}w_{i}[\delta(y_i < t) + 0.5\delta(y_i = t)] }
is the smooth (mid-point) distribution function, and \eqn{\delta(\cdot)} is the indicator variable that takes the value 1 when its argument is true, and the value 0 otherwise. It can be seen that \eqn{\widehat{G}_{w2} = \widehat{G}_{w3}}.


\code{method = 4} (Berger and Gedik-Balay, 2020)

\deqn{\widehat{G}_{w4} = 1 - \displaystyle \frac{\overline{z}_{w}}{\overline{y}_{w}},}
where \eqn{\overline{z}_{w}=\widehat{N}^{-1}\sum_{i \in S}w_{i}z_{i}} and
\deqn{z_{i} = \displaystyle \frac{1}{\widehat{N} - w_{i}}\sum_{ \substack{j \in S\\ j\neq i}}\min(y_{i},y_{j}).}


\code{method = 5} (Lerman and Yitzhaki, 1989)

\deqn{\widehat{G}_{w5} = \displaystyle \frac{2}{\widehat{N}\overline{y}_{w}} \sum_{i \in S} w_{i}[y_{i} - \overline{y}_{w}]\left[ \widehat{F}_{w}^{LY}(y_{i}) - \overline{F}_{w}^{LY} \right], }
where
\deqn{\widehat{F}_{w}^{LY}(y_{i}) = \displaystyle \frac{1}{\widehat{N}}\left(\widehat{N}_{(i-1)} + \frac{w_{(i)}^{\ast}}{2} \right) }
and \eqn{\overline{F}_{w}^{LY}=\widehat{N}^{-1}\sum_{i \in S}w_{i}\widehat{F}_{w}^{LY}(y_{i})}.

}
\examples{
# Income and weights (region "Burgenland") from the 2006 Austrian EU-SILC (Package 'laeken').
data(eusilc, package="laeken")
y <- eusilc$eqIncome[eusilc$db040 == "Burgenland"]
w <- eusilc$rb050[eusilc$db040 == "Burgenland"]

#Comparing the computation time for the various estimation methods and using R
microbenchmark::microbenchmark(
fginindex(y, w, method = 1L,  useRcpp = FALSE),
fginindex(y, w, method = 2L,  useRcpp = FALSE),
fginindex(y, w, method = 3L,  useRcpp = FALSE),
fginindex(y, w, method = 4L,  useRcpp = FALSE),
fginindex(y, w, method = 5L,  useRcpp = FALSE)
)

# Comparing the computation time for the various estimation methods and using Rcpp
microbenchmark::microbenchmark(
fginindex(y, w, method = 1L),
fginindex(y, w, method = 2L),
fginindex(y, w, method = 3L),
fginindex(y, w, method = 4L),
fginindex(y, w, method = 5L)
)



# Estimation of the Gini index using 'method = 4'.
y <- c(30428.83, 14976.54, 18094.09, 29476.79, 20381.93, 6876.17,
       10360.96, 8239.82, 29476.79, 32230.71)
w <- c(357.86, 480.99, 480.99, 476.01, 498.58, 498.58, 476, 498.58, 476.01, 476.01)
fginindex(y, w, method = 4L)

}
\references{
Alfons, A., and Templ, M. (2012). Estimation of social exclusion indicators from complex surveys: The R package laeken. KU Leuven, Faculty of Business and Economics Working Paper.

Berger, Y. G. (2008). A note on the asymptotic equivalence of jackknife and linearization variance estimation for the Gini Coefficient. Journal of Official Statistics, 24(4), 541-555.

Berger, Y. G., and Gedik-Balay, İ.  (2020). Confidence intervals of Gini coefficient under unequal probability sampling. Journal of official statistics, 36(2), 237-249.

Langel, M., and Tillé, Y. (2013). Variance estimation of the Gini index: revisiting a result several times published. Journal of the Royal Statistical Society: Series A (Statistics in Society), 176(2), 521-540.

Lerman, R. I., and Yitzhaki, S. (1989). Improving the accuracy of estimates of Gini coefficients. Journal of econometrics, 42(1), 43-47.

Muñoz, J. F., Moya-Fernández, P. J., and Álvarez-Verdejo, E. (2023). Exploring and Correcting the Bias in the Estimation of the Gini Measure of Inequality. Sociological Methods & Research. https://doi.org/10.1177/00491241231176847
}
\seealso{
\code{\link{fgini}}, \code{\link{fcompareCI}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
