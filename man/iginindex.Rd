\name{iginindex}
\alias{iginindex}
\title{Gini index for infinite populations and different estimation methods.}
\usage{
iginindex(
  y,
  method = 5L,
  bias.correction = TRUE,
  cum.sums = NULL,
  na.rm = TRUE,
  useRcpp = TRUE
)
}
\arguments{
\item{y}{A vector with the non-negative real numbers to be used for estimating the Gini index. This argument can be missing if argument \code{cum.sums} is provided.}

\item{method}{An integer between 1 and 10 selecting one of the 10 methods detailed below for estimating the Gini index in infinite populations. The default method is \code{method = 5L}.}

\item{bias.correction}{A `TRUE/FALSE` logical value indicating whether the bias correction should be applied to the estimation of the Gini index. The default value is \code{bias.correction = TRUE}.}

\item{cum.sums}{A vector with the non-negative real numbers specifying the cumulative sums of the variable used to estimate the Gini index. This argument can be \code{NULL} if argument \code{y} is provided. The default value is \code{cum.sums = NULL}.}

\item{na.rm}{A `TRUE/FALSE` logical value indicating whether \code{NA}'s should be removed before the computation proceeds. The default value is \code{na.rm = TRUE}.}

\item{useRcpp}{A `TRUE/FALSE` logical value indicating whether \code{Rcpp} (\code{useRcpp = TRUE}) or \code{R} (\code{useRcpp = FALSE}) is used for computation. The default value is \code{UseRcpp = TRUE}.}
}
\value{
A single numeric value between 0 and 1 containing the estimation of the Gini index based on the vector \code{y} or the vector \code{cum.sums}.
}
\description{
Estimates the Gini index in infinite populations, using different methods.
}
\details{
For a sample \eqn{S}, with size \eqn{n}, derived from an infinite population, different formulations of the Gini index have  been proposed in the literature, but they only provide two different outputs.

This function estimates the Gini index using the various formulations, and both \code{R} and \verb{C++} codes are implemented. This can be useful for research purposes, and speed comparisons can be made. The argument \code{cum.sums} does not require that the cumulative sums are based on the non-decreasing order of the variable \code{y}.

The different methods for estimating the Gini index are (see Wang et al., 2016; Giorgi and Gigliarano, 2017; Mukhopadhyay and Sengupta, 2021; Muñoz et al., 2023):


\code{method = 1}
\deqn{\widehat{G}_1 = \displaystyle \frac{1}{2\overline{y}n^{2}}\sum_{i \in S}\sum_{j\in S} |y_i-y_j|;}
\deqn{\widehat{G}_{1}^{bc} = \displaystyle \frac{1}{2\overline{y}n(n-1)}\sum_{i \in S} \sum_{j \in S} |y_i-y_j|,}
where \eqn{\overline{y} = n^{-1}\sum_{i \in S}y_i} is the sample mean and the label \eqn{bc} indicates that the bias correction is applied to the estimation of the Gini index.

\code{method = 2}
\deqn{\widehat{G}_{2}      = \displaystyle \frac{n-1}{n}\frac{\sum_{i=1}^{n-1}(p_i-q_i)}{\sum_{i=1}^{n-1}pi};}
\deqn{\widehat{G}_{2}^{bc} = \displaystyle \frac{\sum_{i=1}^{n-1}(p_i-q_i)}{\sum_{i=1}^{n-1}pi},}
where
\deqn{ p_i= \displaystyle \frac{i}{n}; \quad q_i= \frac{y_{i}^{+}}{y_{n}^{+}}, }
and \eqn{y_{i}^{+}=\sum_{j=1}^{i}y_{(j)}}, with \eqn{i=\{1,\ldots,n\}}, are the cumulative sums
of the ordered values \eqn{y_{(i)}} (in non-decreasing order) of the variable of interest \eqn{y}.

\code{method = 3}
\deqn{\widehat{G}_{3}      = \displaystyle \frac{n-1}{n} - \frac{2}{n}\sum_{i=1}^{n-1}q_i; }
\deqn{\widehat{G}_{3}^{bc} = 1  - \displaystyle \frac{2}{n-1}\sum_{i=1}^{n-1}q_i. }

\code{method = 4}
\deqn{\widehat{G}_{4}      = 1 - \displaystyle \sum_{i=0}^{n-1}(q_{i+1} + q_i)(p_{i+1} - p_i); }
\deqn{\widehat{G}_{4}^{bc} = \displaystyle \frac{n}{n-1}\left[1 - \sum_{i=0}^{n-1}(q_{i+1} + q_i)(p_{i+1} - p_i)\right], }
where \eqn{p_0=q_0=0.}

\code{method = 5}
\deqn{\widehat{G}_{5}      = \displaystyle \frac{2}{\overline{y}n^{2}}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n}; }
\deqn{\widehat{G}_{5}^{bc} = \displaystyle \frac{2}{\overline{y}n(n-1)}\sum_{i \in S}iy_{(i)} - \frac{n+1}{n-1}. }

\code{method = 6}
\deqn{\widehat{G}_{6}      = \displaystyle \frac{2}{\overline{y}n}cov(i,y_{(i)}); }
\deqn{\widehat{G}_{6}^{bc} = \displaystyle \frac{2}{\overline{y}(n-1)}cov(i,y_{(i)}). }

\code{method = 7}
\deqn{\widehat{G}_{7}      = \displaystyle \frac{1}{\overline{y}n^2}\sum_{i \in S}\sum_{j\in S}|y_i-y_j|\cdot |\widehat{F}_{n}^{\ast}(y_{i})-\widehat{F}_{n}^{\ast}(y_{j})|; }
\deqn{\widehat{G}_{7}^{bc} = \displaystyle \frac{1}{\overline{y}n(n-1)}\sum_{i\in S}\sum_{j \in S}|y_i-y_j|\cdot |\widehat{F}_{n}^{\ast}(y_{i})-\widehat{F}_{n}^{\ast}(y_{j})|, }
where
\deqn{\widehat{F}_{n}^{\ast}(t)= \displaystyle \frac{1}{n}\sum_{i \in S}[\delta(y_i < t) + 0.5\delta(y_i = t)]}
is the smooth (mid-point) distribution function.

\code{method = 8}
\deqn{\widehat{G}_{8}      = 1 - \displaystyle \frac{1}{\overline{y}n^2}\sum_{i \in S}\sum_{j \in S}min(y_i,y_j); }
\deqn{\widehat{G}_{8}^{bc} = 1 - \displaystyle \frac{1}{\overline{y}n(n-1)}\sum_{i \in S}\sum_{\substack{j \in S\\ j\neq i} }min(y_i,y_j).}

\code{method = 9}
\deqn{\widehat{G}_{9}      = \displaystyle \frac{2}{\overline{y}n}\sum_{i \in S}y_{i}\widehat{F}_{n}^{\ast}(y_{i}) - 1;}
\deqn{\widehat{G}_{9}^{bc} = \displaystyle \frac{2}{\overline{y}(n-1)}\sum_{i \in S}y_{i}\widehat{F}_{n}^{\ast}(y_{i}) - \frac{n}{n-1}.}

\code{method = 10}
\deqn{\widehat{G}_{10}      = \displaystyle \frac{n-1}{2\overline{y}n}\binom{n}{2}^{-1}\sum_{i \leq i_{1} < i_{2} \leq n}|y_{i_{1}}-y_{i_{2}}|;}
\deqn{\widehat{G}_{10}^{bc} = \displaystyle \frac{1}{2\overline{y}}\binom{n}{2}^{-1}\sum_{i \leq i_{1} < i_{2} \leq n}|y_{i_{1}}-y_{i_{2}}|.}
}
\examples{
# Sample, with size 50, from a Lognormal distribution. The true Gini index is 0.5.
set.seed(123)
y <- gsample(n = 50, gini = 0.5, meanlog = 5)

# Estimation of the Gini index using the method = 5, bias correction, and Rcpp.
iginindex(y)

# Estimation of the Gini index using the method = 5, bias correction, and R.
iginindex(y, useRcpp = FALSE)

#Comparing the computation time for the various estimation methods and using R
microbenchmark::microbenchmark(
iginindex(y, method = 1,  useRcpp = FALSE),
iginindex(y, method = 2,  useRcpp = FALSE),
iginindex(y, method = 3,  useRcpp = FALSE),
iginindex(y, method = 4,  useRcpp = FALSE),
iginindex(y, method = 5,  useRcpp = FALSE),
iginindex(y, method = 6,  useRcpp = FALSE),
iginindex(y, method = 7,  useRcpp = FALSE),
iginindex(y, method = 8,  useRcpp = FALSE),
iginindex(y, method = 9,  useRcpp = FALSE),
iginindex(y, method = 10, useRcpp = FALSE)
)

# Comparing the computation time for the various estimation methods and using Rcpp
microbenchmark::microbenchmark(
iginindex(y, method = 1),
iginindex(y, method = 2),
iginindex(y, method = 3),
iginindex(y, method = 4),
iginindex(y, method = 5),
iginindex(y, method = 6),
iginindex(y, method = 7),
iginindex(y, method = 8),
iginindex(y, method = 9),
iginindex(y, method = 10) )

}

\references{

Giorgi, G. M., and Gigliarano, C. (2017). The Gini concentration index: a review of the inference literature. Journal of Economic Surveys, 31(4), 1130-1148.

Mukhopadhyay, N., and Sengupta, P. P. (Eds.). (2021). Gini inequality index: Methods and applications. CRC press.

Muñoz, J. F., Moya-Fernández, P. J., and Álvarez-Verdejo, E. (2023). Exploring and Correcting the Bias in the Estimation of the Gini Measure of Inequality. Sociological Methods & Research. https://doi.org/10.1177/00491241231176847

Wang, D., Zhao, Y., and Gilmore, D. W. (2016). Jackknife empirical likelihood confidence interval for the Gini index. Statistics & Probability Letters, 110, 289-295.

}
\seealso{
\code{\link{igini}}, \code{\link{icompareCI}}
}
\author{
Juan F Munoz \email{jfmunoz@ugr.es}

Jose M Pavia \email{pavia@uv.es}

Encarnacion Alvarez \email{encarniav@ugr.es}
}
