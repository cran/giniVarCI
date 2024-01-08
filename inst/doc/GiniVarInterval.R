## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(giniVarCI)
set.seed(123)
y <- gsample(n = 100, gini = 0.5, distribution = "lognormal")
igini(y)


## -----------------------------------------------------------------------------


#Comparing the computation time for the various estimation methods using R 
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


# Comparing the computation time for the various estimation methods using Rcpp
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




## -----------------------------------------------------------------------------


# Comparing the computation time for estimates of the Gini index in various R packages.

microbenchmark::microbenchmark(
igini(y),
laeken::gini(y),
DescTools::Gini(y),
ineq::Gini(y),
REAT::gini(y))



## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'pbootstrap',

igini(y, interval = "pbootstrap")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'Bca'.

igini(y, interval = "BCa")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'zjackknife'.

igini(y, interval = "zjackknife")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'tjackknife'.

igini(y, interval = "tjackknife")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'zalinearization'.

igini(y, interval = "zalinearization")



# Gini index estimation and confidence interval using 'zblinearization'.

igini(y, interval = "zblinearization")



## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'talinearization'.

igini(y, interval = "talinearization")


# Gini index estimation and confidence interval using 'tblinearization'.

igini(y, interval = "tblinearization")



## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'ELchisq'.

igini(y, interval = "ELchisq")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'ELboot'.

igini(y, interval = "ELboot")


## -----------------------------------------------------------------------------

# Comparisons of variance estimators and confidence intervals.

icompareCI(y, plotCI = FALSE)


## -----------------------------------------------------------------------------
data(eusilc, package="laeken")
y <- eusilc$eqIncome[eusilc$db040 == "Burgenland"]
w <- eusilc$rb050[eusilc$db040 == "Burgenland"]
fgini(y, w)


## -----------------------------------------------------------------------------

#Comparing the computation time for the various estimation methods and using R
microbenchmark::microbenchmark(
fginindex(y, w, method = 1,  useRcpp = FALSE),
fginindex(y, w, method = 2,  useRcpp = FALSE),
fginindex(y, w, method = 3,  useRcpp = FALSE),
fginindex(y, w, method = 4,  useRcpp = FALSE),
fginindex(y, w, method = 5,  useRcpp = FALSE)
)


# Comparing the computation time for the various estimation methods and using Rcpp
microbenchmark::microbenchmark(
fginindex(y, w, method = 1),
fginindex(y, w, method = 2),
fginindex(y, w, method = 3),
fginindex(y, w, method = 4),
fginindex(y, w, method = 5)
)


## -----------------------------------------------------------------------------

# Comparing the computation time for estimates of the Gini index in various R packages.

# Comparing 'method = 2', used also by the laeken package. 

microbenchmark::microbenchmark(
fgini(y,w),
laeken::gini(y,w)
)


# Comparing 'method = 5', used also by the DescTools and REAT packages. 

microbenchmark::microbenchmark(
fgini(y,w, method = 5),
DescTools::Gini(y,w),
REAT::gini(y, weighting = w)
)


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'pbootstrap'.

fgini(y, w, interval = "pbootstrap")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using 'zjackknife'.

fgini(y, w, interval = "zjackknife")


## -----------------------------------------------------------------------------

# Gini index estimation and confidence interval using:
 ## a: The method 2 for point estimation. 
 ## b: The method 'zalinearization' for variance estimation. 
 ## c: The Sen-Yates-Grundy type variance estimator. 
 ## d: The Hàjek approximation for the joint inclusion probabilities. 
fgini(y, w, interval = "zalinearization")

# Gini index estimation and confidence interval using:
 ## a: The method 3 for point estimation. 
 ## b: The method 'zblinearization' for variance estimation. 
 ## c: The Sen-Yates-Grundy type variance estimator. 
 ## d: The Hàjek approximation for the joint inclusion probabilities. 
fgini(y, w, method = 3, interval = "zblinearization")


