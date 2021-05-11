## List of functions required by "1_main.R"

# ker, KDE ########################################################################

# ker computes the mixture kernel function (Gaussian kernel and Cauchy kernel) 
# for point x. w governs the relative weight of the two kernel functions.
# KDE computes the value of the kernel density estimation at point x using the
# observations in the vector data. It uses the kernel function ker with a weight 
# w and the bandwidth specified by h.

ker <- function(x, w) {
  k_heavy <- dcauchy(x)
  k_light <- dnorm(x)
  k <- w * k_heavy + (1 - w) * k_light
  return(k)
}

KDE <- function(x, data, h, w) {
  n <- length(data)
  f <- sapply(x, function(y) 1/ (n) * sum((1 / (h)) * ker(((y - data) / (h)), w)))
  return(f)
}

## lcv_obj  ###################################################################

# Objective function of the likelihood cross validation. 
# Par is the vector of parameters, data is the vector of observations.
# The function returns the value of the objective to be minimized 
# for the given parameter set.

lcv_obj <- function(par, data) {
  n <- length(data)
  h <- par[1]
  w <- par[2]
  index <- seq(1, n, 1)
  ll <- mapply(function(x, y) log((1 / ((n - 1) * h)) * sum(ker(((x - data[-y]) / h), w))),  ## sistemo 
               x = data, y = index)
  J <- sum(ll) / n
  return(-J)
}

## check ######################################################################

# Check if the matrix m is invertible or not. Return 1 if it is, zero otherwise

check <- function(m) "matrix" %in% class(try(solve(m), silent = TRUE))

## sample_gen #################################################################

# Generate synthetic samples of length n from the Gamma, Cauchy or uniform 
# distribution. Distribution is a list that contains the name of the 
# distribution and the parameters.

sample_gen<-function(distribution,n) {
    if (distribution$name == "Gamma"){
    shape <- distribution$par1
    rate <- distribution$par2
    data_set <- rgamma(n, shape = shape, rate = rate)
    } else if (distribution$name == "Cauchy"){
    location <- distribution$par1
    scale <- distribution$par2
    data_set <- rcauchy(n, location = location, scale = scale)
    } else if  (distribution$name == "Uniform") {
      lower <- distribution$par1
    upper <- distribution$par2
    data_set <- runif(n, min = lower, max = upper)
    }
  return(data_set)
}

## quantile_gen ################################################################

# Compute quatiles for the probability levels specified in the vector u for the
# Gamma, Cauchy or uniform distribution. Distribution is a list that contains 
# the name of the distribution and the parameters.

quantile_gen<-function(distribution,u) {
    if (distribution$name == "Gamma"){
    shape <- distribution$par1
    rate <- distribution$par2
    Y <- qgamma(u, shape=shape, rate=rate)
    } else if (distribution$name == "Cauchy"){
    location <- distribution$par1
    scale <- distribution$par2
    Y <- qcauchy(u, location = location, scale = scale)
    } else if  (distribution$name == "Uniform") {
      lower <- distribution$par1
    upper <- distribution$par2
    Y <- qunif(u, min = lower, max = upper)
    }
  return(Y)
}
    