# Quantile Rejection Method

# h: bandwidth
# x0: local point
# X: Samples from P_X 
# Kernel: Assume that the kernel attains maximum at 0
# Output: boolean vector indicating the indices of local samples from rejection sampling
local.sample <- function(x0, h, X, ker = dnorm) {
  # n <- length(X)
  # runif(n) < ker((x0-X)/h)/ker(0)
  weights <-apply(ker((as.numeric(x0) - t(X))/as.numeric(h))/ker(0), 2, prod)
  runif(length(weights)) < weights
}


# input: input samples
# q: quantile
# alpha: 1-alpha is desired coverage
# Output: CI from quantile rejection method 
finiteSampleMedian <- function(input, q, alpha = 0.1) {
  n <- length(input)
  sorted <- sort(input)
  #if(n <= 4) return(c(-Inf, Inf))
  k1 <- qbinom(alpha/2, size = n, prob = q)
  k2 <- qbinom(1-alpha/2, size = n, prob = q)
  lb <- ifelse(k1 == 0, -Inf, sorted[k1])
  ub <- ifelse(k2 == n, Inf, sorted[k2+1])
  c(lb, ub)
}

# Usage
# finiteSampleMedian(Y[local.sample(x0 = 0.1, h = 0.5, X = X)])