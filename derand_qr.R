# Derandomized Quantile Rejection

# ordredY : vector of Y increasing
# theta0 : null hypothesis value of theta being tested
# Output : corresponding p-value
# q: Quantile
pval <- function(orderedY, theta0, q) {
  n = length(orderedY) 
  idx = findInterval(theta0, orderedY)
  if(idx <= n/2) {
    return(2*pbinom(idx, n, q))
  } else {
    return(2*pbinom(n-idx, n, q))
  }
}

# x0: point of interest (assumed 1d for now)
# h: bandwidth
# X: Covariate matrix (assumed n by 1 for now)
# Y: Response vector
# ker: Kernel function
# theta0 : list of null hypothesis value of theta being tested
# K: number of p-values to merge
# q: quantile
# Output: K by length(theta0) matrix of p-values
getPval <- function(x0, h, X, Y, ker, theta0, K = 100, 
                    q = 0.5) {
  n <- length(Y)
  pvalRes <- array(0, c(K, length(theta0)))
  Us <- matrix(runif(n*K), nrow = K)
  for(k in 1:K) {
    y.local <- sort(Y[Us[k,] < ker((x0-X)/h)/ker(0)])
    for(i in 1:length(theta0)){
      pvalRes[k,i] <- pval(y.local, theta0[i])
    }
  }
  pvalRes
}

### mean functions
geo.mean <- function (x) exp(mean(log(x)))
harm.mean <- function (x) 1/(mean(1/x)) 



# pvalRes: K by length(theta0) matrix of p-values
# alpha: 1-alpha is target coverage
# output: 3 CIs from p-value aggregating using AM, GM, HM.
aggregatePval <- function(pvalRes, alpha = 0.1) {
  res <- matrix(nrow = 3, ncol = 2)
  rownames(res) <- c("AM", "GM", "HM")
  colnames(res) <- c("Lower", "Upper")
  ## arithmetic
  a.mean <- 2 *apply(pvalRes, 2, mean)
  res[1,] <- range(theta0[a.mean > alpha])
  ## geometric
  g.mean <- exp(1) * apply(pvalRes, 2, geo.mean)
  res[2,] <- range(theta0[g.mean > alpha])
  ## harmonic
  h.mean <- exp(1) * log(K) * apply(pvalRes, 2, harm.mean)
  res[3,] <- range(theta0[h.mean > alpha]) 
  
  res
}