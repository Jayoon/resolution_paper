# Weighted Quantile Method
list.of.packages <- c("spatstat.geom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)
library(spatstat.geom)

# x0: point of interest (d-dim)
# h: bandwidth (d-dim)
# X: Covariate matrix (assumed n by d)
# Y: Response vector
# ker: Kernel function R -> R (assume we apply same kernel to each covariate)
# both: if TRUE: output both corrected and naive CIs, else: only output corrected CI
# q: quantile
# alpha: 1-alpha is desired coverage
# mcReps: number of Monte Carlo simulations to approximate conditional quantiles
# Output: CI from Weighted Quantile method
weightedQuantileCI.corr <- function(x0, h, X, Y, ker, both = FALSE,
                                    q = 0.5, alpha = 0.1, mcReps = 10000) {
  #weights = ker((x0-X)/h)
  weights = apply(ker((as.numeric(x0) - t(X))/as.numeric(h)), 2, prod)
  res <- numeric(length = mcReps)
  weights.eff <- weights[weights!=0]
  Y.eff <- Y[weights!=0]
  n <- length(Y.eff)
  if(n >= 1) {
    for(i in 1:mcReps){
      ber <- rbinom(n, 1, q)
      res[i] <- sum(ber*weights.eff)
    }
    # Estimate correction factor
    q.tilde <- ewcdf(Y.eff, weights.eff)
    theta.hat <- quantile.ewcdf(q.tilde, q)
    num <- mean(weights^2*ifelse(Y <= theta.hat, 1, 0))
    den <- mean(weights^2)
    corr.val <- (1-2*q)/(q*(1-q))*num/den + q/(1-q)
    alpha.tilde <- 1-pnorm(sqrt(corr.val)*qnorm(1-alpha/2))
    if(both) {
      quantile.ewcdf(q.tilde, 
                     quantile(res, c(alpha/2, 1-alpha/2, alpha.tilde, 1-alpha.tilde))/sum(weights.eff))
    } else {
      pn.tilde <- quantile(res, c(alpha.tilde, 1-alpha.tilde))/sum(weights.eff)
      if(pn.tilde[1] == 0) {
        c(-Inf, quantile.ewcdf(q.tilde, pn.tilde[2]))
      } else {
        quantile.ewcdf(q.tilde, pn.tilde)
      }
    }
  } else{
    if(both) {
      c(-Inf, Inf, -Inf, Inf)
    } else {
      c(-Inf, Inf)
    }
  }
}


