## Compute QR and WQ confidence intervals for different settings
list.of.packages <- c("pbapply")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(pbapply)
source('../R/weighted_quantile.R')
source('../R/quantile_rejection.R')
source('utility_functions.R')

########################################
## Output directory
########################################
out_dir <- "../results/cis"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

signal_info <- list(
  step = list(Y.fn = fn.step, x0.grid = c(1/6, 1/3, 1/2, 2/3, 5/6), max.val = 0.8), 
  blip = list(Y.fn = fn.blip, x0.grid = c(0.2, 0.4, 0.6, 0.8, 0.9), max.val = 0.8),
  bumps = list(Y.fn = fn.bumps.v, x0.grid = c(0.1, 0.15, 0.25, 0.35, 0.4, 0.5, 0.6, 0.7, 0.79,0.85), max.val = 5), 
  spikes = list(Y.fn = fn.spikes, x0.grid = c(0.23, 0.28, 0.33, 0.47, 0.6, 0.69, 0.83, 0.9), max.val = 4),
  parabolas = list(Y.fn = fn.parabolas, x0.grid = c(0.1, 0.2, 0.3, 0.35, 0.37, 0.41, 0.43, 0.5, 0.6, 0.7, 0.8, 0.9),
                   max.val = 0.8), 
  angles = list(Y.fn = fn.angles, x0.grid = c(0.15, 0.2, 0.4, 0.5, 0.6, 0.65, 0.7, 0.85, 0.95), max.val = 1)
)

ker_info <- list(tri = ker.tri, biweight = ker.biweight)

h.list <- c(0.1, 0.08, 0.06, 0.04)
ker.list <- c("tri", "biweight")
signal.list <- c("step", "blip", "bumps", "spikes", "parabolas", "angles")
qs <- c(0.2, 0.5, 0.7)
# n: sample size
# setting : noise distribution (numbered 1 to 3 as in the paper)

getCIs <- function(n, h.list, x0.grid, qs, setting, Y.fn, ker.list) {
  # Output
  results <- array(0, dim = c(length(x0.grid), length(h.list), length(ker.list), length(qs), 2, 2))
  dimnames(results) <- list(x0 = x0.grid,
                            h = h.list,
                            ker = ker.list,
                            q = qs,
                            method = c("WQ", "QR"),
                            CI = c("L", "U"))
  # Data generation
  X <- sort(runif(n))
  if(setting == 1) {
    Y <- Y.fn(X) + 0.3*rnorm(n)
  } else if(setting == 2) {
    Y <- Y.fn(X) + 0.3*(X^2+1)*rnorm(n)
  } else if(setting == 3) {
    Y <- Y.fn(X) + 0.3*(X^2-X+5/4)*rnorm(n)
  }
  
  # Compute CI
  for(i in 1:length(x0.grid)) {
    for(j in 1:length(h.list)) {
      h = h.list[j]
      for(k in 1:length(ker.list)) {
        ker = ker.list[k]
        for(l in 1:length(qs)) {
          q = qs[l]
          ## Weighted quantile
          results[i,j,k,l,1,] <- weightedQuantileCI.corr(x0.grid[i], h = h, X = X, Y = Y, 
                                                     ker = ker_info[[ker]], q = q)
          ## Quantile rejection
          results[i,j,k,l,2,] <- finiteSampleMedian(Y[local.sample(x0 = x0.grid[i], h = h, X = X,
                                                                   ker = ker_info[[ker]])], q = q)
        }
      }
    }
  }
  results
}


n <- 200
nsims <- 1000
for(setting in 1:3) {
  for(signal in signal.list) {
    filename <- paste0(out_dir, "/ci",
                       "_signal", signal,
                       "_setting", setting,
                       ".RData")
    if(!file.exists(filename)) {
      res <- pbreplicate(n = nsims,
                         getCIs(n, h.list, signal_info[[signal]]$x0.grid, 
                                qs, setting, signal_info[[signal]]$Y.fn, 
                                ker.list))
      save(res, file = filename)
    }
  }
}


