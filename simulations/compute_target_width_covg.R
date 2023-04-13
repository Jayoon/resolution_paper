# Compute target --------
########################################
### source code
########################################
source('utility_functions.R')

########################################
## Output directory
########################################
out_dir <- "../results/targets"
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
ker.bound <- 1
derr <- TRUE
px <- dunif
## Homoscedastic case -----
perr <- derr.norm.sd(0.3)
for(h in h.list) {
  for(ker in ker.list) {
    for(signal in signal.list) {
      filename <- paste0(out_dir, "/target",
                         "_h", h,
                         "_ker", ker,
                         "_signal", signal,
                         "_setting", 1, ### homoscedastic
                         ".RData")
      if(!file.exists(filename)) {
        res <- localQuantile(qs = qs, 
                             h = h, 
                             ker= ker_info[[ker]], 
                             px= px, 
                             Y.fn = signal_info[[signal]]$Y.fn, 
                             perr= perr, derr=TRUE, 
                             x0.grid= signal_info[[signal]]$x0.grid, 
                             q.grid=seq(-1, signal_info[[signal]]$max.val + 1, by=0.01), 
                             ker.bound = ker.bound)
        save(res, file = filename)
      }
    }
  }
}

## Heteroscedastic case 1 ----
perr=derr.x2p1
for(h in h.list) {
  for(ker in ker.list) {
    for(signal in signal.list) {
      filename <- paste0(out_dir, "/target",
                         "_h", h,
                         "_ker", ker,
                         "_signal", signal,
                         "_setting", 2, ### heteroscedastic 1
                         ".RData")
      if(!file.exists(filename)) {
        res <- localQuantile2(qs = qs, 
                              h = h, 
                              ker = ker_info[[ker]], 
                              px = px, 
                              Y.fn = signal_info[[signal]]$Y.fn, 
                              perr = perr, derr=TRUE, 
                              x0.grid = signal_info[[signal]]$x0.grid, 
                              q.grid =seq(-1, signal_info[[signal]]$max.val + 1, by=0.01), 
                              ker.bound = ker.bound)
        save(res, file = filename)
      }
    }
  }
}

## Heteroscedastic case 2 ------
perr=derr.x2p2
for(h in h.list) {
  for(ker in ker.list) {
    for(signal in signal.list) {
      filename <- paste0(out_dir, "/target",
                         "_h", h,
                         "_ker", ker,
                         "_signal", signal,
                         "_setting", 3, ### heteroscedastic 2
                         ".RData")
      if(!file.exists(filename)) {
        res <- localQuantile2(qs = qs, 
                              h = h, 
                              ker = ker_info[[ker]], 
                              px = px, 
                              Y.fn = signal_info[[signal]]$Y.fn, 
                              perr = perr, derr=TRUE, 
                              x0.grid = signal_info[[signal]]$x0.grid, 
                              q.grid =seq(-1, signal_info[[signal]]$max.val + 1, by=0.01), 
                              ker.bound = ker.bound)
        save(res, file = filename)
      }
    }
  }
}

# Compute average coverage, width for CIs constructed by WQ and QR. ----
########################################
## Output directory
########################################
out_dir <- "../results/covg_width"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

h.list <- c(0.04, 0.06, 0.08, 0.1)
ker.list <- c("tri", "biweight")
signal.list <- c("step", "blip", "bumps", "spikes", "parabolas", "angles")
qs <- c(0.2, 0.5, 0.7)

for(setting in 1:3) {
  for(signal in signal.list) {
    suffix <- paste0("_signal", signal,
                     "_setting", setting,
                     ".RData")
    load(paste0("../results/cis", "/ci",
                "_signal", signal,
                "_setting", setting,
                ".RData"))
    CIs <- res
    covg <- res[,,,,,1,1]
    width <- res[,,,,,1,1] # average width excluding infinity
    infPer <- res[,,,,,1,1] # fraction of intervals with infinite width
    for(h in h.list) {
      for(ker in ker.list) {
        load(paste0("../results/targets/target", "_h", h,
                    "_ker", ker, suffix))
        true.vals <- res
        for(q in qs) {
          target = true.vals[as.character(q),]
          # WQ
          covg[,as.character(h), ker == ker.list, as.character(q), 1] <- 
            apply(CIs[,as.character(h), ker == ker.list, as.character(q), 1, 1,] <= target &
                    CIs[,as.character(h), ker == ker.list, as.character(q), 1, 2,] >= target,
                  1, mean)
          width_vec <- CIs[,as.character(h), ker == ker.list, as.character(q), 1, 2,] - 
            CIs[,as.character(h), ker == ker.list, as.character(q), 1, 1,]
          width[,as.character(h), ker == ker.list, as.character(q), 1] <-
            apply(width_vec, 1, mean.rminf)
          infPer[,as.character(h), ker == ker.list, as.character(q), 1] <-
            apply(width_vec, 1, per.inf)
          # QR
          covg[,as.character(h), ker == ker.list, as.character(q), 2] <-
            apply(CIs[,as.character(h), ker == ker.list, as.character(q), 2, 1,] <= target &
                    CIs[,as.character(h), ker == ker.list, as.character(q), 2, 2,] >= target,
                  1, mean)
          width_vec <- CIs[,as.character(h), ker == ker.list, as.character(q), 2, 2,] - 
            CIs[,as.character(h), ker == ker.list, as.character(q), 2, 1,]
          width[,as.character(h), ker == ker.list, as.character(q), 2] <-
            apply(width_vec, 1, mean.rminf)
          infPer[,as.character(h), ker == ker.list, as.character(q), 2] <-
            apply(width_vec, 1, per.inf)
        }
      }
    }
    saveRDS(covg, file = paste0(out_dir, "/covg", suffix))
    saveRDS(width, file = paste0(out_dir, "/width", suffix))
    saveRDS(infPer, file = paste0(out_dir, "/infPer", suffix))
    
  }
}
