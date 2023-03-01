source('../R/weighted_quantile.R')
source('../R/quantile_rejection.R')
source('../R/derand_qr.R')
source('utility_functions.R')

list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

########################################
## Output directory
########################################
out_dir <- "../results/compliance"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

# Cholesterol data
# Cholestyramine, a proposed cholesterol lowering drug, was administered to 164 men for an average of seven years each.
# -2.25 to 1.97 on the normalized scale
# norm.comp = 0.0422 * comp - 2.25
# The response variable is a man's decrease in cholesterol level over the course of the experiment.
# The single predictor is compliance, the fraction of intended dose actually taken (standardized)
cholesterol <- read.table(url("https://hastie.su.domains/CASI_files/DATA/cholesterol.txt"),
                          header = TRUE)
head(cholesterol)
plot(cholesterol)
library(dplyr)


###########-----------------
#### Compliance real data analysis
head(cholesterol)
Y <- cholesterol$cholesterol.decrease
X <- cholesterol$compliance
set.seed(1127)
x0.grid <- seq(-2, 2, by = 0.05)
h.grid <- c(0.25, 0.5, 1)
q.grid <- c(0.5, 0.7, 0.3, 0.75, 0.8)
alpha <- 0.1


result <- array(dim = c(length(x0.grid), length(h.grid), length(q.grid), 2, 2),
                dimnames = list(x0 = x0.grid, h = h.grid, q = q.grid, method = c("WQ", "QR"), ci = c("L", "U")))
for (i in 1:length(x0.grid)) {
  for (j in 1:length(h.grid)) {
    for(k in 1:length(q.grid)) {
      result[i,j,k,1,] <- weightedQuantileCI.corr(x0.grid[i], h.grid[j], X, Y, ker = ker.tri, q = q.grid[k], alpha = alpha)
      
      result[i,j,k,2,] <- finiteSampleMedian(Y[local.sample(x0 = x0.grid[i], h = h.grid[j], X = X, ker = ker.tri)], q = q.grid[k], alpha = alpha)
    }
  }
}

saveRDS(result, paste0(out_dir, "/cis_tri.rds"))



