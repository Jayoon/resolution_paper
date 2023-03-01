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
out_dir <- "../results/housing"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

library(dplyr)

# EDA and initial analysis
# California housing prices data 1990
# Downloaded from https://www.kaggle.com/datasets/camnugent/california-housing-prices
# 1. longitude: A measure of how far west a house is; a higher value is farther west
# 2. latitude: A measure of how far north a house is; a higher value is farther north
# 3. housingMedianAge: Median age of a house within a block; a lower number is a newer building
# 4. totalRooms: Total number of rooms within a block
# 5. totalBedrooms: Total number of bedrooms within a block
# 6. population: Total number of people residing within a block
# 7. households: Total number of households, a group of people residing within a home unit, for a block
# 8. medianIncome: Median income for households within a block of houses (measured in tens of thousands of US Dollars)
# 9. medianHouseValue: Median house value for households within a block (measured in US Dollars)
# 10. oceanProximity: Location of the house w.r.t ocean/sea

housing <- read.csv('housing.csv')
# remove missing value
housing <- na.omit(housing)

housing$avg_bedrooms = housing$total_bedrooms/housing$households
housing$avg_rooms = housing$total_rooms/housing$households
housing <- subset(housing, select = -c(ocean_proximity, total_rooms, total_bedrooms))

## Analysis with longitude and latitude -----
Y <- housing$median_house_value
X <- subset(housing, select = -c(median_house_value))
head(X)


x0.grid <- data.frame(long = c(-122.1697, -122.2585, -121.7617, -118.4452, -117.2340, -117.8443, -119.8489, -122.0593, -118.2851),
                      lat = c(37.4275, 37.8719, 38.5382, 34.0689, 32.8801, 33.6405, 34.4140, 36.9821, 34.0224),
                      row.names = c("Stanford", "UC Berkeley", "UC Davis", "UCLA", "UCSD", "UC Irvine", "UCSB", "UC Santa Cruz", "USC"))
h.grid <- data.frame(long = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4), lat = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4))
q.grid <- c(0.5, 0.7, 0.3, 0.75, 0.8)
alpha <- 0.1

set.seed(0209)
result <- array(dim = c(nrow(x0.grid), nrow(h.grid), length(q.grid), 2, 2),
                dimnames = list(x0 = row.names(x0.grid), h = h.grid[,1], q = q.grid, method = c("WQ", "QR"), ci = c("L", "U")))
for (i in 1:nrow(x0.grid)) {
  for (j in 1:nrow(h.grid)) {
    for(k in 1:length(q.grid)) {
      result[i,j,k,1,] <- weightedQuantileCI.corr(x0.grid[i,], h.grid[j,], X[,1:2], Y, ker = ker.tri, q = q.grid[k], alpha = alpha)
      result[i,j,k,2,] <- finiteSampleMedian(Y[local.sample(x0 = as.numeric(x0.grid[i,]), h = h.grid[j,], X = X[,1:2], ker = ker.tri)], q = q.grid[k], alpha = alpha)
    }
  }
}

saveRDS(result, paste0(out_dir, "/cis_housingx1x2.rds"))


## Analysis with 4 covariates -----
loc <- data.frame(long = c(-122.1697, -122.2585, -121.7617, -118.4452, -117.2340, -117.8443, -119.8489, -122.0593, -118.2851),
                  lat = c(37.4275, 37.8719, 38.5382, 34.0689, 32.8801, 33.6405, 34.4140, 36.9821, 34.0224),
                  school = c("Stanford", "UC Berkeley", "UC Davis", "UCLA", "UCSD", "UC Irvine", "UCSB", "UC Santa Cruz", "USC"))
x0.grid <- crossing(loc, age = seq(5, 45, by = 5), rooms = 1:8)
h.grid <- data.frame(long = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4), lat = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                     age = rep(5, 7), rooms = rep(1, 7))
q.grid <- c(0.5, 0.7, 0.3, 0.75, 0.8)
alpha <- 0.1

set.seed(0209)
result <- array(dim = c(nrow(x0.grid), nrow(h.grid), length(q.grid), 2, 2),
                dimnames = list(x0 = row.names(x0.grid), h = h.grid[,1], q = q.grid, method = c("WQ", "QR"), ci = c("L", "U")))
for (i in 1:nrow(x0.grid)) {
  for (j in 1:nrow(h.grid)) {
    for(k in 1:length(q.grid)) {
      result[i,j,k,1,] <- weightedQuantileCI.corr(x0.grid[i,c(1,2,4,5)], h.grid[j,], X[,c(1,2,3,8)], Y, ker = ker.tri, q = q.grid[k], alpha = alpha)
      result[i,j,k,2,] <- finiteSampleMedian(Y[local.sample(x0 = as.numeric(x0.grid[i,c(1,2,4,5)]), h = h.grid[j,], X = X[,c(1,2,3,8)], ker = ker.tri)], q = q.grid[k], alpha = alpha)
    }
  }
}

saveRDS(result, paste0(out_dir, "/cis_housingx1x2x3x8_2.rds"))

