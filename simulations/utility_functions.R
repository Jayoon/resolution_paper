# Utility functions for simulation studies
list.of.packages <- c("cubature", "pracma", "pbapply")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

### Compute local quantile ----

library(cubature) #For 2d integration
library(pracma)
library(pbapply)

# Input
# qs : quantiles to evaluate
# h : bandwidth
# ker: Kernel
# px: pdf of X
# Y.fn: Regression function
# perr: pdf of error if derr = FALSE, cdf of error if derr = TRUE
# derr: indicator of whether pdf 
# x0.grid: Grid to evaluate local quantile function
# q.grid: Grid to evaluate qunatiles
# ker.bound: if the support of kernel is [-a, a]: a 
# Output
# For Y = f(Z) + e, e ~ perr, Z ~ px * K_h(x0-.)
# Compute Qunatile_q(Y)
# TODO 
# 1. Possible optimization of computing integrals
# 2. Choosing the range of q.grid when q*norm.const is not attained
localQuantile <- function(qs = 0.5, h, ker, px, Y.fn, perr, derr,
                          x0.grid, q.grid, ker.bound = Inf,
                          x.supp = c(-Inf, Inf)) {
  result <- matrix(nrow = length(qs), ncol = length(x0.grid))
  colnames(result) <- x0.grid
  rownames(result) <- qs
  integRes <- matrix(nrow = length(x0.grid), ncol = length(q.grid))
  for(i in 1:nrow(integRes)) { # loop over x0's
    if(!derr) {
      fun_sim <- function(x) { #2-d argument (x, y)
        perr(x[2]-Y.fn(x[1])) * px(x[1]) * ker((x0.grid[i]-x[1])/h)
      }
    }
    fun_sim2 <- function(u) {
      px(u) * ker((x0.grid[i]-u)/h)
    }
    norm.const <- cubintegrate(fun_sim2, 
                               lower = max(x0.grid[i]-h*ker.bound, x.supp[1]), 
                               upper = min(x0.grid[i]+h*ker.bound, x.supp[2]),
                               method = "hcubature")$integral
    print(paste("Computing", x0.grid[i]))
    for(j in 1:ncol(integRes)) { #loop over q's
      # Will using the distribution function be faster? (when cdf of error is known)
      if(derr) {
        fun_sim_d <- function(x) {
          perr(q.grid[j]-Y.fn(x)) * px(x) * ker((x0.grid[i]-x)/h)
        }
        if(j %% 100 == 0) print(paste("Computing", q.grid[j]))
        integRes[i,j] <- cubintegrate(f = fun_sim_d,
                                      method = "hcubature",
                                      lower = c(max(x0.grid[i]-h*ker.bound, x.supp[1])),
                                      upper = c(min(x0.grid[i]+h*ker.bound, x.supp[2])))$integral
      } else {
        integRes[i,j] <- cubintegrate(f = fun_sim,
                                      method = "hcubature",
                                      lower = c(max(x0.grid[i]-h*ker.bound, x.supp[1]), -Inf),
                                      upper = c(min(x0.grid[i]+h*ker.bound, x.supp[2]), q.grid[j]))$integral
        # if(integRes[i,j] > q * norm.const) break
      }
    }
    # print(integRes[i,])
    ## how can i do this smartly lol
    for(k in 1:length(qs)) {
      result[k,i] <- q.grid[which.max(integRes[i,] >= qs[k]*norm.const)]
    }
  }
  result
}

### Able to deal with case where p_err depends on x (i.e. heteroscedastic)
### If derr == TRUE, perr(e, x): pdf of e|x
### If derr == FALSE, perr(e, x): cdf of e|x
localQuantile2 <- function(qs = 0.5, h, ker, px, Y.fn, perr, derr,
                           x0.grid, q.grid, ker.bound = Inf,
                           x.supp = c(-Inf, Inf)) {
  result <- matrix(nrow = length(qs), ncol = length(x0.grid))
  colnames(result) <- x0.grid
  rownames(result) <- qs
  integRes <- matrix(nrow = length(x0.grid), ncol = length(q.grid))
  for(i in 1:nrow(integRes)) { # loop over x0's
    if(!derr) {
      fun_sim <- function(x) { #2-d argument (x, y)
        perr(x[2]-Y.fn(x[1]), x[1]) * px(x[1]) * ker((x0.grid[i]-x[1])/h)
      }
    }
    fun_sim2 <- function(u) {
      px(u) * ker((x0.grid[i]-u)/h)
    }
    norm.const <- cubintegrate(fun_sim2, 
                               lower = max(x0.grid[i]-h*ker.bound, x.supp[1]), 
                               upper = min(x0.grid[i]+h*ker.bound, x.supp[2]),
                               method = "hcubature")$integral
    print(paste("Computing", x0.grid[i]))
    for(j in 1:ncol(integRes)) { #loop over q's
      # Will using the distribution function be faster? (when cdf of error is known)
      if(derr) {
        fun_sim_d <- function(x) {
          perr(q.grid[j]-Y.fn(x), x) * px(x) * ker((x0.grid[i]-x)/h)
        }
        if(j %% 100 == 0) print(paste("Computing", q.grid[j]))
        integRes[i,j] <- cubintegrate(f = fun_sim_d,
                                      method = "hcubature",
                                      lower = c(max(x0.grid[i]-h*ker.bound, x.supp[1])),
                                      upper = c(min(x0.grid[i]+h*ker.bound, x.supp[2])))$integral
      } else {
        integRes[i,j] <- cubintegrate(f = fun_sim,
                                      method = "hcubature",
                                      lower = c(max(x0.grid[i]-h*ker.bound, x.supp[1]), -Inf),
                                      upper = c(min(x0.grid[i]+h*ker.bound, x.supp[2]), q.grid[j]))$integral
        # if(integRes[i,j] > q * norm.const) break
      }
    }
    # print(integRes[i,])
    ## how can i do this smartly lol
    for(k in 1:length(qs)) {
      result[k,i] <- q.grid[which.max(integRes[i,] >= qs[k]*norm.const)]
    }
  }
  result
}

# Input
# h : bandwidth
# ker: Kernel
# px: pdf of X
# Y.fn: Regression function
# perr: pdf of error if derr = FALSE, cdf of error if derr = TRUE
# derr: indicator of whether pdf or cdf
# x0.grid: value of x0 (1d)
# q.grid: Grid to evaluate cumulative probabilities
# ker.bound: if the support of kernel is [-a, a]: a 
# x.supp: support of X
# Output
# For Y = f(Z) + e, e ~ perr, Z ~ px * K_h(x0-.)
# Return Compute cdf of Y over q.grid
localCDF <- function(h, ker, px, Y.fn, perr, derr, x0.grid,
                     q.grid, ker.bound = Inf,
                     x.supp = c(-Inf, Inf)) {
  integRes <- vector(mode = "numeric", length(q.grid))
  if(!derr) {
    fun_sim <- function(x) { #2-d argument (x, y)
      perr(x[2]-Y.fn(x[1])) * px(x[1]) * ker((x0.grid-x[1])/h)
    }
  }
  fun_sim2 <- function(u) {
    px(u) * ker((x0.grid-u)/h)
  }
  norm.const <- cubintegrate(fun_sim2, 
                             lower = max(x0.grid-h*ker.bound, x.supp[1]), 
                             upper = min(x0.grid+h*ker.bound, x.supp[2]),
                             method = "hcubature")$integral
  print(norm.const)
  
  for(j in 1:length(integRes)) { #loop over q's
    # Will using the distribution function be faster? (when cdf of error is known)
    if(derr) {
      fun_sim_d <- function(x) {
        perr(q.grid[j]-Y.fn(x)) * px(x) * ker((x0.grid-x)/h)
      }
      if(j %% 100 == 0) print(paste("Computing", q.grid[j]))
      integRes[j] <- cubintegrate(f = fun_sim_d,
                                  method = "hcubature",
                                  lower = c(max(x0.grid-h*ker.bound, x.supp[1])),
                                  upper = c(min(x0.grid+h*ker.bound, x.supp[2])))$integral
    } else {
      integRes[j] <- cubintegrate(f = fun_sim,
                                  method = "hcubature",
                                  lower = c(max(x0.grid-h*ker.bound, x.supp[1]), -Inf),
                                  upper = c(min(x0.grid+h*ker.bound, x.supp[2]), q.grid[j]))$integral
      # if(integRes[i,j] > q * norm.const) break
    }
  }
  integRes/norm.const
}

#### ---- Simulation --------------
#### 1. Signal function -------
fn.step <- function(x) {
  0.2 + 0.6 * ifelse(x > 1/3 & x < 2/3, 1, 0)
}
fn.blip <- function(x) {
  ifelse(x <= 0.8, 0.32 + 0.6*x + 0.3*exp(-100*(x-0.3)^2),
         -0.28 + 0.6*x + 0.3*exp(-100*(x-1.3)^2))
}

tjs <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hjs <- c(4, 5, 3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
wjs <- c(0.005,0.005, 0.006,0.01,0.01,0.03,0.01,0.01,0.005,0.008,0.005)
#ker.quad <- 
fn.bumps <- function(x) {
  sum(hjs*(1+abs((x-tjs)/wjs))^(-4))
}
fn.bumps.v <- Vectorize(fn.bumps)

hjs.2 <- c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
fn.blocks <- function(x) {
  sum(hjs.2*(1+sign(x-tjs))/2)
}
fn.blocks.v <- Vectorize(fn.blocks)

fn.heavisine <- function(x) {
  4*sin(4*pi*x) - sign(x-0.3) - sign(0.72-x)
}
#plot(seq(0,1,by=0.01), fn.heavisine(seq(0,1,by=0.01)),type = "l")


fn.doppler <- function(x) {
  (x*(1-x))^(1/2)*sin(2*pi*(1+0.05)/(x+0.05))
}
# This function is piecewise linear, and continuous, 
# but has big jumps in its first derivatives.
fn.angles <- function(x) {
  (2*x+0.5)*ifelse(x>=0 & x <= 0.15, 1, 0) +
    (-12*(x-0.15)+0.8)*ifelse(x>0.15 & x <= 0.2, 1, 0) +
    0.2*ifelse(x>0.2 & x <= 0.5, 1, 0) +
    (6*(x-0.5)+0.2)*ifelse(x>0.5 & x <= 0.6, 1, 0)+
    (-10*(x-0.6)+0.8)*ifelse(x>0.6 & x <= 0.65, 1, 0) +
    (-5*(x-0.65)+0.3)*ifelse(x>0.65 & x <= 0.85, 1, 0) +
    (2*(x-0.85)+0.2)*ifelse(x>0.85 & x <= 1, 1, 0)
}
# plot(seq(0,1,by=0.01), fn.angles(seq(0,1,by=0.01)),type = "l")

# This function is piecewise parabolic. The function and its first derivative are
# continuous, but there are big jumps in its second derivative.
fn.parabolas <- function(x) {
  0.8-30*(x-0.1)^2*ifelse(x>0.1 & x <= 1, 1, 0) +
    60*(x-0.2)^2*ifelse(x>0.2 & x <= 1, 1, 0) -
    30*(x-0.3)^2*ifelse(x>0.3 & x <= 1, 1, 0) + 
    500*(x-0.35)^2*ifelse(x>0.35 & x <= 1, 1, 0) -
    1000*(x-0.37)^2*ifelse(x>0.37 & x <= 1, 1, 0) + 
    1000*(x-0.41)^2*ifelse(x>0.41 & x <= 1, 1, 0) -
    500*(x-0.43)^2*ifelse(x>0.43 & x <= 1, 1, 0) +
    7.5*(x-0.5)^2*ifelse(x>0.5 & x <= 1, 1, 0) -
    15*(x-0.7)^2*ifelse(x>0.7 & x <= 1, 1, 0) +
    7.5*(x-0.9)^2*ifelse(x>0.9 & x <= 1, 1, 0) 
}
# plot(seq(0,1,by=0.01), fn.parabolas(seq(0,1,by=0.01)),type = "l")


# fn.timeshiftedsine <- function(x) {
#   0.3*sin(3*pi(x + (1-cos(pi*x))/2))+0.5
# }

fn.spikes <- function(x) {
  exp(-500*(x-0.23)^2) + 2*exp(-2000*(x-0.33)^2) +
    4*exp(-8000*(x-0.47)^2) + 3*exp(-16000*(x-0.69)^2) +
    exp(-32000*(x-0.83)^2)
}


#### 2. P_X and P_{XY} -----------------
px.norm <- function(x) { dnorm (x)}
px.unif <- function(x) { ifelse(abs(x/3) <= 1, 1/6, 0)} # Unif[-3,3]
px.mixed <- function(x) { 0.5*dnorm(x, mean = -1) + 0.5*dnorm(x, mean = 1.5, sd = 0.7)}

derr.norm.sd <- function(sd) {
  function(x) { pnorm(x, sd = sd)}
}
derr.norm <- function(x) { pnorm(x, sd = 0.5)}
derr.t1 <- function(x) { pt(2*x, df=1)} # 0.5t(1)
perr.norm <- function(x) { dnorm (x, sd = 0.5)}
perr.t1 <- function(x) {2*dt(2*x, df = 1)} # 0.5t(1)
perr.cauchy <- function(x) {dcauchy(x)}

## heteroscedastic
## 1. sd proportional to given fn
perr.prop <- function(x, fn, sigma = 0.5) {
  dnorm(x, fn(x)^2*sigma^2)
}

# sd = 0.3(x^2+1)
derr.x2p1 <- function(e, x) {
  pnorm(e, sd = 0.3*(x^2+1))
}

# sd = 0.3(x^2-x+5/4)
derr.x2p2 <- function(e, x) {
  pnorm(e, sd = 0.3*(x^2-x+5/4))
}


#### 3. Kernel functions ------------------
ker.gauss <- function(x) {dnorm(x)}
ker.unif <- function(x) {ifelse(abs(x) <= 1, 0.5, 0)}
ker.tri <- function(x) {ifelse(abs(x) <= 1, 1-abs(x), 0)}
ker.para <- function(x) {ifelse(abs(x) <= 1, 0.75*(1-x^2), 0)}
ker.biweight <- function(x) {ifelse(abs(x) <= 1, (15/16)*(1-x^2)^2, 0)}


##### Simulation functions -------------------



###### etc ---------
### percentage of infinity in a vector
per.inf <- function(x) {
  mean(is.infinite(x))
}

## mean after removing inf values in a vector
mean.rminf <- function(x) {
  mean(x[!is.infinite(x)])
}


