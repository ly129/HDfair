setwd("/Users/YLIAN/Library/CloudStorage/Box-Box/FairReg/eesimulation/")
# set.seed(20250121)

library(glmnet)

# data generation
## no. groups
A <- 3L
## no.sites
M <- 2L
## dimension
p <- 20L
## nonzero dimension
p.nz <- 10L
## sample size per site
n.M <- 1000L * seq(M)
## gaussian noise
sd.yi <- sqrt(1)
## proportions of groups
props <- c(0.7, 0.2, 0.1)
## generate groups
pts <- vector(mode = "list", length = M)
for (m in 1:M) {
  pts[[m]] <- props
}
id.cnt <- mapply(FUN = function(x, p) {
  round( x * p )
}, n.M, pts, SIMPLIFY = FALSE)
group <- lapply(id.cnt, FUN = function(x) {sample(rep(seq(A), times = x))})
## true heterogeneity
heterogeneity <- "low"

# pre-allocation
# no. replications
n.rep <- 1
## store data
train.list <- test.list <- vector(mode = "list", length = n.rep)
## store fairness
fair.check.eta <- vector(mode = "list", length = n.rep)
## store variance of outcome for signal-noise ratio
var.y <- var.y.test <- numeric(n.rep)


# generation
## true theta
th <- vector("list", M)
th.base <- rnorm(p.nz, 0, 1)
for (i in seq(M)) {
  th[[i]] <- matrix(0, nrow = p, ncol = A)
  if (heterogeneity == "no") {
    # no heterogeneity
    th[[i]][1:p.nz, 1] <- 0.7
    th[[i]][1:p.nz, 2] <- 0.7
    th[[i]][1:p.nz, 3] <- 0.7
  } else if (heterogeneity == "low") {
    # low heterogeneity
    th[[i]][1:p.nz, 1] <- rep(c(0.6, 0.7, 0.8), length.out = p.nz)
    th[[i]][1:p.nz, 2] <- rep(c(0.7, 0.8, 0.6), length.out = p.nz)
    th[[i]][1:p.nz, 3] <- rep(c(0.8, 0.6, 0.7), length.out = p.nz)
  } else if (heterogeneity == "medium") {
    # medium heterogeneity
    th[[i]][1:p.nz, 1] <- rep(c(0.5, 0.7, 0.9), length.out = p.nz) # 0.5
    th[[i]][1:p.nz, 2] <- rep(c(0.7, 0.9, 0.5), length.out = p.nz) # 0.2
    th[[i]][1:p.nz, 3] <- rep(c(0.9, 0.5, 0.7), length.out = p.nz) # 0.8
  } else if (heterogeneity == "high") {
    # high heterogeneity
    th[[i]][1:p.nz, 1] <- rep(c(0.4, 0.7, 1.0), length.out = p.nz) # 0.5
    th[[i]][1:p.nz, 2] <- rep(c(0.7, 1.0, 0.4), length.out = p.nz) # 0.2
    th[[i]][1:p.nz, 3] <- rep(c(1.0, 0.4, 0.7), length.out = p.nz) # 0.8
  }
}
head(th[[1]], 20)

for (r in 1:n.rep) {
  ## true X
  X <- sapply(n.M,
              function(x) {
                matrix(
                  rnorm(n = p * x, mean = 0, sd = 1),
                  nrow = x, ncol = p
                )
              }, simplify = FALSE
  )

  X.test <- sapply(n.M,
                   function(x) {
                     matrix(
                       rnorm(n = p * x, mean = 0, sd = 1),
                       nrow = x, ncol = p
                     )
                   }, simplify = FALSE
  )

  ## true y
  y <- y.test <- vector(mode = "list", length = M)
  for (m in seq(M)) {
    nm <- n.M[m]
    y[[m]] <- numeric(nm)
    y.test[[m]] <- numeric(nm)
    for (i in seq(nm)) {
      a.tmp <- group[[m]][i]
      X.tmp <- X[[m]][i, ]
      X.test.tmp <- X.test[[m]][i, ]
      th.tmp <- th[[m]][, a.tmp]
      y[[m]][i] <- c(crossprod(X.tmp, th.tmp)) + rnorm(1, mean = 0, sd = sd.yi)
      y.test[[m]][i] <- c(crossprod(X.test.tmp, th.tmp)) + rnorm(1, mean = 0, sd = sd.yi)
    }
    ymean <- mean(c(y[[m]], y.test[[m]]))
    y[[m]] <- y[[m]] - ymean
    y.test[[m]] <- y.test[[m]] - ymean
  }
  var.y[r] <- var(y[[m]])
  var.y.test[r] <- var(y.test[[m]])
  cat("Var(y) =", var.y[r], "| Var(y.test) =", var.y.test[r], "\n")

  train.list[[r]]$x <- X; train.list[[r]]$y <- y
  test.list[[r]]$x <- X.test; test.list[[r]]$y <- y.test
}

# # check data generation with glmnet
# glmnetcv <- cv.glmnet(train.list[[1]]$x[[1]], train.list[[1]]$y[[1]])
# plot(glmnetcv)
# plot(glmnetcv$glmnet.fit, xvar = "lambda")


# Algorithm
## parameters
eta_fr <- 0.1
eta_ee <- 2e-1  # it seems like 1e-1 is a sweet spot
# tuning
eps <- 1e-6
rho_fr <- 1.1
rho_ee <- rho_fr
stepsize <- 0.03
maxit <- 1e4
verbose <- TRUE
# loose eta_ee similar to heavier penalization
# smaller eta requires smaller rho to converge
# however larger rho helps reach etas
# fix eta_ee and rho_ee is a good idea
# needs an algorithm to reduce rho?

# theta starts from all zeros
# use smaller rho/stepsize to keep all fairness constraints satisfied
## Not really feasible, there are methods to keep the solution within the feasibility region
## e.g. barrier methods


# Latest tuning experience - at least for simulation
## First find a suitable eta_ee, which generates a rather promising sparsity pattern
## It is probably OK to use the same rho for the two sets of constraints
## In the algorithm, the slack variables are set to eta * 0.999999, really helps convergence

custom_esti_func <- function(x, y, theta) {
  1/nrow(x) * crossprod(x, (x %*% theta - y) )
}

custom_esti_func_grad <- function(x, y, theta) {
  1/nrow(x) * crossprod(x)
}


X <- train.list[[1]]$x
y <- train.list[[1]]$y

source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/multifairEE_base.R")

res <- multifairee_base(
  X,
  y,
  group = group,
  eta_fr = eta_fr,
  eta_ee = eta_ee,
  stepsize = stepsize,
  custom_esti_func = custom_esti_func,
  custom_esti_func_grad = custom_esti_func_grad,
  rho_fr = rho_fr,
  rho_ee = rho_ee,
  maxit = maxit,
  eps = 1e-6,
  verbose = verbose
)


( theta <- res$Estimates )
theta_bar <- apply(theta, MARGIN = 1:2, FUN = mean)
theta_ctrd <- apply(theta, MARGIN = 3, FUN = "-", theta_bar, simplify = FALSE)
theta_ctrd <- array(unlist(theta_ctrd), dim = c(p, M, A))

fair <- apply(theta_ctrd^2, 3, sum)/2; fair

ee <- array(dim = c(p, M, A))

for (m in 1:M) {
  xm <- X[[m]]
  ym <- y[[m]]
  grp <- group[[m]]

  for (a in 1:A) {
    a.id <- which(grp == a)
    xma <- xm[a.id, ]
    yma <- ym[a.id]
    thetama <- theta[, m, a]
    ee[, m, a] <- custom_esti_func(xma, yma, thetama)
  }
}
apply(ee, MARGIN = 2:3, function(xxx) {xxx[which.max(abs(xxx))]})
