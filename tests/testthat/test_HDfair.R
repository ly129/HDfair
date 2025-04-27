rm(list = ls())

source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/HDfair_base.R")
source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/HDfair_sp_lambda.R")
source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/HDfair_sp_eta.R")
source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/HDfair_cv_lambda.R")
source("~/Library/CloudStorage/Box-Box/FairReg/HDfair/R/HDfair_cv_eta.R")


A <- 3L
M <- 2L
p <- 1000L
p.nz <- 10L
n.M <- 500L * seq(M)
N <- sum(n.M)
props <- c(0.7, 0.2, 0.1)

heterogeneity <- "high"

# generate true theta

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

th[[1]][1:p.nz, ] <- th[[1]][1:p.nz, ] + 1.4
th[[2]][1:p.nz, ] <- th[[2]][1:p.nz, ] + 1.2

th.true <- array(dim = c(p,A,M))

for (m in 1:M) {
  th.true[, , m] <- th[[m]]
}

th.true[1:20,,]

# gaussian noise
sd.yi <- sqrt(1)

ma = matrix(0,N,2)
ma[,1] = rep(1:M, n.M)

for (m in 1:M) {
  ma[ma[, 1] == m, 2] <- rep(1:A, n.M[m] * props)
}


X <- matrix(NA, N, p)
y <- numeric(N)

# x, y and group (indicator) are all lists of length M
for (m in 1:M) {
  for (a in 1:A) {
    nma <- n.M[m] * props[a]
    maids <- ma[, 1] == m & ma[, 2] == a
    X[maids, ] <- matrix(rnorm(nma * p, mean = 0, sd = 1), nrow = nma, ncol = p)
    y[maids] <- X[maids, ] %*% th.true[, a, m] + rnorm(nma, mean = 0, sd = sd.yi)
  }
}


eta_init <- 100
# lambda <- 1e-1
rho <- 1
adj <- 1
eps <- 1e-6
maxiter <- 1e4
# th_init <- NULL
# delta_init <- NULL
verbose <- TRUE

lambda_length <- eta_length <- 10
lambda_ratio <- eta_ratio <- 1e-2

# th_init <- th.true
#
# fit <- HDfair(
#   X = X,
#   y = y,
#   ma = ma,
#   lambda = lambda,
#   eta = eta,
#   rho = rho,
#   th_init=th_init*10,
#   delta_init=delta_init,
#   adj=adj,
#   eps=eps,
#   maxiter=maxiter,
#   verbose = verbose
# )
#
# fit$th[1:10, , ]; th.true[1:10, , ]
# fit$iter
#
#
#
# sp_lambda <- HDfair_sp_lambda(
#   X = X,
#   y = y,
#   ma = ma,
#   lambda_length = lambda_length,
#   lambda_ratio = lambda_ratio,
#   eta = eta,
#   rho = rho,
#   adj=adj,
#   eps=eps,
#   maxiter=maxiter,
#   verbose = verbose
# )
# plot_sp_lambda(sp_lambda)
# sp_lambda$iterations
#
#
#
#
# sp_eta <- HDfair_sp_eta(
#   X = X,
#   y = y,
#   ma = ma,
#   lambda = lambda,
#   eta_length = eta_length,
#   eta_ratio = eta_ratio,
#   rho = rho,
#   adj=adj,
#   eps=eps,
#   maxiter=maxiter,
#   verbose = verbose
# )
# plot_sp_eta(sp_eta)



cv_lambda <- HDfair_cv_lambda(
  X = X,
  y = y,
  ma = ma,
  lambda_length = lambda_length,
  lambda_ratio = lambda_ratio,
  nfold = 5,
  eta = eta_init,
  rho = rho,
  adj=adj,
  eps=eps,
  maxiter=maxiter,
  verbose = verbose
)
plot_cv_lambda(cv_lambda)
plot_sp_lambda(cv_lambda$sp)


cv_eta <- HDfair_cv_eta(
  X = X,
  y = y,
  ma = ma,
  lambda = cv_lambda$lambda.min,
  eta_length = eta_length,
  eta_ratio = eta_ratio,
  nfold = 5,
  rho = rho,
  adj=adj,
  eps=eps,
  maxiter=maxiter,
  verbose = verbose
)
plot_cv_eta(cv_eta)
plot_sp_eta(cv_eta$sp)

cv_lambda$sp$estimates[,,,cv_lambda$index[1]]
cv_eta$sp$estimates[,,,cv_eta$index[1]]
