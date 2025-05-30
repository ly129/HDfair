A <- 3L
M <- 1L
p <- 100L
n.M <- rep(100, M)
N <- sum(n.M)
fixed.th <- FALSE
# weighted true: loss function is sum of per-task MSE
# weighted false: loss function is sum of per-task sum of squared errors
weighted <- FALSE

# # AR-1 correlation in X
# rho <- 0
# Sigma <- toeplitz(rho^(0:(p-1)))


# between site heterogeneity
sd.th.m <- 0
if (M == 1) sd.th.m <- 0
# between group heterogeneity
sd.th.a <- 0.2

# used to generate true theta.
## true.theta = th.base + group heterogeneity + site heterogeneity
th.base <- rep(0.3, 10)
# th.base <- 1:10 * 0.1
p.nz <- length(th.base)

# proportions of groups
props <- matrix(c(0.7, 0.2, 0.1), nrow = M)
# props <- matrix(c(1/3, 1/3, 1/3), nrow = M)
# props <- matrix(c(0.6, 0.1, 0.1, 0.1, 0.1), nrow = M)
if (M > 1) {
  props <- matrix(props, nrow = M, ncol = A, byrow = TRUE)
  props[2, ] <- props[2, c(2, 1, 3)]
}
# group indicator
ma = matrix(0,N,2)
ma[,1] = rep(1:M, n.M)
for (m in 1:M) {
  ma[ma[, 1] == m, 2] <- rep(1:A, n.M[m] * props[m, ])
}

table(ma[, 1], ma[, 2])

# gaussian noise
sd.yi <- 2

th.true <- array(0, dim = c(p,A,M))
for (m in seq(M)) {
  epsilon.m <- rep(0, p.nz)
  if (M > 1) {
    # epsilon.m <- runif(n = p.nz, min = 0, max = sd.th.m)
    epsilon.m <- rnorm(n = p.nz, mean = 0, sd = sd.th.m)
  }
  for (a in seq(A)) {
    # epsilon.a <- runif(n = p.nz, min = 0, max = sd.th.a)
    epsilon.a <- rnorm(n = p.nz, mean = 0, sd = sd.th.a)
    th.true[1:p.nz, a, m] <- th.base + epsilon.m + epsilon.a
  }
}

th.true[1:min(15, p.nz), ,]


# generate X
X <- X.test <- matrix(NA, N, p)
y <- y.test <- numeric(N)

# x, y and group (indicator) are all lists of length M
for (m in 1:M) {
  for (a in 1:A) {
    nma <- n.M[m] * props[m, a]
    maids <- ma[, 1] == m & ma[, 2] == a

    # X[maids, ] <- mvrnorm(n = nma, mu = rep(0, p), Sigma = Sigma)
    X[maids, ] <- matrix(rnorm(nma * p, mean = 0, sd = 1), nrow = nma, ncol = p)
    y[maids] <- X[maids, ] %*% th.true[, a, m] + rnorm(nma, mean = 0, sd = sd.yi)

    # X.test[maids, ] <- mvrnorm(n = nma, mu = rep(0, p), Sigma = Sigma)
    X.test[maids, ] <- matrix(rnorm(nma * p, mean = 0, sd = 1), nrow = nma, ncol = p)
    y.test[maids] <- X.test[maids, ] %*% th.true[, a, m] + rnorm(nma, mean = 0, sd = sd.yi)
  }
}




rho <- 1
adj <- 1
eps <- 1e-6
maxiter <- 1e4
verbose <- FALSE
lambda_length <- eta_length <- 20
lambda_ratio <- 1e-2
eta_ratio <- 1e-2
nfold <- 5
