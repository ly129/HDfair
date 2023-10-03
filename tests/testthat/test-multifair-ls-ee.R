# set.seed(20230305)

A <- 3L
M <- 2L
p <- 20L
p.nz <- floor(0.25 * p)
n.M <- 800L * seq(M)

# generate X
X <- sapply(n.M,
            function(x) {
              matrix(
                rnorm(n = p * x, mean = 0, sd = 1),
                nrow = x, ncol = p, byrow = TRUE
              )
            }, simplify = FALSE
)

# generate true theta
th <- vector("list", M)
th.base <- 1
for (i in seq(M)) {
  th[[i]] <- matrix(0, nrow = p, ncol = A)
  epsilon.m <- rnorm(n = p.nz, mean = 0, sd = 0.3)
  for (j in seq(A)) {
    epsilon.a <- rnorm(n = p.nz, mean = 0, sd = 0.4)
    th[[i]][1:p.nz, j] <- th.base - 0.1 * i + 0.2 * j#+ epsilon.m + epsilon.a
  }
}

# generate groups
pt.maj <- seq(0.6, 0.8, length.out = M)
pts <- vector(mode = "list", length = M)
for (m in 1:M) {
  pts[[m]] <- c(pt.maj[m], rep((1 - pt.maj[m])/(A - 1), times = A - 1))
}

id.cnt <- mapply(FUN = function(x, p) {
  round( x * p )
}, n.M, pts, SIMPLIFY = FALSE)
id.grp <- lapply(id.cnt, FUN = function(x) {sample(rep(seq(A), times = x))})

# generate y
y <- vector(mode = "list", length = M)
for (m in seq(M)) {
  nm <- n.M[m]
  y[[m]] <- numeric(nm)
  for (i in seq(nm)) {
    a.tmp <- id.grp[[m]][i]
    X.tmp <- X[[m]][i, ]
    th.tmp <- th[[m]][, a.tmp]
    y[[m]][i] <- c(crossprod(X.tmp, th.tmp)) + rnorm(1, mean = 0, sd = 1)
  }
}


# estimating function
custom_esti_func <- ufunc <- function(X_ma, y_ma, th_ma) {
  n_ma <- length(y_ma)
  resid <- y_ma - X_ma %*% th_ma
  ufunc <- -c(crossprod(X_ma, resid)/n_ma)
  return(ufunc)
}

custom_esti_func_gr <- ugrad <- function(X_ma, y_ma, th_ma) {
  n_ma <- length(y_ma)
  loss_gr <- crossprod(X_ma)/n_ma
  return(loss_gr)
}

# fit
# lam <- 10^(seq(2, 7, 0.1))
# lam <- 1
intercept <- FALSE
# crit <- "BGL"
# eta <- 2
crit <- "metric"
eta_fr <- 0.1
eta_ee <- 10^(seq(2, -4, -0.5))
neta_ee <- length(eta_ee)
rho_fr <- 1
rho_ee <- 1
reg <- "group-lasso"
group <- id.grp
stepsize <- 1e-4

fit.m <- multifairee(
  X,
  y,
  group = group,
  eta_fr,
  eta_ee,
  stepsize = stepsize,
  intercept = FALSE,
  custom_esti_func = ufunc,
  custom_esti_func_gr = ugrad,
  reg = "group-lasso",
  crit = crit,
  rho_fr = rho_fr,
  rho_ee = rho_ee,
  maxit = 1e4,
  eps = 1e-10,
  verbose = TRUE
)

# solution path
par(mfrow = c(M, A))
for (m in 1:M) {
  for (a in 1:A) {
    matplot(x = eta_ee, y = t(fit.m$Estimates[,m,a,]),
            main = paste(c(m, a)), type = "l", log = "x")
  }
}

# check metric fairness
if (crit == "metric") {
  metricfair <- apply(fit.m$Estimates, MARGIN = 4, FUN = metric_check, p, M, A)
  View(t(metricfair))
}

fit.m$Iterations

# check ee
ee.max <- array(dim = c(M, A, neta_ee))
for (l in 1:neta_ee) {
  for (m in 1:M) {
    for (a in 1:A) {
      Xma <- X[[m]][id.grp[[m]] == a, ]
      yma <- y[[m]][id.grp[[m]] == a]
      thmal <- fit.m$Estimates[,m,a,l]
      ee.max[m, a, l] <- max(abs(ufunc(Xma, yma, thmal)))
    }
  }
}
ee.max





# compare to FREE-based algorithm
lam <- 10^(seq(0, -5, -0.1))
nlam <- length(lam)
stepsize.free <- 1
fit.m.free <- multifair(X,
                        y,
                        group = group,
                        lam = lam,
                        eta = eta_fr,
                        stepsize = stepsize.free,
                        intercept = FALSE,
                        type = "custom",
                        custom_esti_func = ufunc,
                        reg = "group-lasso",
                        crit = "metric",
                        rho = 1,
                        maxit = 1e3,
                        eps = 1e-20,
                        verbose = TRUE)

par(mfrow = c(M, A))
for (m in 1:M) {
  for (a in 1:A) {
    matplot(x = lam, y = t(fit.m.free$Estimates[,m,a,]),
            main = paste(c(m, a)), type = "l", log = "x")
  }
}

# check metric fairness
if (crit == "metric") {
  metricfair <- apply(fit.m.free$Estimates, MARGIN = 4, FUN = metric_check, p, M, A)
  View(t(metricfair))
}


# check ee
ee.max.free <- array(dim = c(M, A, nlam))
for (l in 1:nlam) {
  for (m in 1:M) {
    for (a in 1:A) {
      Xma <- X[[m]][id.grp[[m]] == a, ]
      yma <- y[[m]][id.grp[[m]] == a]
      thmal <- fit.m.free$Estimates[,m,a,l]
      ee.max.free[m, a, l] <- mean(abs(ufunc(Xma, yma, thmal)))
    }
  }
}
ee.max.free

fit.m.free$Iterations
