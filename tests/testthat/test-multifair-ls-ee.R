# set.seed(20230305)

A <- 3L
M <- 2L
p <- 8L
p.nz <- floor(0.25 * p)
n.M <- 200L * seq(M)

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
    th[[i]][1:p.nz, j] <- th.base #+ epsilon.m + epsilon.a
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
lam <- 1
nlam <- length(lam)
intercept <- FALSE
# crit <- "BGL"
# eta <- 2
crit <- "metric"
eta_fr <- 1e-1
eta_ee <- 1e-3
rho_fr <- 1
rho_ee <- 5e1
reg <- "group-lasso"
group <- id.grp
stepsize <- 1e-4

fit.m <- multifairee(
  X,
  y,
  group = group,
  lam,
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
  eps = 1e-20,
  verbose = TRUE
)

# solution path
par(mfrow = c(M, A))
for (m in 1:M) {
  for (a in 1:A) {
    matplot(x = lam, y = t(fit.m$Estimates[,m,a,]),
            main = paste(c(m, a)), type = "l", log = "x")
  }
}

# check bgl fairness
if (crit == "BGL") {
  bgl.check.m <- matrix(0, nrow = nlam, ncol = A)
  for (l in 1:nlam) {
    for (a in 1:A) {
      for (m in 1:M) {
        id.ma <- id.grp[[m]] == a
        Xma <- X[[m]][id.ma,]
        yma <- y[[m]][id.ma]
        nma <- sum(id.ma)
        thma <- fit.m$Estimates[, m, a, l]
        bglma <- fair_bgl(X_ma = Xma, y_ma = yma, th_ma = thma,
                          type = type)
        bgl.check.m[l, a] <- bgl.check.m[l, a] + bglma$fair
      }
    }
  }
  View(bgl.check.m/M)
}


# check metric fairness
if (crit == "metric") {
  metricfair <- apply(fit.m$Estimates, MARGIN = 4, FUN = metric_check, p, M, A)
  View(t(metricfair))
}

fit.m$Iterations

# check ee
ee_max <- array(dim = c(M, A, nlam))
for (l in 1:nlam) {
  for (m in 1:M) {
    for (a in 1:A) {
      Xma <- X[[m]][id.grp[[m]] == a, ]
      yma <- y[[m]][id.grp[[m]] == a]
      thmal <- fit.m$Estimates[,m,a,l]
      ee_max[m, a, l] <- max(abs(ufunc(Xma, yma, thmal)))
    }
  }
}
ee_max
