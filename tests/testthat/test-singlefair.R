set.seed(20230305)

A <- 3L
M <- 2L
p <- 20L
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

# fit
lam <- 10^(seq(-3.5, 0.5, 0.1))
nlam <- length(lam)
eta <- 1
type <- "continuous"
intercept <- FALSE
crit <- "BGL"
rho <- 1
reg <- "group-lasso"

fit.s <- singlefair(
  X = X,
  y = y,
  group = id.grp,
  lam = lam,
  eta = eta,
  stepsize = 0.1,
  type = type,
  reg = reg,
  crit = crit,
  rho = rho,
  maxit = 1e4,
  eps = 1e-20,
  verbose = TRUE
)

# solution path
par(mfrow = c(1, M))
for (m in 1:M) {
  matplot(x = lam, y = t(fit.s$Estimates[,m,]),
          main = paste(m), type = "l", log = "x")
}

# check fairness
bgl.check.s <- matrix(0, nrow = nlam, ncol = A)
for (l in 1:nlam) {
  for (a in 1:A) {
    for (m in 1:M) {
      id.ma <- id.grp[[m]] == a
      Xma <- X[[m]][id.ma,]
      yma <- y[[m]][id.ma]
      nma <- sum(id.ma)
      thma <- fit.s$Estimates[, m, l]
      bglma <- fair_bgl(X_ma = Xma, y_ma = yma,
                        n_ma = nma, th_ma = thma,
                        type = "continuous")
      bgl.check.s[l, ] <- bgl.check.s[l, ] + bglma$fair
    }
  }
}
View(bgl.check.s/M/A)
fit.s$Iterations
