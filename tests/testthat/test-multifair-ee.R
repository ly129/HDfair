# set.seed(20230305)

A <- 3L
M <- 2L
p <- 20L
p.nz <- floor(0.25 * p)
n.M <- 200L * seq(M)

library(survival)
data_gen <- function(X, sigma, beta_true, censor, c = 0) {
  p <- length(beta_true)
  noise <- rnorm(n = 1, mean = 0, sd = sigma)
  TT <- X %*% beta_true + noise
  CC <- NULL
  YY <- TT
  Delta <- NULL
  if (censor) {
    if (c >= 0) {
      CC <- runif(1, min = 0, max = c)
    } else {
      CC <- runif(1, min = c, max = 0)
    }
    Delta <- as.integer((TT < CC))
    YY <- pmin(TT, CC)
  }
  return(c(YY, Delta))
}

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
yy <- delta <- vector(mode = "list", length = M)
for (m in seq(M)) {
  nm <- n.M[m]
  yy[[m]] <- numeric(nm)
  delta[[m]] <- integer(nm)
  for (i in seq(nm)) {
    a.tmp <- id.grp[[m]][i]
    X.tmp <- X[[m]][i, ]
    th.tmp <- th[[m]][, a.tmp]
    dat.tmp <- data_gen(X = X.tmp, sigma = 1,
                        beta_true = th.tmp,
                        censor = TRUE, c = 5)
    yy[[m]][i] <- dat.tmp[1]
    delta[[m]][i] <- dat.tmp[2]
  }
}

y <- mapply(FUN = Surv, time = yy, event = delta)

# estimating function
ufunc <- function(x, y, b) {
  n <- nrow(x)
  p <- ncol(x)
  yy <- y[, 1]
  status <- y[, 2]
  timeorig <- yy
  dummystrat <- factor(rep(1, n))
  ypred <- x %*% b
  ehat <- timeorig - ypred
  state <- status
  state[ehat == max(ehat)] <- 1
  S <- structure(cbind(ehat, state),
                 class = "Surv", type = "right")
  KM.ehat <- survfitKM(dummystrat, S,
                       conf.type = "none",
                       se.fit = FALSE)
  n.risk <- KM.ehat$n.risk
  surv <- KM.ehat$surv
  repeats <- c(diff( - n.risk),
               n.risk[length(n.risk)])
  surv <- rep(surv, repeats)
  w <-  - diff(c(1, surv))
  m <- order(ehat,  - status)
  bla <- cumsum((w * ehat[m]))
  bla <- (bla[length(bla)] - bla)/(surv + state[m])
  bl <- bla
  bl[(1 : n)[m]] <- bla
  yhat <- if (p == 0) bl else x %*% b + bl
  yy[state == 0] <- yhat[state == 0]
  return( c(- crossprod(x, (yy - x %*% b))/n) )
}

# fit
lam <- 10^(seq(-3.5, 0.5, 0.1))
nlam <- length(lam)
type <- "custom"
intercept <- FALSE
crit <- "metric"
eta <- 0.1
rho <- 1
reg <- "group-lasso"

fit.m <- multifair(
  X = X,
  y = y,
  group = id.grp,
  lam = lam,
  eta = eta,
  stepsize = 0.1,
  type = type,
  custom_esti_func = ufunc,
  reg = reg,
  crit = crit,
  rho = rho,
  maxit = 1e3,
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
