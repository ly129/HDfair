A <- 3L
M <- 1L
p <- 20L
p.nz <- floor(0.25 * p)
n.M <- 1000L * seq(M)

# generate true theta
th <- vector("list", M)
th.base <- 0.5
for (i in seq(M)) {
  th[[i]] <- matrix(0, nrow = p, ncol = A)
  # epsilon.m <- rnorm(n = p.nz, mean = 0, sd = 0.3)
  # for (j in seq(A)) {
  #   epsilon.a <- rnorm(n = p.nz, mean = 0, sd = 0.4)
  #   th[[i]][1:p.nz, j] <- th.base + (0.1 * i - 0.03 * j) * rep(c(1, -1), length.out = p.nz) * ((1:p.nz) - mean(p.nz)) + epsilon.m + epsilon.a
  # }
  th[[i]][1:p.nz, 1] <- rep(c(0.2, 0.5, 0.8), length.out = p.nz) # 0.5
  th[[i]][1:p.nz, 2] <-rep(c(0.5, 0.8, 0.2), length.out = p.nz) # 0.2
  th[[i]][1:p.nz, 3] <- rep(c(0.8, 0.2, 0.5), length.out = p.nz) # 0.8 # rep(c(0.6, 0.7, 0.8, 0.9, 1.0), length.out = p.nz)
}
head(th[[1]], p.nz)

# generate groups
pt.maj <- rep(0.90, M)
pts <- vector(mode = "list", length = M)

for (m in 1:M) {
  pts[[m]] <- c(pt.maj[m], 0.06, 0.04)
}

id.cnt <- mapply(FUN = function(x, p) {
  round( x * p )
}, n.M, pts, SIMPLIFY = FALSE)
id.grp <- lapply(id.cnt, FUN = function(x) {sample(rep(seq(A), times = x))})


x <- sapply(n.M,
            function(x) {
              matrix(
                rnorm(n = p * x, mean = 0, sd = 1),
                nrow = x, ncol = p, byrow = TRUE
              )
            }, simplify = FALSE
)

x.test <- sapply(n.M,
                 function(x) {
                   matrix(
                     rnorm(n = p * x, mean = 0, sd = 1),
                     nrow = x, ncol = p, byrow = TRUE
                   )
                 }, simplify = FALSE
)

y <- y.test <- vector(mode = "list", length = M)
for (m in seq(M)) {
  nm <- n.M[m]
  y[[m]] <- numeric(nm)
  y.test[[m]] <- numeric(nm)
  for (i in seq(nm)) {
    a.tmp <- id.grp[[m]][i]
    x.tmp <- x[[m]][i, ]
    x.test.tmp <- x.test[[m]][i, ]
    th.tmp <- th[[m]][, a.tmp]
    y[[m]][i] <- c(crossprod(x.tmp, th.tmp)) + rnorm(1, mean = 0, sd = 1)
    y.test[[m]][i] <- c(crossprod(x.test.tmp, th.tmp)) + rnorm(1, mean = 0, sd = 1)
  }
}


reg <- "group-lasso"


results <- multifair_base(
  x = x,
  y = y,
  group = id.grp,
  lambda = 1e-6,
  eta = 0.2,
  stepsize = 0.1,
  intercept = FALSE,
  # theta_init = array(0.5, dim = c(p, M, A)),
  type = "continuous",
  reg = reg,
  crit = "metric",
  rho = 0.1,
  maxit = 1e4,
  tol = 1e-8,
  verbose = FALSE
)
results$Iterations
(metric.fair <- metric_check(results$Estimates))


sp.results.lam <- multifair_sp_lambda(
  x = x,
  y = y,
  group = id.grp,
  lambda_seq = 10^seq(1, -5, -1),
  eta = 1,
  stepsize = 0.1,
  intercept = FALSE,
  # theta_init = array(0.5, dim = c(p, M, A)),
  type = "continuous",
  reg = reg,
  crit = "metric",
  rho = 0.1,
  maxit = 1e3,
  tol = 1e-8,
  verbose = FALSE
)

sp.results.lam$Iterations
plot_multifair_sp_lam(sp.results.lam)

sp.results.eta <- multifair_sp_eta(
  x = x,
  y = y,
  group = id.grp,
  lambda = 1e-1,
  eta_seq = c(2, seq(1, 0.1, -0.1), 0.01, 0.001),
  stepsize = 0.1,
  intercept = FALSE,
  # theta_init = array(0.5, dim = c(p, M, A)),
  type = "continuous",
  reg = reg,
  crit = "metric",
  rho = 6,
  maxit = 1e4,
  tol = 1e-8,
  verbose = FALSE
)

sp.results.eta$Iterations
plot_multifair_sp_eta(sp.results.eta, n = seq(p.nz))

################################################################################
# weights <- ifelse(id.grp[[1]] == 1, 1/0.9,
#                   ifelse(id.grp[[1]] == 2,
#                          1/0.06,
#                          1/0.04))
# weights <- weights/sum(weights) * n.M[1]
#
# weighted.glm <- glmnet::glmnet(
#   x = x[[1]],
#   y = y[[1]],
#   lambda = 0,
#   intercept = FALSE,
#   weights = weights,
#   thresh = 1e-20
# )
# # plot(weighted.glm)
#
# th.glmnet <- as.vector(coef(weighted.glm)[-1])
# th.glmnet
#
# th.init <- array(dim = c(p, M, A))
# th.init[, 1:M, 1:A] <- th.glmnet
#
# sp.results.eta <- multifair_sp_eta(
#   x = x,
#   y = y,
#   group = id.grp,
#   lambda = 0,
#   eta_seq = c(2, seq(1, 0.1, -0.1), 0.01, 0.001, 0.0001, 0.00001, 0.000001),
#   stepsize = 0.1,
#   intercept = FALSE,
#   # theta_init = th.init,
#   type = "continuous",
#   reg = reg,
#   crit = "metric",
#   rho = 5,
#   maxit = 1e4,
#   tol = 1e-8,
#   verbose = FALSE
# )
# sp.results.eta$Estimates[,1,1,1]
#
# sp.results.eta$Iterations
# plot_multifair_sp_eta(sp.results.eta, n = seq(p.nz))
#
# apply(sp.results.eta$Estimates, 4, FUN = metric_check)
#
# matplot(x = sp.results.eta$Etas,
#         y = t(sp.results.eta$Estimates[,1,1,]),
#         type = "l", log = "x", lwd = 2)
# abline(h = th.glmnet[1:5], lty = 1:5, col = 1:6, lwd = 2)
#
# matplot(x = sp.results.eta$Etas,
#         y = t(sp.results.eta$Estimates[,1,3,]),
#         type = "l", log = "x", lwd = 2)
# abline(h = th.glmnet[1:5], lty = 1:5, col = 1:6, lwd = 2)
################################################################################


n.alternate <- 10

etasss <- lamsss <- rep(NA, n.alternate)
etasss[1] <- 10

for (nn in 1:n.alternate) {
  cv.results.lam <- multifair_cv_lambda(
    x = x,
    y = y,
    group = id.grp,
    lambda_seq = 10^seq(0, -2, -0.1),
    eta = etasss[nn],
    nfolds = 5L,
    foldid = NULL,
    stepsize = 0.1,
    intercept = FALSE,
    theta_init = NULL,
    type = "continuous",
    reg = reg,
    crit = "metric",
    rho = 1,
    maxit = 1e3,
    tol = 1e-6,
    warmstart = TRUE,
    verbose = FALSE
  )
  (lamsss[nn] <- cv.results.lam$lambda.1se)
  plot_multifair_cv_lambda(cv.results.lam)


  cv.results.eta <- multifair_cv_eta(
    x = x,
    y = y,
    group = id.grp,
    lambda = lamsss[nn],
    eta_seq = c(seq(0.5, 0.1, -0.05), 0.05, 0.01),
    nfolds = 5L,
    foldid = NULL,
    stepsize = 0.1,
    intercept = FALSE,
    theta_init = NULL,
    type = "continuous",
    reg = reg,
    crit = "metric",
    rho = 1,
    maxit = 1e3,
    tol = 1e-6,
    warmstart = TRUE,
    verbose = FALSE
  )
  (etasss[nn + 1] <- cv.results.eta$eta.min)
  plot_multifair_cv_eta(cv.results.eta)
}

lamsss; etasss

# results2 <- multifair(
#   x = x,
#   y = y,
#   group = id.grp,
#   lam = 1e-6,
#   eta = 1,
#   stepsize = 0.1,
#   intercept = FALSE,
#   # theta_init = rep(1, p + 1),
#   type = "continuous",
#   reg = reg,
#   crit = "metric",
#   rho = 0.1,
#   maxit = 1e3,
#   eps = 1e-6,
#   verbose = FALSE
# )
# results2$Iterations
#
# identical(results$Estimates, results2$Estimates[,,,1])
# (metric.fair2 <- metric_check(results2$Estimates[,,,1], p, M, A))
#
# microbenchmark::microbenchmark(results <- HDfair::multifair_base(
#   x = x,
#   y = y,
#   group = id.grp,
#   lam = 1e-6,
#   eta = 1,
#   stepsize = 0.1,
#   intercept = FALSE,
#   # theta_init = array(0.5, dim = c(p, M, A)),
#   type = "continuous",
#   reg = reg,
#   crit = "metric",
#   rho = 0.1,
#   maxit = 1e3,
#   eps = 1e-6,
#   verbose = FALSE
# ),
# results2 <- multifair(
#   x = x,
#   y = y,
#   group = id.grp,
#   lam = 1e-6,
#   eta = 1,
#   stepsize = 0.1,
#   intercept = FALSE,
#   # theta_init = rep(1, p + 1),
#   type = "continuous",
#   reg = reg,
#   crit = "metric",
#   rho = 0.1,
#   maxit = 1e3,
#   eps = 1e-6,
#   verbose = FALSE
# ))
