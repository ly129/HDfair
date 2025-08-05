
rm(list = ls())

library(glmnet)
library(RMTL)
library(HDfair)
library(reticulate)
sklearn       <- import("sklearn.linear_model")
MultiTaskLassoCV <- sklearn$MultiTaskLassoCV
MultiTaskLasso <- sklearn$MultiTaskLasso

simulate_multitask_data <- function(
    n = 200,     # sample sizes for 3 groups
    p = 100,                    # number of predictors
    s = 10,                     # number of nonzero features (shared support)
    n_task = 3,
    beta_center = 0.6,          # base signal for nonzero coefficients
    hetero_sd = 0.3,            # coefficient heterogeneity across tasks
    sigma = 0.5,                # residual SD (same across tasks)
    seed = 1
) {
  set.seed(seed)
  support_idx <- 1:s   # shared support

  beta_mat <- matrix(0, nrow = p, ncol = n_task)

  for (j in 1:n_task) {
    beta_mat[support_idx, j] <- rep(beta_center, s) + rnorm(s, mean = 0, sd = hetero_sd)
  }

  X_list <- list()
  y_list <- list()

  X_t <- matrix(rnorm(n * p), nrow = n, ncol = p)

  for (t in 1:n_task) {
    y_t <- X_t %*% beta_mat[, t] + rnorm(n, sd = sigma)

    X_list[[t]] <- X_t
    y_list[[t]] <- y_t
  }

  return(list(
    X_list = X_list,
    y_list = y_list,
    beta_mat = beta_mat,
    support = support_idx
  ))
}



n_sims   <- 200
nfolds   <- 5
n <- 500
p <- 100
s <- 10
beta_center <- 0.3
hetero_sd <- 0.1
sigma <- 2
seed <- sample.int(1e8, size = 1)
n_task <- 3


data <- simulate_multitask_data(
  p           = p,
  s           = s,
  n = n,
  n_task = n_task,
  beta_center = beta_center,
  hetero_sd   = hetero_sd,
  sigma       = sigma,
  seed        = seed
)
X_list <- data$X_list
y_list <- data$y_list
true_support <- data$support

# 2) POOLED LASSO
X_stack <- do.call(rbind, X_list)
y_stack <- unlist(y_list)
# weights = 1/n_t
w_vec   <- unlist(mapply(function(n) rep(1/n, n), n = rep(n, n_task), SIMPLIFY = FALSE))



ma <- matrix(nrow = n * n_task, ncol = 2)
ma[, 1] <- 1
ma[, 2] <- rep(1:n_task, each = n)

sp_lambda <- HDfair_sp_lambda(
  X = X_stack,
  y = y_stack,
  ma = ma,
  lambda_length = 50,
  lambda_ratio = 1e-2,
  eta = p * 1000, # large eta s.t. fairness constraint is inactive
  rho = 1,
  weighted = TRUE,
  adj=1,
  eps=1e-6,
  maxiter=1e3,
  verbose = FALSE
)
par(mfrow = c(1, n_task))
for (tt in 1:n_task) {
  hdfair.tt <- sp_lambda$estimates[, tt, 1,  ]
  matplot(x = sp_lambda$lambda,
          y = t(hdfair.tt),
          type = "l",
          log = "x",
          ylim = c(-0.2, 1))
}


sp.rmtl <- array(dim = c(p, n_task, 50))
for (ii in 1:50) {
  lam1 <- sp_lambda$lambda[ii]
  fit.rmtl <- MTL(
    X             = X_list,
    Y             = y_list,
    type          = "Regression",
    Regularization = "L21",
    Lam1          = lam1,
    Lam2          = 0,
    opts          = list(init = 0, tol = 1e-10, maxIter = 1000)
  )
  sp.rmtl[,,ii] <- fit.rmtl$W
}
par(mfrow = c(1, n_task))
for (tt in 1:n_task) {
  rmtl.tt <- sp.rmtl[, tt, ]
  matplot(x = sp_lambda$lambda,
          y = t(rmtl.tt),
          type = "l",
          log = "x",
          ylim = c(-0.2, 1))
}

sp.sklearn <- array(dim = c(p, n_task, 50))
for (ii in 1:50) {
  lam <- sp_lambda$lambda[ii]*n_task
  sklearn.mtl <- MultiTaskLasso(alpha = lam, fit_intercept = FALSE)
  sklearn.mtl$fit(X_list[[1]], matrix(unlist(y_list), ncol = n_task))
  sp.sklearn[,,ii] <- t(sklearn.mtl$coef_)
}
par(mfrow = c(1, n_task))
for (tt in 1:n_task) {
  sklearn.tt <- sp.sklearn[, tt, ]
  matplot(x = sp_lambda$lambda*n_task,
          y = t(sklearn.tt),
          type = "l",
          log = "x",
          ylim = c(-0.2, 1))
}

mtl_cv <- MultiTaskLassoCV(fit_intercept = FALSE, alphas = sp_lambda$lambda)
mtl_cv$fit(X_list[[1]], matrix(unlist(y_list), ncol = n_task))

# 1) Pull alphas and mse_path from Python
alphas_py   <- mtl_cv$alphas_       # numpy array of shape (n_alphas,)
mse_path_py <- mtl_cv$mse_path_     # numpy array of shape (n_alphas, n_folds)

# 2) Convert to R
alphas   <- py_to_r(alphas_py)
mse_path <- py_to_r(mse_path_py)

# 3) Compute mean and standard error across folds
mean_mse <- rowMeans(mse_path)
se_mse   <- apply(mse_path, 1, sd) / sqrt(ncol(mse_path))

# 4) Plot CV‐curve
par(mfrow = c(1, 1))
plot(log(alphas), mean_mse,
     type = "b",                     # lines + points
     # log  = "x",                     # log‐scale on x
     xlim = range(log(alphas)),      # reverse x‐axis so large→small
     xlab = expression(alpha),
     ylab = "Mean CV MSE",
     main = "MultiTaskLassoCV: CV error vs. alpha")

# 5) Add error bars for ±1 SE
arrows(log(alphas), mean_mse + se_mse,
       log(alphas), mean_mse - se_mse,
       angle = 90, code = 3, length = 0.02)

# (Optionally) highlight the best alpha
abline(v = log(mtl_cv$alpha_), lty = 2, col = "red")




cv_lambda <- HDfair_cv_lambda(
  X = X_stack,
  y = y_stack,
  ma = ma,
  lambda_length = 50,
  lambda_ratio = 1e-2,
  nfold = 5,
  eta = p * 1000, # large eta s.t. fairness constraint is inactive
  rho = 1,
  weighted = TRUE,
  adj=1,
  eps=1e-6,
  maxiter=1e3,
  verbose = FALSE
)
par(mfrow = c(1, 1))
plot_cv_lambda(cv_lambda)



cv_lasso.p <- cv.glmnet(
  x = X_stack,
  y = y_stack,
  foldid = cv_lambda$foldid,
  lambda.min.ratio = 1e-2,
  intercept = FALSE)
par(mfrow = c(1, 1))
plot(cv_lasso.p)

par(mfrow = c(1, 1))
plot(cv_lasso.p$glmnet.fit, xvar = "lambda")

# sp_lambda2 <- HDfair_sp_lambda(
#   X = X_stack,
#   y = y_stack,
#   ma = ma,
#   lambda_length = 20,
#   lambda_ratio = 1e-2,
#   eta = 1e-6,
#   rho = 1,
#   weighted = FALSE,
#   adj=1,
#   eps=1e-6,
#   maxiter=1e4,
#   verbose = TRUE
# )
# par(mfrow = c(1, 1))
# for (tt in 1:n_task) {
#   hdfair.tt <- sp_lambda2$estimates[, tt, 1,  ]
#   matplot(x = sp_lambda2$lambda,
#           y = t(hdfair.tt),
#           type = "l",
#           log = "x",
#           ylim = c(-0.2, 0.6))
# }





sp_lambda <- HDfair_sp_lambda(
  X = X_stack,
  y = y_stack,
  ma = ma,
  eta = 1e-3, # large eta s.t. fairness constraint is inactive
  rho = 1,
  weighted = TRUE,
  adj=1,
  eps=1e-6,
  maxiter=1e4,
  verbose = TRUE
)
plot_sp_lambda(sp_lambda)

lasso.p <- glmnet(x = X_stack,
                  y = y_stack,
                  lambda = sp_lambda$lambdas * sqrt(n_task),
                  weights = w_vec,
                  intercept = FALSE)
plot(lasso.p, xvar = "lambda")








sp_lambda <- HDfair_sp_lambda(
  X = X_stack,
  y = y_stack,
  ma = ma,
  eta = p * 1e6, # large eta s.t. fairness constraint is inactive
  rho = 1,
  weighted = FALSE,
  adj=1,
  eps=1e-6,
  maxiter=1e4,
  verbose = TRUE
)

index <- 5
sum(sp_lambda$estimates[,1,1,index] != 0)


sp_eta <- HDfair::HDfair_sp_eta(
  X = X_stack,
  y = y_stack,
  ma = ma,
  eta_length = 20,
  eta_ratio = 1e-2,
  lambda = sp_lambda$lambdas[index],
  rho = 1,
  weighted = FALSE,
  adj=1,
  eps=1e-6,
  maxiter=1e4,
  verbose = TRUE
)

colSums(sp_eta$estimates[,1,1,] != 0)
