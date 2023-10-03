multifair_sp_lambda <- function(
    x,
    y,
    group,
    lambda_seq,
    eta,
    stepsize,
    intercept = FALSE,
    theta_init = NULL,
    type = "continuous",
    reg = "group-lasso",
    crit = "metric",
    rho = 1,
    maxit = 1e3,
    tol = 1e-6,
    warmstart = TRUE,
    verbose = FALSE
) {
  # Checks
  # M
  M <- length(x)
  # p
  p <- sapply(x, ncol, simplify = TRUE)
  if (length(unique(p)) == 1) {
    p <- unique(p)
  } else {
    Stop("Different datasets have different number of predictors.")
  }
  # n
  n_M <- sapply(x, nrow, simplify = TRUE)
  # A
  A <- as.integer(max(unique(unlist(group))))

  # intercept
  if (intercept) p <- p + 1

  nlam <- length(lambda_seq)
  lambda_seq <- sort(lambda_seq, decreasing = TRUE)

  Theta <- array(NA, dim = c(p, M, A, nlam))
  Iterations <- integer(nlam)

  if (is.null(theta_init)) theta_init <- array(0, dim = c(p, M, A))

  fit <- multifair_base(
    x = x,
    y = y,
    group = group,
    lambda = lambda_seq[1],
    eta = eta,
    stepsize = stepsize,
    intercept = intercept,
    theta_init = theta_init,
    type = type,
    reg = reg,
    crit = crit,
    rho = rho,
    maxit = maxit,
    tol = tol,
    verbose = verbose
  )

  Theta[, , , 1] <- fit$Estimates
  Iterations[1] <- fit$Iterations

  for (l in 2:nlam) {
    lam_tmp <- lambda_seq[l]

    if (warmstart) {
      th_init <- fit$Estimates
    } else {
      th_init <- theta_init
    }

    fit <- multifair_base(
      x = x,
      y = y,
      group = group,
      lambda = lam_tmp,
      eta = eta,
      stepsize = stepsize,
      intercept = intercept,
      theta_init = th_init,
      type = type,
      reg = reg,
      crit = crit,
      rho = rho,
      maxit = maxit,
      tol = tol,
      verbose = verbose
    )

    Theta[, , , l] <- fit$Estimates
    Iterations[l] <- fit$Iterations
  }
  return(list(Estimates = Theta,
              Iterations = Iterations,
              Lambdas = lambda_seq))
}


plot_multifair_sp <- function(multifair_sp, type = "l", log = "x", ...) {
  dims <- dim(multifair_sp$Estimates)
  M <- dims[2]
  A <- dims[3]

  if (!is.null(multifair_sp$Etas)) {
    sp_type = "eta"
  } else if (!is.null(multifair_sp$Lambdas)) {
    sp_type = "lambda"
  }

  par(mfrow = c(M, A))

  for (m in 1:M) {
    for (a in 1:A) {
      ylab <- paste0("Data = ", m, ", Group = ", a)
      matplot(x = multifair_sp[[3]],
              y = t(multifair_sp$Estimates[, m, a, ]),
              type = type,
              log = log,
              xlab = sp_type,
              ylab = ylab,
              ...)
    }
  }
  par(mfrow = c(1, 1))
}
