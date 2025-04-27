#' @export
multifair_sp_eta <- function(
    x,
    y,
    group,
    lambda,
    eta_seq,
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

  neta <- length(eta_seq)
  eta_seq <- sort(eta_seq, decreasing = FALSE)

  Theta <- array(NA, dim = c(p, M, A, neta))
  Iterations <- integer(neta)

  if (is.null(theta_init)) theta_init <- array(0, dim = c(p, M, A))

  fit <- multifair_base(
    x = x,
    y = y,
    group = group,
    lambda = lambda,
    eta = eta_seq[1],
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

  for (l in 2:neta) {
    eta_tmp <- eta_seq[l]

    if (warmstart) {
      th_init <- fit$Estimates
    } else {
      th_init <- theta_init
    }

    fit <- multifair_base(
      x = x,
      y = y,
      group = group,
      lambda = lambda,
      eta = eta_tmp,
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
              Etas = eta_seq))
}





plot_multifair_sp_eta <- function(multifair_sp, n = 1:6, type = "l", log = "x", ...) {
  dims <- dim(multifair_sp$Estimates)
  M <- dims[2]
  A <- dims[3]

  par(mfrow = c(1, length(n)))

  for (i in n) {
      ylab <- paste0("Variable No. ", i)
      yy <- multifair_sp$Estimates[i, , , ]
      if (!is.matrix(yy)) {
        yy <- apply(yy, MARGIN = 3, FUN = as.vector)
      }
      matplot(x = multifair_sp[[3]],
              y = t(yy),
              type = type,
              log = log,
              xlab = expression(paste(eta)),
              ylab = ylab,
              ...)
  }
  par(mfrow = c(1, 1))
}
