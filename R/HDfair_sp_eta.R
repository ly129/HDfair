#' @export
HDfair_sp_eta <- function(
    X,
    y,
    ma,
    lambda,
    eta_length = 10,
    eta_ratio = 1e-2,
    eta_seq = NULL,
    rho,
    th_init = NULL,
    delta_init = NULL,
    adj=1,
    eps=1e-6,
    maxiter=1e4,
    verbose = FALSE
) {
  N = nrow(X)
  p = ncol(X)
  M = max(ma[,1])
  A = max(ma[,2])

  if (!is.null(eta_seq)) eta_length <- length(eta_seq)

  thetas <- array(dim = c(p, A, M, eta_length))
  g <- matrix(nrow = A, ncol = eta_length)
  iters <- integer(eta_length)

  if (is.null(eta_seq)) {
    # determine eta max that the constraints begin to have an effect
    eta_init <- 1e3
    while (TRUE) {
      fit = HDfair(X = X,
                   y = y,
                   ma = ma,
                   lambda = lambda,
                   eta = eta_init,
                   rho = rho,
                   th_init = th_init,
                   delta_init = NULL,
                   adj = adj,
                   eps = eps,
                   maxiter = maxiter)
      eta_max <- max(fit$g)
      if (eta_max < eta_init) {
        break
      } else {
        eta_init <- eta_init * 10
      }
    }
    etas <- exp(seq(log(eta_max), log(eta_max * eta_ratio), length.out = eta_length))
  } else {
    etas <- eta_seq
  }

  theta <- array(0, dim = c(p, A, M))
  delta <- numeric(A)

  for (i in 1:eta_length) {
    eta.tmp <- etas[i]

    fit = HDfair(X = X,
                 y = y,
                 ma = ma,
                 lambda = lambda,
                 eta = eta.tmp,
                 rho = rho,
                 th_init = theta,
                 delta_init = delta,
                 adj = adj,
                 eps = eps,
                 maxiter = maxiter,
                 verbose = verbose)

    iters[i] <- fit$iter
    theta <- fit$th
    delta <- fit$delta
    thetas[, , , i] <- theta
    g[, i] <- fit$g
  }

  return.obj <- list(estimates = thetas,
                     etas = etas,
                     iterations = iters,
                     g = g)
  class(return.obj) <- "HDfair_sp"
  return(return.obj)
}



#' @export

plot_sp_eta <- function(sp_eta, eta = NULL, vars = 1:5) {
  thetas <- sp_eta$estimates
  etas <- sp_eta$etas

  M <- dim(thetas)[3]
  A <- dim(thetas)[2]

  for (m in 1:M) {
    for (v in vars) {
      matplot(etas, t(thetas[v, , m, ]), type = "l", log = "x")
      if (!is.null(eta)) abline(v = eta, lty = 2)
    }
  }
}
