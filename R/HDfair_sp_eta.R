#' HDfair_sp_eta: Compute HDfair Solution Path over Eta
#'
#' Computes HDfair coefficient estimates across a sequence of \eqn{\eta} values for a fixed \code{lambda}, using warm starts for efficiency. Returns an object of class "HDfair_sp" containing estimates, eta sequence, iteration counts, and fairness penalty values for each \eqn{\eta}.
#'
#' @param X Numeric matrix of dimension \eqn{n \times p}; the design matrix.
#' @param y Numeric vector of length \eqn{n}; the response variable.
#' @param ma Integer matrix of dimension \eqn{n \times 2}; the first column is the source indicator (1, 2, ...), and the second column is the group indicator (1, 2, ...).
#' @param lambda Numeric scalar; tuning parameter for the group lasso penalty (fixed for the path).
#' @param eta_length Integer; number of eta values to evaluate if \code{eta_seq} is not provided.
#' @param eta_ratio Numeric scalar; ratio between minimum and maximum eta when generating the sequence (on the log scale).
#' @param eta_seq Numeric vector; optional user-specified sequence of \eqn{\eta} values. Overrides \code{eta_length} and \code{eta_ratio}.
#' @param rho Numeric scalar; augmented Lagrangian parameter for the optimization.
#' @param weighted Logical; if \code{TRUE}, computes the loss as the sum of mean losses over all sources and groups.
#' @param adj Numeric scalar; optional factor applied to the Lagrangian multiplier updates.
#' @param eps Numeric scalar; convergence tolerance based on the L2 norm of the coefficient updates.
#' @param maxiter Integer; maximum number of augmented Lagrangian iterations.
#' @param verbose Logical; if \code{TRUE}, prints iteration details.
#'
#' @return An object of class "HDfair_sp" containing:
#' \describe{
#'   \item{estimates}{Array of dimension \eqn{p \times A \times M \times L}, where \eqn{L} is the number of \eqn{\eta} values, with estimated coefficients for each group, source, and eta.}
#'   \item{etas}{Numeric vector of length \code{eta_length} with the sequence of eta values evaluated.}
#'   \item{iterations}{Integer vector of length \code{eta_length} with the number of iterations for each fit.}
#'   \item{g}{Matrix of dimension \eqn{A \times L} with fairness penalty values for each group and eta.}
#' }
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
    weighted = FALSE,
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
    eta_init <- p * 1e6
    while (TRUE) {
      fit <- HDfair(X = X,
                    y = y,
                    ma = ma,
                    lambda = lambda,
                    eta = eta_init,
                    rho = rho,
                    weighted = weighted,
                    th_init = NULL,
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

  # if (weighted) {
  #   wts <- rep(1/t(table(ma[, 1], ma[, 2])), table(ma[, 1], ma[, 2]))
  # } else {
  #   wts <- rep(1, N)
  # }
  # pooled.lasso <- glmnet::glmnet(
  #   x = X,
  #   y = y,
  #   family = "gaussian",
  #   lambda = lambda*sqrt(M * A),
  #   weights = wts,
  #   intercept = FALSE
  # )

  # theta <- array(as.vector(pooled.lasso$beta), dim = c(p, A, M))
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
                 weighted = weighted,
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


#' plot_sp_eta: Plot HDfair Coefficient Paths over Eta
#'
#' Plots the estimated coefficient trajectories for selected variables across \eqn{\eta} values on a log-scale axis.
#'
#' @param sp_eta An object of class "HDfair_sp" returned by \code{HDfair_sp_eta}.
#' @param eta Numeric scalar; optional value of eta to highlight with a vertical dashed line.
#' @param vars Integer vector; indices of predictor variables to plot (default: first five).
#' @export

plot_sp_eta <- function(sp_eta, eta = NULL, vars = 1:5) {
  thetas <- sp_eta$estimates
  etas <- sp_eta$etas

  M <- dim(thetas)[3]
  A <- dim(thetas)[2]


  if (M > 1) {
    for (v in vars) {
      yminmax <- c(min(thetas[v, , , ]), max(thetas[v, , , ]))

      matplot(etas, t(thetas[v, , 1, ]), type = "l", log = "x", lty = 1, ylim = yminmax)
      if (!is.null(eta)) abline(v = eta, lty = 2)

      for (m in 2:M) {
        matplot(etas, t(thetas[v, , m, ]), type = "l", log = "x", lty = m, add = TRUE)
      }
    }

  } else {
    for (v in vars) {
      matplot(etas, t(thetas[v, , 1, ]), type = "l", log = "x")
      if (!is.null(eta)) abline(v = eta, lty = 2)
    }
  }


}
