#' HDfair_sp_lambda: Compute HDfair Solution Path over Lambda
#'
#' Computes HDfair coefficient estimates across a sequence of \eqn{\lambda} values for a fixed \code{eta}, using warm starts to improve computational efficiency.
#' Returns an object of class "HDfair_sp" containing estimates, lambda sequence, iteration counts, and fairness penalty values for each \eqn{\lambda}.
#'
#' @param X Numeric matrix of dimension \eqn{n \times p}; the design matrix.
#' @param y Numeric vector of length \eqn{n}; the response variable.
#' @param ma Integer matrix of dimension \eqn{n \times 2}; the first column is the source indicator and the second column is the group indicator.
#' @param lambda_length Integer; number of \eqn{\lambda} values in the default sequence (ignored if \code{lambda_seq} is provided).
#' @param lambda_ratio Numeric; ratio between the smallest and largest \eqn{\lambda} in the default sequence (ignored if \code{lambda_seq} is provided).
#' @param lambda_seq Numeric vector; optional user-specified sequence of \eqn{\lambda} values. If provided, overrides \code{lambda_length} and \code{lambda_ratio}.
#' @param eta Numeric scalar; fairness-penalty parameter (constraint threshold).
#' @param rho Numeric scalar; augmented Lagrangian parameter.
#' @param weighted Logical; if \code{TRUE}, computes the loss as the sum of mean losses over sources and groups; if \code{FALSE}, computes loss as the sum over all observations.
#' @param adj Numeric scalar; optional factor applied to the Lagrangian multiplier updates.
#' @param eps Numeric scalar; convergence tolerance based on the L2 norm of the coefficient updates.
#' @param maxiter Integer; maximum number of augmented Lagrangian iterations per \eqn{\lambda}.
#' @param verbose Logical; if \code{TRUE}, prints iteration details (e.g., parameter updates and multiplier values) for each \eqn{\lambda}.
#'
#' @return An object of class "HDfair_sp" with components:
#' \describe{
#'   \item{estimates}{Array of dimension \eqn{p \times A \times M \times L}, where \eqn{L} is the number of \eqn{\lambda} values; estimated coefficients for each source/group and \eqn{\lambda}.}
#'   \item{lambdas}{Numeric vector of \eqn{\lambda} values used.}
#'   \item{iterations}{Integer vector of length \eqn{L}; iterations taken to converge for each \eqn{\lambda}.}
#'   \item{g}{Matrix of dimension \eqn{A \times L}; final fairness penalty values for each group and \eqn{\lambda}.}
#' }
#'
#' @seealso \code{\link{HDfair}}, \code{\link{plot_sp_lambda}}
#' @export
HDfair_sp_lambda <- function(
    X,
    y,
    ma,
    lambda_length = 10,
    lambda_ratio = 1e-2,
    lambda_seq = NULL,
    eta,
    rho,
    weighted = FALSE,
    adj=1,
    eps=1e-6,
    maxiter=1e4,
    verbose=FALSE
) {
  N = nrow(X)
  p = ncol(X)
  M = max(ma[,1])
  A = max(ma[,2])

  if (!is.null(lambda_seq)) lambda_length <- length(lambda_seq)

  thetas <- array(dim = c(p, A, M, lambda_length))
  g <- matrix(nrow = A, ncol = lambda_length)
  iters <- integer(lambda_length)

  if (is.null(lambda_seq)) {
    # determine lambda_max that set all thetas to zero
    Xy <- matrix(NA, p, M*A)
    for (m in 1:M) {
      for (a in 1:A) {
        maids <- ma[, 1] == m & ma[, 2] == a
        Xma <- X[maids, ]
        yma <- y[maids]

        if (weighted) {
          Xy[, (m-1)*A+a] <- crossprod(Xma, yma)/sum(maids)/M/A
        } else {
          Xy[, (m-1)*A+a] <- crossprod(Xma, yma)/N
        }
      }
    }
    row_norms <- sqrt(rowSums(Xy^2))
    lambda_max <- max(row_norms)

    lambdas <- exp(seq(log(lambda_max), log(lambda_max * lambda_ratio),
                       length.out = lambda_length))
  } else {
    lambdas <- lambda_seq
  }

  delta <- rep(0, A)
  theta <- array(0, dim = c(p,A,M))

  for (i in 1:lambda_length) {
    lambda.tmp <- lambdas[i]

    fit = HDfair(X = X,
                 y = y,
                 ma = ma,
                 lambda = lambda.tmp,
                 eta = eta,
                 rho = rho,
                 weighted = weighted,
                 th_init = theta,
                 delta_init = delta,
                 adj = adj,
                 eps = eps,
                 maxiter = maxiter,
                 verbose = verbose)

    theta <- fit$th
    delta <- fit$delta
    thetas[, , , i] <- theta
    g[, i] <- fit$g
    iters[i] <- fit$iter
  }

  return.obj <- list(estimates = thetas,
                     lambdas = lambdas,
                     iterations = iters,
                     g = g)
  class(return.obj) <- "HDfair_sp"

  return(return.obj)
}

#' plot_sp_lambda: Plot HDfair Solution Path over Lambda
#'
#' Plots the HDfair coefficient paths for each source and group against \eqn{\lambda} on a log scale.
#'
#' @param sp_lambda An object of class "HDfair_sp", as returned by \code{HDfair_sp_lambda}.
#' @param lam Numeric scalar; optional \eqn{\lambda} value at which to draw a vertical reference line.
#'
#' @return NULL; this function is called for its side effect of generating plots.
#' @export

plot_sp_lambda <- function(sp_lambda, lam = NULL) {
  thetas <- sp_lambda$estimates
  lambdas <- sp_lambda$lambdas

  M <- dim(thetas)[3]
  A <- dim(thetas)[2]

  for (m in 1:M) {
    for (a in 1:A) {
      matplot(lambdas, t(thetas[, a, m, ]), type = "l", log = "x",
              main = paste0("Site ", m, ", Group ", a))
      if (!is.null(lam)) abline(v = lam, lty = 2)
    }
  }
}
