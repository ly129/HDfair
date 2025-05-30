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
          Xy[, (m-1)*A+a] <- crossprod(Xma, yma)/sum(maids)
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
