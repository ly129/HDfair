#' HDfair_cv_eta: Cross-Validation for Optimal Eta
#'
#' Performs K-fold cross-validation over a sequence of \eqn{\eta} values to select an optimal
#' fairness-penalty parameter for HDfair, using warm-started solution paths from \code{HDfair_sp_eta}.
#'
#' @param X Numeric matrix of dimension \eqn{n \times p}; the design matrix.
#' @param y Numeric vector of length \eqn{n}; the response variable.
#' @param ma Integer matrix of dimension \eqn{n \times 2}; first column is source indicator (1,...,M), second column is group indicator (1,...,A).
#' @param lambda Numeric scalar; fixed group-lasso penalty parameter.
#' @param eta_length Integer; number of \eqn{\eta} values to consider if \code{eta_seq} is NULL.
#' @param eta_ratio Numeric scalar; ratio of minimum to maximum \eqn{\eta} when constructing a default sequence.
#' @param nfold Integer; number of cross-validation folds.
#' @param foldid Optional integer vector of length \eqn{n} specifying fold assignments; default generates stratified folds by group-source.
#' @param rho Numeric scalar; augmented Lagrangian parameter for optimization.
#' @param weighted Logical; if \code{TRUE}, computes validation loss as mean over sources and groups; otherwise uses pooled sum over all observations.
#' @param adj Numeric scalar; adjustment factor applied to multiplier updates.
#' @param eps Numeric scalar; convergence tolerance for coefficient updates.
#' @param maxiter Integer; maximum number of augmented Lagrangian iterations.
#' @param verbose Logical; if \code{TRUE}, prints iteration details.
#'
#' @return An object of class \code{"HDfair_cv"} containing:
#' \describe{
#'   \item{eta}{Vector of \eqn{\eta} values considered.}
#'   \item{cvm}{Mean cross-validation error for each \eqn{\eta}.}
#'   \item{cvsd}{Standard error of CV error across folds.}
#'   \item{cvup}{Upper bound = \code{cvm + cvsd}.}
#'   \item{cvlo}{Lower bound = \code{cvm - cvsd}.}
#'   \item{nzero}{Array of dimension \eqn{M \times A \times L}, where \eqn{L} is the length of \eqn{\eta}, giving number of nonzeros in each fit.}
#'   \item{eta.min}{Selected \eqn{\eta} with minimum CV error.}
#'   \item{eta.1se}{Largest \eqn{\eta} within one standard error of minimum CV error.}
#'   \item{index}{Named indices \code{eta.min} and \code{eta.1se}.}
#'   \item{foldid}{Fold assignments used.}
#'   \item{sp}{Solution path object returned by \code{HDfair_sp_eta}.}
#' }
#'
#' @seealso \code{\link{HDfair_sp_eta}}, \code{\link{plot_cv_eta}}
#' #' @export
HDfair_cv_eta <- function(
    X,
    y,
    ma,
    lambda,
    eta_length,
    eta_ratio,
    nfold,
    foldid = NULL,
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

  thetas <- array(dim = c(p, A, M, eta_length))
  g <- matrix(nrow = A, ncol = eta_length)
  iters <- integer(eta_length)

  ### Initial solution path
  sp <- HDfair_sp_eta(
    X = X,
    y = y,
    ma = ma,
    lambda = lambda,
    eta_length = eta_length,
    eta_ratio = eta_ratio,
    rho = rho,
    weighted = weighted,
    adj=adj,
    eps=eps,
    maxiter=maxiter,
    verbose = verbose
  )
  nzero <- array(dim = c(M, A, eta_length))
  for (m in 1:M) {
    for (a in 1:A) {
      nzero[m, a, ] <- rowSums(apply(sp$estimates[, a, m, ], MARGIN = 1, "!=", 0))
    }
  }

  eta.seq <- sp$etas

  ### train-test split and pre-allocate result storage
  if (is.null(foldid)) {
    grp <- with(as.data.frame(ma), paste(V1, V2, sep = "_"))
    fold_id <- integer(nrow(ma))

    for (g in unique(grp)) {
      idx <- which(grp == g)
      n   <- length(idx)
      # assign 1:K in a recycled way, then shuffle
      fold_id[idx] <- sample(rep(seq_len(nfold), length.out = n))
    }
  } else {
    fold_id <- foldid
  }

  loss_mat <- matrix(0, nrow = nfold, ncol = eta_length)

  ### cv
  for (i in 1:nfold) {
    which_train <- fold_id != i
    which_val <- fold_id == i

    train_x <- X[which_train, ]
    train_y <- y[which_train]
    train_ma <- ma[which_train, ]
    val_x <- X[which_val, ]
    val_y <- y[which_val]
    val_ma <- ma[which_val, ]

    sp_fold <- HDfair_sp_eta(
      X = train_x,
      y = train_y,
      ma = train_ma,
      lambda = lambda,
      eta_seq = eta.seq,
      rho = rho,
      weighted = weighted,
      adj=adj,
      eps=eps,
      maxiter=maxiter,
      verbose = verbose
    )

    for (m in 1:M) {
      for (a in 1:A) {
        index_ma <- val_ma[, 1] == m & val_ma[, 2] == a
        val_xma <- val_x[index_ma, ]
        val_yma <- val_y[index_ma]
        for (l in 1:eta_length) {
          th_mal <- sp_fold$estimates[, a, m, l]
          loss_mat[i, l] <- loss_mat[i, l] + sum((val_yma - val_xma %*% th_mal)^2)/sum(which_val)
        }
      }
    }
  }
  if (weighted) loss_mat <- loss_mat/M/A

  ### cv evaluation
  cvm <- apply(loss_mat, MARGIN = 2, FUN = mean)
  cvsd <- apply(loss_mat, MARGIN = 2, FUN = sd)/sqrt(nfold)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  eta.min <- which.min(cvm)
  eta.1se <- sum(cvm[seq(eta.min)] > cvup[eta.min]) + 1



  # ylims <- range(c(cvup, cvlo))
  # plot(x = log(lam.seq), y = cvm, ylim = ylims, col = 0,
  #      ylab = "CV Error",
  #      xlab = expression(paste("Log(", lambda, ")")))
  # for (i in 1:nlam) {
  #   segments(x0 = log(lam.seq)[i], y0 = cvup[i],
  #            x1 = log(lam.seq)[i], y1 = cvlo[i],
  #            col = 'grey60')
  # }
  # points(x = log(lam.seq), y = cvup, pch = 95, cex = 1, col = 'grey60')
  # points(x = log(lam.seq), y = cvlo, pch = 95, cex = 1, col = 'grey60')
  # points(x = log(lam.seq), y = cvm, pch = 20, cex = 1, col = 'red')
  # abline(v = log(lam.seq[lam.min]), lty = 3)
  # abline(v = log(lam.seq[lam.1se]), lty = 3)
  # axis(side = 3, at = log(lam.seq), labels = nzero, tick = FALSE)

  return.obj <- list(eta = eta.seq,
                     cvm = cvm,
                     cvsd = cvsd,
                     cvup = cvup,
                     cvlo = cvlo,
                     nzero = nzero,
                     eta.min = eta.seq[eta.min],
                     eta.1se = eta.seq[eta.1se],
                     index = c(eta.min = eta.min, eta.1se = eta.1se),
                     foldid = fold_id,
                     sp = sp)
  class(return.obj) <- "HDfair_cv"
  return(return.obj)
}



#' plot_cv_eta: Plot Cross-Validation Results over Eta
#'
#' Plots mean cross-validation error and error bars against log(eta), marking
#' the selected \code{eta.min} and \code{eta.1se} values.
#'
#' @param cv_eta An object of class \code{"HDfair_cv"} returned by \code{HDfair_cv_eta}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @export

plot_cv_eta <- function(cv_eta) {
  cvup <- cv_eta$cvup
  cvlo <- cv_eta$cvlo
  cvm <- cv_eta$cvm

  eta.seq <- cv_eta$eta
  eta.min <- cv_eta$eta.min
  eta.1se <- cv_eta$eta.1se

  ylims <- range(c(cvup, cvlo))
  plot(x = log(eta.seq), y = cvm, ylim = ylims, col = 0,
       ylab = "CV Error",
       xlab = expression(paste("Log(", eta, ")")))
  for (i in 1:length(cvm)) {
    segments(x0 = log(eta.seq)[i], y0 = cvup[i],
             x1 = log(eta.seq)[i], y1 = cvlo[i],
             col = 'grey60')
  }
  points(x = log(eta.seq), y = cvup, pch = 95, cex = 1, col = 'grey60')
  points(x = log(eta.seq), y = cvlo, pch = 95, cex = 1, col = 'grey60')
  points(x = log(eta.seq), y = cvm, pch = 20, cex = 1, col = 'red')
  abline(v = log(eta.min), lty = 3)
  abline(v = log(eta.1se), lty = 3)

  axis(side = 3,
       at = log(eta.seq),
       labels = apply(cv_eta$nzero, MARGIN = 3, FUN = mean),
       tick = FALSE)
}
