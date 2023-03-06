#' Main solver of \code{FairReg}
#' @param data A list containing the datasets from \eqn{M} sources.
#' @param outcome Character string or integer indicating the name of the outcome variable or the column index of the outcome variable.
#' @param group Character string or integer indicating the name of the grouping variable or the column index of the grouping variable.
multifair <- function(X,
                      y,
                      group,
                      lam,
                      eta,
                      stepsize,
                      intercept = FALSE,
                      type = "continuous",
                      reg = "group-lasso",
                      crit = "Metric",
                      rho = 1,
                      maxit = 1e3,
                      eps = 1e-6,
                      verbose = FALSE)
{
  # M
  M <- length(X)

  # p
  p <- sapply(X, ncol, simplify = TRUE)
  if (length(unique(p)) == 1) {
    p <- unique(p)
  } else {
    Stop("Different datasets have different number of predictors.")
  }

  # n
  n_M <- sapply(X, nrow, simplify = TRUE)
  # A
  A <- as.integer(max(unique(unlist(group))))

  nlam <- length(lam)

  Th <- array(0, dim = c(p, M, A, nlam))
  th <- array(0, dim = c(p, M, A))

  iters <- integer(length = nlam)

  MA <- M * A
  Xma <- yma <- vector(mode = "list", length = MA)
  nma <- matrix(NA, nrow = M, ncol = A)
  for (m in seq(M)) {
    Xm <- X[[m]]
    ym <- y[[m]]
    for (a in seq(A)) {
      grpma <- group[[m]] == a
      tmp_list_id <- (m - 1) * A + a
      nma[m, a] <- sum(grpma)
      if (intercept) {
        Xma[[tmp_list_id]] <- cbind(1, Xm[grpma, ])
      } else {
        Xma[[tmp_list_id]] <- Xm[grpma, ]
      }
      yma[[tmp_list_id]] <- ym[grpma]
    }
  }

  # initial fairness
  U <- rep(0, A)
  delta <- rep(0, A)
  iters <- integer(length = nlam)

  loss_full <- matrix(nrow = M, ncol = A)
  loss_grad_full <- array(dim = c(p, M, A))

  for (l in seq(nlam)) {
    lam_tmp <- lam[l]

    for (it in seq(maxit)) {
      if (crit == "BGL") {
        for (a in seq(A)) {
          for (m in seq(M)) {
            th_ma <- th[, m, a]
            ls_ma <- loss_cts(
              X_ma = Xma[[(m - 1) * A + a]],
              y_ma = yma[[(m - 1) * A + a]],
              n_ma = nma[m, a],
              th_ma = th_ma
            )
            loss_full[m, a] <- ls_ma$loss
            loss_grad_full[, m, a] <- ls_ma$loss_gr
          }
        }

        fr <- apply(loss_full, MARGIN = 2, FUN = mean)
        fr_grad <- loss_grad_full/M
        ls <- mean(loss_full)
        ls_grad <- loss_grad_full/M/A

        tmp <- delta + rho * fr - U
        grad_update <- ls_grad + fr_grad * rep(tmp, each = p * M)
      } else if (crit == "metric") {
        for (a in seq(A)) {
          for (m in seq(M)) {
            th_ma <- th[, m, a]
            ls_ma <- loss_cts(
              X_ma = Xma[[(m - 1) * A + a]],
              y_ma = yma[[(m - 1) * A + a]],
              n_ma = nma[m, a],
              th_ma = th_ma
            )
            loss_full[m, a] <- ls_ma$loss
            loss_grad_full[, m, a] <- ls_ma$loss_gr
          }
        }
        ls <- mean(loss_full)
        ls_grad <- loss_grad_full/M/A

        ff <- fair_metric(th = th, p = p, M = M, A = A)
        fr <- ff$fair
        tmp <- matrix(ff$fair_gr, ncol = A)%*% (delta + rho * (fr - U))
        grad_update <- ls_grad + array(tmp, c(p, M, A))
      }

      th_new <- prox(th - stepsize * grad_update, lam_tmp * stepsize, reg)

      # step 2b
      U <- pmax(-eta, pmin(eta, fr + delta/rho))

      # step 2c
      delta <- delta + rho * c(fr - U)

      if (it == maxit) message("Maximum iteration reached at lambda = ", lam_tmp, "\n")

      # convergence check
      th_update_norm <- sum((th_new - th)^2)/p/M/A

      if (verbose) {
        cat("Lambda = ", lam_tmp, ", Iteration ", it, ",\n Theta update = ", th_update_norm, "\n", sep = "")
      }

      if (crit == "BGL") {
        if (sum(fr < eta) == A && th_update_norm < eps) break
      } else {
        if ( th_update_norm < eps ) break
      }
      th <- th_new
    }
    Th[, , , l] <- th
    iters[l] <- it
  }

  return(list(Estimates = Th,
              Iterations = iters))
}
