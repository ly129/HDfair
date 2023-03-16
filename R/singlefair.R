#' Main solver of \code{FairReg}
#' @param X Predictors.
singlefair <- function(X,
                       y,
                       group,
                       lam,
                       eta,
                       stepsize,
                       intercept = FALSE,
                       type = "continuous",
                       esti_func = NULL,
                       reg = "lasso",
                       crit = "metric",
                       rho = 1,
                       maxit = 1e3,
                       eps = 1e-6,
                       verbose = FALSE)
{
  # Checks
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

  if (type == "custom") {
    if ( is.null(custom_esti_func) || !is.function(custom_esti_func )) {
      stop("An R function needs to be provided to `custom_esti_func`")
    }
  }

  # Pre-allocations
  nlam <- length(lam)

  Th <- array(0, dim = c(p, M, nlam))
  th <- array(0, dim = c(p, M))

  iters <- integer(length = nlam)

  MA <- M * A
  Xma <- yma <- vector(mode = "list", length = MA)
  nma <- matrix(NA, nrow = M, ncol = A)
  for (m in seq(M)) {
    for (a in seq(A)) {
      grpma <- group[[m]] == a
      tmp_list_id <- (m - 1) * A + a
      nma[m, a] <- sum(grpma)
      if (intercept) {
        Xma[[tmp_list_id]] <- cbind(1, X[[m]][grpma, ])
      } else {
        Xma[[tmp_list_id]] <- X[[m]][grpma, ]
      }
      yma[[tmp_list_id]] <- y[[m]][grpma]
    }
  }

  wts <- nma/n_M * A

  # initial fairness
  U <- rep(0, A)
  delta <- rep(0, A)
  iters <- integer(length = nlam)

  # loss_full <-
  fair_full <- matrix(nrow = M, ncol = A)
  loss_grad_full <- fair_grad_full <- array(dim = c(p, M, A))

  for (l in seq(nlam)) {
    lam_tmp <- lam[l]

    for (it in seq(maxit)) {
      if (crit == "BGL") {
        for (a in seq(A)) {
          for (m in seq(M)) {
            th_m <- th[, m]
            tmp_list_id <- (m - 1) * A + a
            ls_ma <- loss_cts(
              X_ma = Xma[[tmp_list_id]],
              y_ma = yma[[tmp_list_id]],
              # n_ma = nma[m, a],
              th_ma = th_m
            )
            # loss_full[m, a] <- ls_ma$loss * wts[m, a]
            loss_grad_full[, m, a] <- ls_ma$loss_gr * wts[m, a]
            fair_full[m, a] <- ls_ma$loss
            fair_grad_full[, m, a] <- ls_ma$loss_gr
          }
        }

        fr <- apply(fair_full, MARGIN = 2, FUN = mean)
        fr_grad <- fair_grad_full/M

        # ls <- mean(loss_full)
        ls_grad <- loss_grad_full/M/A

        tmp <- delta + rho * fr - U
        grad_update <- apply(ls_grad, MARGIN = 1:2, FUN = mean) +
          apply(fr_grad * rep(tmp, each = p * M), MARGIN = 1:2, FUN = mean)

        # for (a in seq(A)) {
        #   tmp <- delta[a] + rho * fr[a] - U[a]
        #   fr_grad_a <- fr_grad[, , a] * tmp
        #   ls_grad_a <- ls_grad[, , a]
        #   th[, , a] <- th[, , a] - stepsize * (fr_grad_a + ls_grad_a)
        # }
      }

      th_new <- prox(th - stepsize * grad_update, lam_tmp * stepsize, reg)

      # step 2b
      U <- pmax(-eta, pmin(eta, fr + delta/rho))

      # step 2c
      delta <- delta + rho * c(fr - U)

      if (it == maxit) message("Maximum iteration reached at lambda = ", lam_tmp, "\n")

      # convergence check
      th_update_norm <- sum((th_new - th)^2)/p/M

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
    Th[, , l] <- th
    iters[l] <- it
  }

  return(list(Estimates = Th,
              Iterations = iters))
}
