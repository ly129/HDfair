#' Main solver of \code{HDfair}
#' @param data A list containing the datasets from \eqn{M} sources.
#' @param outcome Character string or integer indicating the name of the outcome variable or the column index of the outcome variable.
#' @param group Character string or integer indicating the name of the grouping variable or the column index of the grouping variable.
multifair_base <- function(
    x,
    y,
    group,
    lambda,
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

  # Pre-allocations
  if (is.null(theta_init)) {
    th <- array(0, dim = c(p, M, A))
  } else {
    th <- theta_init
  }

  xma <- yma <- vector(mode = "list", length = M * A)
  for (m in seq(M)) {
    for (a in seq(A)) {
      grpma <- group[[m]] == a
      tmp_list_id <- (m - 1) * A + a
      if (intercept) {
        xma[[tmp_list_id]] <- cbind(1, x[[m]][grpma, ])
      } else {
        xma[[tmp_list_id]] <- x[[m]][grpma, ]
      }
      yma[[tmp_list_id]] <- y[[m]][grpma]
    }
  }

  # initial fairness
  U <- delta <- rep(0, A)
  # loss_full <- matrix(nrow = M, ncol = A)
  loss_grad_full <- array(dim = c(p, M, A))

  for (it in seq(maxit)) {
    for (m in seq(M)) {
      for (a in seq(A)) {
        tmp_list_id <- (m - 1) * A + a
        loss_grad_full[, m, a] <- loss_cts_grad(xma[[tmp_list_id]],
                                                yma[[tmp_list_id]],
                                                th[, m, a])/sum(n_M)
        # thma[[tmp_list_id]] <- th[, m, a]
      }
    }

    # loss_tmp <- mapply(loss_cts_loss, xma, yma, thma, SIMPLIFY = TRUE)
    # # loss_full <- loss_tmp[1, ]
    # for (m in seq(M)) {
    #   for (a in seq(A)) {
    #     loss_grad_full[, m, a] <- loss_tmp[, (m - 1) * A + a]
    #   }
    # }

    if (crit == "metric") {
      ls_grad <- loss_grad_full/M/A

      ff <- fair_metric(th = th, p = p, M = M, A = A)
      fr <- ff$fair
      tmp <- matrix(ff$fair_gr, ncol = A) %*% (delta + rho * (fr - U))
      grad_update <- ls_grad + array(tmp, c(p, M, A))
    }

    th_new <- prox(th - stepsize * grad_update, lambda * stepsize, reg)

    # step 2b
    U <- pmax(-eta, pmin(eta, fr + delta/rho))

    # step 2c
    delta <- delta + rho * c(fr - U)

    # if (it == maxit) message("Maximum iteration reached at lambda = ", lambda, ", eta = ", eta, "\n")

    # convergence check
    th_update_norm <- norm((th_new - th), "2")/p/M/A

    if (verbose) {
      cat("Lambda = ", lambda, ", Eta = ", eta, ", Iteration ", it, ", Theta update = ", th_update_norm, "\n", sep = "")
    }

    if (crit == "metric") {
      if (sum(fr < eta) == A && th_update_norm < tol) break
    } else {
      if ( th_update_norm < tol ) break
    }
    th <- th_new
  }

  return(list(Estimates = th,
              Iterations = it))
}
