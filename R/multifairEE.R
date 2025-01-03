multifairee <- function(
    X,
    y,
    group,
    eta_fr,
    eta_ee,
    stepsize = 1,
    intercept = FALSE,
    custom_esti_func = NULL,
    custom_esti_func_gr = NULL,
    reg = "group-lasso",
    crit = "metric",
    rho_fr = 1,
    rho_ee = 1,
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

  # Pre-allocations
  neta_ee <- length(eta_ee)

  Th <- array(0, dim = c(p, M, A, neta_ee))
  th <- array(0, dim = c(p, M, A))

  iters <- integer(length = neta_ee)

  MA <- M * A

  # initial fairness
  U <- delta <- rep(0, A)
  # initial EE
  # A list of length MA, each item is a vector 0_p
  V <- xi <- vector(mode = "list", length = MA)

  Xma <- yma <- thma <- vector(mode = "list", length = MA)
  for (m in seq(M)) {
    for (a in seq(A)) {
      grpma <- group[[m]] == a
      tmp_list_id <- (a - 1) * M + m
      if (intercept) {
        Xma[[tmp_list_id]] <- cbind(1, X[[m]][grpma, ])
      } else {
        Xma[[tmp_list_id]] <- X[[m]][grpma, ]
      }
      yma[[tmp_list_id]] <- y[[m]][grpma]
      V[[tmp_list_id]] <- xi[[tmp_list_id]] <- rep(0, p)
    }
  }

  iters <- integer(length = neta_ee)

  for (l in seq(neta_ee)) {
    eta_ee_tmp <- eta_ee[l]

    for (it in seq(maxit)) {
      thma <- lapply(seq(MA) - 1, FUN = function(x) th[, x %% M + 1, x %/% M + 1])
      ufunc <- mapply(custom_esti_func, Xma, yma, thma, SIMPLIFY = FALSE)
      ugrad <- mapply(custom_esti_func_gr, Xma, yma, thma, SIMPLIFY = FALSE)
      ufunc_V <- mapply("-", ufunc, V, SIMPLIFY = FALSE)
      rho_ufunc_V <- lapply(ufunc_V, "*", rho_ee)
      xi_rho_ufunc_V <- mapply("+", rho_ufunc_V, xi, SIMPLIFY = FALSE)
      tmp_ee <- mapply("%*%", ugrad, xi_rho_ufunc_V)

      ff <- fair_metric(th = th, p = p, M = M, A = A)
      fr <- ff$fair
      tmp_ff <- matrix(ff$fair_gr, ncol = A) %*% (delta + rho_fr * (fr - U))

      grad_update <- array(tmp_ee, c(p, M, A)) + array(tmp_ff, c(p, M, A))

      th_new <- prox(th - stepsize * grad_update, stepsize, reg)

      # step 2b
      U <- pmax(-eta_fr, pmin(eta_fr, fr + delta/rho_fr))
      V <- mapply(function(x, y){
        pmax(-eta_ee_tmp, pmin(eta_ee_tmp, x + y/rho_ee))
      }, ufunc, xi, SIMPLIFY = FALSE)

      # step 2c
      delta <- delta + rho_fr * c(fr - U)
      xi <- mapply(function(x, y, z) {
        x + rho_ee * (y - z)
      }, xi, ufunc, V, SIMPLIFY = FALSE)

      if (it == maxit) message("Maximum iteration reached at eta_ee = ", eta_ee_tmp, "\n")

      # convergence check
      th_update_norm <- sum((th_new - th)^2)/p/M/A

      if (verbose) {
        cat("eta_ee = ", eta_ee_tmp, ", Iteration ", it, ", Theta update = ", th_update_norm, "\n", sep = "")
      }

      if (crit == "BGL") {
        if (sum(fr < eta) == A && th_update_norm < eps) break
      } else if (crit == "metric") {
        cat("etafr", sum(fr < eta_fr), "maxee", max(abs(unlist(ufunc))), "\n")
        if (sum(fr < eta_fr) == A && max(abs(unlist(ufunc))) < eta_ee_tmp && th_update_norm < eps) break
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
