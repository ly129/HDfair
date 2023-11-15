multifair_cv_eta <- function(
    x,
    y,
    group,
    lambda,
    eta_seq,
    nfolds = 5L,
    foldid = NULL,
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
  eta_seq <- sort(eta_seq, decreasing = TRUE)

  ### Initial solution path
  sp <- multifair_sp_eta(
    x,
    y,
    group,
    lambda,
    eta_seq,
    stepsize,
    intercept,
    theta_init,
    type,
    reg,
    crit,
    rho,
    maxit,
    tol,
    warmstart,
    verbose
  )
  nzero <- array(dim = c(M, A, neta))
  for (m in 1:M) {
    for (a in 1:A) {
      nzero[m, a, ] <- rowSums(apply(sp$Estimates[, m, a, ], MARGIN = 1, "!=", 0))
    }
  }

  ### train-test split and pre-allocate result storage
  n <- lapply(x, FUN = nrow)
  # if (is.null(foldid)) {
  foldid <- lapply(n, integer)
  for (m in 1:M) {
    id_by_grp <- split(1:n[[m]], group[[m]])
    foldid_by_grp <- unlist(lapply(id_by_grp, FUN = function(x) sample(rep(seq(nfolds), length.out = length(x)))))
    foldid[[m]][unlist(id_by_grp)] <- foldid_by_grp
  }
  # }

  loss_mat <- matrix(0, nrow = nfolds, ncol = neta)

  ### cv
  for (i in 1:nfolds) {
    which_train <- lapply(foldid, FUN = "!=", i)
    which_val <- lapply(foldid, FUN = "==", i)

    train_x <- mapply(subset, x, which_train, SIMPLIFY = FALSE)
    train_y <- mapply(subset, y, which_train, SIMPLIFY = FALSE)
    train_grp <- mapply(subset, group, which_train, SIMPLIFY = FALSE)
    val_x <- mapply(subset, x, which_val, SIMPLIFY = FALSE)
    val_y <- mapply(subset, y, which_val, SIMPLIFY = FALSE)
    val_grp <- mapply(subset, group, which_val, SIMPLIFY = FALSE)

    sp_fold <- multifair_sp_eta(
      train_x,
      train_y,
      train_grp,
      lambda,
      eta_seq,
      stepsize,
      intercept,
      theta_init,
      type,
      reg,
      crit,
      rho,
      maxit,
      tol,
      warmstart,
      verbose
    )

    for (m in 1:M) {
      for (a in 1:A) {
        index_ma <- val_grp[[m]] == a
        val_xma <- val_x[[m]][index_ma, ]
        val_yma <- val_y[[m]][index_ma]
        for (l in 1:neta) {
          th_mal <- sp_fold$Estimates[, m, a, l]
          loss_mat[i, l] <- loss_mat[i, l] + loss_cts(val_xma, val_yma, th_mal)
        }
      }
    }
  }
  # loss_mat

  ### cv evaluation
  cvm <- apply(loss_mat, MARGIN = 2, FUN = mean)
  cvsd <- apply(loss_mat, MARGIN = 2, FUN = sd)/sqrt(nfolds)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  eta.min <- which.min(cvm)
  # lam.1se <- sum(cvm[seq(lam.min)] > cvup[lam.min]) + 1



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

  return.obj <- list(eta = eta_seq,
                     cvm = cvm,
                     cvsd = cvsd,
                     cvup = cvup,
                     cvlo = cvlo,
                     nzero = nzero,
                     eta.min = eta_seq[eta.min],
                     index = eta.min,
                     # lambda.1se = lambda_seq[lam.1se],
                     Regularization = reg,
                     sp = sp)
  return(return.obj)
}





plot_multifair_cv_eta <- function(multifair_cv_eta_obj) {
  cvup <- multifair_cv_eta_obj$cvup
  cvlo <- multifair_cv_eta_obj$cvlo
  cvm <- multifair_cv_eta_obj$cvm

  eta.seq <- multifair_cv_eta_obj$eta
  eta.min <- multifair_cv_eta_obj$eta.min
  # eta.1se <- multifair_cv_eta_obj$eta.1se

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
  # abline(v = log(lambda.1se), lty = 3)

  axis(side = 3,
       at = log(eta.seq),
       labels = apply(multifair_cv_eta_obj$nzero, MARGIN = 3, FUN = mean),
       tick = FALSE)
}
