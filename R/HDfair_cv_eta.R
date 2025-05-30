#' @export
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
