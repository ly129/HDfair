#' @export
HDfair_cv_lambda <- function(
    X,
    y,
    ma,
    lambda_length,
    lambda_ratio,
    nfold,
    foldid = NULL,
    eta,
    rho,
    adj=1,
    eps=1e-6,
    maxiter=1e4,
    verbose = FALSE
) {
  N = nrow(X)
  p = ncol(X)
  M = max(ma[,1])
  A = max(ma[,2])

  thetas <- array(dim = c(p, A, M, lambda_length))
  g <- matrix(nrow = A, ncol = lambda_length)
  iters <- integer(lambda_length)

  ### Initial solution path
  sp <- HDfair_sp_lambda(
    X = X,
    y = y,
    ma = ma,
    lambda_length = lambda_length,
    lambda_ratio = lambda_ratio,
    eta = eta,
    rho = rho,
    adj=adj,
    eps=eps,
    maxiter=maxiter,
    verbose = verbose
  )
  nzero <- array(dim = c(M, A, lambda_length))
  for (m in 1:M) {
    for (a in 1:A) {
      nzero[m, a, ] <- rowSums(apply(sp$estimates[, a, m, ], MARGIN = 1, "!=", 0))
    }
  }

  lam.seq <- sp$lambdas

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

  loss_mat <- matrix(0, nrow = nfold, ncol = lambda_length)

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

    sp_fold <- HDfair_sp_lambda(
      X = train_x,
      y = train_y,
      ma = train_ma,
      lambda_seq = lam.seq,
      eta = eta,
      rho = rho,
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
        for (l in 1:lambda_length) {
          th_mal <- sp_fold$estimates[, a, m, l]
          loss_mat[i, l] <- loss_mat[i, l] + sum((val_yma - val_xma %*% th_mal)^2)/sum(which_val)
        }
      }
    }
  }
  # loss_mat

  ### cv evaluation
  cvm <- apply(loss_mat, MARGIN = 2, FUN = mean)
  cvsd <- apply(loss_mat, MARGIN = 2, FUN = sd)/sqrt(nfold)
  cvup <- cvm + cvsd
  cvlo <- cvm - cvsd
  lam.min <- which.min(cvm)
  lam.1se <- sum(cvm[seq(lam.min)] > cvup[lam.min]) + 1



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

  return.obj <- list(lambda = lam.seq,
                     cvm = cvm,
                     cvsd = cvsd,
                     cvup = cvup,
                     cvlo = cvlo,
                     nzero = nzero,
                     lambda.min = lam.seq[lam.min],
                     lambda.1se = lam.seq[lam.1se],
                     index = c(lam.min = lam.min, lam.1se = lam.1se),
                     foldid = fold_id,
                     sp = sp)
  class(return.obj) <- "HDfair_cv"
  return(return.obj)
}




#' @export

plot_cv_lambda <- function(cv_lambda) {
  cvup <- cv_lambda$cvup
  cvlo <- cv_lambda$cvlo
  cvm <- cv_lambda$cvm

  lam.seq <- cv_lambda$lambda
  lambda.min <- cv_lambda$lambda.min
  lambda.1se <- cv_lambda$lambda.1se

  ylims <- range(c(cvup, cvlo))
  plot(x = log(lam.seq), y = cvm, ylim = ylims, col = 0,
       ylab = "CV Error",
       xlab = expression(paste("Log(", lambda, ")")))
  for (i in 1:length(cvm)) {
    segments(x0 = log(lam.seq)[i], y0 = cvup[i],
             x1 = log(lam.seq)[i], y1 = cvlo[i],
             col = 'grey60')
  }
  points(x = log(lam.seq), y = cvup, pch = 95, cex = 1, col = 'grey60')
  points(x = log(lam.seq), y = cvlo, pch = 95, cex = 1, col = 'grey60')
  points(x = log(lam.seq), y = cvm, pch = 20, cex = 1, col = 'red')
  abline(v = log(lambda.min), lty = 3)
  abline(v = log(lambda.1se), lty = 3)

  axis(side = 3,
       at = log(lam.seq),
       labels = apply(cv_lambda$nzero, MARGIN = 3, FUN = mean),
       tick = FALSE)
}
