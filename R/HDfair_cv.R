HDfair_cv <- function(
    X,
    y,
    ma,
    lambda_length,
    lambda_ratio,
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

  thetas <- array(dim = c(p, A, M, lambda_length, eta_length))
  g <- array(dim = c(A, lambda_length, eta_length))
  iters <- matrix(nrow = lambda_length, ncol = eta_length)

  ### Initial solution path
  sp <- HDfair_sp_lambda(
    X = X,
    y = y,
    ma = ma,
    lambda_length = lambda_length,
    lambda_ratio = lambda_ratio,
    eta = p * 1e6,
    rho = rho,
    weighted = weighted,
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

  loss_mat <- array(0, dim = c(nfold, lambda_length, eta_length))

  for (ll in 2:lambda_length) {
    lam.tmp <- sp$lambdas[ll]

    cv_eta <- HDfair_cv_eta(
      X = X,
      y = y,
      ma = ma,
      lambda = lam.tmp,
      eta_length = eta_length,
      eta_ratio = eta_ratio,
      nfold = nfold,
      foldid = fold_id,
      rho = rho,
      adj=adj,
      eps=eps,
      maxiter=maxiter,
      verbose = TRUE
    )
  }
}
