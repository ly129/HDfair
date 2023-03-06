# least squares loss functions for a single dataset-group pair
# loss = 1/2 * || y - x %*% th ||_2^2
loss_cts <- function(X_ma, y_ma, n_ma, th_ma) {
  resid <- y_ma - X_ma %*% th_ma
  loss <- sum(resid^2)/2/n_ma
  loss_gr <- -c(crossprod(X_ma, resid)/n_ma)
  return(list(loss = loss, loss_gr = loss_gr))
}

# BGL fair function
fair_bgl <- function(X_ma, y_ma, n_ma, th_ma, type) {
  if (type == "continuous") {
    bgl <- loss_cts(X_ma, y_ma, n_ma, th_ma)
    return(list(fair = bgl$loss,
                fair_gr = bgl$loss_gr))
  }
}

# metric fair function
fair_metric <- function(th, p, M, A) {
  th_bar <- apply(th, 1:2, mean)
  th_ctrd <- apply(X = th, MARGIN = 3,
                   FUN = "-", th_bar, simplify = FALSE)
  th_ctrd <- array(unlist(th_ctrd), dim = c(p, M, A))
  fair <- apply(th_ctrd^2, 3, sum)/2

  gr = array(0, c(p, M, A, A))
  for ( a in seq(A) ) {
    gr[,,,a] = -th_ctrd[,,a]/A
    gr[,,a,a] = (A-1)*th_ctrd[,,a]/A
  }
  return(list(fair = fair, fair_gr = gr))
}

# check metric
metric_check <- function(th, p, M, A) {
  return(fair_metric(th, p, M, A)$fair)
}

