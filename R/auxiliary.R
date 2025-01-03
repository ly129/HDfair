# least squares loss functions for a single dataset-group pair
# loss = 1/2 * || y - x %*% th ||_2^2
loss_cts <- function(X_ma, y_ma, th_ma) {
  resid <- c(y_ma - X_ma %*% th_ma)
  return(sum(resid * resid)/2)
}

loss_cts_grad <- function(X_ma, y_ma, th_ma) {
  resid <- c(y_ma - X_ma %*% th_ma)
  return(-c(crossprod(X_ma, resid)))
}

# BGL fair function
fair_bgl <- function(X_ma, y_ma, th_ma, type) {
  if (type == "continuous") {
    bgl <- loss_cts(X_ma, y_ma, th_ma)
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

#' @export
# check metric
metric_check <- function(th) {
  dims <- dim(th)
  p <- dims[1]
  M <- dims[2]
  A <- dims[3]
  return(fair_metric(th, p, M, A)$fair)
}


# fair metric with th of list object
fair_metric_list <- function(th_ma, p, M, A) {
  th <- array(unlist(th_ma), dim = c(p, M, A))
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
