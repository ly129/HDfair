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



