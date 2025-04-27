prox <- function(th, lam, reg) {
  if ( reg == "group-lasso" ) {
    type <- ifelse(length(dim(th)) == 3L, "F", "2")
    l2 <- apply(th, 1, norm, type = type)
    a <- pmax(0, 1 - lam/l2)
    a[l2 == 0] <- 0
    return(th * a)
  } else if ( reg == "lasso" ) {
      a <- pmax(0, 1 - lam/abs(th))
      a[th == 0] <- 0
      return(th * a)
  }
}
