# GEE

# cd4 <- read.delim("~/Desktop/Research/FairReg/cd4.txt")
tv90 <- read.csv("~/Desktop/Research/FairReg/TV1990.csv")

n.unique <- length(unique(tv90$subject))

V <- matrix(nrow = 4, ncol = 4)
for (i in 1:4) {
  for (j in 1:4) {
    V[i, j] <- 0.6^abs(i-j)
  }
}

X <- tv90[, c("trt", "base", "age", "period")]
X$trt <- ifelse(X$trt == "placebo", 0, 1)
Y <- tv90$y

beta <- rep(1, 4)
V.inv <- solve(V)

gee.linear <- function(x, y, beta) {
  x_i <- split(x, f = tv90$subject)
  y_i <- split(y, f = tv90$subject)
  gee_i <- mapply(
    FUN = function(xx, yy) {
      partial_mu <- as.matrix(xx)
      y_mu <- (yy - as.matrix(xx) %*% beta)
      return(t(partial_mu) %*% V.inv %*% y_mu)
    }, x_i, y_i
  )
  return(rowMeans(gee_i))
}

gee.linear.grad <- function(x, y, beta) {
  x_i <- split(x, f = tv90$subject)
  x_i <- lapply(x_i, as.matrix)
  gee_grad_i <- lapply(x_i, crossprod)
  return(Reduce("+", gee_grad_i)/length(x_i))
}

gee.linear(x = X, y = Y, beta = rep(1, 4))
gee.linear.grad(x = X, y = Y, beta = rep(1, 4))
