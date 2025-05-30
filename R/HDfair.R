#' @export
HDfair <- function(
    X,
    y,
    ma,
    lambda,
    eta,
    rho,
    weighted = FALSE,
    th_init=NULL,
    delta_init=NULL,
    adj=1,
    eps=1e-6,
    maxiter=1e4,
    verbose = FALSE)
{
  N = length(y)
  p = ncol(X)
  M = max(ma[,1])
  A = max(ma[,2])
  n <- table(ma[, 1], ma[, 2])

  if (is.null(th_init)) {
    th = th0 = array(0,c(p,A,M))
  } else {
    th = th_init
    mth = apply(th,c(1,3),mean)
    th0 = th - as.vector(mth[,rep(1:M,each=A)])
  }

  if (is.null(delta_init)) {
    delta = rep(0,A)
  } else {
    delta <- delta_init
  }

  # fairness constraint
  g = apply(th0^2,2,sum)/2
  ggr = array(0,c(A,p,A,M))

  # loss function
  h <- 0
  hgr <- array(0, c(p, A, M))
  for (m in 1:M) {
    for (a in 1:A) {
      maids <- ma[, 1] == m & ma[, 2] == a
      Xma <- X[maids, ]
      yma <- y[maids]
      thma <- th[, a, m]

      Xth <- Xma %*% thma
      y_Xth <- yma - Xth

      # sum of mean at each m and a (weighted)
      # vs
      # sum of all obs at each m and a (unweighted): basically upweighs minority

      # also needs to modify the loss function in cv_lambda and cv_eta
      # also modify lambda max in sp_lambda

      if (weighted) {
        h <- h + 0.5 * mean(y_Xth^2)
        hgr[, a, m] <- - t(Xma) %*% y_Xth/n[m,a]
      } else {
        h <- h + 0.5 * sum(y_Xth^2)/N
        hgr[, a, m] <- - t(Xma) %*% y_Xth/N
      }
    }
  }

  iter = 0
  step = 1
  while ( iter < maxiter )
  {
    iter = iter + 1
    pth = th

    # Update u
    u = pmax(pmin(g+delta/rho/adj,eta),0)

    # Update delta
    delta = delta + rho*adj*(g-u)
    delta[u<eta] = 0

    # gradient of g
    for ( a in 1:A )
    {
      ggr[a,,,] = -th0[,rep(a,A),]/A
      ggr[a,,a,] = ggr[a,,a,]*(1-A)
    }

    # objective function and gradient
    L = sum((g-u)^2)*rho*adj^2/2 + sum(delta*(g-u))*adj + h
    #    f = L + sum(apply(th,1,function(x) norm(x,"F")))
    Lgr = apply(ggr*(delta+rho*adj*(g-u))*adj,2:4,sum) + hgr

    # Proximal gradient descent with backtracking line search
    while (TRUE)
    {
      nth = prox(th-step*Lgr,step*lambda)
      mth = apply(nth,c(1,3),mean)
      th0 = nth - as.vector(mth[,rep(1:M,each=A)])
      ng = apply(th0^2,2,sum)/2
      nh <- 0
      nhgr <- array(0, c(p, A, M))
      for (m in 1:M) {
        for (a in 1:A) {
          maids <- ma[, 1] == m & ma[, 2] == a
          Xma <- X[maids, ]
          yma <- y[maids]
          thma <- nth[, a, m]

          Xth <- Xma %*% thma
          y_Xth <- yma - Xth

          if (weighted) {
            nh <- nh + 0.5 * mean(y_Xth^2)
            nhgr[, a, m] <- - t(Xma) %*% y_Xth/n[m,a]
          } else {
            nh <- nh + 0.5 * sum(y_Xth^2)/N
            nhgr[, a, m] <- - t(Xma) %*% y_Xth/N
          }
        }
      }

      nL = sum((ng-u)^2)*rho*adj^2/2 + sum(delta*(ng-u))*adj + nh
      G = (th-nth)/step
      if ( nL <= L - step*sum(Lgr*G) + step*sum(G^2)/2 )
      {
        th = nth
        g = ng
        h = nh
        hgr <- nhgr
        break
      }
      step = step/2
    }
    step = step*1.1

    rho = min(rho*1.001,1000)

    if (verbose) {
      print(sprintf("lambda=%f,eta=%f,max(g)=%f,max(abs(g-u))=%f,max(dth)=%f,step=%f,rho=%f",lambda,eta,max(g),max(abs(g-u)),max(abs(th-pth)),step,rho))

    }
    if ( max(abs(g-u)) < eps & max(abs(th-pth)) < eps )
      break
  }

  list(th=th,g=g,delta=delta,h=h,iter=iter,rho=rho)
}

prox <- function(gr,t)
{
  grF = apply(gr,1,function(x) norm(x,"F"))
  gr*pmax(1-t/grF,0)
}


