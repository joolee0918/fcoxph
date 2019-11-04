

TWeiRandom.f <- function(tt,  lam, alp, tau) {
  v <- runif(1, min = 0, max = 1)

  nomin <- - log(1 - v)
  denom <- lam
  term <- (nomin/denom)
  term <- ifelse(is.na(term), tau, term)
  return(term)
}


data.generator <- function(nSimu, N, lam, alp, gamma1, gamma2, rangeval, probC, tau)
{
  nbasis=50+5-2
  data.basis <- create.bspline.basis(rangeval=c(0, 1),norder=5,nbasis=nbasis)

  getdata.f <- function(id,  tau, lam, alp, gamma1, gamma2, W1, W2, Xbeta) {

    lam <- lam * exp(Xbeta) #gamma1 * W1 + gamma2*W2 +

    #cur.t <- TWeiRandom.f(tt = 0,  lam = lam, alp = alp, tau = tau)
    cur.t <- rexp(1, rate=lam)
    if (cur.t >= tau) {
      estart <- 0
      estop <- tau
      estatus <- 0
    } else{
      estart <- 0
      estop <- cur.t
      estatus <- 1
    }
    tmp <- data.frame(id = id, estart = estart, estop = estop, estatus = estatus,
                      tau = tau, W1 = W1, W2 = W2)
    return(tmp)
  }



  G1 <- matrix(0,nbasis,1)
  for(j in 1:nbasis) G1[j] <- inner.prod(beta.func,data.basis,j)

  for(i in 1:nSimu)
  {

    if (probC == 0) {
      CC <- rep(tau, N)
    } else {
      CC <- rexp(N, rate = ((-1) * log(1 - probC)))
      CC <- ifelse(CC > tau, tau, CC)
    }

    cMat1 <- matrix(rnorm(N*nbasis),N,nbasis)

    W1 <- rbinom(N, 1, 0.5)
    W2 <- rnorm(N, 0, 1)
    Xbeta <- cMat1 %*% G1

    data.basis <- create.bspline.basis(rangeval=c(0, 1),norder=5,nbasis=nbasis)
    knots <- c(0,data.basis$params,1)
    X <- cMat1%*%t(eval.basis(knots, data.basis))
    event <- lapply(1:N, function(i) getdata.f(id = i, W1 = W1[i], W2 = W2[i], Xbeta = Xbeta[i],
                                               tau = CC[i], lam = lam, alp = alp, gamma1 = gamma1, gamma2 = gamma2))
    data <- do.call(rbind, event)

  }

  return(list(data, X=X))
}


inner.prod <- function(f,basis,j)
{
  rng <- getbasisrange(basis)
  knots <- c(rng[1],basis$params,1)
  nbasis <- basis$nbasis
  norder <- basis$nbasis - length(knots) + 2

  a <- rng[1]
  if(j-norder > 0) a <- knots[j-norder+1]

  b <- rng[2]
  if (j <= nbasis-norder) b <- knots[j+1]

  bfun <- function(t)
  {
    mat <- eval.basis(t,basis)
    z <- t(mat[,j])
    return(z)
  }

  y <- integrate(function(t) {f(t)*bfun(t)},a,b)
  y$value
}


beta.func <- function(t)
{
  t1 <- t[t<=0.5]
  # t3 <- t[t>=0.7]
  y <- matrix(0,1,length(t))
  y[t<=0.5] <- 2*(1-t1)*sin(2*pi*(t1+0.2))
  # y[t>=0.7] <- 2*t3*sin(2*pi*(t3-0.2))
  y
}


require('fda')
nSimu <- 1000
n <- 500
lam <- 1.5
alp <- 1
tau <- 1
probC <- 0.2
gamma1 <- log(0.8)
gamma2 <- log(1.2)
rangeval <- c(0,1)


#data <- data.generator(nSimu, n, lam, alp, gamma1, gamma2, rangeval, probC, tau)


#data1 <- cbind(data$data, X=data$X)
#m1 <- fcoxph(Surv(estop, estatus)~lf(X, k=30, bs="ps", integration="riemann"),  data=data1, sparse ="none")

#data2 <- cbind(data$data, m1$pcox$smoothdata[[1]])
#m2 <-  gam(estop~s(X.smat,by=X.LX, k=30, bs="cr"),
#                family=cox.ph(),data=data2,weights=estatus)

#m30 <- fcoxph(Surv(estop, estatus)~lf(X, k=30, bs="ps", integration="riemann"),  data=data1, penalty="Lasso", tuning.method = "GCV", sparse ="tail", nlambda=20)
#m31 <- fcoxph(Surv(estop, estatus)~lf(X, k=30, bs="ps", integration="riemann"),  data=data1, penalty="SCAD", tuning.method = "GCV", sparse ="tail", nlambda=20)
#m32 <- fcoxph(Surv(estop, estatus)~lf(X, k=30, bs="ps", integration="riemann"),  data=data1, penalty="MCP", tuning.method = "GCV", sparse ="tail", nlambda=20)





