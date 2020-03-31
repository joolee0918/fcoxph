
#' @importFrom fda eval.basis
TWeiRandom.f <- function(tt,  lam, alp, tau) {
  v <- runif(1, min = 0, max = 1)

  nomin <- - log(1 - v)
  denom <- lam
  term <- (nomin/denom)
  term <- ifelse(is.na(term), tau, term)
  return(term)
}

#' @export
data.generator <- function(N, lam, alp, gamma1, gamma2, rangeval, probC, tau)
{
  nbasis=50+5-2
  data.basis <- fda::create.bspline.basis(rangeval=c(0, 1),norder=5,nbasis=nbasis)

  getdata.f <- function(id,  tau, lam, alp, gamma1, gamma2, W1, W2, Xbeta) {

    lam <- lam * exp(Xbeta + gamma1 * W1)

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

    data.basis <- fda::create.bspline.basis(rangeval=c(0, 1),norder=5,nbasis=nbasis)
    knots <- c(0,data.basis$params,1)
    X <- cMat1%*%t(fda::eval.basis(knots, data.basis))

    event <- lapply(1:N, function(i) getdata.f(id = i, W1 = W1[i], W2 = W2[i], Xbeta = Xbeta[i],
                                               tau = CC[i], lam = lam, alp = alp, gamma1 = gamma1, gamma2 = gamma2))
    data <- do.call(rbind, event)

    data1 <- list(data=data, X = X)

  return(data1)
}


inner.prod <- function(f,basis,j)
{
  rng <- fda::getbasisrange(basis)
  knots <- c(rng[1],basis$params,1)
  nbasis <- basis$nbasis
  norder <- basis$nbasis - length(knots) + 2

  a <- rng[1]
  if(j-norder > 0) a <- knots[j-norder+1]

  b <- rng[2]
  if (j <= nbasis-norder) b <- knots[j+1]

  bfun <- function(t)
  {
    mat <- fda::eval.basis(t,basis)
    z <- t(mat[,j])
    return(z)
  }

  y <- integrate(function(t) {f(t)*bfun(t)},a,b)
  y$value
}


