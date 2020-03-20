
#' @importFrom fda getbasisrange inprod
#' @importFrom Matrix bandSparse
#'
compute.W <- function(j, basis)
{
  L <- basis$nbasis
  rng <- fda::getbasisrange(basis)
  breaks <- c(rng[1],basis$params,rng[2])
  M <- length(breaks) - 1
  norder <- L-M+1
  #W1 <- fda::inprod(basis,basis,rng=c(breaks[j],rng[2]))
  #W2 <- fda::inprod(basis,basis,rng=c(rng[1], breaks[j]))
  w <- fda::inprod(basis,basis,rng=c(rng[1], rng[2]))
  W1 <- W2 <- matrix(0, L, L)
  W2[1:(j-1), 1:(j-1)]  <- w[1:(j-1), 1:(j-1)]
  W1[j:L, j:L] <- w[j:L, j:L]
  W <- list(W1, W2)
  #fda::inprod(basis,basis,rng=c(rng[1], breaks[j-1]))
  #W[j:ncol(W), j:ncol(W)]
}

positivepart <- function(fx) {
  return(ifelse(fx>=0, fx, 0))
}


scadderiv <- function(ftheta, fa, flambda) {
  return(flambda*(1-(1-apply(as.matrix(fa*flambda-abs(ftheta)), 1, positivepart)/((fa-1)*flambda))*as.numeric(abs(ftheta)>flambda)))
}

mcpderiv <- function(ftheta, fa, flambda) {
  return( ifelse(abs(ftheta) < fa*flambda, (flambda - abs(ftheta)/fa)*sign(ftheta), 0))
}

mu = function(b, gamma, tau, M, d)
{
  mu <- rep(0, M)
  for(j in 1:(M+1)){
    mu[j] <- ((1/gamma-1)/tau)^gamma*(sum(abs(b)[j:(j+d-1)]))^gamma
  }
  return(mu)
}

g.pf = function(mu, gamma, k, M, d)
{

  A <- as.matrix(Matrix::bandSparse(M+1, M+d, rep(list(rep(1, M+1)), d), k=seq(0, d-1)))
  as.vector(t(ifelse(mu==0, 1e10, mu^(1-1/gamma)))%*%A)
}

g.pf2 = function(mu, gamma, k, M, d)
{

  A <- as.matrix(Matrix::bandSparse(M+1, M+d, rep(list(rep(1, M+1)), d), k=seq(0, d-1)))
  as.vector(t(ifelse(mu==0, 0, mu^(1-1/gamma)))%*%A)
}


list2df <- function (l)
{
  nrows <- sapply(l, function(x) nrow(as.matrix(x)))
  print(nrows)
  stopifnot(length(unique(nrows)) == 1)
  ret <- data.frame(rep(NA, nrows[1]))
  for (i in 1:length(l)) ret[[i]] <- l[[i]]
  names(ret) <- names(l)
  return(ret)
}


na.omit_pcox <- function (object, ...)
{
  n <- length(object)
  omit <- logical(nrow(object))
  vars <- seq_len(n)
  for (j in vars) {
    x <- object[[j]]
    if (!is.atomic(x))
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L)
      omit <- omit | x
    else omit <- omit | apply(x, 1, all)
  }
  xx <- object[!omit, , drop = FALSE]
  if (any(omit > 0L)) {
    temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
    attr(temp, "class") <- "omit"
    attr(xx, "na.action") <- temp
  }
  xx
}


