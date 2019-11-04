
#' @importFrom fda getbasisrange inprod
#'
compute.W <- function(j, basis)
{
  L <- basis$nbasis
  rng <- fda::getbasisrange(basis)
  breaks <- c(rng[1],basis$params,rng[2])
  M <- length(breaks) - 1
  norder <- L-M+1
  W <- fda::inprod(basis,basis,rng=c(breaks[j],rng[2]))

  W[j:ncol(W), j:ncol(W)]
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


list2df <- function (l)
{
  nrows <- sapply(l, function(x) nrow(as.matrix(x)))
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


