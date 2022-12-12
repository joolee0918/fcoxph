#' @importFrom mgcv s smoothCon
#' @importFrom fda eval.penalty

#' @export
fs <- function(X, argvals = NULL, xind = NULL, breaks = NULL, integration = c("dlm", "simpson","trapezoidal", "riemann"),
                presmooth = NULL, presmooth.opts = NULL, sparse = c("none", "local"), tuning.method=c("aic", "bic", "gcv"),
                theta = NULL, lambda = NULL, penalty = c("lasso", "MCP", "gBridge"), m = c(3,2),
          ...)
{
  dots <- list(...)

  if (!is.null(xind)) {
    cat("Argument xind is placed by argvals. xind will not be supported in the next\n        version of refund.")
    argvals = xind
  }
  if (class(X) == "fd") {
    if (is.null(argvals))
      argvals <- argvals <- seq(X$basis$rangeval[1], X$basis$rangeval[2],
                                length = length(X$fdnames[[1]]))
    X <- t(eval.fd(argvals, X))
  }
  else if (is.null(argvals))
    argvals <- seq(0, 1, l = ncol(X))
  xrange <- c(argvals[1], argvals[length(argvals)])
  xind = argvals
  n = nrow(X)
  nt = ncol(X)
  integration <- match.arg(integration)

  tindname <- paste(deparse(substitute(X)), ".smat", sep = "")

  LXname <- paste("B", ".smat", sep = "")
  basistype = "s"

  norder <- m[1] + 1
  mm <- m[2]
  if(is.null(breaks))
    nbasis <- dots$k
  else
    nbasis <- norder + length(breaks) - 2

  M <- nbasis-norder
  beta.basis <- fda::create.bspline.basis(xrange, nbasis=nbasis, norder=norder, breaks=breaks)
  beta.basismat = fda::eval.basis(xind, beta.basis)

  smooth <- list()
  smooth$xind <- xind
  smooth$term <- tindname
  smooth$m <- m

  newcall <- list(as.symbol(basistype))


  if (is.null(dim(xind))) {
    xind <- t(xind)
    stopifnot(ncol(xind) == nt)
    if (nrow(xind) == 1) {
      xind <- matrix(as.vector(xind), nrow = n, ncol = nt,
                     byrow = T)
    }
    stopifnot(nrow(xind) == n)
  }
  if (!is.null(presmooth)) {
    prep.func = refund:::create.prep.func(X, argvals = xind[1, ], method = presmooth, options = presmooth.opts)
    X <- prep.func(newX = X)
  }
    L <- switch(integration, dlm={
      matrix(1, nrow=n, ncol=nt)
    }, simpson = {
      ((xind[, nt] - xind[, 1])/nt)/3 * matrix(c(1, rep(c(4,
                                                          2), length = nt - 2), 1), nrow = n, ncol = nt,
                                               byrow = T)
    }, trapezoidal = {
      diffs <- t(apply(xind, 1, diff))
      0.5 * cbind(t(apply(diffs, 1, filter,
                                      filter = c(1, 1)))[, -(nt - 1)], diffs[, (nt -
                                                                                  1)])
    }, riemann = {
      diffs <- t(apply(xind, 1, diff))
      diffs
     })

  LX <- as.matrix(L * X)

 smooth$X <- LX %*%beta.basismat
 smooth$beta_factor <- as.vector(L[1,]%*%beta.basismat)

  ## Penalty
  dmat <- diag(nbasis)

  if(dots$bs!="ps") smooth$S <- fda::eval.penalty(beta.basis, fda::int2Lfd(mm))/M^3
  else {
    smooth$D<- apply(dmat, 2, diff, 1, mm)
    smooth$S <- t(smooth$D)%*%smooth$D/16
  }

  smooth$beta.basis <- beta.basis


  if(sparse == "none") X <- pterm(smooth, theta,  method = tuning.method, eps = 1e-06, n=n)
  else X <- pterm1(smooth, theta, lambda)

  smooth$X <- NULL
  #if(dots$bs!="ps") smooth$S <- fda::eval.penalty(basis,int2Lfd(m))/M^(3)


  names <- paste0(basistype, "(", tindname,  ", ", "by = ", LXname, ")")
  smooth$label <- names
  smooth$plot.me <- TRUE

  res <- list(names=names, X=X, sm = smooth, argvals = argvals, xind = xind[1,], L = L, tindname=tindname,
              LXname=LXname )
  return(res)
}



pterm <- function (sm, theta, method = c("aic", "bic", "gcv", "fixed"),
                   eps = 1e-06, n)
{

  method <- match.arg(method)

  if (!is.null(theta) ) {
    method <- 'fixed'
    if (theta <=0 || theta >=1) stop("Invalid value for theta")
  }

  W <- sm$X
  #D <- as.matrix(sm$S[[1]])
  D <- as.matrix(sm$S)

  pfun.lFunc <- function(coef, theta, nevent, D) {
    lambda <- ifelse(theta <= 0, 0, theta/(1 - theta))
    list(penalty = as.numeric(t(coef) %*% D %*% coef) * lambda/2,
         first = lambda * D %*% coef, second = lambda * D,
         flag = FALSE)
  }
  printfun <- function(coef, var, var2, df, history, cbase) {
    cmat <- matrix(c(NA, NA, NA, NA, df, NA), nrow = 1)
    nn <- nrow(history$history)
    theta <- ifelse(length(nn), history$history[nn, 1], history$theta)
    list(coef = cmat, history = paste("theta:", format(theta)))
  }
  temp <- switch(method, fixed = list(pfun = pfun.lFunc, cfun = function(parms,
                                                                         ...) {
    list(theta = parms$theta, done = TRUE)
  }, diag = FALSE, pparm = D, cparm = list(theta = theta),
  printfun = printfun),
  aic = list(pfun = pfun.lFunc,cfun = control.tuning, cparm = list(eps = eps, init = c(0.5,0.95), lower = 0, upper = 1, type = "aic", n=n),
             diag = FALSE, pparm = D, cargs = c("neff", "df", "plik"), printfun = printfun),

  bic = list(pfun = pfun.lFunc, cfun =control.tuning, cparm =  list(eps = eps, init = c(0.5,0.95), lower = 0, upper = 1, type = "bic", n=n),
             diag = FALSE, pparm = D,  cargs = c("neff", "df","plik"), printfun = printfun),
  gcv = list(pfun = pfun.lFunc, cfun =control.tuning, cparm =  list(eps = eps, init = c(0.5,0.95), lower = 0, upper = 1, type = "gcv", n=n),
             diag = FALSE, pparm = D,  cargs = c("neff", "df","plik"), printfun = printfun))

  class(W) <- "fcoxph.penalty"
  attributes(W) <- c(attributes(W), temp)
  W
}


pterm1 <- function (sm, theta, lambda)
{


  if(is.null(theta)) theta <- rev(c(0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999))
  #theta <- 0
  W <- sm$X
  #D <- sm$S[[1]]
  D <- as.matrix(sm$S)

  pfun.lFunc <- function(coef, theta, lambda,  penalty.f, init,  penalty, D) {

    kappa <- ifelse(theta <= 0, 0, theta/(1 - theta))

    nvar <- length(coef)

    H <- ifelse(coef==0, 0, abs(coef))
    if(lambda == 0) sparse.penalty <- sparse.first <- sparse.second <- lampen <- 0
    else {
      lampen <- switch(penalty,
                       lasso =  lambda*penalty.f,
                 #      alasso = lambda*penalty.f,
                #       scad = scadderiv(H, alpha, lambda*penalty.f),
                       MCP = mcpderiv(H, alpha, lambda*penalty.f))

      sparse.penalty <- switch(penalty,
                               lasso = sum(as.numeric(lampen)*coef^2/H/2),
            #                   alasso = sum(as.numeric(lampen)*coef^2/H/abs(init))/2,
            #                   scad = sum(as.numeric(lampen)*coef^2/H/2),
                               MCP = sum(as.numeric(lampen)*coef^2/H/2))

      sparse.first <- switch(penalty,
                             lasso = as.numeric(lampen)*coef/H,
            #                 alasso = as.numeric(lampen)*(coef/H/abs(init)),
            #                 scad = as.numeric(lampen)*coef/H,
                             MCP = as.numeric(lampen)*coef/H)

      sparse.second <- switch(penalty,
                              lasso = diag(lampen/H, ncol=nvar, nrow=nvar),
          #                    alasso = diag(lampen/H/abs(init), ncol=nvar, nrow=nvar),
          #                    scad = diag(lampen/H,  ncol=nvar, nrow=nvar),
                              MCP = diag(lampen/H,  ncol=nvar, nrow=nvar))
    }


    list(penalty = as.numeric(t(coef) %*% D %*% coef) * kappa/2 + sparse.penalty ,
         first = kappa * D %*% coef + sparse.first , second = kappa * D + sparse.second,
         flag = FALSE)
  }

  temp <-  list(pfun = pfun.lFunc,
                diag = FALSE, cparm = list(theta = theta),
                pparm = D)

  class(W) <- "fcoxph.penalty"
  attributes(W) <- c(attributes(W), temp)
  W
}




f_override <- function (...)
{
  callstring <- deparse(match.call(), width.cutoff = 500L)
  if (length(callstring) > 1)
    callstring <- paste0(callstring)
  get(callstring, parent.frame())
}



control.fcoxph <- function (eps = 1e-06, eps2 = 1e-04, toler.chol = .Machine$double.eps^0.75,
          iter.max = 100, toler.inf = sqrt(eps), outer.max = 5, timefix = TRUE)
{
  if (iter.max < 0)
    stop("Invalid value for iterations")
  if (eps <= 0)
    stop("Invalid convergence criteria")
  if (eps2 <= 0)
    stop("Invalid convergence criteria")

  if (eps <= toler.chol)
    warning("For numerical accuracy, tolerance should be < eps")
  if (toler.inf <= 0)
    stop("The inf.warn setting must be >0")
  if (!is.logical(timefix))
    stop("timefix must be TRUE or FALSE")
  list(eps = eps, eps2 = eps2, toler.chol = toler.chol, iter.max = as.integer(iter.max),
       toler.inf = toler.inf, outer.max = outer.max,
       timefix = timefix)
}





control.tuning<- function (parms, iter, old, n, df, loglik)
{
  if (iter == 0) {
    if (is.null(parms$init))
      theta <- 0.005
    else theta <- parms$init[1]
    return(list(theta = theta, done = FALSE))
  }
  if (length(parms$type))
    type <- parms$type
  else type <- "aic"
  if (n < df + 2)
    dfc <- (df - n) + (df + 1) * df/2 - 1
  else dfc <- -1 + (df + 1)/(1 - ((df + 2)/n))
  if (iter == 1) {
    history <- c(theta = old$theta, loglik = loglik, df = df,
                 aic = loglik - df, bic = loglik - log(parms$n)*df, gcv = -loglik/(parms$n*(1-df/parms$n)^2))
    if (length(parms$init) < 2)
      theta <- 1
    else theta <- parms$init[2]
    temp <- list(theta = theta, done = FALSE, history = history)
    return(temp)
  }
  history <- rbind(old$history, c(old$theta, loglik, df, loglik -
                                    df,  loglik - log(parms$n)*df,  -loglik/(parms$n*(1-df/parms$n)^2)))
  if (is.null(parms$trace))
    trace <- FALSE
  else trace <- parms$trace
  if (iter == 2) {
    theta <- mean(history[, 1])
    return(list(theta = theta, done = FALSE, history = history,
                tst = 4))
  }
  else if (type == "bic"){
    aic <- history[, 5]
  }
  else if (type == "gcv"){
    aic <- history[, 6]
  } else{
    aic <- history[, 4]
  }
  done <- (abs(1 - aic[iter]/aic[iter - 1]) < parms$eps)
  x <- history[, 1]
  if (x[iter] == max(aic) && x[iter] == max(x))
    newtheta <- 2 * max(x)
  else newtheta <- survival:::frailty.brent(x, aic, lower = parms$lower,
                                            upper = parms$upper)
  if (length(parms$trace) && parms$trace) {
    print(history)
    cat("    new theta=", format(newtheta), "\n\n")
  }
  list(theta = newtheta, done = done, history = history, tst = 4)
}



tuning.summary <- function(tuning.method, loglik, df, n){
c(aic = loglik -  df, bic =loglik - log(n)*df,  gcv = -loglik/(n*(1-df/n)^2))
}
