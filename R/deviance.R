cox.deviance <- function(x, y, strata, offset, init, control,
                      weights, method)
{
  n <-  nrow(y)
  if (is.matrix(x)) nvar <- ncol(x)
  else {
    if (length(x)==0) nvar <-0
    else nvar <-1
  }
  time <- y[,1]
  status <- y[,2]

  # Sort the data (or rather, get a list of sorted indices)
  if (length(strata)==0) {
    sorted <- order(time)
    strata <- NULL
    newstrat <- as.integer(rep(0,n))
  }
  else {
    sorted <- order(strata, time)
    strata <- strata[sorted]
    newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
  }
  if (missing(offset) || is.null(offset)) offset <- rep(0,n)
  if (missing(weights)|| is.null(weights))weights<- rep(1,n)
  else {
    if (any(weights<=0)) stop("Invalid weights, must be >0")
    weights <- weights[sorted]
  }
  stime <- as.double(time[sorted])
  sstat <- as.integer(status[sorted])

  if (nvar==0) {
    # A special case: Null model.
    #  (This is why I need the rownames arg- can't use x' names)
    # Set things up for 0 iterations on a dummy variable
    x <- as.matrix(rep(1.0, n))
    nullmodel <- TRUE
    nvar <- 1
    init <- 0
    maxiter <- 0
  }
  else {
    nullmodel <- FALSE
  }

  storage.mode(weights) <-  "double"
  loglik <- rep(0, ncol(init))
  for(i in 1:ncol(init)){
    beta <- init[,i]
    storage.mode(beta) <- "double"
    loglik[i] <- .Call(survival:::Ccoxfit6,
                 as.integer(0),
                 stime,
                 sstat,
                 x[sorted,],
                 as.double(offset[sorted]),
                 weights,
                 newstrat,
                 as.integer(method=="efron"),
                 as.double(control$eps),
                 as.double(control$toler.chol),
                 as.vector(beta),
                 as.integer(1))$loglik[1]  # internally rescale

  # loglik<- loglik_coxfit(
  #                stime,
  #                sstat,
  #                x[sorted,],
  #                as.double(offset[sorted]),
  #                weights,
  #                newstrat,
  #                as.integer(method=="efron"),
  #                as.vector(init),
  #                as.integer(1))  # internally rescale

  }
  return(-2*loglik)
}
