score <-
  function(object, x, y, strat = NULL, collapse=FALSE, weighted=FALSE, ...)
  {

    n <- length(object$residuals)
    rr <- object$residuals
    vv <- drop(object$naive.var)
    if (is.null(vv)) vv <- drop(object$var)
    weights <- object$weights
    if (is.null(weights)) weights <- rep(1,n)
     method <- object$method

     ny <- ncol(y)
     status <- y[,ny,drop=TRUE]

     nstrat <- as.numeric(strat)
     nvar <- ncol(x)
      if (is.null(strat)) {
          ord <- order(y[,ny-1], -status)
          newstrat <- rep(0,n)
        }
        else {
          ord <- order(nstrat, y[,ny-1], -status)
          newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
        }
        newstrat[n] <- 1

        # sort the data
        x <- x[ord,]
        y <- y[ord,]
        score <- exp(object$linear.predictors)[ord]



   if (ny==2) {
        resid <- .C(survival:::Ccoxscore, as.integer(n),
                    as.integer(nvar),
                    as.double(y),
                    x=as.double(x),
                    as.integer(newstrat),
                    as.double(score),
                    as.double(weights[ord]),
                    as.integer(method=='efron'),
                    resid= double(n*nvar),
                    double(2*nvar))$resid
      }
      else {
        resid<- .C(survival:::Cagscore,
                   as.integer(n),
                   as.integer(nvar),
                   as.double(y),
                   as.double(x),
                   as.integer(newstrat),
                   as.double(score),
                   as.double(weights[ord]),
                   as.integer(method=='efron'),
                   resid=double(n*nvar),
                   double(nvar*6))$resid
      }
      if (nvar >1) {
        rr <- matrix(0, n, nvar)
        rr[ord,] <- matrix(resid, ncol=nvar)
        dimnames(rr) <- list(names(object$residuals),
                             names(object$coefficients))
      }
      else rr[ord] <- resid


    #
    # Multiply up by case weights (which will be 1 for unweighted)
    #
    if (weighted) rr <- rr * weights

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
      rr <- naresid(object$na.action, rr)
      if (is.matrix(rr)) n <- nrow(rr)
      else               n <- length(rr)
    }


    # Collapse if desired
    if (!missing(collapse)) {
      if (length(collapse) !=n) stop("Wrong length for 'collapse'")
      rr <- drop(rowsum(rr, collapse))
    }

   return(rr)
   }

