#
# General penalized likelihood
#

#' @importFrom gglasso gglasso
#' @importFrom fda create.bspline.basis
#' @import grpreg
#' @import glmnet glmnet
fcoxpenal.fit <- function(x, y, strata, offset, init, control,
                          weights, method,
                          pcols, pattr, assign, npcols, tuning.method, sm, alpha, gamma, theta, lambda, lambda.min.ratio, nlambda = NULL,
                          penalty, sparse.what, argvals, group.multiplier,
                          cv.fit = FALSE)
{
  eps <- control$eps
  n <-  nrow(y)
  if (is.matrix(x)) nvar <- ncol(x)
  else if (length(x)==0) stop("Must have an X variable")
  else nvar <-1

  if (missing(offset) || is.null(offset)) offset <- rep(0.0,n)
  if (missing(weights)|| is.null(weights))weights<- rep(1.0,n)
  else {
    if (any(weights<=0)) stop("Invalid weights, must be >0")
  }

  # Get the list of sort indices, but don't sort the data itself
  if (ncol(y) ==3) {
    if (length(strata) ==0) {
      sort.end  <- order(-y[,2]) -1L
      sort.start<- order(-y[,1]) -1L

      sort <- cbind(order(-y[,2], y[,3]),
                    order(-y[,1])) -1L

      newstrat <- as.integer(n)
    }else {
      sort.end  <- order(strata, -y[,2]) -1L
      sort.start<- order(strata, -y[,1]) -1L

      sort <- cbind(order(strata, -y[,2], y[,3]),
                    order(strata, -y[,1])) -1L

      newstrat  <- as.integer(cumsum(table(strata)))
    }
    status <- y[,3]
    andersen <- TRUE
  }else {
    if (length(strata) ==0) {
      sort <- order(-y[,1], y[,2]) -1L
      newstrat <- as.integer(n)

      sorted <- order(y[,1])
      strata <- NULL
      cox.newstrat <- as.integer(rep(0,n))
    }else {
      sort <- order(strata, -y[,1], y[,2]) -1L
      newstrat <-  as.integer(cumsum(table(strata)))

      sorted <- order(strata, y[,1])
      strata <- strata[sorted]
      cox.newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
    }
    time <- y[,1]
    status <- y[,2]
    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])
    andersen <- FALSE
  }



  if (is.null(nvar) || nvar==0) {
    # A special case: Null model.  Just return obvious stuff
    #  To keep the C code to a small set, we call the usual routines, but
    #  with a dummy X matrix and 0 iterations
    nvar <- 1
    x <- matrix(as.double(1:n), ncol=1)  #keep the .C call happy
    nullmodel <- TRUE
    if (length(init) !=0) stop("Wrong length for inital values")
    init <- 0.0  #dummy value to keep a .C call happy (doesn't like 0 length)
  } else {
    nullmodel <- FALSE

    if (is.null(init)) init <- rep(0., nvar)
    if (length(init) != nvar) stop("Wrong length for inital values")
  }


  #
  # are there any sparse frailty terms?
  #
  npenal <- length(pattr)
  if (npenal == 0 || length(pcols) != npenal)
    stop("Invalid pcols or pattr arg")
  sparse <- sapply(pattr, function(x) !is.null(x$sparse) &&  x$sparse)
  if (sum(sparse) >0) stop("sparse frailty penalty term are not allowed with functionl covariates")

  ## Here we do not allow sparse frailty penalty

  ##Remove 'intercept'
  pterms <- rep(0, length(assign))
  names(pterms) <- names(assign)
  pindex <- rep(0, npenal)
  for (i in 1:npenal) {
    temp <- unlist(lapply(assign, function(x,y) (length(x) == length(y) &&
                                                   all(x==y)), pcols[[i]]))
    if (sparse[i]) pterms[temp] <- 2
    else pterms[temp] <- 1
    pindex[i] <- (seq(along.with=temp))[temp]
  }
  if ((sum(pterms==2) != sum(sparse)) || (sum(pterms>0) != npenal))
    stop("pcols and assign arguments disagree")
  if (any(pindex != sort(pindex))) {
    temp <- order(pindex)
    pindex <- pindex[temp]
    pcols <- pcols[temp]
    pattr <- pattr[temp]
  }

  # ptype= 1 or 3 if a sparse term exists, 2 or 3 if a non-sparse exists
  ptype <- any(sparse) + 2*(any(!sparse))

  pen <- as.integer(switch(penalty, lasso = 1, MCP = 2, gBridge = 3)[1])


  ## Make sure these get defined <TSL>
  f.expr1<-function(coef) NULL
  f.expr2<-function(coef) NULL

  ## normalize xx
  #sd.x <- sqrt(apply(x,2,var)*(n-1))
  #xx <- apply(x, 2, normalize)
  xx <- x
  frailx <- 0
  nfrail <- 0

  # Now the non-sparse penalties
  if (sum(!sparse) >0) {
    full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
    ipenal <- (1:length(pattr))[!sparse]   #index for non-sparse terms
    f.expr2 <- function(coef){
      coxlist2$coef<-coef ##<TSL>
      pentot <- 0
      for (i in ipenal) {
        pen.col <- pcols[[i]]
        coef <- coxlist2$coef[pen.col]

        if (is.null(extralist[[i]]))
          temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]], 1, init, penalty)
        else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]], 1, init, penalty, extralist[[i]])
        if (!is.null(temp$recenter))
          coxlist2$coef[pen.col] <- coxlist2$coef[pen.col]-
          temp$recenter
        if (temp$flag) coxlist2$flag[pen.col] <- TRUE
        else {
          coxlist2$flag[pen.col] <- FALSE
          coxlist2$first[pen.col] <- -temp$first
          if (full.imat) {
            tmat <- matrix(coxlist2$second, nvar, nvar)
            tmat[pen.col,pen.col] <- temp$second
            coxlist2$second <- c(tmat)
          }
          else coxlist2$second[pen.col] <- temp$second
        }
        pentot <- pentot - temp$penalty
      }
      coxlist2$penalty <- as.double(pentot)
      if (any(sapply(coxlist2, length) != length2))
        stop("Length error in coxlist2")
      coxlist2
    }
    if (!is.null(getOption("survdebug")))
      debug(f.expr2)
    if (full.imat) {
      coxlist2 <- list(coef=double(nvar), first=double(nvar),
                       second= double(nvar*nvar), penalty=0.0, flag=rep(FALSE,nvar))
      length2 <- c(nvar, nvar, nvar*nvar, 1, nvar)
    }
    else {
      coxlist2 <- list(coef=double(nvar), first=double(nvar),
                       second=double(nvar), penalty= 0.0, flag=rep(FALSE,nvar))
      length2 <- c(nvar, nvar, nvar, 1, nvar)
    }
    ## in R, f.expr2 is passed as an argument later
    ##.C("init_coxcall2", as.integer(sys.nframe()), expr2)
  }else full.imat <- FALSE

  #
  # Set up initial values for the coefficients
  #  If there are no sparse terms, finit is set to a vector of length 1
  #  rather than length 0, just to stop some "zero length" errors for
  #  later statements where fcoef is saved (but not used)
  #

  finit <- 0
  if (!missing(init) && !is.null(init)) {
    if (length(init) != nvar) {
      stop("Wrong length for inital values")
    }
  }else init <- double(nvar)

  #
  # "Unpack" the passed in paramter list,
  #   and make the initial call to each of the external routines
  #

  m <- length(pcols)

  cfun <- lapply(pattr, function(x) x$cfun)
  parmlist <- lapply(pattr, function(x,eps) c(x$cparm, eps2=eps), sqrt(eps))


  beta.basis <- lapply(1:m, function(i) fda::create.bspline.basis(rangeval=c(argvals[[i]][1], argvals[[i]][length(argvals[[i]])]), nbasis=sm[[i]]$bs.dim))
  B <- as.vector(sapply(1:m, function(i) inprod(beta.basis[[i]])))

  keep.extra <- extralist <- lapply(pattr, function(x) x$pparm)
  iterlist <- thetalist <- lambdalist <-  D <- vector('list', length(cfun))

  if(!is.null(theta)) {
    L <- length(theta)
    Theta <- theta
  }else {
    L <- control$outer.max
    Theta <- parmlist[[1]]$theta
  }


  for (i in 1:m) {
    if(sparse.what == "global") {
      W[[i]] <- compute.W(1, beta.basis[[i]])
    }
    thetalist[[i]] <- Theta[1]
    lambdalist[[i]] <- 0
  }

  penalty.where <- as.numeric(unlist(pcols))
  n.nonpar <- length(penalty.where)


  #
  # Last of the setup: create the vector of variable names
  #
  varnames <- dimnames(xx)[[2]]
  for (i in 1:length(cfun)) {
    if (!is.null(pattr[[i]]$varname))
      varnames[pcols[[i]]] <- pattr[[i]]$varname
  }

  ## need the current environment for callbacks
  rho<-environment()

  ### Score, Hessian
  storage.mode(y) <- storage.mode(weights) <-  "double"
  storage.mode(xx) <- storage.mode(offset) <- "double"
  storage.mode(newstrat) <- "integer"

  ## calculate lambda

  if (nvar==n.nonpar) {
    if (andersen) {coxfit <-  fagfit_init(y, xx, newstrat, weights,
                                          offset, as.vector(rep(0, nvar)),
                                          sort.start, sort.end,
                                          as.integer(method=="efron"),
                                          as.double(control$eps))

    }else{ coxfit <- fcoxfit_init(stime,   sstat, xx[sorted,],
                                  as.double(offset[sorted]), weights[sorted],
                                  as.integer(cox.newstrat), as.double(control$eps),as.integer(method=="efron"),
                                  as.vector(rep(0, nvar)))

  }
    }else{
    if (andersen) {coxfit0 <- .Call(survival:::Cagfit4,
                                   y, xx[, -penalty.where], newstrat, weights,
                                   offset,
                                   as.double(init[-penalty.where]),
                                   sort.start, sort.end,
                                   as.integer(method=="efron"),
                                   as.integer(control$iter.max),
                                   as.double(control$eps),
                                   as.double(control$toler.chol),
                                   as.integer(1))

    }else{ coxfit0 <- .Call(survival:::Ccoxfit6,
                           as.integer(control$iter.max),
                           stime,
                           sstat,
                           xx[sorted, -penalty.where],
                           as.double(offset[sorted]),
                           weights[sorted],
                           cox.newstrat,
                           as.integer(method=="efron"),
                           as.double(control$eps),
                           as.double(control$toler.chol),
                           as.vector(init[-penalty.where]),
                           as.integer(1))
    }

    tmpinit <- rep(0, nvar)
    tmpinit[-penalty.where] <- coxfit0$coef

    if (andersen) {coxfit <-  fagfit_init(y, xx, newstrat, weights,
                                          offset, as.vector(tmpinit),
                                          sort.start, sort.end,
                                          as.integer(method=="efron"),
                                          as.double(control$eps))

    }else{ coxfit <- fcoxfit_init(stime,   sstat, xx[sorted,],
                                  as.double(offset[sorted]), weights[sorted],
                                  as.integer(cox.newstrat), as.double(control$eps),as.integer(method=="efron"),
                                  tmpinit)
    }
    }

  S <-coxfit$u

  if(is.null(lambda)) {
    lambda.max <- max(S[penalty.where])/n
    p.lambda <-  exp(seq(log(lambda.max),log(lambda.min.ratio*lambda.max),len=nlambda))
    if(penalty=="gBridge") p.lambda <- p.lambda*30
  }else {
    p.lambda <- lambda
  }

  nlambda <- length(p.lambda)


  ## Fit without sparse penalty for initial parameter

  if (andersen) {
    coxfit <- .C(survival:::Cagfit5a,
                 as.integer(n),
                 as.integer(nvar),
                 y,
                 xx ,
                 offset,
                 weights,
                 newstrat,
                 sort,
                 means= double(nvar),
                 coef= as.double(init),
                 u = double(nvar),
                 loglik=double(1),
                 as.integer(method=='efron'),
                 as.integer(ptype),
                 as.integer(full.imat),
                 as.integer(nfrail),
                 as.integer(frailx),
                 #R callback additions
                 f.expr1,f.expr2,rho)
  }else{
    coxfit <- .C(survival:::Ccoxfit5a,as.integer(n),
                 as.integer(nvar),
                 y,
                 xx,
                 offset,
                 weights,
                 newstrat,
                 sort,
                 means= double(nvar),
                 coef= as.double(init),
                 u = double(nvar),
                 loglik=double(1),
                 as.integer(method=='efron'),
                 as.integer(ptype),
                 as.integer(full.imat),
                 as.integer(nfrail),
                 as.integer(frailx),
                 f.expr1,f.expr2,rho)
  }
  loglik0 <- coxfit$loglik


  ## Broup band matrix
  d <- 4
  M <- sm[[1]]$bs.dim - d
  H <- as.matrix(Matrix::bandSparse(M+1, M+d, rep(list(rep(1, M+1)), d), k=seq(0, d-1)))


  ## Fitting
  var <- A <- matrix(0, ncol = L*nlambda, nrow=nvar*nvar)
  df <- loglik <-  penalty <- rep(0, L*nlambda)
  coef <-  u <-  matrix(0, ncol = L*nlambda, nrow=nvar)
  for(iter in 1:L){
    for(i in 1:m){
      thetalist[[i]] <- Theta[iter]
      lambdalist[[i]] <- 0
      G <- sm[[1]]$S[[1]]*thetalist[[1]]/(1-thetalist[[1]])

      if(thetalist[[i]]== 0){
        D[[i]] <- NULL
        Dstar <- NULL
        Dncol <- 0
        Dnrow <- 0
        nystar <- nvar
      } else{
        if(is.null(sm[[i]]$D)) {
          eig = eigen(sm[[i]]$S[[1]])
          eig$values[eig$values < 0] = 0
          D[[i]]   = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
        }else {
          D[[i]] <- sm[[i]]$D/4
        }
        Dncol <- sum(sapply(D, ncol))
        Dnrow <- sum(sapply(D, nrow))
        Dstar <- matrix(0, nrow=Dnrow, ncol=nvar)
        nystar <- nvar + Dnrow
      }

    }

    ### initial values estimated without sparse penalty
    if (andersen) { coxfit0 <- .C(survival:::Cagfit5b,
                                  iter=as.integer(control$iter.max),
                                  as.integer(n),
                                  as.integer(nvar),
                                  as.integer(newstrat),
                                  coef = as.double(init),
                                  u    = double(nvar+nfrail),
                                  hmat = double(nvar*(nvar+nfrail)),
                                  hinv = double(nvar*(nvar+nfrail)),
                                  loglik = double(1),
                                  flag = integer(1),
                                  as.double(control$eps),
                                  as.double(control$toler.chol),
                                  as.integer(method=='efron'),
                                  as.integer(nfrail),
                                  fcoef = as.double(finit),
                                  fdiag = double(nfrail+nvar),
                                  f.expr1,f.expr2,rho)
    }else  {coxfit0 <- .C(survival:::Ccoxfit5b,
                          iter=as.integer(control$iter.max),
                          as.integer(n),
                          as.integer(nvar),
                          as.integer(newstrat),
                          coef = as.double(init),
                          u    = double(nvar+nfrail),
                          hmat = double(nvar*(nvar+nfrail)),
                          hinv = double(nvar*(nvar+nfrail)),
                          loglik = double(1),
                          flag = integer(1),
                          as.double(control$eps),
                          as.double(control$toler.chol),
                          as.integer(method=='efron'),
                          as.integer(nfrail),
                          fcoef = as.double(finit),
                          fdiag = double(nfrail+nvar),
                          f.expr1,f.expr2,rho)
    }
    init <- coxfit0$coef

    fit.beta <- matrix(0, nvar, nlambda)
    if(!is.null(Dstar)) Dstar[,penalty.where] <- as.matrix(bdiag(lapply(1:m, function(i) D[[i]]*sqrt(thetalist[[i]]/(1-thetalist[[i]])) )))


    if(andersen){
      fit <- fagfit_cpp(y, xx, newstrat, weights,
                         offset, as.double(init),
                         sort.start, sort.end,
                         as.integer(method=="efron"),
                         as.integer(control$iter.max),
                         as.double(control$eps),
                         H, Dstar, G,  p.lambda, alpha,
                         gamma,  M, d, n.nonpar,  Dnrow, pen, penalty.where, chol, df.f)

    } else{

     fit <- fcoxfit_cpp(stime,   sstat, xx[sorted,],
                       as.double(offset[sorted]), weights[sorted],
                       as.integer(cox.newstrat), as.integer(control$iter.max), as.double(control$eps),
                       H, Dstar, G, as.integer(method=="efron"), init,  p.lambda, alpha,
                       gamma,  M, d, n.nonpar,  Dnrow, pen, penalty.where, chol, df.f)
    }

    df[((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)]  <- fit$df
    loglik[((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- fit$loglik
    var[, ((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- fit$var
    A[, ((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- fit$A
    coef[, ((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- fit$beta
    u[, ((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- fit$u



   for(i in 1:nlambda)  penalty[((iter-1)*nlambda+1): ((iter-1)*nlambda + nlambda)] <- as.numeric(t(fit$beta[penalty.where,i])%*%G%*%fit$beta[penalty.where,i]/2)+ sum(as.numeric(n*p.lambda[i]*sqrt(H%*%abs(fit$beta[penalty.where,i]))))


    }

    if(cv.fit) {
      list(beta=coef)
    }else{
      list(beta = coef,
           df=df,
           var = var,
           A = A,
           u = u,
           loglik0 = loglik0,
           loglik = loglik,
           penalty = penalty,
           lambda = p.lambda,
           theta = Theta,
           H = H,
           pterms = pterms, assign2=assign,
           class = c('fcoxph.penal','fcoxph'))
    }


  }




