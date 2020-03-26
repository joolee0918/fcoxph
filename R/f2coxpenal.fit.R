#
# General penalized likelihood
#

#' @importFrom gglasso gglasso
#' @importFrom fda create.bspline.basis
#' @import grpreg
#' @import glmnet glmnet
f22coxpenal.fit <- function(x, y, strata, offset, init, control,
                          weights, method,
                          pcols, pattr, assign, npcols, tuning.method, sm, alpha, gamma, theta, lambda, nlambda = NULL,
                          penalty, sparse.what, argvals, group.multiplier,
                          cv.fit = FALSE)
{
  eps <- control$eps
  n <-  nrow(y)
  if (is.matrix(x)) nvar <- ncol(x)
  else if (length(x)==0) stop("Must have an X variable")
  else nvar <-1

  if (missing(offset) || is.null(offset)) offset <- rep(0,n)
  if (missing(weights)|| is.null(weights))weights<- rep(1,n)
  else {
    if (any(weights<=0)) stop("Invalid weights, must be >0")
  }

  # Get the list of sort indices, but don't sort the data itself
  if (ncol(y) ==3) {
    if (length(strata) ==0) {
      sort.end  <- order(strata, -y[,2]) -1L
      sort.start<- order(strata, -y[,1]) -1L

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
          temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]], 1, init, nystar, penalty)  #init, nystar, penalty)
        else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]], 1, init, nystar, penalty, extralist[[i]])# init, nystar, penalty, extralist[[i]])
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
  npenalty.where <- as.numeric(unlist(npcols))
  n.nonpar <- length(penalty.where)
  n.par <- length(npenalty.where)



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

  ## calculate lambda

  if(is.null(lambda)) {
    nlambda <- ifelse(is.null(nlambda), 20, nlambda)
    p.lambda <- glmnet(xx, y, family="cox",  nlambda=nlambda, standardize=FALSE, thresh=1)$lambda*10
  }else {
    p.lambda <- lambda
  }
  nlambda <- length(p.lambda)

  ## Fit without sparse penalty

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
  #means   <- coxfit$means

  #
  #  Now for the actual fit
  #

  # for (i in 1:m) {
  #    D[[i]] <- sm[[i]]$D/4
  #  }
  #print(D)

  ## Broup band matrix
  d <- 4
  M <- sm[[1]]$bs.dim - d
  H <- as.matrix(Matrix::bandSparse(M+1, M+d, rep(list(rep(1, M+1)), d), k=seq(0, d-1)))


  ## Fitting
  var <- vector('list', L*nlambda)
  df <- loglik <- p.loglik <- rep(0, L*nlambda)
  coef <- matrix(0, ncol = L*nlambda, nrow=nvar)
  for(iter in 1:L){
    for(i in 1:m){
      thetalist[[i]] <- Theta[iter]
      lambdalist[[i]] <- 0

      if(thetalist[[i]]== 0){
        D[[i]] <- NULL
        Dstar <- NULL
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
        Dnrow <- sum(sapply(D, nrow))
        Dstar <- matrix(0, nrow=Dnrow, ncol=nvar)
        nystar <- nvar + Dnrow
      }

    }

#    BB <- rep(1, nvar)
#    BB[penalty.where] <- B


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
    loglik0 <- coxfit$loglik[1]

    n.penalty <- ifelse(penalty == "gBridge", "lasso", penalty)

    fit.beta <- matrix(0, nvar, nlambda)
    if(!is.null(Dstar)) Dstar[,penalty.where] <- as.matrix(bdiag(lapply(1:m, function(i) D[[i]]*sqrt(thetalist[[i]]/(1-thetalist[[i]])) )))
    G <- sm[[1]]$S[[1]]*thetalist[[1]]/(1-thetalist[[1]])



    for(i in 1:nlambda){
      oldbeta <- init
      error <- 1
      iter2 <- 1

      while(iter2 < control$iter.max){
        if (andersen) {coxfit <- .Call(survival:::Cagfit4,
                                       y, xx, cox.newstrat, weights,
                                       offset,
                                       as.double(oldbeta),
                                       sort.start, sort.end,
                                       as.integer(method=="efron"),
                                       as.integer(0),
                                       as.double(control$eps),
                                       as.double(control$toler.chol),
                                       as.integer(1))

        }else{ coxfit <- .Call(survival:::Ccoxfit6,
                               as.integer(0),
                               stime,
                               sstat,
                               xx[sorted,],
                               as.double(offset[sorted]),
                               weights[sorted],
                               cox.newstrat,
                               as.integer(method=="efron"),
                               as.double(control$eps),
                               as.double(control$toler.chol),
                               as.vector(oldbeta),
                               as.integer(1))
        }
        means   <- coxfit$means
        S <- coxfit$u
        I <- solve(matrix(coxfit$imat, nvar, nvar))
        V <- chol(I)
        Y <- solve(t(V))%*%(I%*%oldbeta + S)

        penalty.f <- rep(0, nvar)
        if(penalty == "alasso"){
          penalty.f[penalty.where] <- 1/abs(oldbeta)
        }else if(penalty == "gBridge") {
          tau0 <- (p.lambda[i])^(1/(1-gamma))*gamma^(gamma/(1-gamma))*(1-gamma)
          M <- sm[[1]]$bs.dim - 4
          mu0 <- mu(oldbeta, gamma, tau0, M, 4)
          penalty.f[penalty.where] <- g.pf(mu0, gamma, M, 4)
        }else {
          penalty.f[penalty.where] <- 1
        }



        Ystar <- c(Y, rep(0, Dnrow))
        Vstar <- rbind(V, Dstar) %*% diag(1/penalty.f)

       # newbeta <- wshoot1(nvar, Vstar, Ystar, oldbeta, 1, penalty.f, control$iter.max, control$eps, n)
        #ifelse(penalty=="alasso", "lasso", penalty)
        #p.fit1 <- ncpen1(Ystar, Vstar, family = "gaussian", penalty=ifelse(n.penalty=="alasso", "lasso", n.penalty),  lambda =ifelse(penalty=="gBridge", n/nystar, p.lambda[i]*n/nystar) , x.standardize = FALSE,  intercept=FALSE)

        p.fit1 <- glmnet(Vstar, Ystar, family = "gaussian",  lambda = n/nystar, standardize = FALSE,  intercept=FALSE)
        newbeta <- as.vector(p.fit1$beta)/penalty.f
        lamf <- p.fit1$lambda
        error.r = rep(0, length(newbeta))
        idx0 = which((newbeta- oldbeta) == 0)
        if(length(idx0) == 0){
          error.r <- (newbeta - oldbeta)/oldbeta
        }else {
          error.r[-idx0] <- (newbeta[-idx0] - oldbeta[-idx0])/oldbeta[-idx0]
        }
        error <- max(abs(error.r))
        if(error < control$eps2) break

        oldbeta <- newbeta
        iter2 <- iter2+1
      }


      fit.beta[,i] <- newbeta
      if(cv.fit) {
        coef[, (iter-1)*nlambda + i] <- fit.beta[,i]
      } else{

    #    p.beta.nonzero <- nonzero <- vector('list', m)
    #    n.nonzero <- rep(0, m)
    #    for(j in 1:m){
    #      nonzero[[j]] <- fit.beta[,i][pcols[[j]]] !=0
    #      n.nonzero[j] <- sum(nonzero[[j]])
    #      tmp <- fit.beta[pcols[[j]], i]
    #      p.beta.nonzero[[j]] <- tmp[tmp!=0]
    #    }

        anonzero <- fit.beta[,i] !=0
        an.nonzero <- sum(anonzero)

        if(an.nonzero == 0) {
          var[[(iter-1)*nlambda + i]]  <- NA
          df[(iter-1)*nlambda + i] <- 0
          loglik[(iter-1)*nlambda + i] <- coxfit$loglik[2]
          p.loglik[(iter-1)*nlambda + i] <- NA
          coef[, (iter-1)*nlambda + i] <- fit.beta[,i]
        } else {
      #    newxx <- as.matrix(xx[, anonzero])
      #    beta.nonzero <- fit.beta[anonzero,i]


          #lambda <- p.lambda[i]*n/nystar

          if (andersen) {coxfit <- .Call(survival:::Cagfit4,
                                         y, xx, cox.newstrat, weights,
                                         offset,
                                         as.double(newbeta),
                                         sort.start, sort.end,
                                         as.integer(method=="efron"),
                                         as.integer(0),
                                         as.double(control$eps),
                                         as.double(control$toler.chol),
                                         as.integer(1))

          }else{ coxfit <- .Call(survival:::Ccoxfit6,
                                 as.integer(0),
                                 stime,
                                 sstat,
                                 xx[sorted,],
                                 as.double(offset[sorted]),
                                 weights[sorted],
                                 cox.newstrat,
                                 as.integer(method=="efron"),
                                 as.double(control$eps),
                                 as.double(control$toler.chol),
                                 as.vector(newbeta),
                                 as.integer(1))
          }
          I <- solve(matrix(coxfit$imat, nvar, nvar))

     #     I.nonzero <- matrix(0, nvar, nvar)
    #      I.nonzero[anonzero, anonzero] <- I0 <- solve(matrix(coxfit$imat, an.nonzero, an.nonzero))
    #      pen.tot <- 0

          w <- rep(0, nvar)
          if(penalty == "alasso"){
            w[penalty.where] <- 1/abs(newbeta)
          }else if(penalty == "gBridge") {
            tau0 <- p.lambda[i]^(1/(1-gamma))*gamma^(gamma/(1-gamma))*(1-gamma)
            M <- sm[[1]]$bs.dim - 4
            mu0 <- mu(newbeta, gamma, tau0, M, 4)
            w[penalty.where] <- g.pf2(mu0, gamma, M, 4)
          }else {
            w[penalty.where] <- 1
          }

          A <- numeric(nvar)
          G <- matrix(0, nvar, nvar)
          A[penalty.where & abs(newbeta) > 0] <- w[abs(newbeta)>0]/abs(newbeta[abs(newbeta) > 0])
          A[penalty.where & newbeta == 0]  <- 0.0
          G[penalty.where, penalty.where] <- sm[[1]]$S[[1]]*thetalist[[1]]/(1-thetalist[[1]])
          H <- I + n*diag(A) + G
          H2 <- I +  G
          var[[(iter-1)*nlambda + i]] <- (solve(H)%*%I%*%solve(H))[anonzero, anonzero]

    #        keep.extra[[j]]  <- extralist[[j]][nonzero[[j]], nonzero[[j]]]
    #        temp <- ((pattr[[j]])$pfun)(p.beta.nonzero[[j]], thetalist[[j]], p.lambda[i]*n/nystar, penalty.f[pcols[[j]]][nonzero[[j]]], init[pcols[[j]]][nonzero[[j]]], nystar, penalty, keep.extra[[j]])
    #        I.nonzero[pcols[[j]], pcols[[j]]][nonzero[[j]], nonzero[[j]]] <- I.nonzero[pcols[[j]], pcols[[j]]][nonzero[[j]], nonzero[[j]]] + temp$second
    #        pen.tot <- pen.tot - temp$penalty
    #      }
    #      I.nonzero <- I.nonzero[anonzero, anonzero]
    #      var[[(iter-1)*nlambda + i]] <- solve(I.nonzero)
    #      df[(iter-1)*nlambda + i] <- sum(diag((solve(I.nonzero)%*%I0)))

          df[(iter-1)*nlambda + i]  <- sum( diag((solve(H2[anonzero, anonzero])%*%I[anonzero, anonzero])))
          loglik[(iter-1)*nlambda + i] <- coxfit$loglik[2]
    #      p.loglik[(iter-1)*nlambda + i] <- coxfit$loglik[2] + pen.tot
          coef[, (iter-1)*nlambda + i] <- fit.beta[,i]
        }
      }
    }
    iter <- iter+1
  }

  if(cv.fit) {
    list(beta=coef)
  }else{
    list(beta = coef,
         var = var,
         df=df,
         loglik0 = loglik0,
         loglik = loglik,
         p.loglik = p.loglik,
         lambda = p.lambda,
         theta = Theta,
         pterms = pterms, assign2=assign,
         class = c('fcoxph.penal','fcoxph'))
  }


}




