#
# General penalized likelihood
#

#' @importFrom gglasso gglasso
#' @importFrom fda create.bspline.basis
#' @import grpreg
ffcoxpenal.fit <- function(x, y, strata, offset, init, control,
                          weights, method, rownames,
                          pcols, pattr, assign, npcols, tuning.method, sm, l, alpha, theta, lambda, lambda.min, nlambda, penalty, sparse.what, argvals, group.multiplier) {
  eps <- control$eps
  n <-  nrow(y)
  if (is.matrix(x)) nvar <- ncol(x)
  else  if (length(x)==0) stop("Must have an X variable")
  else nvar <-1

  if (missing(offset) || is.null(offset)) offset <- rep(0,n)
  if (missing(weights)|| is.null(weights))weights<- rep(1,n)
  else {
    if (any(weights<=0)) stop("Invalid weights, must be >0")
  }

  # Get the list of sort indices, but don't sort the data itself
  if (ncol(y) ==3) {
    if (length(strata) ==0) {
      sorted <- cbind(order(-y[,2], y[,3]),
                      order(-y[,1])) -1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- cbind(order(strata, -y[,2], y[,3]),
                      order(strata, -y[,1])) -1L
      newstrat  <- as.integer(cumsum(table(strata)))
    }
    status <- y[,3]
    andersen <- TRUE
  }
  else {
    if (length(strata) ==0) {
      sorted <- order(-y[,1], y[,2]) -1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- order(strata, -y[,1], y[,2]) -1L
      newstrat <-  as.integer(cumsum(table(strata)))
    }
    status <- y[,2]
    andersen <- FALSE
  }

  n.eff <- sum(y[,ncol(y)])  #effective n for a Cox model is #events
  n.coef <- ncol(x)
  #
  # are there any sparse frailty terms?
  #
  npenal <- length(pattr)
  if (npenal == 0 || length(pcols) != npenal)
    stop("Invalid pcols or pattr arg")
  sparse <- sapply(pattr, function(x) !is.null(x$sparse) &&  x$sparse)
  if (sum(sparse) >1) stop("Only one sparse penalty term allowed")

  #
  # Create a marking vector for the terms, the same length as assign
  #    with pterms == 0=ordinary term, 1=penalized, 2=sparse,
  #    pindex = length of pcols = position in pterms
  #
  # Make sure that pcols is a strict subset of assign, so that the
  #   df computation (and printing) can unambiguously decide which cols of
  #   X are penalized and which are not when doing "terms" like actions.
  # To make some downstream things easier, order pcols and pattr to be
  #   in the same relative order as the terms in 'assign'
  #
  ## can't compute assign attribute in R without terms
  ## if (missing(assign)) assign <- attr(x, 'assign')[-1]
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


  if (any(sparse)) {
    sparse.attr <- (pattr[sparse])[[1]]  #can't use [[sparse]] directly
    # if 'sparse' is a T/F vector
    fcol <- unlist(pcols[sparse])
    if (length(fcol) > 1) stop("Sparse term must be single column")

    # Remove the sparse term from the X matrix
    xx <- x[, -fcol, drop=FALSE]
    for (i in 1:length(assign)){
      j <- assign[[i]]
      if (j[1] > fcol) assign[[i]] <- j-1
    }
    for (i in 1:npenal) {
      j <- pcols[[i]]
      if (j[1] > fcol) pcols[[i]] <- j-1
    }

    frailx <- x[, fcol]
    frailx <- match(frailx, sort(unique(frailx)))
    nfrail <- max(frailx)
    nvar <- nvar - 1

    #Set up the callback for the sparse frailty term
    pfun1 <- sparse.attr$pfun
    ### In R we use a function and eval() it, not an expression
    f.expr1 <- function(coef){
      coxlist1$coef <- coef
      if (is.null(extra1)) temp <- pfun1(coef, theta1, n.eff)
      else  temp <- pfun1(coef, theta1, n.eff, extra1)

      if (!is.null(temp$recenter))
        coxlist1$coef <- coxlist1$coef - as.double(temp$recenter)
      if (!temp$flag) {
        coxlist1$first <- -as.double(temp$first)
        coxlist1$second <- as.double(temp$second)
      }
      coxlist1$penalty <- -as.double(temp$penalty)
      coxlist1$flag   <- as.logical(temp$flag)
      if (any(sapply(coxlist1, length) != c(rep(nfrail,3), 1, 1)))
        stop("Incorrect length in coxlist1")
      coxlist1
    }
    if (!is.null(getOption("survdebug"))) debug(f.expr1)

    coxlist1 <- list(coef=double(nfrail), first=double(nfrail),
                     second=double(nfrail), penalty=0.0, flag=FALSE)
    ## we pass f.expr1 in as an argument in R
    ##.C("init_coxcall1", as.integer(sys.nframe()), expr1)
  }
  else {
    xx <- x
    frailx <- 0
    nfrail <- 0
  }



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
          temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]], W[[i]])
        else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambdalist[[i]],
                                            W[[i]], extralist[[i]][[1]], n, extralist[[i]][[2]])
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
  }
  else full.imat <- FALSE

  #
  # Set up initial values for the coefficients
  #  If there are no sparse terms, finit is set to a vector of length 1
  #  rather than length 0, just to stop some "zero length" errors for
  #  later statements where fcoef is saved (but not used)
  #
  if (nfrail >0) finit <- rep(0,nfrail)
  else finit <- 0
  if (!missing(init) && !is.null(init)) {
    if (length(init) != nvar) {
      if (length(init) == (nvar+nfrail)) {
        finit <- init[-(1:nvar)]
        init  <- init[1:nvar]
      }
      else stop("Wrong length for inital values")
    }
  }
  else init <- double(nvar)

  #
  # "Unpack" the passed in paramter list,
  #   and make the initial call to each of the external routines
  #

  m <- length(pcols)

  cfun <- lapply(pattr, function(x) x$cfun)
  parmlist <- lapply(pattr, function(x,eps) c(x$cparm, eps2=eps), sqrt(eps))


  beta.basis <- lapply(1:m, function(i) fda::create.bspline.basis(rangeval=c(argvals[[i]][1], argvals[[i]][length(argvals[[i]])]), nbasis=sm[[i]]$bs.dim))
  extralist<- lapply(pattr, function(x) x$pparm)
  iterlist <- thetalist <- lambdalist <- W <- cutoff <- sparse.where <- W0 <-  vector('list', length(cfun))
  printfun  <- lapply(pattr, function(x) x$printfun)


  for (i in 1:m) {

    if(!is.null(l)) cutoff[[i]] <- l[i]
    else cutoff[[i]] <-  length(beta.basis[[i]]$params) + 1


    if(sparse.what == "global") {
      W[[i]] <- compute.W(1, beta.basis[[i]])
    } else if (sparse.what == "local"){
      W[[i]] <- compute.W(cutoff[[i]], beta.basis[[i]])
      W[[i]] <- as.matrix(Matrix::bdiag(matrix(0, cutoff[[i]]-1, cutoff[[i]]-1), W[[i]]))
    } else if (sparse.what == "tail"){
      W[[i]] <- compute.W(cutoff[[i]], beta.basis[[i]])
      W0[[i]] <- as.matrix(bdiag(diag(1, cutoff[[i]]-1), solve(chol(W[[i]]))))
      W[[i]] <- as.matrix(Matrix::bdiag(matrix(0, cutoff[[i]]-1, cutoff[[i]]-1), W[[i]]))
    }
  }



  n.nonp <- 0
  if(m>0) for(i in 1:m){
    n.nonp <- n.nonp + length(pcols[[i]])
    sparse.where[[i]] <-  seq(pcols[[i]][1] + cutoff[[i]] -1, pcols[[i]][n.nonp])
  }


  penalty.where <- as.numeric(unlist(pcols))
  npenalty.where <- as.numeric(unlist(npcols))
  n.par <- length(npenalty.where)

  group <- seq(1:n.coef)
  group[npenalty.where] <- 0
  for(i in 1:m) group[pcols[[i]]] <- i



  ## Set up lambda using grpsurv library function

  tempW0 <- list()
  for (i in 1:m) {
    tempW0[[i]] <- compute.W(1, beta.basis[[i]])
  }
  Wb <- chol(as.matrix(Matrix::bdiag(tempW0)))

  ind <- order(y[,1])
  XG <- grpreg:::newXG(xx[ind,], group, group.multiplier, 1, FALSE)
  K <- table(group)
  K1 <- as.integer(if (min(group)==0) cumsum(K) else c(0, cumsum(K)))

  nullFit <- survival:::coxph.fit(XG$X[, group==0, drop=FALSE], y[ind,], strata[ind], offset[ind], init=rep(0, length(group[group==0])),
                                  control, weights[ind], method, rownames=row.names(XG$X))

  s <- weights*residuals(nullFit, type="martingale" )
  lambda.max <- .Call(grpreg:::maxgrad, XG$X%*%solve(Wb), s, K1, as.double(XG$m)) / n

  if (lambda.min==0) lambda <- c(exp(seq(log(.001*lambda.max),log(lambda.max), len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.min*lambda.max),log(lambda.max), len=nlambda))

  for (i in 1:m) {
    temp <- (cfun[[i]])(parmlist[[i]], iter=0)
    if (sparse[i]) {
      theta1 <- temp$theta
      extra1 <- extralist[[i]]
    }
    thetalist[[i]] <- temp$theta
    lambdalist[[i]] <- lambda[1]*XG$m[i]
    iterlist[[i]] <- temp
  }

  #
  # Manufacture the list of calls to cfun, with appropriate arguments
  #
  ## Amazingly, all this works in R, so I don't need to understand it.
  ##
  temp1 <- c('x', 'coef', 'plik', 'loglik', 'status', 'neff', 'df', 'trH')
  temp2 <- c('frailx', 'coxfit$fcoef', 'loglik1',  'coxfit$loglik', 'status',
             'n.eff')
  temp3 <- c('xx[,pen.col]', 'coxfit$coef[pen.col]','loglik1',
             'coxfit$loglik', 'status', 'n.eff')
  calls <- vector('expression', length(cfun))
  cargs <- lapply(pattr, function(x) x$cargs)
  for (i in 1:length(cfun)) {
    tempchar <- paste("(cfun[[", i, "]])(parmlist[[", i, "]], iter,",
                      "iterlist[[", i, "]]")
    temp2b <- c(temp2, paste('pdf[', i, ']'), paste('trH[', i, ']'))
    temp3b <- c(temp3, paste('pdf[', i, ']'), paste('trH[', i, ']'))
    if (length(cargs[[i]])==0)
      calls[i] <- parse(text=paste(tempchar, ")"))
    else {
      temp <- match(cargs[[i]], temp1)
      if (any(is.na(temp))) stop(paste((cargs[[i]])[is.na(temp)],
                                       "not matched"))
      if (sparse[i]) temp4 <- paste(temp2b[temp], collapse=',')
      else           temp4 <- paste(temp3b[temp], collapse=',')

      calls[i] <- parse(text=paste(paste(tempchar,temp4,sep=','),')'))
    }
  }
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

    #
    # Have C store the data, and get the loglik for beta=initial, frailty=0
    #
    storage.mode(y) <- storage.mode(weights) <-  "double"
    storage.mode(xx) <- storage.mode(offset) <- "double"
    if (andersen) coxfit <- .C(survival:::Cagfit5a,
                               as.integer(n),
                               as.integer(nvar),
                               y,
                               xx ,
                               offset,
                               weights,
                               newstrat,
                               sorted,
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
    else       coxfit <- .C(survival:::Ccoxfit5a,as.integer(n),
                            as.integer(nvar),
                            y,
                            xx,
                            offset,
                            weights,
                            newstrat,
                            sorted,
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

    loglik0 <- coxfit$loglik
    means   <- coxfit$means

    #
    #  Now for the actual fit
    #
    # for (outer1 in 1:control$outer.max1) {
    #    betazero <- FALSE
    #    for(outer2 in 1:control$outer.max2) {
    #iter <- (outer1-1)*control$outer.max2 + outer2

    for (i in 1:m) {
      lambdalist[[i]] <- 0
    }

    iter <- 0
    for(i in 1:control$outer.max){
      if (andersen)  coxfit <- .C(survival:::Cagfit5b,
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
      else   coxfit <- .C(survival:::Ccoxfit5b,
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

      iter <- i
      print(iter)
      print(coxfit$coef)
      if(iter==1) coef0 <- coxfit$coef
      else coef0 <- cbind(coef0, coxfit$coef)

      temp <- rep(FALSE, nvar+nfrail)
      if (nfrail>0) temp[1:nfrail] <- coxlist1$flag
      if (ptype >1) temp[nfrail+ 1:nvar] <- coxlist2$flag
      fdiag <- ifelse(temp, 0, coxfit$fdiag)


      if (nfrail>0) temp1 <- coxlist1$second
      else 	  temp1 <- 0
      if (ptype>1)  temp2 <- coxlist2$second
      else          temp2 <- 0

      dftemp <- survival:::coxpenal.df(matrix(coxfit$hmat,ncol=nvar),
                                       matrix(coxfit$hinv,ncol=nvar),  fdiag,
                                       assign, ptype, nvar,
                                       temp1, temp2, pindex[sparse])
      df <- dftemp$df
      pdf <- df[pterms>0]	          # df's for penalized terms
      trH <- dftemp$trH[pterms>0]   # trace H

      if (nfrail >0)  penalty <- -coxlist1$penalty
      else            penalty <- 0
      if (ptype >1) penalty <- penalty - coxlist2$penalty
      loglik1 <- coxfit$loglik + penalty  #C code returns PL - penalty
      if (iter==1) penalty0 <- penalty


      done <- TRUE
      for (i in 1:m) {
        temp <- eval(calls[i])
        if (sparse[i]) theta1 <- temp$theta
        done <- done & temp$done
        thetalist[[i]] <- temp$theta
        iterlist[[i]] <- temp
      }
      if(done) break

    }


    theta.init <- vector('list', m)
    for (i in 1:m)
      theta.init[[i]] <- iterlist[[i]]$history[,1]


  print(coef0)
  print(theta.init)

    for (i in 1:m) {
      temp <- (cfun[[i]])(parmlist[[i]], iter=0)
      if (sparse[i]) {
        theta1 <- temp$theta
        extra1 <- extralist[[i]]
      }
      thetalist[[i]] <- temp$theta
      lambdalist[[i]] <- lambda[1]*XG$m[i]
      iterlist[[i]] <- temp
    }

   stime <- y[ind, 1]
   sstat <- y[ind, 2]



   sparse.all <- unlist(sparse.where)

   group <- rep(1, n.coef)
   for(i in 1:m) group[pcols[[i]]] <- (i+1)
   if(length(sparse.all)==0) group <- rep(1, n.coef)
   else group[-sparse.all] <- 1
   print(group)
   ngroup <- length(table(group))


        kappa <- D <- vector('list', m)

        for(i in 1:m){
          kappa[[i]] <- theta.init[[i]]/(1-theta.init[[i]])
          D[[i]] <- sm[[i]]$D*0.25
        }

         K <- length(kappa[[1]])
        tuning.par <- expand.grid(lambda, kappa[[1]])


        if(sparse == "global" | !is.null(l)) TT <- 1
        else TT <-  length(beta.basis[[i]]$params)+1

        iter3 <- 1

        for(t in TT:1){


        newbeta <- matrix(0, nrow = K*nlambda, ncol = n.coef)
        gcv <- aic <- bic <- df <- loglik <- rep(0, K*nlambda)



        for(k in 1:K){
          init <- coef0[,k]
          print(init)
          lp <- c(x %*% init) + offset - sum(init*colMeans(x))
          score <- exp(lp[ind])
          resid <- Score(n, as.integer(method=='efron'), stime, sstat, newstrat,   score, weights)
          ss <- ii <- double(n)
          ss[ind] <- resid$score
          ii[ind] <- resid$infor


          DD <- lapply(1:m, function(l) sqrt(kappa[[l]][k])*D[[l]])
          DD <- as.matrix(bdiag(diag(0, n.par), as.matrix(bdiag(DD))))

          newz <- c(sqrt(ii)*(x%*%init + 1/ii*ss), rep(0, nrow(DD)))
          newx <-  rbind(sqrt(ii)*x%*%as.matrix(bdiag(W0)), -DD)

          if(ngroup == 1) nfactor <- 1
          else nfactor <- c(0, rep(1, ngroup-1))

          WW <- as.matrix(bdiag(W))
          H <- as.numeric(sqrt(t(init)%*%WW%*%init))
          print(H)
          for(j in 1:nlambda){
            print(scadderiv(H, alpha, lambda[j]))

            if(scadderiv(H, alpha, lambda[j])>0){
              xa <- newx[,-sparse.all]
              xb <- newx[,sparse.all]

              #Vb <- VV[,k:30]
              #Va <- VV[,1:(k-1)]


              xb <- xb*lambda[j]/scadderiv(H, alpha, lambda[j])
              Ha <- xa%*%solve(t(xa)%*%xa)%*%t(xa)

              zstar <- as.numeric(newz-Ha%*%newz)
              xstarb <- xb - Ha%*%xb

              print(dim(W0[[1]][sparse.all, sparse.all]))
              print(dim(xstarb))
              x2starb <- xstarb%*%W0[[1]][sparse.all, sparse.all]

              xb <- xb%*%W0[[1]][sparse.all, sparse.all]
              lasso <- gglasso:::gglasso(x=x2starb, y=zstar, group=rep(1,ncol(x2starb)), intercept=FALSE, lambda=lambda[j])# /sqrt(1+kappa))

              newbetab <- lasso$beta
              print(newbetab)
              newbetaa <- solve(t(xa)%*%xa)%*%t(xa)%*%(newz - xb%*%newbetab)
              newbeta[(k-1)*(nlambda) + j,] <- as.numeric(c(newbetaa, W0[[1]][sparse.all, sparse.all]%*%newbetab*lambda[j]/scadderiv(H, alpha, lambda[j]) ))#*sqrt(1+kappa)))
              coxLik  <-  survival:::coxph.fit(x[ind,], y[ind,], strata, offset, newbeta[(k-1)*(nlambda) + j,] , list(iter.max=0), weights=weights,
                                    method="breslow", row.names(x[ind,]))

              Omega <- t(DD)%*%DD/2

              H <- as.numeric(sqrt(t(newbeta[(k-1)*(nlambda) + j, ] )%*%WW%*%newbeta[(k-1)*(nlambda) + j, ] ))

              if(H == 0)  Sigma <- matrix(0, n.coef, n.coef)
              else Sigma <- n*WW/H/2
              nonzero <- newbeta[(k-1)*(nlambda) + j,]  != 0

              d2l <- solve(-coxLik$var)
              df[(k-1)*(nlambda) + j] <- sum(diag( solve(d2l[nonzero, nonzero] - Omega[nonzero, nonzero] - Sigma[nonzero, nonzero])%*%d2l[nonzero, nonzero] ) )
              loglik[(k-1)*(nlambda) + j] <- coxLik$loglik[1]

              if(tuning.method == "GCV") gcv[(k-1)*(nlambda) + j] <- -coxLik$loglik[1]/(n*(1-df[(k-1)*(nlambda) + j]/n)^2)
              else if (tuning.method == "AIC") aic[(k-1)*(nlambda) + j] <- -coxLik$loglik[1] + df[(k-1)*(nlambda) + j]
              else if (tuning.method == "BIC") bic[(k-1)*(nlambda) + j] <- -coxLik$loglik[1] + log(n)*df[(k-1)*(nlambda) + j]
            } else{
              newbeta[(k-1)*(nlambda) + j,] <- init
              coxLik  <-  survival:::coxph.fit(xx[ind,], y[ind,], strata, offset, init, list(iter.max=0), weights=weights,
                                    method="breslow", row.names(x[ind,]))

              Omega <- t(DD)%*%DD/2

              d2l <- solve(-coxLik$var)
              df[(k-1)*(nlambda) + j] <- sum(diag( solve(d2l - Omega)%*%d2l ) )
              loglik[(k-1)*(nlambda) + j] <- coxLik$loglik[1]

              if(tuning.method == "GCV") gcv[(k-1)*(nlambda) + j]  <- -coxLik$loglik[1]/n*(1-df[(k-1)*(nlambda) + j] /n)^2
              else if (tuning.method == "AIC") aic[(k-1)*(nlambda) + j]  <- -coxLik$loglik[1] + df[(k-1)*(nlambda) + j]
              else if (tuning.method == "BIC") bic[(k-1)*(nlambda) + j]  <- -coxLik$loglik[1] + log(n)*df[(k-1)*(nlambda) + j]
            }


            betaTRUE <- FALSE
            for (i in 1:m) {
              if(all(newbeta[(k-1)*(nlambda) + j, sparse.all] == 0)){
                betaTRUE <- TRUE
              }
            }

            if(betaTRUE) break
          }
        }

        niter <- nrow(aic)
        if(tuning.method == "AIC") which <- min((1:niter)[aic==min(aic)])
        else if (tuning.method == "BIC") which <- min((1:niter)[bic==min(bic)])
        else if (tuning.method == "GCV") which <- min((1:niter)[gcv==min(gcv)])


        coef <- newbeta[which,]
        ftheta <- tuning.par[which, 2]
        flambda <- tuning.par[which, 1]

        if(!any(coef==0)) break

        lbeta <- coef
        lftheta <- ftheta
        lflambda <- flambda
        lW <- W
        lsparse.where <- sparse.where


    if(t==1) break

    for (i in 1:m) {
      cutoff[[i]] <- cutoff[[i]] - 1

      if (sparse.what == "tail"){
        W[[i]] <- compute.W(cutoff[[i]], beta.basis[[i]])
        W0[[i]] <- as.matrix(bdiag(diag(1, cutoff[[i]]-1), solve(chol(W[[i]]))))
        W[[i]] <- as.matrix(Matrix::bdiag(matrix(0, cutoff[[i]]-1, cutoff[[i]]-1), W[[i]]))
      }
      sparse.where[[i]] <-  seq(pcols[[i]][1] + cutoff[[i]] -1, pcols[[i]][n.nonp])
    }

  }


  ### at convergence
  if(iter3==1) {
    lbeta <- beta
   #lfbeta <- fbeta
    lftheta <- ftheta
    lflambda <- flambda
   # lhistory <- history
    lW <- W
    lsparse.where <- sparse.where
  }

  for (i in 1:length(cfun)) {
    thetalist[[i]] <- lftheta
    lambdalist[[i]] <- lflambda
  }
  W <- lW

  if (andersen)  coxfit <- .C(survival:::Cagfit5b,
                              iter=as.integer(0),
                              as.integer(n),
                              as.integer(nvar),
                              as.integer(newstrat),
                              coef = as.double(lbeta),
                              u    = double(nvar+nfrail),
                              hmat = double(nvar*(nvar+nfrail)),
                              hinv = double(nvar*(nvar+nfrail)),
                              loglik = double(1),
                              flag = integer(1),
                              as.double(control$eps),
                              as.double(control$toler.chol),
                              as.integer(method=='efron'),
                              as.integer(nfrail),
                              fcoef = as.double(lfbeta),
                              fdiag = double(nfrail+nvar),
                              f.expr1,f.expr2,rho)
  else   coxfit <- .C(survival:::Ccoxfit5b,
                      iter=as.integer(0),
                      as.integer(n),
                      as.integer(nvar),
                      as.integer(newstrat),
                      coef = as.double(lbeta),
                      u    = double(nvar+nfrail),
                      hmat = double(nvar*(nvar+nfrail)),
                      hinv = double(nvar*(nvar+nfrail)),
                      loglik = double(1),
                      flag = integer(1),
                      as.double(control$eps),
                      as.double(control$toler.chol),
                      as.integer(method=='efron'),
                      as.integer(nfrail),
                      fcoef = as.double(lfbeta),
                      fdiag = double(nfrail+nvar),
                      f.expr1,f.expr2,rho)


  fit <- survival:::coxph.fit(xx, y, strata, offset, init=coxfit$coef, control=list(iter.max=0), weights, method, row.names(xx))

  if (nfrail >0) {
    lp <- offset + coxfit$fcoef[frailx]
    if (nvar >0)    #sparse frailties and covariates
      lp <- lp + x[,-fcol,drop=FALSE] %*%coxfit$coef -
        sum(means*coxfit1$coef)
  }
  else  lp <- offset + as.vector(x%*%coxfit$coef) - sum(means*coxfit$coef)

  # release the memory
  if (andersen) {
    .C(survival:::Cagfit5c, as.integer(nvar)) #release the memory
    resid <- .Call(survival:::Cagmart3,
                   y, exp(lp),
                   weights,
                   newstrat,
                   sorted,
                   as.integer(method=='efron'))
  }
  else  {
    expect <- .C(survival:::Ccoxfit5c, as.integer(n),
                 as.integer(nvar),
                 as.integer(newstrat),
                 as.integer(method=='efron'),
                 expect= double(n))$expect
    resid <- status - expect
  }
  names(resid) <- rownames

  #get the penalty portion of the second derive matrix
  if (nfrail>0) temp1 <- coxlist1$second
  else 	      temp1 <- 0
  if (ptype>1)  temp2 <- coxlist2$second
  else          temp2 <- 0

  ## Calculate df
  coef <- coxfit$coef
  names(coef) <- varnames
  nonzero <- coef != 0


  I0 <- solve(fit$var)
  S0 <- as.matrix(weights*score(fit, xx, y, strata))
  ES0 <- as.matrix(colMeans(S0))
  B <- crossprod(S0)/n -   ES0%*%t(ES0)

  hinv.full <- solve(-I0[nonzero, nonzero] - matrix(coxlist2$second,n.coef, n.coef)[nonzero, nonzero])
  B <- B[nonzero, nonzero]

  var <- n*hinv.full%*%B%*%hinv.full
  var0 <- -hinv.full

  var2 <- hinv.full%*%(I0[nonzero, nonzero])%*%hinv.full

  if (length(assign)==1){
    df=sum(diag( hinv.full%*%(-I0[nonzero, nonzero])) )
    trH=sum(diag(hinv.full))
  }else {
    df <- trH <- NULL
    d2 <- diag(hinv.full)
    for (i in assign) {
      temp <- coxph.wtest(hinv.full[i,i], var0[i,i])$solve
      if (is.matrix(temp)) df <- c(df, sum(diag(temp)))
      else                 df <- c(df, sum(temp))
      trH<- c(trH, sum(d2[i]))
    }
  }

  trH <- trH[pterms>0]   # trace H



  if (control$iter.max >1 && length(iterfail)>0)
    warning(paste("Inner loop failed to coverge for iterations",
                  paste(iterfail, collapse=' ')))



  names(iterlist) <- names(pterms[pterms>0])
  if (nfrail >0) {
    if (nvar >0) {   #sparse frailties and covariates
      list(coefficients  = coef,
           var    = var,
           var2   = var2,
           loglik = c(loglik0, loglik1),
           iter   = c(iter, iter2),
           linear.predictors = as.vector(lp),
           residuals = resid,
           means = means,
           method = method,
           class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
           df = df,
           penalty= c(penalty0, penalty),
           tuning.par = c( lftheta, lflambda),
           pterms = pterms, assign2=assign,
           history = lhistory,
           coxlist1=coxlist1,
           printfun=printfun)
    }
    else {  #sparse frailties only
      list( loglik = c(loglik0, loglik1),
            iter   = c(iter, iter2),
            linear.predictors = as.vector(lp),
            residuals = resid,
            means = means,
            method = method,
            class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
            df = df,
            penalty = c(penalty0, penalty),
            tuning.par = c( lftheta, lflambda),
            pterms = pterms, assign2=assign,
            history = lhistory,
            printfun=printfun)
    }
  }
  else {  #no sparse terms
    list(coefficients  = coef,
         var    = var,
         var2   = var2,
         loglik = c(loglik0, loglik1),
         iter   = c(iter, iter2),
         linear.predictors = lp,
         residuals = resid,
         means = means,
         method = method,
         class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
         df = df,
         penalty= c(penalty0, penalty),
         tuning.par = c( lftheta, lflambda),
         pterms = pterms, assign2=assign,
         history = lhistory,
         coxlist2=coxlist2,
         printfun= printfun)
  }
}


