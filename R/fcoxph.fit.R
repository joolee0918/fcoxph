#' @importFrom stats model.frame

fcoxph.fit <- function(formula, data, weights, subset, na.action,
                  init, control, ties= c("efron", "breslow", "exact"),
                  singular.ok =TRUE, robust=FALSE,
                  model=FALSE, x=FALSE, y=TRUE,  tt, method=ties,
                  tuning.method, nfolds, foldid, sm, alpha, gamma, theta, lambda, lambda.min.ratio, nlambda, penalty, L2penalty, sparse, argvals, group.multiplier,
                  parallel, ncluster,
                  ...) {

  ties <- match.arg(ties)
  Call <- match.call()


  ## We want to pass any ... args to coxph.control, but not pass things
  ##  like "dats=mydata" where someone just made a typo.  The use of ...
  ##  is simply to allow things like "eps=1e6" with easier typing
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(control.fcoxph)) #legal arg names
    indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched",
                    names(extraArgs)[indx==0L]), domain = NA)
  }
  if (missing(control)) control <- control.fcoxph(...)
  # create a call to model.frame() that contains the formula (required)
  #  and any other of the relevant optional arguments
  # then evaluate it in the proper frame
  indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                names(Call), nomatch=0)
  if (indx[1] ==0) stop("A formula argument is required")
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame)  # change the function called

  special <- c("strata", "cluster", "tt")
  temp$formula <- if(missing(data)) terms(formula, special)
  else              terms(formula, special, data=data)
  # Make "tt" visible for coxph formulas, without making it visible elsewhere
  if (!is.null(attr(temp$formula, "specials")$tt)) {
    coxenv <- new.env(parent= environment(formula))
    assign("tt", function(x) x, env=coxenv)
    environment(temp$formula) <- coxenv
  }


  mf <- eval(temp, parent.frame())
  if (nrow(mf) ==0) stop("No (non-missing) observations")
  Terms <- terms(mf)


  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type!='right' && type!='counting')
    stop(paste("Cox model doesn't support \"", type,
               "\" survival data", sep=''))
  data.n <- nrow(Y)   #remember this before any time transforms

  if (control$timefix) Y <- aeqSurv(Y)
  if (ncol(Y)==2) {
    case.n = sum(Y[,2])
  } else case.n = sum(Y[,3])

  if (length(attr(Terms, 'variables')) > 2) { # a ~1 formula has length 2
    ytemp <- survival:::terms.inner(formula[1:2])
    xtemp <- survival:::terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp))))
      warning("a variable appears on both the left and right sides of the formula")
  }

  # The time transform will expand the data frame mf.  To do this
  #  it needs Y and the strata.  Everything else (cluster, offset, weights)
  #  should be extracted after the transform
  #
  strats <- attr(Terms, "specials")$strata
  if (length(strats)) {
    stemp <- untangle.specials(Terms, 'strata', 1)
    if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
    strats <- as.numeric(strata.keep)
  }

  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt)) tt <- NULL
  if (length(timetrans)) {
    timetrans <- untangle.specials(Terms, 'tt')
    ntrans <- length(timetrans$terms)

    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights){ #default to O'Brien's logit rank
        obrien <- function(x) {
          r <- rank(x)
          (r-.5)/(.5+length(r)-r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt)) tt <- list(tt)  #single function becomes a list

    if (is.list(tt)) {
      if (any(!sapply(tt, is.function)))
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) ==1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")

    if (ncol(Y)==2) {

      if (length(strats)==0) {
        sorted <- order(-Y[,1], Y[,2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(strats, -Y[,1], Y[,2])
        #newstrat marks the first obs of each strata
        newstrat <-  as.integer(c(1, 1*(diff(strats[sorted])!=0)))
      }
      if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
      counts <- .Call(survival:::Ccoxcount1, Y[sorted,],
                      as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats)==0) {
        sort.end  <- order(-Y[,2], Y[,3])
        sort.start<- order(-Y[,1])
        newstrat  <- c(1L, rep(0, nrow(Y) -1))
      }
      else {
        sort.end  <- order(strats, -Y[,2], Y[,3])
        sort.start<- order(strats, -Y[,1])
        newstrat  <- c(1L, as.integer(diff(strats[sort.end])!=0))
      }
      if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
      counts <- .Call(survival:::Ccoxcount2, Y,
                      as.integer(sort.start -1L),
                      as.integer(sort.end -1L),
                      as.integer(newstrat))
      tindex <- counts$index
    }
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- 'right'  # new Y is right censored, even if the old was (start, stop]

    mf <- mf[tindex,]
    strats <- rep(1:length(counts$nrisk), counts$nrisk)
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights)))
      stop("weights must be finite")

    tcall <- attr(Terms, 'variables')[timetrans$terms+2]
    pvars <- attr(Terms, 'predvars')
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[,1], strats, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) { # It has a makepredictcall method
        dummy <- as.call(list(as.name(class(newtt)[1]), tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i]+2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }

  cluster<- attr(Terms, "specials")$cluster
  if(ncol(Y)==3 & length(cluster)==0) stop("subject ID must be given as cluster for time-dependent covariate")
  if (length(cluster)) {
    robust <- TRUE  #flag to later compute a robust variance
    tempc <- untangle.specials(Terms, 'cluster', 1:10)
    ord <- attr(Terms, 'order')[tempc$terms]
    if (any(ord>1)) stop ("Cluster can not be used in an interaction")
    cluster <- strata(mf[,tempc$vars], shortlabel=TRUE)  #allow multiples
    dropterms <- tempc$terms
    # Save away xlevels after removing cluster (we don't want to save upteen
    #  levels of that variable, which we will never need).
    xlevels <- .getXlevels(Terms[-tempc$terms], mf)
  }
  else {
    if (missing(robust)) robust <- FALSE
    xlevels <- .getXlevels(Terms, mf)
    dropterms <- NULL
  }

  contrast.arg <- NULL  #due to shared code with model.matrix.coxph
  attr(Terms, "intercept") <- 1
  stemp <- untangle.specials(Terms, 'strata', 1)
  hasinteractions <- FALSE
  if (length(stemp$vars) > 0) {  #if there is a strata statement
    for (i in stemp$vars) {  #multiple strata terms are allowed
      # The factors attr has one row for each variable in the frame, one
      #   col for each term in the model.  Pick rows for each strata
      #   var, and find if it participates in any interactions.
      if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
        hasinteractions <- TRUE
    }
    if (!hasinteractions) dropterms <- c(dropterms, stemp$terms)
  }

  if (length(dropterms)) {
    Terms2 <- Terms[ -dropterms]
    X <- model.matrix(Terms2, mf, constrasts=contrast.arg)
    # we want to number the terms wrt the original model matrix
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    temp <- temp + 1*(shift[1] <= temp)
    if (length(shift)==2) temp + 1*(shift[2] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts=contrast.arg)

  # infinite covariates are not screened out by the na.omit routines
  if (!all(is.finite(X)))
    stop("data contains an infinite predictor")

  # drop the intercept after the fact, and also drop strata if necessary
  Xatt <- attributes(X)
  if (hasinteractions) adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset==0)) offset <- rep(0., nrow(mf))
  else if (any(!is.finite(offset))) stop("offsets must be finite")

  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights)))
    stop("weights must be finite")

  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (missing(init)) init <- NULL
  else {
    if (length(init) != ncol(X)) stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init)
    if (any(temp < .Machine$double.min.exp | temp > .Machine$double.max.exp))
      stop("initial values lead to overflow or underflow of the exp function")
  }


  ncoef <- ncol(X)

  if (sum(Y[, ncol(Y)]) == 0) {
    # No events in the data!
    ncoef <- ncol(X)
    ctemp <- rep (NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance= c(concordant=0, discordant=0, tied.x=0, tied.y=0, tied.xy=0,
                   concordance=NA, std=NA)
    rval <- list(coefficients= ctemp,
                 var = matrix(0.0, ncoef, ncoef),
                 loglik=c(0,0),
                 score =0,
                 iter =0,
                 linear.predictors = offset,
                 residuals = rep(0.0, data.n),
                 means = colMeans(X), method=method,
                 n = data.n, nevent=0, terms=Terms, assign=assign,
                 concordance=concordance,
                 y = Y, call=Call)
    class(rval) <- c("fcoxph", "coxph")
    return(rval)
  }


  pterms <- sapply(mf, inherits, 'fcoxph.penalty')
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    #
    # Check the order of any penalty terms
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 'term.labels'))]
    if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
    pcols <- assign[match(pname, names(assign))]

    penalty.where <- as.numeric(unlist(pcols))


 if(sparse == "none"){

    fit <- survival:::coxpenal.fit(X, Y, strats, offset, init=init,
                         control,
                         weights=weights, method=method,
                         row.names(mf), pcols, pattr, assign)
  }else if(sparse != "none"){

    fit0 <- fcoxpenal.fit(x = X, y =Y, strata = strats,  offset = offset, init=init,
                           control = control, weights=weights, method=method,
                           pcols = pcols, pattr = pattr, assign = assign, npcols = npcols, tuning.method = tuning.method,
                           sm = sm,  gamma = gamma, alpha = alpha, theta = theta, lambda = lambda, lambda.min.ratio = lambda.min.ratio, nlambda = nlambda, penalty = penalty,
                          L2penalty = L2penalty, sparse.what = sparse, argvals = argvals, group.multiplier = group.multiplier, cv.fit=FALSE)


    lambda <- fit0$lambda
    theta <- fit0$theta

    minv <- sel <- minv1 <- minv2 <- minv3 <- sel1 <- sel2 <- sel3 <- rep(0, length(theta))
    fmin.all <- fsel.all <- rep(0, 3)
    for(i in 1:length(theta)){

    if (tuning.method == "aic") {
      minv[i] <- min(-2*fit0$loglik[[i]] + 2*fit0$df[[i]])
      sel[i] <- which.min(-2*fit0$loglik[[i]] + 2*fit0$df[[i]]) # choose lambda given theta
    }else if (tuning.method == "bic") {
      minv[i] <- min(-2*fit0$loglik[[i]]+ log(data.n)*fit0$df[[i]])
      sel[i] <- which.min(-2*fit0$loglik[[i]]+ log(data.n)*fit0$df[[i]]) #+ 0.5*fit0$df*log(length(penalty.where)))
    }else if(tuning.method == "gcv") {
      minv[i] <- min(-fit0$loglik[[i]]/(data.n*(1-fit0$df[[i]]/data.n)^2) )
      sel[i] <- which.min(-fit0$loglik[[i]]/(data.n*(1-fit0$df[[i]]/data.n)^2) )
    }else {
      minv1[i] <- .min(-2*fit0$loglik[[i]] + 2*fit0$df[[i]])
      minv2[i] <- min(-2*fit0$loglik[[i]]+ log(data.n)*fit0$df[[i]])
      minv3[i] <- min(-fit0$loglik[[i]]/(data.n*(1-fit0$df[[i]]/data.n)^2) )
      sel1[i] <- which.min(-2*fit0$loglik[[i]] + 2*fit0$df[[i]])
      sel2[i] <- which.min(-2*fit0$loglik[[i]]+ log(data.n)*fit0$df[[i]]) #+ 0.5*fit0$df*log(length(penalty.where)))
      sel3[i] <-  which.min(-fit0$loglik[[i]]/(data.n*(1-fit0$df[[i]]/data.n)^2) )
    }
    }

    if(is.null(tuning.method)){
      fminv.all[1] <- min(minv1)
      fminv.all[2] <- min(minv2)
      fminv.all[3] <- min(minv3)
      fsel.all[1] <- which.min(minv1)
      fsel.all[2] <- which.min(minv2)
      fsel.all[3] <- which.min(minv3)
    } else{
    fminv <- min(minv)
    fsel <- which.min(minv) # choose theta
    }

    if(!is.null(tuning.method)){

    fit <- list()

    fit$coefficients <- fit0$beta[[fsel]][, sel[fsel]]
    names(fit$coefficients) <- fit0$varnames
    fit$history <- fit0$beta
    nvar <- length(fit$coefficients)
    fit$var <- matrix(fit0$var[[fsel]][,sel[fsel]], nvar, nvar)
    fit$A <- matrix(fit0$A[[fsel]][,sel[fsel]], nvar, nvar)
    fit$u <- fit0$u[[fsel]][, sel[fsel]]
    zero <- penalty.where[fit$coefficients[penalty.where]==0]

    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {

      fit$naive.var <- fit$var

      if(length(zero) == nvar) {
        fit$var = rep(0, nvar*nvar)
      }else{
        ny <- ncol(Y)
        nstrat <- as.numeric(strats)

        status <- Y[,ny,drop=TRUE]

        if (is.null(strats)) {
          ord <- order(Y[,ny-1], -status)
          newstrat <- rep(0,data.n)
        }else {
          ord <- order(nstrat, Y[,ny-1], -status)
          newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
        }
        newstrat[data.n] <- 1

        # sort the data
        xx <- X[ord,]
        yy <- Y[ord,]

        if (is.null(weights)) weights <- rep(1, data.n)
        temp2 = sum(weights)
        means = sapply(1:nvar, function(i) sum(weights*X[, i])/temp2)
        score <- exp(c(as.matrix(X) %*% fit$coefficients) + offset - sum(fit$coefficients*means))[ord]
        if (ny==2) {
          resid <- .C(survival:::Ccoxscore, as.integer(data.n),
                      as.integer(nvar),
                      as.double(yy),
                      x=as.double(xx),
                      as.integer(newstrat),
                      as.double(score),
                      as.double(weights[ord]),
                      as.integer(method=='efron'),
                      resid= double(n*nvar),
                      double(2*nvar))$resid
        }
        else {
          resid<- .C(survival:::Cagscore,
                     as.integer(data.n),
                     as.integer(nvar),
                     as.double(yy),
                     as.double(xx),
                     as.integer(newstrat),
                     as.double(score),
                     as.double(weights[ord]),
                     as.integer(method=='efron'),
                     resid=double(n*nvar),
                     double(nvar*6))$resid
        }

        if (nvar >1) {
          rr <- matrix(0, data.n, nvar)
          rr[ord,] <- matrix(resid, ncol=nvar)
        }else rr[ord] <- resid

        if (!missing(cluster)) {
          if (length(cluster) !=data.n) stop("Wrong length for 'cluster'")
          rr <- drop(rowsum(rr, cluster))
        }

        rr <- rr * weights

        A <- matrix(fit0$A[[fsel]][,sel[fsel]], nvar, nvar)
        B <- t(rr)%*%rr
        fit$var <- matrix(0, nvar, nvar)

        if(length(zero) ==0) fit$var <- solve(A)%*%B%*%solve(A)
        else fit$var[-zero, -zero] <- (solve(A)%*%B%*%solve(A))[-zero, -zero]
      }
    }

    fit$penalty <- fit$penalty
    fit$loglik <- c(fit0$loglik0, fit0$loglik[sel])

    fit$aic <-  -2*fit0$loglik[[fsel]][sel[fsel]] + 2*fit0$df[[fsel]][sel[fsel]]
    fit$bic <- -2*fit0$loglik[[fsel]][sel[fsel]] + log(data.n)*fit0$df[[fsel]][sel[fsel]]

    if(sel%%length(lambda) ==0) lambda.where <- length(lambda)
    else lambda.where <- sel%%length(lambda)
    fit$penalty.par <- c(theta = theta[ceiling(sel/length(lambda))], lambda = lambda[lambda.where])
    fit$lambda <- fit0$lambda
    fit$theta <- fit0$theta
    fit$pterms <- pterms

    }
  } else{


    fit.all <- lapply(1:3, function(k){
    fminv <- fmin.all[k]
    fsel <- fsel.all[k]
    sel <- get(paste0("sel"),k)

    fit <- list()

    fit$coefficients <- fit0$beta[[fsel]][, sel[fsel]]
    names(fit$coefficients) <- fit0$varnames
    fit$history <- fit0$beta
    nvar <- length(fit$coefficients)
    fit$var <- matrix(fit0$var[[fsel]][,sel[fsel]], nvar, nvar)
    fit$A <- matrix(fit0$A[[fsel]][,sel[fsel]], nvar, nvar)
    fit$u <- fit0$u[[fsel]][, sel[fsel]]
    zero <- penalty.where[fit$coefficients[penalty.where]==0]

    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {

      fit$naive.var <- fit$var

      if(length(zero) == nvar) {
        fit$var = rep(0, nvar*nvar)
      }else{
        ny <- ncol(Y)
        nstrat <- as.numeric(strats)

        status <- Y[,ny,drop=TRUE]

        if (is.null(strats)) {
          ord <- order(Y[,ny-1], -status)
          newstrat <- rep(0,data.n)
        }else {
          ord <- order(nstrat, Y[,ny-1], -status)
          newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
        }
        newstrat[data.n] <- 1

        # sort the data
        xx <- X[ord,]
        yy <- Y[ord,]

        if (is.null(weights)) weights <- rep(1, data.n)
        temp2 = sum(weights)
        means = sapply(1:nvar, function(i) sum(weights*X[, i])/temp2)
        score <- exp(c(as.matrix(X) %*% fit$coefficients) + offset - sum(fit$coefficients*means))[ord]
        if (ny==2) {
          resid <- .C(survival:::Ccoxscore, as.integer(data.n),
                      as.integer(nvar),
                      as.double(yy),
                      x=as.double(xx),
                      as.integer(newstrat),
                      as.double(score),
                      as.double(weights[ord]),
                      as.integer(method=='efron'),
                      resid= double(n*nvar),
                      double(2*nvar))$resid
        }
        else {
          resid<- .C(survival:::Cagscore,
                     as.integer(data.n),
                     as.integer(nvar),
                     as.double(yy),
                     as.double(xx),
                     as.integer(newstrat),
                     as.double(score),
                     as.double(weights[ord]),
                     as.integer(method=='efron'),
                     resid=double(n*nvar),
                     double(nvar*6))$resid
        }

        if (nvar >1) {
          rr <- matrix(0, data.n, nvar)
          rr[ord,] <- matrix(resid, ncol=nvar)
        }else rr[ord] <- resid

        if (!missing(cluster)) {
          if (length(cluster) !=data.n) stop("Wrong length for 'cluster'")
          rr <- drop(rowsum(rr, cluster))
        }

        rr <- rr * weights

        A <- matrix(fit0$A[[fsel]][,sel[fsel]], nvar, nvar)
        B <- t(rr)%*%rr
        fit$var <- matrix(0, nvar, nvar)

        if(length(zero) ==0) fit$var <- solve(A)%*%B%*%solve(A)
        else fit$var[-zero, -zero] <- (solve(A)%*%B%*%solve(A))[-zero, -zero]
      }
    }

    fit$penalty <- fit$penalty
    fit$loglik <- c(fit0$loglik0, fit0$loglik[sel])

    fit$aic <-  -2*fit0$loglik[[fsel]][sel[fsel]] + 2*fit0$df[[fsel]][sel[fsel]]
    fit$bic <- -2*fit0$loglik[[fsel]][sel[fsel]] + log(data.n)*fit0$df[[fsel]][sel[fsel]]

    if(sel%%length(lambda) ==0) lambda.where <- length(lambda)
    else lambda.where <- sel%%length(lambda)
    fit$penalty.par <- c(theta = theta[ceiling(sel/length(lambda))], lambda = lambda[lambda.where])
    fit$lambda <- fit0$lambda
    fit$theta <- fit0$theta
    fit$pterms <- pterms
    })
    fit <- fit.all
    names(fit) <- c("aic", "bic", "gcv")
  }
  }


    for (i in 1:length(sm)) {
      tmpid <- pcols[[i]]
      start <- 1
      idx <- tmpid[start:(start + length(pcols[[i]]) -1)]
      names(fit$coefficients)[idx] <- paste(sm[[i]]$label,1:length(idx), sep = ".")
      sm[[i]]$first.para <- min(idx)
      sm[[i]]$last.para <- max(idx)
      start <- start + length(idx)
    }

  }


  fit$n <- data.n
  fit$nevent <- sum(Y[,ncol(Y)])
  fit$terms <- Terms
  fit$assign <- assign
  fit$assign2 <- pcols
  fit$smooth <- sm



  if (model) {
    if (length(timetrans)) {
      # Fix up the model frame -- still in the thinking stage
      mf[[".surv."]]   <- Y
      mf[[".strata."]] <- strats
      stop("'model=TRUE' not supported for models with tt terms")
    }
    fit$model <- mf
  }
  if (x)  {
    fit$x <- X
    if (length(strats)) {
      if (length(timetrans)) fit$strata <- strats
      else     fit$strata <- strata.keep
    }
  }
  if (y)  fit$y <- Y

  class(fit) <- c('fcoxph', 'fcoxph.penal')

  return(fit)

}

