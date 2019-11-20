#' @importFrom stats model.frame

fcoxph.fit <- function(formula, data, weights, subset, na.action,
                  init, control, ties= c("efron", "breslow", "exact"),
                  singular.ok =TRUE, robust=FALSE,
                  model=FALSE, x=FALSE, y=TRUE,  tt, method=ties,
                  tuning.method, sm, cutoff, alpha, theta, lambda, lambda.min, nlambda, penalty, sparse, argvals, group.multiplier,
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


  pterms <- sapply(mf, inherits, 'coxph.penalty')
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    #
    # Check the order of any penalty terms
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 'term.labels'))]
    if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
    pcols <- assign[match(pname, names(assign))]
    npcols <- assign[!match(pname, names(assign))]

   if(sparse == "none"){
    fit <- survival:::coxpenal.fit(X, Y, strats, offset, init=init,
                         control,
                         weights=weights, method=method,
                         row.names(mf), pcols, pattr, assign)



  } else if(sparse != "none"){

    fit <- fcoxpenal.fit(X, Y, strats, offset, init=init,
                        control,
                        weights=weights, method=method,
                       row.names(mf), pcols, pattr, assign, npcols, tuning.method, sm, cutoff, alpha, theta, lambda, lambda.min, nlambda, penalty, sparse, argvals, group.multiplier)


  }

  }

  fit$n <- data.n
  fit$nevent <- sum(Y[,ncol(Y)])
  fit$terms <- Terms
  fit$assign <- assign
  class(fit) <- fit$class
  fit$class <- NULL



  return(fit)
}

