#' @importFrom parallel makeCluster stop
#'@importFrom doparallel registerDoParallel
#' @importFrom foreach foreach
cv.fcoxph <- function (fitobj, x, y, strats, cluster, weights, offset = NULL, control, init, lambda, nfolds, foldid,
            parallel = FALSE,  pcols, pattr, assign, npcols = npcols, tuning.method,
            sm, gamma, alpha, theta, nlambda, penalty, method,
            sparse.what, argvals, group.multiplier, ncluster)
  {
   ###Next line is commented out so each call generates its own lambda sequence
    # lambda=fitobj$lambda
    if(is.null(cluster)) {
      N <- nrow(x)
    }else {
      N <- length(unique(cluster))
    }

    if (is.null(foldid)){
      foldid = sample(rep(seq(nfolds), length = N))
    }else{ nfolds = max(foldid)
    }
    if(is.null(offset)) offset <- rep(0, N)
    if(is.null(weights)) weights <- rep(1, N)
    if (nfolds < 3)
      stop("nfolds must be bigger than 3; nfolds=10 recommended")

    lambda <- fitobj$lambda
    theta <- fitobj$theta
    cvraw = matrix(NA, nfolds, length(theta)*length(lambda))

    if (parallel) {
      cl <- parallel::makeCluster(ncluster)
      doParallel::registerDoParallel(cl)

      #  if (parallel && require(foreach)) {
      cvraw = foreach::foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
      {
        if(length(cluster)) {
          which = (1:N)[foldid==i]
          which = cluster %in% which
        } else which = foldid==i

        #      if (is.matrix(y))
        y_sub = y[!which, ]
        offset_sub =  offset[!which]
        weights_sub = weights[!which]

        if (length(strats)){
          strats_sub = strats[!which]
        }else {
          strats_sub = strats
        }

        out = fcoxpenal.fit(x = x[!which, , drop=FALSE], y =y_sub, strata = strats_sub, offset = offset_sub, init=init,
                                     control = control, weights=weights_sub, method=method,
                                     pcols = pcols, pattr = pattr, assign = assign, npcols = npcols, tuning.method = tuning.method,
                                     sm = sm,  gamma = gamma, alpha = alpha, theta = theta, lambda = lambda, penalty = penalty,
                                     sparse.what = sparse.what, argvals = argvals, group.multiplier = group.multiplier, cv.fit=TRUE)

        coefmat = out$beta
        if(ncol(y)==2){
          plfull = cox.deviance(x = x, y = y, offset = offset, strata = strats,
                                weights = weights, init = coefmat, method = method, control=control)
          plminusk = cox.deviance(x = x[!which, ], y = y_sub, strata = strats_sub, offset = offset_sub, weights = weights_sub,
                                  init = coefmat, method = method, control=control)
          plfull - plminusk

        } else if(ncol(y)==3){
          plfull = ag.deviance(x = x, y = y, offset = offset, strata = strats,
                               weights = weights, beta = coefmat, method = method)
          plminusk = ag.deviance(x = x[!which, ], y = y_sub, strata = strats_sub, offset = offset_sub, weights = weights_sub,
                                 beta = coefmat, method = method)
          plfull - plminusk

        }
      }
      parallel::stopCluster(cl)
      cvraw <- do.call('rbind', cvraw)
    }
    else {
      for (i in seq(nfolds)) {
        if(length(cluster)) {
          which = (1:N)[foldid==i]
          which = cluster %in% which
        } else which = foldid==i

       y_sub = y[!which, ]
       offset_sub =  offset[!which]
       weights_sub = weights[!which]

       if (length(strats)){
        strats_sub = strats[!which]
       }else {
         strats_sub = strats
       }

        out = fcoxpenal.fit(x = x[!which, , drop=FALSE], y =y_sub, strata = strats_sub, offset = offset_sub, init=init,
                                     control = control, weights=weights_sub, method=method,
                                     pcols = pcols, pattr = pattr, assign = assign, npcols = npcols, tuning.method = tuning.method,
                                     sm = sm,  gamma = gamma, alpha = alpha, theta = theta, lambda = lambda, penalty = penalty,
                                     sparse.what = sparse.what, argvals = argvals, group.multiplier = group.multiplier, cv.fit=TRUE)

        coefmat = out$beta
        if(ncol(y)==2){
          plfull = cox.deviance(x = x, y = y, offset = offset, strata = strats,
                                weights = weights, init = coefmat, method = method, control=control)
          plminusk = cox.deviance(x = x[!which, ], y = y_sub, strata = strats_sub, offset = offset_sub, weights = weights_sub,
                                  init = coefmat, method = method, control=control)
          cvraw[i, seq(along = plfull)] = plfull - plminusk

        } else if(ncol(y)==3){
          plfull = ag.deviance(x = x, y = y, offset = offset, strata = strats,
                                weights = weights, beta = coefmat, method = method)
          plminusk = ag.deviance(x = x[!which, ], y = y_sub, strata = strats_sub, offset = offset_sub, weights = weights_sub,
                                  beta = coefmat, method = method)
          cvraw[i, seq(along = plfull)] = plfull - plminusk

        }
    }
    }

   cvm = apply(cvraw, 2, mean, na.rm = TRUE)
   cvmin = which.min(cvm)

   return(cvmin)

  }

