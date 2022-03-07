#' @importFrom parallel makeCluster stopCluster
#'@importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
cv.fcoxph <- function (fitobj, x, y, strats, cluster, weights, offset = NULL, control, init, lambda, lambda.min.ratio, nfolds, foldid,
            parallel = FALSE,  pcols, pattr, assign, npcols = npcols, tuning.method,
            sm, gamma, alpha, theta, nlambda, penalty, method,
            L2penalty, sparse.what, argvals, group.multiplier, ncluster)
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
    nlambda <- length(lambda)
    ntheta <- length(theta)
    cvraw = matrix(NA, nfolds, ntheta*nlambda)
    plk <- matrix(rep(0, ntheta*nlambda), nrow=ntheta)

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
                                     sm = sm,  gamma = gamma, alpha = alpha, theta = theta, lambda = lambda, lambda.min.ratio = lambda.min.ratio, penalty = penalty,
                                     L2penalty = L2penalty, sparse.what = sparse.what, argvals = argvals, group.multiplier = group.multiplier, cv.fit=TRUE)

        coefmat = out$beta
        theta <- out$theta
        lambda <- out$lambda
        if(ncol(y)==2){
          for(k in 1:length(theta)){
            plk[k,] = cox.deviance(x = x[which, ], y = y[which,], strata = strats[which], offset = offset[which], weights = weights[which],
                                   init = coefmat[[k]], method = method, control=control)
          }
            
            cvraw[i, ] = as.vector(plk)
            
            
          } else if(ncol(y)==3){
            for(k in 1:length(theta)){
            plk[k,] = ag.deviance(x = x[which, ], y = y[which,], strata = strats[which], offset = offset[which], weights = weights[which],
                                  init = coefmat[[k]], method = method, control=control)
            }
            cvraw[i, ] = as.vector(plk)
            
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
                                     sm = sm,  gamma = gamma, alpha = alpha, theta = theta, lambda = lambda,lambda.min.ratio = lambda.min.ratio, penalty = penalty,
                                     L2penalty = L2penalty, sparse.what = sparse.what, argvals = argvals, group.multiplier = group.multiplier, cv.fit=TRUE)

        coefmat = out$beta
        
        if(ncol(y)==2){
          
          for(k in 1:length(theta)){
          plk[k,] = cox.deviance(x = x[which, ], y = y[which,], strata = strats[which], offset = offset[which], weights = weights[which],
                                  init = coefmat[[k]], method = method, control=control)
          }
          
          cvraw[i, ] = as.vector(plk)

        } else if(ncol(y)==3){
          for(k in 1:length(theta)){
          plk[k,] = ag.deviance(x = x[which, ], y = y[which,], strata = strats[which], offset = offset[which], weights = weights[which],
                                  init = coefmat[[k]], method = method, control=control)
          }
          cvraw[i, ] = as.vector(plk)

      
    }
      }
    }

   
   cvm = apply(cvraw, 2, mean, na.rm = TRUE)
   cvmin = which.min(cvm)
   which.lambda <- ceiling(cvmin / ntheta)
   which.theta <- cvmin %% ntheta
   if(which.theta==0) which.theta <- ntheta
   
   return(list(cvmin=cvmin, which.lambda=which.lambda, which.theta=which.theta, opt.lambda=lambda[which.lambda], opt.theta =theta[which.theta]))

  }

