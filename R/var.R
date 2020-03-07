
## df / Variance
if(var = TRUE){

  for(i in 1:nlambda){
    nonzero <- fit.beta[,i] !=0
    n.nonzero <- sum(nonzero)
    newxx <- xx[, nonzero]
    beta.nonzero <- fit.beta[nonzero,i]


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
            temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambda, init, nystar, D)  #init, nystar, penalty)
          else    temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], lambda, init, nystar, extralist[[i]])# init, nystar, penalty, extralist[[i]])
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
    # "Unpack" the passed in paramter list,
    #   and make the initial call to each of the external routines
    #

    rho<-environment()

    if (andersen){ coxfit <- .C(survival:::Cagfit5a,
                                as.integer(n),
                                as.integer(n.nonzero),
                                y,
                                newxx ,
                                offset,
                                weights,
                                newstrat,
                                sort,
                                means= double(n.nonzero),
                                coef= as.double(beta.nonzero),
                                u = double(n.nonzero),
                                loglik=double(1),
                                as.integer(method=='efron'),
                                as.integer(ptype),
                                as.integer(full.imat),
                                as.integer(nfrail),
                                as.integer(frailx),
                                #R callback additions
                                f.expr1,f.expr2,rho)

    }else      { coxfit <- .C(survival:::Ccoxfit5a,as.integer(0),
                              as.integer(n.nonzero),
                              y,
                              newxx,
                              offset,
                              weights,
                              newstrat,
                              sort,
                              means= double(n.nonzero),
                              coef= as.double(beta.nonzero),
                              u = double(n.nonzero),
                              loglik=double(1),
                              as.integer(method=='efron'),
                              as.integer(ptype),
                              as.integer(full.imat),
                              as.integer(nfrail),
                              as.integer(frailx),
                              f.expr1,f.expr2,rho)

    }


    if (andersen) { coxfit0 <- .C(survival:::Cagfit5b,
                                  iter=as.integer(0),
                                  as.integer(n),
                                  as.integer(n.nonzero),
                                  as.integer(newstrat),
                                  coef = as.double(p.betar[nonzero,i]),
                                  u    = double(n.nonzero+nfrail),
                                  hmat = double(n.nonzero*(n.nonzero+nfrail)),
                                  hinv = double(n.nonzero*(n.nonzero+nfrail)),
                                  loglik = double(1),
                                  flag = integer(1),
                                  as.double(control$eps),
                                  as.double(control$toler.chol),
                                  as.integer(method=='efron'),
                                  as.integer(nfrail),
                                  fcoef = as.double(finit),
                                  fdiag = double(nfrail+n.nonzero),
                                  f.expr1,f.expr2,rho)
    }else  {coxfit0 <- .C(survival:::Ccoxfit5b,
                          iter=as.integer(0),
                          as.integer(n),
                          as.integer(n.nonzero),
                          as.integer(newstrat),
                          coef = as.double(p.betar[,i]),
                          u    = double(n.nonzero+nfrail),
                          hmat = double(n.nonzero*(n.nonzero+nfrail)),
                          hinv = double(n.nonzero*(n.nonzero+nfrail)),
                          loglik = double(1),
                          flag = integer(1),
                          as.double(control$eps),
                          as.double(control$toler.chol),
                          as.integer(method=='efron'),
                          as.integer(nfrail),
                          fcoef = as.double(finit),
                          fdiag = double(nfrail+n.nonzero),
                          f.expr1,f.expr2,rho)
    }


    I <- solve(matrix(coxfit$imat, nvar, nvar))

    if (nfrail>0) temp1 <- coxlist1$second
    else 	      temp1 <- 0
    if (ptype>1)  temp2 <- coxlist2$second
    else          temp2 <- 0

    dftemp <-coxpenal.df(matrix(coxfit$hmat,ncol=nvar),
                         matrix(coxfit$hinv,ncol=nvar),  fdiag,
                         assign, ptype, nvar,
                         temp1, temp2, pindex[sparse])
    df <- dftemp$df
    trH <- dftemp$trH
    var <- dftemp$var
    var2  <- dftemp$var2

    nonzero <- p.fit$beta[,i]!=0
    if(sum(nonzero) != 0){
      eval.penalty <- penalty.v(init[nonzero], p.fit$beta[nonzero,i], thetalist[[1]], p.lambda[i], alpha,
                                nystar, extralist[[1]][nonzero, nonzero], penalty, I[nonzero, nonzero])
      eval.tuning <- tuning.summary(tuning.method, loglik0[2], eval.penalty$df, n)
      tuning.s[i,] <- c(p.lambda[i], eval.penalty$df, eval.tuning)
    } else{
      tuning.s[i,] <-c(rep(NA, 5))
    }
  }



  names(iterlist) <- names(pterms[pterms>0])
  if (nfrail >0) {
    if (nvar >0) {   #sparse frailties and covariates
      list(coefficients  = coef,
           var    = var,
           var2   = var2,
           loglik = c(loglik0, loglik1),
           iter   = c(iter, iter3),
           linear.predictors = as.vector(lp),
           residuals = resid,
           means = means,
           method = method,
           class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
           df = df,
           penalty= c(penalty0, penalty),
           tuning.par = c( lftheta, lflambda),
           pterms = pterms, assign2=assign,
           #history = lhistory,
           coxlist1=coxlist1,
           printfun=printfun)
    }
    else {  #sparse frailties only
      list( loglik = c(loglik0, loglik1),
            iter   = c(iter, iter3),
            linear.predictors = as.vector(lp),
            residuals = resid,
            means = means,
            method = method,
            class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
            df = df,
            penalty = c(penalty0, penalty),
            tuning.par = c( lftheta, lflambda),
            pterms = pterms, assign2=assign,
            #history = lhistory,
            printfun=printfun)
    }
  }
  else {  #no sparse terms
    list(coefficients  = coef,
         var    = var,
         var2   = var2,
         loglik = c(loglik0, loglik1),
         iter   = c(iter, iter3),
         linear.predictors = lp,
         residuals = resid,
         means = means,
         method = method,
         class = c('fcoxph.penal', 'coxph.penal', 'coxph'),
         df = df,
         penalty= c(penalty0, penalty),
         tuning.par = c( lftheta, lflambda),
         pterms = pterms, assign2=assign,
         #history = lhistory,
         coxlist2=coxlist2,
         printfun= printfun)
  }
}
}
