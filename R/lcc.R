wshoot1 <- function(p,x,y,init,weight,lambda,maxiter,tol,n)
{
  Q = t(x)%*%x
  B = t(x)%*%y
  i=0
  status = 0
  lams =lambda*weight
  oldbeta <- init
  tmpbeta <- oldbeta

  while (i<maxiter && status==0){
    for (j in 1:p){
      s<-ss(j,tmpbeta,Q,B,n)
      if (s > lams[j])
        tmpbeta[j]<-(lams[j]-s)/(1/n*Q[j,j])
      else if (s < (-lams[j]) )
        tmpbeta[j]<-(-lams[j]-s)/(1/n*Q[j,j])
      else
        tmpbeta[j]<- 0.0
    }
    dx<-max(abs(tmpbeta-oldbeta))
    oldbeta <- tmpbeta
    if (dx<=tol)
      status <- 1
    i <- i+1
  }
  tmpbeta
}

ss <- function(j,tmpb,Q,B, n)
{
  a <- sum(tmpb*Q[,j])-tmpb[j]*Q[j,j]
  s <- 1/n*(a-B[j])
  return(s)
}


df.f <- function(beta, penalty.where, dA, G, I){
  nvar <- length(beta)
  H <- I
  H[penalty.where, penalty.where] <- H[penalty.where, penalty.where] + G

  A <- I + dA
  A[penalty.where, penalty.where] <- A[penalty.where, penalty.where] + G


  zero <- penalty.where[beta[penalty.where]==0]
  var <- matrix(0, nvar, nvar)
  if(sum(nzero) == 0) {
    df = 0
    var = rep(0, nvar*nvar)
  }else{
    df  <- sum( diag((solve(H[-zero, -zero])%*%I[-zero, -zero])))
    var[-zero, -zero] <-  (solve(A[-zero, -zero])%*%I[-zero, -zero]%*%solve(A[-zero, -zero]))
  }
   res <- list(df=df, var=as.vector(var), A = as.vector(A))
  return(res)
}
