wshoot <- function(p,x,y,init,weight,lambda,maxiter,tol,n)
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
      s<-ss2(j,tmpbeta,Q,B,n)
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

ss2 <- function(j,tmpb,Q,B, n)
{
  a <- sum(tmpb*Q[,j])-tmpb[j]*Q[j,j]
  s <- 1/n*(a-B[j])
  return(s)
}

