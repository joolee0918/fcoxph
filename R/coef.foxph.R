
#' @importFrom mgcv PredictMat
#' @importFrom fda fd eval.fd


#' @export
T1.hat = function(coef.res)
{
  knots <- coef.res[,1]
  beta <- coef.res[,2]
  T1n = length(beta)
  t1hatfuntemp=knots[T1n]
  if(sum(abs(beta)[1:T1n])==0){t1hatfuntemp=knots[1]}else{
    for(t1i in 1:(T1n-1))
    {
      if(beta[t1i]!=0 && sum(abs(beta)[(t1i+1):T1n])==0)
      {t1hatfuntemp = knots[t1i+1]}
    }}
  return(t1hatfuntemp)
}


#' @export
coef.fcoxph <-  function (x, n=NULL){

   m <- length(x$smooth)
   fit <- vector(mode = "list", length=m)

   for(i in 1:m){
     if(is.null(n)) n <-  x$smooth[[i]]$beta.basis$nbasis - x$smooth[[i]]$m[1] + 1
     xx <- seq(x$smooth[[i]]$beta.basis$rangeval[1], x$smooth[[i]]$beta.basis$rangeval[2], length = n)

  first <- x$smooth[[i]]$first.para
  last <- x$smooth[[i]]$last.para
  p <- x$coefficients[first:last]
  beta.basis <- x$smooth[[i]]$beta.basis
  X <- fda::eval.basis(xx, beta.basis)
  betafd <- fda::fd(coef=p,basisobj=beta.basis)

  fit[[i]]$s <- xx
  fit[[i]]$value <- fda::eval.fd(xx, betafd)
  fit[[i]]$se <- sqrt(pmax(0, rowSums((X %*% x$var[first:last, first:last, drop = FALSE]) *
                                   X)))
  if(!is.null(x$naive.var))  fit[[i]]$naive.se <- sqrt(pmax(0, rowSums((X %*% x$naive.var[first:last, first:last, drop = FALSE]) *
                                                         X)))

  fit[[i]] <- as.data.frame(fit[[i]])
  colnames(fit[[i]]) <- c("s", "value", "se")
  if(!is.null(x$naive.var)) colnames(fit[[i]])[4] <- c("naive.se")

   }

   names(fit) <- sapply(1:m, function(i) x$smooth[[i]]$term)

   return(fit)
}
