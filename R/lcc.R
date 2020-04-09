## Calculation of degress of freedom and variacne

df.f <- function(beta, penalty.where, dA, G, I, Dnrow){
  nvar <- length(beta)
  H <- I
  A <- I + diag(dA)

  if(Dnrow !=0){
    H[penalty.where, penalty.where] <- H[penalty.where, penalty.where] + G
    A[penalty.where, penalty.where] <- A[penalty.where, penalty.where] + G
  }


  zero <- penalty.where[beta[penalty.where]==0]

  var <- matrix(0, nvar, nvar)
  if(length(zero) == nvar) {
    df = 0
    var = rep(0, nvar*nvar)
  }else if(length(zero) ==0){
    if(Dnrow!=0) {
      df  <- sum( diag((solve(H)%*%I)))
    } else df <- nvar
    var <-  (solve(A)%*%I%*%solve(A))
  } else{
    if(Dnrow!=0){
      df  <- sum( diag((solve(H[-zero, -zero])%*%I[-zero, -zero])))
    } else df <- nvar - length(zero)
    var[-zero, -zero] <-  (solve(A)%*%I%*%solve(A))[-zero, -zero]
  }
   res <- list(df=df, var=as.vector(var), A = as.vector(A))
  return(res)
}
