
#' @importFrom mgcv PredictMat
#' @export
coef.fcoxph <-  function (x,  n=100){

   m <- length(x$smooth)
   fit <- vector(mode = "list", length=m)
   for(i in 1:m){
   raw <- x$fcoxph$fs[[i]]$xind
   xx <- seq(min(raw), max(raw), length = n)

  first <- x$smooth[[i]]$first.para
  last <- x$smooth[[i]]$last.para
  p <- x$coefficients[first:last]
  beta.basis <- x$smooth[[i]]$beta.basis
  X <- fda::eval.basis(xind, beta.basis)


  fit[[i]]$s <- xx
  fit[[i]]$value <- eval.fd(xx, p)
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
