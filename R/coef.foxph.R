
#' @importFrom mgcv PredictMat
#' @export
coef.fcoxph <-  function (x,  n){

   m <- length(x$smooth)
   fit <- vector(mode = "list", length=m)
   for(i in 1:m){
   raw <- x$fcoxph$fs[[i]]$xind
   xx <- seq(min(raw), max(raw), length = n)

   if (x$smooth[[i]]$by != "NA") {
  by <- rep(1, n)
  dat <- data.frame(x = xx, by = by)
  names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)
  }else {
  dat <- data.frame(x = xx)
  names(dat) <- x$smooth[[i]]$term
}
  X <- mgcv::PredictMat(x$smooth[[i]], dat)

  first <- x$smooth[[i]]$first.para
  last <- x$smooth[[i]]$last.para
  p <- x$coefficients[first:last]
  fit[[i]]$s <- dat[x$smooth[[i]]$term]
  fit[[i]]$value <- X %*% p
  fit[[i]]$se <- sqrt(pmax(0, rowSums((X %*% x$var[first:last, first:last, drop = FALSE]) *
                                   X)))
  if(x$naive.var)  fit[[i]]$naive.se <- sqrt(pmax(0, rowSums((X %*% x$naive.var[first:last, first:last, drop = FALSE]) *
                                                         X)))

  fit[[i]] <- as.data.frame(fit[[i]])

  names(fit[[i]]) <- x$smooth[[i]]$term
   }

   return(fit)
}
