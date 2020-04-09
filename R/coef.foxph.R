
#' @export
coef.fcoxph <-  function (x,  n){

   m <- length(x$smooth)
   fit <- sd.fit <- vector(mode = "list", length=m)
   for(i in 1:m){
   raw <- x$fcoxph$fs[[i]]$xind
   xx <- seq(min(raw), max(raw), length = n)

   if (x$by != "NA") {
  by <- rep(1, n)
  dat <- data.frame(x = xx, by = by)
  names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)
  }
else {
  dat <- data.frame(x = xx)
  names(dat) <- x$smooth[[i]]$term
}
  X <- PredictMat(x, dat)

  first <- x$smooth[[i]]$first.para
  last <- x$smooth[[i]]$last.para
  p <- x$coefficients[first:last]
  fit[[i]] <- X %*% p
  se.fit[[i]] <- sqrt(pmax(0, rowSums((X %*% x$var[first:last, first:last, drop = FALSE]) *
                                   X)))

  names(p) <- x$smooth[[i]]$term
   }

   return(p)
}
