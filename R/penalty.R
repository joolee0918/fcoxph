
#' @importFrom fda getbasisrange inprod
#' @importFrom Matrix bandSparse
#'

#positivepart <- function(fx) {
#  return(ifelse(fx>=0, fx, 0))
#}


#scadderiv <- function(ftheta, fa, flambda) {
#  return(flambda*(1-(1-apply(as.matrix(fa*flambda-abs(ftheta)), 1, positivepart)/((fa-1)*flambda))*as.numeric(abs(ftheta)>flambda)))
#}

#mcpderiv <- function(ftheta, fa, flambda) {
#  return( ifelse(abs(ftheta) < fa*flambda, (flambda - abs(ftheta)/fa)*sign(ftheta), 0))
#}


na.omit_pcox <- function (object, ...)
{
  n <- length(object)
  omit <- logical(nrow(object))
  vars <- seq_len(n)
  for (j in vars) {
    x <- object[[j]]
    if (!is.atomic(x))
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L)
      omit <- omit | x
    else omit <- omit | apply(x, 1, all)
  }
  xx <- object[!omit, , drop = FALSE]
  if (any(omit > 0L)) {
    temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
    attr(temp, "class") <- "omit"
    attr(xx, "na.action") <- temp
  }
  xx
}


