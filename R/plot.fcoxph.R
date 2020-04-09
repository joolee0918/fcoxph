plot.fcoxph <- function (x, rug = TRUE, se = TRUE, pages = 0,
                          scale = -1, n = 100,
                          jit = FALSE, xlab = NULL, ylab = NULL,
                          main = NULL, ylim = NULL, xlim = NULL,
                          shade = FALSE, shade.col = "gray80", shift = 0, scheme) {
  sub.edf <- function(lab, edf) {
    pos <- regexpr(":", lab)[1]
    if (pos < 0) {
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab, start = 1, stop = pos),
                   ",", round(edf, digits = 2), ")", sep = "")
    }
    else {
      lab1 <- substr(lab, start = 1, stop = pos - 2)
      lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
      lab <- paste(lab1, ",", round(edf, digits = 2), lab2,
                   sep = "")
    }
    lab
  }

  m <- length(x$smooth)

  order <- if (is.list(x$pterms))
    unlist(lapply(x$pterms, attr, "order"))
  else attr(x$pterms, "order")

  if (se) {
    if (is.numeric(se))
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0)
      se1.mult <- 0
    if (se2.mult < 0)
      se2.mult <- 0
  }
  else se1.mult <- se2.mult <- 1
  if (se && x$var[1, 1] < 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  pd <- list()
  i <- 1
  if (m > 0)
    for (i in 1:m) {
      first <- x$smooth[[i]]$first.para
      last <- x$smooth[[i]]$last.para
      attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
      P <- refund:::plot.mgcv.smooth(x$smooth[[i]], P = NULL, data = x$fcoxph$ft[[i]]$data,
                se = se,
                scale = scale, n = n, jit = jit, xlab = xlab, ylab = ylab,
                main = main, label = term.lab, ylim = ylim, xlim = xlim,
                shade = shade, shade.col = shade.col, shift=shift, scheme=scheme[i],
                se1.mult = se1.mult, se2.mult = se2.mult,...)
      if (is.null(P))
        pd[[i]] <- list(plot.me = FALSE)
      else if (is.null(P$fit)) {
        p <- x$coefficients[first:last]
        offset <- attr(P$X, "offset")
        if (is.null(offset))
          P$fit <- fit <- P$X %*% p
        else P$fit <- P$X %*% p + offset
        if (!is.null(P$exclude))
          P$fit[P$exclude] <- NA
        if (se && P$se) {
         se.fit <- sqrt(pmax(0, rowSums((P$X %*%x$var[first:last, first:last, drop = FALSE]) * P$X)))
          if (!is.null(P$exclude))
            P$se.fit[P$exclude] <- NA
        }
       if (se && P$se)
          P$se <- se.fit * P$se.mult
        P$X <- NULL
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
      else {
         P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
    }
  n.plots <- n.para
  if (m > 0)
    for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me)
  if (n.plots == 0)
    stop("No terms to plot - nothing for plot.gam() to do.")
  if (pages > n.plots)
    pages <- n.plots
  if (pages < 0)
    pages <- 0
  if (pages != 0) {
    ppp <- n.plots%/%pages
    if (n.plots%%pages != 0) {
      ppp <- ppp + 1
      while (ppp * (pages - 1) >= n.plots) pages <- pages -
          1
    }
    c <- r <- trunc(sqrt(ppp))
    if (c < 1)
      r <- c <- 1
    if (c * r < ppp)
      c <- c + 1
    if (c * r < ppp)
      r <- r + 1
    oldpar <- par(mfrow = c(r, c))
  }
  else {
    ppp <- 1
    oldpar <- par()
  }
  if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
      pages > 1 && dev.interactive())
    ask <- TRUE
  else ask <- FALSE
  if (!is.null(select)) {
    ask <- FALSE
  }
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (scale == -1 && is.null(ylim)) {
    k <- 0
    if (m > 0)
      for (i in 1:m) if (pd[[i]]$plot.me && pd[[i]]$scale) {
        if (se && length(pd[[i]]$se) > 1) {
          ul <- pd[[i]]$fit + pd[[i]]$se
          ll <- pd[[i]]$fit - pd[[i]]$se
          if (k == 0) {
            ylim <- c(min(ll, na.rm = TRUE), max(ul,
                                                 na.rm = TRUE))
            k <- 1
          }
          else {
            if (min(ll, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(ll, na.rm = TRUE)
            if (max(ul, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(ul, na.rm = TRUE)
          }
        }
        else {
          if (k == 0) {
            ylim <- range(pd[[i]]$fit, na.rm = TRUE)
            k <- 1
          }
          else {
            if (min(pd[[i]]$fit, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(pd[[i]]$fit, na.rm = TRUE)
            if (max(pd[[i]]$fit, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(pd[[i]]$fit, na.rm = TRUE)
          }
        }
      }
  }
  if (m > 0)
    for (i in 1:m) if (pd[[i]]$plot.me && (is.null(select) ||
                                           i == select)) {
      refund:::plot.mgcv.smooth(x$smooth[[i]], P = pd[[i]], partial.resids = NULL,
           rug = rug, se = se, scale = scale, n = n,
           jit = jit,
           xlab = xlab, ylab = ylab, main = main, ylim = ylim,
           xlim = xlim, shade = shade,
           shade.col = shade.col, shift = shift,
           scheme = scheme[i], ...)
    }
   if (pages > 0)
    par(oldpar)
  invisible(pd)

  return(list(beta.fit = fit, beta.sd = se.fit))
}

