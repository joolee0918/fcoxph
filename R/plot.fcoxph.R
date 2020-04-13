
#' @export
plot.fcoxph <- function (x, rug = TRUE, se = TRUE, pages = 0, select = NULL,
                          scale = -1, n = 100,
                          jit = FALSE, xlab = NULL, ylab = NULL,
                          main = NULL, ylim = NULL, xlim = NULL,
                          shade = FALSE, shade.col = "gray80", shift = 0, trans = I, scheme =0, ...) {
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
      P <- plot.coef.smooth(x$smooth[[i]], P = NULL,
                se = se,
                scale = scale, n = n, jit = jit, xlab = xlab, ylab = ylab,
                main = main, label = x$smooth[[i]]$label, ylim = ylim, xlim = xlim,
                shade = shade, shade.col = shade.col, shift=shift, trans = trans, scheme=scheme[i],
                se1.mult = se1.mult, se2.mult = se2.mult, ...)
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
  n.plots <-  0
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
     plot.coef.smooth(x$smooth[[i]], P = pd[[i]],
           rug = rug, se = se, scale = scale, n = n,
           jit = jit,
           xlab = xlab, ylab = ylab, main = main, ylim = ylim,
           xlim = xlim, shade = shade,
           shade.col = shade.col, shift = shift, trans = trans,
           scheme = scheme[i], ...)
    }
   if (pages > 0)
    par(oldpar)
  invisible(pd)

}


plot.coef.smooth <- function (x, P = NULL, label = "", se1.mult = 1,
                              se2.mult = 2, rug = TRUE, se = TRUE,
                              scale = -1, n = 100,
                              jit = FALSE, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL,
                              xlim = NULL,  shade = FALSE, shade.col = "gray80",
                              shift = 0, trans = I, scheme = 0, ...) {
  sp.contour <- function(x, y, z, zse, xlab = "", ylab = "",
                         zlab = "", titleOnly = FALSE, se.plot = TRUE, se.mult = 1,
                         trans = I, shift = 0, ...) {
    gap <- median(zse, na.rm = TRUE)
    zr <- max(trans(z + zse + shift), na.rm = TRUE) - min(trans(z -
                                                                  zse + shift), na.rm = TRUE)
    n <- 10
    while (n > 1 && zr/n < 2.5 * gap) n <- n - 1
    zrange <- c(min(trans(z - zse + shift), na.rm = TRUE),
                max(trans(z + zse + shift), na.rm = TRUE))
    zlev <- pretty(zrange, n)
    yrange <- range(y)
    yr <- yrange[2] - yrange[1]
    xrange <- range(x)
    xr <- xrange[2] - xrange[1]
    ypos <- yrange[2] + yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$type = "n"
    args$xlab <- args$ylab <- ""
    args$axes <- FALSE
    do.call("plot", args)
    cs <- (yr/10)/strheight(zlab)
    if (cs > 1)
      cs <- 1
    tl <- strwidth(zlab)
    if (tl * cs > 3 * xr/10)
      cs <- (3 * xr/10)/tl
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z + shift)
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$z <- substitute(zz)
    if (!"levels" %in% n.args)
      args$levels <- substitute(zlev)
    if (!"lwd" %in% n.args)
      args$lwd <- 2
    if (!"labcex" %in% n.args)
      args$labcex <- cs * 0.65
    if (!"axes" %in% n.args)
      args$axes <- FALSE
    if (!"add" %in% n.args)
      args$add <- TRUE
    do.call("contour", args)
    if (is.null(args$cex.main))
      cm <- 1
    else cm <- args$cex.main
    if (titleOnly)
      title(zlab, cex.main = cm)
    else {
      xpos <- xrange[1] + 3 * xr/10
      xl <- c(xpos, xpos + xr/10)
      yl <- c(ypos, ypos)
      lines(xl, yl, xpd = TRUE, lwd = args$lwd)
      text(xpos + xr/10, ypos, zlab, xpd = TRUE, pos = 4,
           cex = cs * cm, off = 0.5 * cs * cm)
    }
    if (is.null(args$cex.axis))
      cma <- 1
    else cma <- args$cex.axis
    axis(1, cex.axis = cs * cma)
    axis(2, cex.axis = cs * cma)
    box()
    if (is.null(args$cex.lab))
      cma <- 1
    else cma <- args$cex.lab
    mtext(xlab, 1, 2.5, cex = cs * cma)
    mtext(ylab, 2, 2.5, cex = cs * cma)
    if (!"lwd" %in% n.args)
      args$lwd <- 1
    if (!"lty" %in% n.args)
      args$lty <- 2
    if (!"col" %in% n.args)
      args$col <- 2
    if (!"labcex" %in% n.args)
      args$labcex <- cs * 0.5
    zz <- trans(z + zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[1]
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("-", round(se.mult),
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs *
             cm, off = 0.5 * cs * cm)
    }
    if (!"lty" %in% n.args)
      args$lty <- 3
    if (!"col" %in% n.args)
      args$col <- 3
    zz <- trans(z - zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[2] - xr/5
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("+", round(se.mult),
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs *
             cm, off = 0.5 * cs * cm)
    }
  }
  if (is.null(P)) {
    if (!x$plot.me || x$dim > 2)
      return(NULL)

      raw <- x$xind
      if (is.null(xlim))
        xx <- seq(x$beta.basis$rangeval[1], x$beta.basis$rangeval[2], length = n)
      else xx <- seq(xlim[1], xlim[2], length = n)
      beta.basis <- x$beta.basis
      X <- fda::eval.basis(xx, beta.basis)

       if (is.null(xlab))
        xlabel <- x$term
      else xlabel <- xlab
      if (is.null(ylab))
        ylabel <- label
      else ylabel <- ylab
      if (is.null(xlim))
        xlim <- range(xx)
      return(list(X = X, x = xx, scale = TRUE, se = TRUE,
                  raw = raw, xlab = xlabel, ylab = ylabel, main = main,
                  se.mult = se1.mult, xlim = xlim))
  }
  else {
    if (se) {
        if (scheme == 1)
          shade <- TRUE
        ul <- P$fit + P$se
        ll <- P$fit - P$se
        if (scale == 0 && is.null(ylim)) {
          ylimit <- c(min(ll), max(ul))
           }
        if (!is.null(ylim))
          ylimit <- ylim
        if (shade) {
          plot(P$x, trans(P$fit + shift), type = "n",
               xlab = P$xlab, ylim = trans(ylimit + shift),
               xlim = P$xlim, ylab = P$ylab, main = P$main,
               ...)
          polygon(c(P$x, P$x[n:1], P$x[1]), trans(c(ul,
                                                    ll[n:1], ul[1]) + shift), col = shade.col,
                  border = NA)
          lines(P$x, trans(P$fit + shift), ...)
        }
        else {
          plot(P$x, trans(P$fit + shift), type = "l",
               xlab = P$xlab, ylim = trans(ylimit + shift),
               xlim = P$xlim, ylab = P$ylab, main = P$main,
               ...)
          if (is.null(list(...)[["lty"]])) {
            lines(P$x, trans(ul + shift), lty = 2, ...)
            lines(P$x, trans(ll + shift), lty = 2, ...)
          }
          else {
            lines(P$x, trans(ul + shift), ...)
            lines(P$x, trans(ll + shift), ...)
          }
        }
         if (rug) {
          if (jit)
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }

      else {
        warning("no automatic plotting for smooths of more than two variables")
      }
    }
    else {
        if (scale == 0 && is.null(ylim)) {
          ylimit <- range(P$fit)
        }
        if (!is.null(ylim))
          ylimit <- ylim
        plot(P$x, trans(P$fit + shift), type = "l", xlab = P$xlab,
             ylab = P$ylab, ylim = trans(ylimit + shift),
             xlim = P$xlim, main = P$main, ...)
        if (rug) {
          if (jit)
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }
       else {
        warning("no automatic plotting for smooths of more than one variable")
      }
    }
  }
}





