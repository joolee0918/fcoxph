

#' @import survival refund
#' @importFrom pryr modify_call
#' @importFrom Matrix bdiag

#' @export
fcoxph <- function (formula, data, weights, na.action, init, control, knots = NULL, penalty = c("lasso", "MCP", "gBridge"),
                    L2penalty = c("none", "smooth"), theta = switch(L2penalty, none = 0, ridge = NULL, smooth = NULL), gamma = 0.5, lambda = NULL, lambda.min.ratio = 0.001, nlambda =50,
                    alpha = switch(penalty, lasso = 0, SCAD = 3.7, MCP = 3, gBridge=0),
                    sparse = c("none", "local"),  group.multiplier=NULL, tuning.method = c("aic", "bic", "gcv", "all", "cv"), nfolds = 10, foldid = NULL,
                    method = "breslow", parallel = FALSE, ncluster = NULL, x = FALSE, model = FALSE, y = TRUE, ...)
{

  call <- match.call()
  dots <- list(...)
  penalty <- match.arg(penalty)
  allVars <- all.vars(formula)

  tf <- terms.formula(formula, specials = c("s", "fs", "fpc"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]],
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.s <- specials$s - 1
  where.fs <- specials$fs - 1
  where.fp <- specials$fpc - 1

  where.all <- c( where.s, where.fs, where.fp)
  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% where.all))
  }else {
    where.par <- numeric(0)
  }
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(refund:::safeDeparse(responsename), "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call))
    eval.parent(call$data)
  else NULL

  surv <- eval(responsename, envir = evalenv)
  nobs <- nrow(surv)
  if (is.call(responsename)) {
    sapply(as.list(responsename[-1]), function(x) {
      assign(x = deparse(x), value = eval(x, envir = evalenv,
                                          enclos = frmlenv), envir = newfrmlenv)
      invisible(NULL)
    })
  }
  else {
    assign(x = deparse(responsename), value = eval(responsename,
                                                   envir = evalenv, enclos = frmlenv), envir = newfrmlenv)
  }

  newtrmstrings <- attr(tf, "term.labels")

  assign("s", f_override, envir = parent.env(newfrmlenv))
#  assign("te", f_override, envir = parent.env(newfrmlenv))
#  assign("t2", f_override, envir = parent.env(newfrmlenv))

  where.refund <- c(where.fs, where.fp)
   if (length(where.refund)) {
       fterms <- lapply(terms[where.refund], function(x) {
         newx <- match.call(fs, call = x)
         newx$sparse <- sparse
         newx$theta <- theta
         newx$lambda <- lambda
         newx$penalty <- penalty
         newx$tuning.method <- tuning.method

      eval(newx, envir = evalenv)
    })

    newtrmstrings[where.refund] <- sapply(fterms, function(x) {
      x$names
    })

    nm <- smooth <- smooth.terms <- argvals <- theta <- lambda <- vector(mode = "list", length = length(fterms))
    for(i in 1:length(fterms)){
      nm[[i]] <- fterms[[i]]$names
      smooth[[i]] <- fterms[[i]]$sm
      argvals[[i]] <- fterms[[i]]$argvals
    }



    lapply(fterms, function(x) {
        assign(x = x$names, x$X, envir = newfrmlenv)
        invisible(NULL)
    })
  }
  else fterms <- NULL


  where.mgcv <- c(where.s)
  if (length(where.mgcv)) stop("s term has not been implemented yet")
#  if (length(where.mgcv)) {
#    if ("data" %in% names(call))
#      frmlenv <- list2env(eval(call$data), frmlenv)
#    lapply(terms[where.mgcv], function(x) {
#      nms <- if (!is.null(names(x))) {
#        all.vars(x[names(x) == ""])
#      }
#      else all.vars(x)
#      sapply(nms, function(nm) {
#        stopifnot(length(get(nm, envir = frmlenv)) ==
#                    nobs)
#        assign(x = nm, value = get(nm, envir = frmlenv),
#               envir = newfrmlenv)
#
#        invisible(NULL)
#      })
#      invisible(NULL)
#    })
#  }



  varmap <- vector(mode = "list",   length = length(trmstrings))
  if (length(where.par)) {
  for (i in where.par) {
     term.i <- terms[i]
    nms <- all.vars(term.i[[1]])
    varmap[[i]] <- nms
    sapply(nms, function(nm) {
      v <- eval(as.name(nm), envir = evalenv, enclos = frmlenv)
      stopifnot(length(v) == nobs)
      assign(x = nm, value = v, envir = newfrmlenv)
      invisible(NULL)
    })
  }
}
suppressMessages(trace(coxph, at = which(sapply(as.list(body(coxph)),
                                                function(x) any(grepl(x, pattern = "mf[[timetrans$var[i]]]",
                                                                      fixed = TRUE)))) + 1, print = FALSE, tracer = quote({
                                                                        tmp <- attr(mf, "na.action")
                                                                        attr(mf, "na.action") <- NULL
                                                                        mf <- na.action(mf)
                                                                        omit <- attr(mf, "na.action")
                                                                        if (!is.null(omit)) {
                                                                          Y <- Y[-omit, ]
                                                                          strats <- strats[-omit]
                                                                        }
                                                                      })))
on.exit({
  suppressMessages(try(untrace(coxph), silent = TRUE))
})

  newfrml <- formula(paste(newfrml, paste(newtrmstrings, collapse = "+")))
  environment(newfrml) <- newfrmlenv
  fcoxphdata <- refund:::list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), function(x) {
    if (is.numeric(x) | is.logical(x)) {
      mean(x)
    }
    else if (is.factor(x)) {
      names(which.max(table(x)))
    }
    else NA
  }, simplify = FALSE)

  newcall <- refund:::expand.call(fcoxph, call)
  newcall$formula <- newfrml
  newcall$fitter <- NULL
  newcall <- pryr::modify_call(newcall, dots)
  newcall$data <- quote(fcoxphdata)
  newcall$na.action <- na.omit_pcox


  if(length(trmstrings)==0) {
    newcall$eps <-  newcall$knots <- newcall$argvals <- newcall$penalty <- newcall$nfolds <- newcall$foldid <-  newcall$L2penalty <- NULL
    newcall$theta <- newcall$lambda <- newcall$lambda.min <- newcall$nlambda <- newcall$alpha <- NULL
    newcall$sparse <- newcall$group.multiplier <- newcall$tuning.method  <- newcall$parallel <- newcall$ncluster  <- NULL
    newcall[[1]] <- as.symbol("coxph")

    }else{
    newcall$eps <-  newcall$knots <- NULL
    newcall$argvals <- argvals
    newcall$sm <- smooth
    newcall[[1]] <- as.symbol("fcoxph.fit")
}
  res <- eval(newcall)

  termtype <- rep("par", length(terms))
  for (i in 1:length(specials))
    termtype[specials[[i]]-1] <- names(specials)[i]


  fcoxterm <- list(formula = formula,
              responsename = responsename, nobs = nobs,
              termnames = names(terms),
              termtype = termtype, datameans=datameans)

  res$fcoxph <- fcoxterm
  return(res)
  class(fit) <- res$class

}

