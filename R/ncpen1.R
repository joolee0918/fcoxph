ncpen1 = function(y.vec,x.mat,
                 family=c("gaussian","linear","binomial","logit","poisson","multinomial","cox"),
                 penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                 x.standardize=TRUE,intercept=TRUE,
                 lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                 df.max=50,cf.max=100,proj.min=10,add.max=10,niter.max=30,qiter.max=10,aiter.max=100,
                 b.eps=1e-7,k.eps=1e-4,c.eps=1e-6,cut=TRUE,local=FALSE,local.initial=NULL){
  if(family[1] == "linear") {
    family[1] = "gaussian";
  } else if(family[1] == "logit") {
    family[1] = "binomial";
  }

  # Data type conversion ------------------------
  if(!is.vector(y.vec)) {
    y.vec = as.vector(y.vec);
  }

  if(!is.matrix(x.mat)) {
    x.mat = as.matrix(x.mat);
  }
  # ---------------------------------------------
  family = match.arg(family); penalty = match.arg(penalty)
  tun = control.ncpen1(y.vec,x.mat,family,penalty,x.standardize,intercept,
                      lambda,n.lambda,r.lambda,w.lambda,gamma,tau,alpha,aiter.max,b.eps)
  if(local==TRUE){
    if(is.null(local.initial)){ stop(" supply 'local.initial' \n") }
    warning("       'local==TRUE' option may take a long time \n")
  } else { local.initial = rep(0,length(tun$w.lambda)) }
  fit = ncpen:::native_cpp_ncpen_fun_(tun$y.vec,tun$x.mat,
                              tun$w.lambda,tun$lambda,tun$gamma,tun$tau,tun$alpha,
                              df.max,niter.max,qiter.max,aiter.max,
                              b.eps,k.eps,proj.min,cut,c.eps,add.max,family,penalty,local,local.initial,cf.max)
  if(x.standardize==TRUE){ fit$beta = fit$beta/tun$std  }
  # 201800906 -------------------------
  # coefficient names
  coef.names = NULL;
  n.coef = ncol(x.mat);
  if(family == "cox") {
    n.coef = n.coef - 1;
  }

  if(!is.null(colnames(x.mat))) { # if matrix has colnames
    if(intercept == TRUE && family != "cox") {
      coef.names = c("intercept", colnames(x.mat)[1:n.coef]);
    } else {
      coef.names = colnames(x.mat)[1:n.coef];
    }
  } else {
    if(intercept == TRUE && family != "cox") {
      coef.names = c("intercept", paste("x", 1:n.coef, sep = ""));
    } else {
      coef.names = paste("x", 1:n.coef, sep = "");
    }
  }

  if(family == "multinomial") {
    coef.names = rep(coef.names, nrow(fit$beta)/length(coef.names));
  }
  rownames(fit$beta) = coef.names;
  # -----------------------------------
  ret = list(y.vec=tun$y.vec,x.mat=tun$x.mat,
             family=family,penalty=penalty,
             x.standardize=x.standardize,intercept=tun$intercept,std=tun$std,
             lambda=drop(fit$lambda),w.lambda=tun$w.lambda,gamma=tun$gamma,tau=tun$tau,alpha=tun$alpha,
             local.initial=local.initial,
             beta=fit$beta,df=drop(fit$df))

  class(ret) = "ncpen";
  return(ret);
}



control.ncpen1 = function(y.vec,x.mat,
                         family=c("gaussian","binomial","poisson","multinomial","cox"),
                         penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                         x.standardize=TRUE,intercept=TRUE,
                         lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                         aiter.max=1e+2,b.eps = 1e-7){
  family = match.arg(family); penalty = match.arg(penalty);
  if(is.null(alpha)){ alpha = switch(penalty,scad=1,mcp=1,tlp=1,lasso=1,classo=1,ridge=0,sridge=1,mbridge=1,mlog=1)
  } else {
    if(penalty=="ridge"){ if(alpha!=0){ warning("alpha is set to 0 \n") }; alpha = 0;
    } else if(penalty=="sridge"){ if(alpha!=1){ warning("alpha is set to 1 \n") }; alpha = 1;
    } else { if(alpha<0){ warning("alpha is set to 0 \n"); alpha = 0 }; if(alpha>1){ warning("alpha is set to 1 \n"); alpha = 1 } }
  }
  if(is.null(tau)){ tau = switch(penalty,scad=3.7,mcp=2.1,tlp=0.001,lasso=2,classo=2.1,ridge=2,sridge=2.1,mbridge=0.001,mlog=0.001)
  } else {
    if(penalty=="scad"){ if(tau<=2){ tau = 2.001; warning("tau is set to ",tau,"\n") }
    } else if(penalty=="mcp"){ if(tau<=1){ tau = 1.001; warning("tau is set to ",tau,"\n") }
    } else if((penalty=="classo")|(penalty=="sridge")){ if(tau<=1){ tau = 1.001; warning("tau is set to ",tau,"\n") }
    } else { if(tau<=0){ tau = 0.001; warning("tau is set to ",tau,"\n") } }
  }
  n = dim(x.mat)[1]; p = dim(x.mat)[2]; if(family=="cox"){ p = p-1 }
  if(is.null(w.lambda)){ w.lambda = rep(1,p)
  } else {
    if(length(w.lambda)!=p){ stop("the number of elements in w.lambda should be the number of input variables")
    } else if(sum(w.lambda==0)>=n){ stop("the number of zero elements in w.lambda should be less than the sample size")
    } else if(min(w.lambda)<0){ stop("elements in w.lambda should be non-negative")
    } else { w.lambda = (p-sum(w.lambda==0))*w.lambda/sum(w.lambda) }
  }
  if(x.standardize==FALSE){ std = rep(1,p)
  } else {
    if(family=="cox"){ std = sqrt(colSums(x.mat[,-(p+1)]^2)/n); std[std==0] = 1; x.mat[,-(p+1)] = sweep(x.mat[,-(p+1)],2,std,"/")
    } else { std = sqrt(colSums(x.mat^2)/n); std[std==0] = 1; x.mat = sweep(x.mat,2,std,"/") }
  }
  if(family=="cox") intercept = FALSE;
  if(intercept==TRUE){ x.mat = cbind(1,x.mat); std = c(1,std); w.lambda = c(0,w.lambda); p = p+1; }
  if(family=="multinomial"){
    k = max(y.vec)
    if(length(unique(y.vec))!=k){ stop("label must be denoted by 1,2,...: ") }
    x.mat = kronecker(diag(k-1),x.mat); w.lambda = rep(w.lambda,k-1); p = p*(k-1); std = rep(std,k-1);
  }
  if(is.null(r.lambda)){
    if((family=="cox")|(family=="multinomial")){
      r.lambda = ifelse(n<p,0.1,0.01);
    } else {
      r.lambda = ifelse(n<p,0.01,0.001);
    }
  } else {
    if((r.lambda<=0)|(r.lambda>=1)) stop("r.lambda should be between 0 and 1");
  }
  if(is.null(n.lambda)){ n.lambda = 100;
  } else { if((n.lambda<1e+1)|(n.lambda>1e+3)) stop("n.lambda should be between 10 and 1000") }
  if(sum(w.lambda==0)==0){ b.vec = rep(0,p);
  } else { b.vec = rep(0,p);
  if(family!="cox"){ a.set = c(1:p)[w.lambda==0]; } else { a.set = c(c(1:p)[w.lambda==0],p+1); }
  ax.mat = x.mat[,a.set,drop=F]; b.vec[a.set] = drop(native_cpp_nr_fun_(family,y.vec,ax.mat,aiter.max,b.eps));
  }
  g.vec = abs(native_cpp_obj_grad_fun_(family,y.vec,x.mat,b.vec))
  if(alpha!=0){
    lam.max = max(g.vec[w.lambda!=0]/w.lambda[w.lambda!=0])/alpha; lam.max = lam.max+lam.max/10
    lam.vec = exp(seq(log(lam.max),log(lam.max*r.lambda),length.out=n.lambda))
    if(is.null(lambda)){ lambda = lam.vec
    } else {
      if(min(lambda)<=0){ stop("elements in lambda should be positive");
      } else { lambda = sort(lambda,decreasing=TRUE);
      #if(max(lambda)<lam.max){ lambda = c(lam.vec[lam.vec>lambda[1]],lambda); warning("lambda is extended up to ",lam.max,"\n") }
      }
    }
  }
  if(alpha==0) lambda = exp(seq(log(1e-2),log(1e-5),length.out=n.lambda))
  if(is.null(gamma)){
    gamma = switch(penalty,scad=0,mcp=0,tlp=0,lasso=0,classo=min(lambda)/2,ridge=0,sridge=min(lambda)/2,mbridge=0,mlog=0)
  } else {
    if((penalty!="classo")&(penalty!="sridge")){ gamma = 0
    } else {
      if(gamma<0){ gamma = min(lambda)/2; warning("gamma is negative and set to ",gamma,"\n") }
      if(gamma>=max(lambda)){ gamma = max(lambda)/2; warning("gamma is larger than lambda and set to ",gamma,"\n") }
      if(gamma>min(lambda)){ lambda = lambda[lambda>=gamma]; warning("lambda is set larger than gamma=",gamma,"\n") }
    }
  }
  ret = list(y.vec=y.vec,x.mat=x.mat,
             family=family,penalty=penalty,
             x.standardize=x.standardize,intercept=intercept,std=std,
             lambda=lambda,n.lambda=n.lambda,r.lambda=r.lambda,w.lambda=w.lambda,gamma=gamma,tau=tau,alpha=alpha)
  class(ret) = "control.ncpen";
  return(ret)
}
