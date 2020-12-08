//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <float.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

double ss2(int j, NumericVector tmpb, arma::mat Q, arma::vec B, int n)
{
  double a, a1, s;
  a1 = 0;
  for(int i=0; i<Q.n_rows; i++) a1 += tmpb[i]*Q(i, j);
  a = a1 - tmpb[j]*Q(j,j);
  s = (a-B(j))/n;

  return(s);
}

//[[Rcpp::export()]]
NumericVector wshoot1 (arma::mat x, arma::vec y, NumericVector init, int pen, NumericVector weight, NumericVector wbeta, double lambda, double alpha, int maxiter, double tol, int n)
{
  int nrow = x.n_rows;
  int ncol = x.n_cols;

  arma::vec B(nrow);
  int i, j;
  double s, eps;
  NumericVector oldbeta(ncol), tmpbeta(ncol), lams(ncol);

  arma::mat Q = x.t()*x;

  for(j=0; j<ncol; j++){
    B(j)=0;
    for(i=0; i<nrow; i++) B(j) += x(i, j)*y(i);
  }


  for(i=0; i<ncol; i++) {
    lams[i] = lambda*weight[i]*wbeta[i];
    oldbeta[i] = init[i];
    tmpbeta[i] = oldbeta[i];
  }

  for (i=0; i<maxiter; i++){
    for (j=0; j<ncol; j++){
      s = ss2(j,tmpbeta,Q,B,n);
      if(fabs(s) <= lams[j]) tmpbeta[j] = 0.0;
      else {
        if(pen==2){
          if(1/fabs(tmpbeta[j]) > 1/(lams[j]*alpha)) tmpbeta[j] = (-s)/(Q(j,j)/n + lams[j]/fabs(tmpbeta[j]) - 1/alpha); // MCP
          else tmpbeta[j] = (-s)/(Q(j,j)/n);
        }else tmpbeta[j] = (-s)/(Q(j,j)/n + lams[j]/fabs(tmpbeta[j])); //LASSO + gBridge
      }
    }
    eps = max(abs(tmpbeta-oldbeta));
    for(i=0; i<ncol; i++) oldbeta[i] =  tmpbeta[i];
    if (eps <=tol) break;

    i = i+1;
  }
  return(tmpbeta);
}


NumericVector wshoot (arma::mat x, arma::vec y, NumericVector init, NumericVector weight, double lambda, int maxiter, double tol, int n)
{
  int nrow = x.n_rows;
  int ncol = x.n_cols;

  arma::vec B(nrow);
  int i, j;
  double s, eps;
  NumericVector oldbeta(ncol), tmpbeta(ncol), lams(ncol);

  arma::mat Q = x.t()*x;

  for(j=0; j<ncol; j++){
    B(j)=0;
    for(i=0; i<nrow; i++) B(j) += x(i, j)*y(i);
  }


  for(i=0; i<ncol; i++) {
    lams[i] = lambda*weight[i];
    oldbeta[i] = init[i];
    tmpbeta[i] = oldbeta[i];
  }

  for (i=0; i<maxiter; i++){
    for (j=0; j<ncol; j++){
      s = ss2(j,tmpbeta,Q,B,n);
      if (s > lams[j])
        tmpbeta[j] = (lams[j]-s)/(Q(j,j)/n);
      else if (s < (-lams[j]) )
        tmpbeta[j] = (-lams[j]-s)/(Q(j,j)/n);
      else
        tmpbeta[j] = 0.0;
    }

    eps = max(abs(tmpbeta-oldbeta));
    for(i=0; i<ncol; i++) oldbeta[i] =  tmpbeta[i];
    if (eps <=tol) break;

    i = i+1;
  }
  return(tmpbeta);
}


NumericVector muf(NumericVector b, int M, int d)
{
  double tmpb;
  NumericVector mu(M+1);

  for(int j=0; j<(M+1); j++){
    tmpb = 0;
    for(int k=0; k<d; k++) tmpb += fabs(b[j+k]);

    mu[j] = tmpb;
  }

  return(mu);
}

