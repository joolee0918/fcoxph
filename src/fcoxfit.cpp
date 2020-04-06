
//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <float.h>
#include <iostream>
using namespace Rcpp;
#include "optimize.h"


//[[Rcpp::export()]]
List fcoxfit_cpp(NumericVector time,   IntegerVector status,
                 NumericMatrix covar2,    NumericVector offset, NumericVector weights,
                 IntegerVector strata2, int maxiter, double eps,
                 NumericMatrix H, NumericMatrix Dstar, NumericMatrix G, int method, NumericVector ibeta,  NumericVector lambda, double alpha,
                 double gamma, int M, int d, int n_npvar,  int Dnrow, int pen, IntegerVector penalty_where, Function f, Function df_f) {


  double tdeath, temp, temp2, zbeta, risk;
  double denom, dtime, deadwt, denom2, wtave;
  double loglik;
  double error;


  int i, j, k, person, ilam, iter;
  int nused, nvar;
  int ndead, nrisk;

  int nlambda = lambda.size();
  nused = offset.size();
  nvar  = covar2.ncol();
  int n_pvar = nvar - n_npvar;

  NumericVector beta(nvar), newbeta(nvar), means(nvar),  pbeta(n_npvar);
  NumericVector penalty_f(nvar);
  NumericVector a(nvar), a2(nvar), u(nvar), u2(nvar);
  NumericMatrix imat(nvar, nvar), cmat(nvar, nvar), cmat2(nvar, nvar), V(nvar, nvar);
  arma::vec yy(nvar);
  arma::vec Ystar(nvar + Dnrow);
  arma::mat Vstar(nvar + Dnrow, nvar);

  NumericVector mu(M+1), theta(M+1);
  NumericMatrix fit_beta(nvar, nlambda);
  NumericVector df(nlambda), logl(nlambda);
  NumericMatrix var(nvar*nvar, nlambda), A(nvar*nvar, nlambda);
  NumericVector dA(nvar);
  NumericMatrix covar = clone(covar2);
  IntegerVector strata = clone(strata2);


  /*
  ** Subtract the mean from each covar, as this makes the regression
  **  much more stable.
  */
  tdeath=0; temp2=0;
  for (i=0; i<nused; i++) {
    temp2 += weights[i];
    tdeath += weights[i] * status[i];
  }
  for (i=0; i<nvar; i++) {
    temp=0;
    for (person=0; person<nused; person++)
      temp += weights[person] * covar(person, i);
    temp /= temp2;
    means[i] = temp;
    for (person=0; person<nused; person++) covar(person, i) -=temp;
  }


  strata[nused-1] =1;



  for(ilam=0; ilam<nlambda; ilam++){

    error = 1;
    for(i=0; i<nvar; i++) beta[i] = ibeta[i];

    for(iter=1; iter <= maxiter; iter++){
      loglik =0;

      for (i=0; i<nvar; i++) {
        u[i] =0;
        for (j=0; j<nvar; j++)
          imat(i, j) =0;
      }

      for (person=nused-1; person>=0; ) {
        if (strata[person] == 1) {
          nrisk =0 ;
          denom = 0;
          for (i=0; i<nvar; i++) {
            a[i] = 0;
            for (j=0; j<nvar; j++) cmat(i, j) = 0;
          }
        }


        dtime = time[person];
        ndead =0; /*number of deaths at this time point */
  deadwt =0;  /* sum of weights for the deaths */
  denom2=0;  /* sum of weighted risks for the deaths */
  while(person >=0 &&time[person]==dtime) {
    /* walk through the this set of tied times */
    nrisk++;
    zbeta = offset[person];    /* form the term beta*z (vector mult) */
    for (i=0; i<nvar; i++)
      zbeta += beta[i]*covar(person, i);
    risk = exp(zbeta) * weights[person];
    if (status[person] ==0) {
      denom += risk;
      for (i=0; i<nvar; i++) {
        a[i] += risk*covar(person, i);
        for (j=0; j<=i; j++)
          cmat(i, j) += risk*covar(person, i)*covar(person, j);
      }
    }
    else {
      ndead++;
      deadwt += weights[person];
      denom2 += risk;
      loglik += weights[person]*zbeta;
      for (i=0; i<nvar; i++) {
        u[i] += weights[person]*covar(person, i);
        a2[i] +=  risk*covar(person, i);
        for (j=0; j<=i; j++)
          cmat2(i, j) += risk*covar(person, i)*covar(person, j);
      }
    }
    person--;
    if (person>=0 && strata[person]==1) break;
  }

  if (ndead >0) {  /* we need to add to the main terms */
    if (method==0) { /* Breslow */
    denom += denom2;
      loglik -= deadwt* log(denom);
      for (i=0; i<nvar; i++) {
        a[i] += a2[i];
        temp2= a[i]/ denom;  /* mean */
    u[i] -=  deadwt* temp2;
    for (j=0; j<=i; j++) {
      cmat(i, j) += cmat2(i, j);
      imat(j, i) += deadwt*(cmat(i, j) - temp2*a[j])/denom;
    }
      }
    }
    else { /* Efron */
    /*
      ** If there are 3 deaths we have 3 terms: in the first the
      **  three deaths are all in, in the second they are 2/3
      **  in the sums, and in the last 1/3 in the sum.  Let k go
      **  1 to ndead: we sequentially add a2/ndead and cmat2/ndead
      **  and efron_wt/ndead to the totals.
      */
    wtave = deadwt/ndead;
      for (k=0; k<ndead; k++) {
        denom += denom2/ndead;
        loglik -= wtave* log(denom);
        for (i=0; i<nvar; i++) {
          a[i] += a2[i]/ndead;
          temp2 = a[i]/denom;
          u[i] -= wtave *temp2;
          for (j=0; j<=i; j++) {
            cmat(i, j) += cmat2(i, j)/ndead;
            imat(j, i) += wtave*(cmat(i, j) - temp2*a[j])/denom;
          }
        }

      }
    }

    for (i=0; i<nvar; i++) {
      a2[i]=0;
      for (j=0; j<nvar; j++) cmat2(i, j)=0;
    }

  }
      }

      for(i=0; i<nvar; i++){
        for (j=0; j<i; j++) {
          imat(i, j) = imat(j, i);
        }
      }
      if( sum(abs(beta))==0) break;
      if(error < eps) break;


      V = f(imat);
      arma::mat V0 = as<arma::mat>(V);

      for(i=0; i<nvar; i++) {
        u2[i]=0;
        for(j=0; j<nvar; j++){
          u2[i] += imat(i, j)*beta[j];
        }
      }



      arma::mat Vt = V0.t();
      arma::vec u0(nvar);

      for(i=0; i<nvar; i++) u0(i) = u2[i] + u[i];
      for(j=0; j<nvar; j++){
        yy = arma::solve(Vt, u0);
      }


      for(i=0; i<n_npvar; i++){
        pbeta[i] = beta(penalty_where[i]-1);
      }

      penalty_f.fill(0);

      if(pen == 3) {
      mu= muf(pbeta, gamma, lambda[ilam], M, d);

      for(j=0; j<(M+1); j++){
        if(mu[j] == 0) theta[j] = 1e10;
        else theta[j] = pow(mu[j], 1-1/gamma);
      }

      for(k=0; k<n_npvar; k++){
        for(j=0; j<(M+1); j++) penalty_f[penalty_where[k]-1] += theta[j]*H(j, k);
      }
      } else{
        for(k=0; k<n_npvar; k++) penalty_f[penalty_where[k]-1] = 1;
      }


      if(Dnrow == 0){
        newbeta = wshoot1(V, Y, beta, pen, penalty_f, lambda[ilam], alpha, maxiter, eps, nused);
      } else{
      Ystar.fill(0);
      for(i=0; i<nvar; i++) Ystar(i) = yy(i);

      for(j=0; j<nvar; j++){
        for(i=0; i<nvar; i++){
          Vstar(i, j) = V(i, j);
        }
        for(i=0; i<Dnrow; i++){
          Vstar(i+nvar, j) = Dstar(i, j);
        }
      }

      newbeta = wshoot1(Vstar, Ystar, beta, pen, penalty_f, lambda[ilam], alpha, maxiter, eps, nused);

      }
      error = max(abs(newbeta - beta));
      for(i=0; i<nvar; i++) beta[i] = newbeta[i];

    }

    for(i=0; i<nvar; i++) fit_beta(i,ilam) = newbeta[i];


    dA.fill(0);
    for(i=0; i<nvar; i++) if(newbeta[i]!=0) {
      if(pen==2) dA[i] = nused*penalty_f[i]/fabs(newbeta[i]) - nused*1/alpha;
      else dA[i] = nused*penalty_f[i]/fabs(newbeta[i]);
    }

    List df_var = df_f(newbeta, penalty_where, dA, G, imat);

    df[ilam] = df_var["df"];
    NumericVector tmpvar = df_var["var"];
    NumericVector tmpA = df_var["A"];
    var(_, ilam) = tmpvar;
    A(_, ilam) = tmpA;
    logl[ilam] = loglik;
  }
  Rcpp::List res = List::create(Named("loglik")= logl,
                                Named("beta") = fit_beta,
                                Named("df")=df, Named("var")=var, Named("A") = A, Named("u") = u);

  return(res);


}
