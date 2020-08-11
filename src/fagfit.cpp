//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <float.h>
#include <iostream>
using namespace Rcpp;

#include "optimize.h"

//[[Rcpp::export()]]
List fagfit_cpp(NumericMatrix surv2,
                 NumericMatrix covar2,   IntegerVector strata2, NumericVector weights, NumericVector offset,
                 NumericVector ibeta, IntegerVector sort1, IntegerVector sort2, int method, int maxiter, double eps,
                 NumericMatrix H, NumericMatrix Dstar, NumericMatrix G,  NumericVector wbeta, NumericVector lambda, double alpha,
                 double gamma, int M, int d, int n_npvar, int Dnrow, int pen,  IntegerVector penalty_where, Function f, Function df_f) {


  double temp, temp2, zbeta, risk;
  double denom, dtime, deadwt, denom2, wtave, meanwt, etasum;
  double loglik;
  double error;

  int i, j, k, person, ilam, iter;
  int indx1, istrat, p, p1;
  int nused, nvar;
  int ndead, nrisk, deaths;
  int rank, rank2;

  int nlambda = lambda.size();
  nused = covar2.nrow();
  nvar  = covar2.ncol();
  int n_pvar = nvar - n_npvar;

  NumericVector eta(nused), keep(nused);
  NumericVector beta(nvar), newbeta(nvar), means(nvar), pbeta(n_npvar);
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


  NumericVector start = surv2(_, 0);
  NumericVector tstop = surv2(_, 1);
  NumericVector event = surv2(_, 2);

  int nstrat = strata.size();
  /*
  ** Subtract the mean from each covar, as this makes the regression
  **  much more stable.
  */
  for (i=0; i<nvar; i++) {
    person=0;
    for (istrat=0; istrat<nstrat; istrat++) {
      temp=0;
      temp2 =0;
      for (k=person; k<strata[istrat]; k++) {
        j = sort2[k];
        temp += weights[j] * covar(j,i);
        temp2 += weights[j];
      }
      temp /= temp2;   /* mean for this covariate, this strata */
  for (; person< strata[istrat]; person++) {
    j = sort2[person];
    covar(j, i) -=temp;
  }
    }
     }


  indx1 =0;
  deaths =0;
  istrat =0;

  for (person=0; person < nused;) {
    if (person == strata[istrat]) {  /* first subject in a new stratum */
  /* finish the work for the prior stratum */
  for (; indx1<person; indx1++) {
    p1 = sort1[indx1];
    keep[p1] += deaths;
  }
  deaths=0;
      istrat++;
    }
    p = sort2[person];
    keep[p] = - deaths;
    if (event[p]==1) {
      dtime = tstop[p];
      for(person =person+1; person < strata[istrat]; person++) {
        /* walk forward over any tied times */
        p = sort2[person];
        if (tstop[p] != dtime) break;
        keep[p] = -deaths;
      }
      for (; indx1<person; indx1++) {
        p1 = sort1[indx1];
        if (start[p1] < dtime) break;
        keep[p1] += deaths;
      }
      deaths++;
    } else person++;
  }
  for (; indx1<nused; indx1++) {  /* finish up the last strata */
        p1 = sort1[indx1];
    keep[p1] += deaths;
  }


  for(ilam=0; ilam<nlambda; ilam++){

    error = 1;
    for(i=0; i<nvar; i++) beta[i] = ibeta[i];

    for(iter=1; iter <= maxiter; iter++){
      loglik =0;

      /* First iteration, which has different ending criteria */
      for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
      for (i=0; i<nvar; i++)
        zbeta += beta[i]*covar(person, i);
      eta[person] = zbeta + offset[person];
      }

      /*
      **  'person' walks through the the data from 1 to n,
      **     sort1[0] points to the largest stop time, sort1[1] the next, ...
      **  'dtime' is a scratch variable holding the time of current interest
      **  'indx1' walks through the start times.
      */
      for (i=0; i<nvar; i++) {
        u[i] =0;
        for (j=0; j<nvar; j++) imat(i, j) =0;
      }
      person =0;
      indx1 =0;
      istrat =0;

      /* this next set is rezeroed at the start of each stratum */
      denom=0;
      nrisk=0;
      etasum =0;
      for (i=0; i<nvar; i++) {
        a[i] =0;
        for (j=0; j<nvar; j++) cmat(i, j) =0;
      }
      /* end of the per-stratum set */

      while (person < nused) {
        /* find the next death time */
        for (k=person; k< nused; k++) {
          if (k == strata[istrat]) {
            /* hit a new stratum; reset temporary sums */
            istrat++;
            denom = 0;
            nrisk = 0;
            etasum =0;
            for (i=0; i<nvar; i++) {
              a[i] =0;
              for (j=0; j<nvar; j++) cmat(i, j) =0;
            }
            person =k;  /* skip to end of stratum */
            indx1  =k;
          }
          p = sort2[k];
          if (event[p] == 1) {
            dtime = tstop[p];
            break;
          }
        }
        if (k == nused) person =k;  /* no more deaths to be processed */
            else {
              /* remove any subjects no longer at risk */
              /*
              ** subtract out the subjects whose start time is to the right
              ** If everyone is removed reset the totals to zero.  (This happens when
              ** the survSplit function is used, so it is worth checking).
              */
              for (; indx1<strata[istrat]; indx1++) {
                p1 = sort1[indx1];
                if (start[p1] < dtime) break;
                if (keep[p1] == 0) continue;  /* skip any never-at-risk rows */
              nrisk--;
              if (nrisk ==0) {
                etasum =0;
                denom =0;
                for (i=0; i<nvar; i++) {
                  a[i] =0;
                  for (j=0; j<=i; j++) cmat(i, j) =0;
                }
              }
              else {
                etasum -= eta[p1];
                risk = exp(eta[p1]) * weights[p1];
                denom -= risk;
                for (i=0; i<nvar; i++) {
                  a[i] -= risk*covar(p1, i);
                  for (j=0; j<=i; j++)
                    cmat(i, j) -= risk*covar(p1, i)*covar(p1, j);
                }
              }
              /*
              ** We must avoid overflow in the exp function (~750 on Intel)
              ** and want to act well before that, but not take action very often.
              ** One of the case-cohort papers suggests an offset of -100 meaning
              ** that etas of 50-100 can occur in "ok" data, so make it larger
              ** than this.
              ** If the range of eta is more then log(1e16) = 37 then the data is
              **  hopeless: some observations will have effectively 0 weight.  Keeping
              **  the mean sensible suffices to keep the max in check for all other
              *   data sets.
              */
              if (fabs(etasum/nrisk) > 200) {
              temp = etasum/nrisk;
              for (i=0; i<nused; i++) eta[i] -= temp;
              temp = exp(-temp);
              denom *= temp;
              for (i=0; i<nvar; i++) {
                a[i] *= temp;
                for (j=0; j<nvar; j++) {
                  cmat(i, j)*= temp;
                }
              }
              etasum =0;
              }
              }

              /*
              ** add any new subjects who are at risk
              ** denom2, a2, cmat2, meanwt and deaths count only the deaths
              */
              denom2= 0;
              meanwt =0;
              deaths=0;
              for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                  cmat2(i, j)=0;
                }
              }

              for (; person<strata[istrat]; person++) {
                p = sort2[person];
                if (tstop[p] < dtime) break; /* no more to add */
              risk = exp(eta[p]) * weights[p];

              if (event[p] ==1 ){
                nrisk++;
                etasum += eta[p];
                deaths++;
                denom2 += risk*event[p];
                meanwt += weights[p];
                loglik += weights[p]* eta[p];
                for (i=0; i<nvar; i++) {
                  u[i] += weights[p] * covar(p, i);
                  a2[i]+= risk*covar(p, i);
                  for (j=0; j<=i; j++)
                    cmat2(i, j) += risk*covar(p, i)*covar(p, j);
                }
              }
              else if (keep[p] >0) {
                nrisk++;
                etasum += eta[p];
                denom += risk;
                for (i=0; i<nvar; i++) {
                  a[i] += risk*covar(p, i);
                  for (j=0; j<=i; j++)
                    cmat(i, j) += risk*covar(p, i)*covar(p, j);
                }
              }
              }
              /*
              ** Add results into u and imat for all events at this time point
              */
              if (method==0 || deaths ==1) { /*Breslow */
              denom += denom2;
                loglik -= meanwt*log(denom);  /* sum of death weights*/
              for (i=0; i<nvar; i++) {
                a[i] += a2[i];
                temp = a[i]/denom;   /*mean covariate at this time */
              u[i] -= meanwt*temp;
              for (j=0; j<=i; j++) {
                cmat(i, j) += cmat2(i, j);
                imat(j, i) += meanwt*((cmat(i, j)- temp*a[j])/denom);
              }
              }
              }
              else {
                meanwt /= deaths;
                for (k=0; k<deaths; k++) {
                  denom += denom2/deaths;
                  loglik -= meanwt*log(denom);
                  for (i=0; i<nvar; i++) {
                    a[i] += a2[i]/deaths;
                    temp = a[i]/denom;
                    u[i] -= meanwt*temp;
                    for (j=0; j<=i; j++) {
                      cmat(i, j) += cmat2(i, j)/deaths;
                      imat(j, i) += meanwt*((cmat(i, j)- temp*a[j])/denom);
                    }
                  }
                }
              }
              /*
              ** We must avoid overflow in the exp function (~750 on Intel)
              ** and want to act well before that, but not take action very often.
              ** One of the case-cohort papers suggests an offset of -100 meaning
              ** that etas of 50-100 can occur in "ok" data, so make it larger
              ** than this.
              ** If the range of eta is more then log(1e16) = 37 then the data is
              **  hopeless: some observations will have effectively 0 weight.  Keeping
              **  the mean sensible suffices to keep the max in check for all other
              *   data sets.
              */
              if (fabs(etasum/nrisk) > 200) {
               temp = etasum/nrisk;
              for (i=0; i<nused; i++) eta[i] -= temp;
              temp = exp(-temp);
              denom *= temp;
              for (i=0; i<nvar; i++) {
                a[i] *= temp;
                for (j=0; j<nvar; j++) {
                  cmat(i, j)*= temp;
                }
              }
              etasum =0;
              }
            }
      }   /* end  of accumulation loop */


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
        mu= muf(pbeta, M, d);

        for(j=0; j<(M+1); j++){
          if(mu[j] == 0) theta[j] = 1e10;
          else theta[j] = gamma*pow(mu[j], gamma-1);
        }

        for(k=0; k<n_npvar; k++){
          for(j=0; j<(M+1); j++) penalty_f[penalty_where[k]-1] += theta[j]*H(j, k);
        }
      } else{
        for(k=0; k<n_npvar; k++) penalty_f[penalty_where[k]-1] = 1;
      }


     if(Dnrow==0){
       newbeta = wshoot1(V0, yy, beta, pen, penalty_f, wbeta, lambda[ilam], alpha, maxiter, eps, nused);
     }else{
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

      newbeta = wshoot1(Vstar, Ystar, beta, pen, penalty_f, wbeta, lambda[ilam], alpha, maxiter, eps, nused);
     }
      error = max(abs(newbeta - beta));

      for(i=0; i<nvar; i++) beta[i] = newbeta[i];
    }

    for(i=0; i<nvar; i++) fit_beta(i,ilam) = newbeta[i];

    for(i=0; i<n_npvar; i++){
      pbeta[i] = newbeta(penalty_where[i]-1);
    }

    penalty_f.fill(0);

    if(pen == 3) {
      mu= muf(pbeta, M, d);

      for(j=0; j<(M+1); j++){
        if(mu[j] == 0) theta[j] = 1e10;
        else theta[j] = gamma*pow(mu[j], gamma-1);
      }

      for(k=0; k<n_npvar; k++){
        for(j=0; j<(M+1); j++) penalty_f[penalty_where[k]-1] += theta[j]*H(j, k);
      }
    } else{
      for(k=0; k<n_npvar; k++) penalty_f[penalty_where[k]-1] = 1;
    }

    dA.fill(0);
    for(i=0; i<nvar; i++) if(newbeta[i]!=0) {
      if(pen==2) dA[i] = nused*penalty_f[i]/fabs(newbeta[i]*wbeta[i]) - nused*1/alpha;
      else dA[i] = nused*penalty_f[i]/fabs(newbeta[i]*wbeta[i]);
    }

    List df_var = df_f(newbeta, penalty_where, dA, G, imat, Dnrow);
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
