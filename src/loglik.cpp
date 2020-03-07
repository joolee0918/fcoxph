#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <float.h>
#include <iostream>

//[[Rcpp::export()]]
using namespace Rcpp;


double loglik_coxfit(NumericVector time,   IntegerVector status,
        NumericMatrix covar,    NumericVector offset, NumericVector weights,
        IntegerVector strata,   int method, NumericVector beta,    int doscale) {


  double tdeath, temp, temp2, zbeta, risk;
  double denom, dtime, deadwt, denom2, newlk;
  double loglik;

  int i, j, person;
  int nused, nvar;
  int ndead, nrisk;
  int method;

  nused = offset.size();
  nvar  = covar.ncol();

  NumericVector beta(nvar), means(nvar), scale(nvar);


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
    if (doscale==1) {  /* and also scale it */
  temp =0;
      for (person=0; person<nused; person++) {
        temp += weights[person] * fabs(covar(person, i));
      }
      if (temp > 0) temp = temp2/temp;   /* scaling */
  else temp=1.0; /* rare case of a constant covariate */
  scale[i] = temp;
  for (person=0; person<nused; person++) {
    covar(person, i) *= temp;
  }
    }
  }

  if (doscale==1) {
    for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */
  }
  else {
    for (i=0; i<nvar; i++) scale[i] = 1.0;
  }

  strata[nused-1] =1;
  loglik =0;

  for (person=nused-1; person>=0; ) {
    if (strata[person] == 1) {
      nrisk =0 ;
      denom = 0;
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
    }
    else {
      ndead++;
      deadwt += weights[person];
      denom2 += risk;
      loglik += weights[person]*zbeta;
    }
    person--;
    if (person>=0 && strata[person]==1) break;
  }

  if (ndead >0) {  /* we need to add to the main terms */
    if (method==0) { /* Breslow */
    denom += denom2;
      loglik -= deadwt* log(denom);
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
        loglik[1] -= wtave* log(denom);
      }
    }
  }
  }


  return(loglik);
}




double loglik_agfit(SEXP surv2,      SEXP covar2,    SEXP strata2,
            SEXP weights2,   SEXP offset2,   SEXP ibeta2,
            SEXP sort12,     SEXP sort22,    SEXP method2,
            SEXP maxiter2,   SEXP  eps2,     SEXP tolerance2,
            SEXP doscale2) {

  int i,j,k, person;
  int indx1, istrat, p, p1;
  int nrisk;
  int nused, nvar;
  int rank, rank2, fail;

  double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *oldbeta;
    double *scale;
    double *a2, **cmat2;
    double *eta;
    double  denom, zbeta, risk;
    double  dtime;
    double  temp, temp2;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    double  tol_chol, eps;
    double  meanwt;
    int deaths;
    double denom2, etasum;
    int *keep;               /* marker for useless obs */

    /* inputs */
    double *start, *tstop, *event;
    double *weights, *offset;
    int *sort1, *sort2, maxiter;
    int *strata, nstrat;
    double method;  /* saving this as double forces some double arithmetic */
    int doscale;

    /* returned objects */
    SEXP imat2, beta2, u2, loglik2;
    double *beta, *u, *loglik;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist;
    static const char *outnames[]={"coef", "u", "imat", "loglik",
                                   "sctest", "flag", "iter", ""};
    int nprotect;  /* number of protect calls I have issued */

    /* get sizes and constants */
    nused = nrows(covar2);
    nvar  = ncols(covar2);
    method= asInteger(method2);
    eps   = asReal(eps2);
    tol_chol = asReal(tolerance2);
    maxiter = asInteger(maxiter2);
    doscale = asInteger(doscale2);
    nstrat = LENGTH(strata2);

    /* input arguments */
    start = REAL(surv2);
    tstop  = start + nused;
    event = tstop + nused;
    weights = REAL(weights2);
    offset = REAL(offset2);
    sort1  = INTEGER(sort12);
    sort2  = INTEGER(sort22);
    strata = INTEGER(strata2);

    /*
    ** scratch space
    **  nvar: a, a2, oldbeta, scale
    **  nvar*nvar: cmat, cmat2
    **  nused:  eta, keep
    */
    eta = (double *) R_alloc(nused + 4*nvar + 2*nvar*nvar, sizeof(double));
    a = eta + nused;
    a2= a + nvar;
    scale  = a2 + nvar;
    oldbeta = scale + nvar;
    keep = (int *) R_alloc(nused, sizeof(int));

    /*
    **  Set up the ragged arrays
    **  covar2 might not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  But check
    */
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
    nprotect =1;
    if (MAYBE_REFERENCED(covar2)) {
      PROTECT(covar2 = duplicate(covar2));
      nprotect++;
    }
    covar= dmatrix(REAL(covar2), nused, nvar);
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    cmat = dmatrix(oldbeta+ nvar,   nvar, nvar);
    cmat2= dmatrix(oldbeta+ nvar + nvar*nvar, nvar, nvar);

    /*
    ** create the output structures
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    nprotect++;
    beta2 = SET_VECTOR_ELT(rlist, 0, duplicate(ibeta2));
    beta  = REAL(beta2);
    u2 =    SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, nvar));
    u = REAL(u2);

    SET_VECTOR_ELT(rlist, 2, imat2);
    loglik2 = SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, 2));
    loglik  = REAL(loglik2);

    sctest2 = SET_VECTOR_ELT(rlist, 4, allocVector(REALSXP, 1));
    sctest =  REAL(sctest2);
    flag2  =  SET_VECTOR_ELT(rlist, 5, allocVector(INTSXP, 3));
    flag   =  INTEGER(flag2);
    for (i=0; i<3; i++) flag[i]=0;

    iter2  =  SET_VECTOR_ELT(rlist, 6, allocVector(INTSXP, 1));
    iter = INTEGER(iter2);

    /*
    ** Subtract the mean from each covar, as this makes the variance
    **  computation much more stable.  The mean is taken per stratum,
    **  the scaling is overall.
    */
    if (nvar==1) doscale =0;  /* scaling has no impact, so skip it */
    for (i=0; i<nvar; i++) {
      person=0;
      for (istrat=0; istrat<nstrat; istrat++) {
        temp=0;
        temp2 =0;
        for (k=person; k<strata[istrat]; k++) {
          j = sort2[k];
          temp += weights[j] * covar[i][j];
          temp2 += weights[j];
        }
        temp /= temp2;   /* mean for this covariate, this strata */
    for (; person< strata[istrat]; person++) {
      j = sort2[person];
      covar[i][j] -=temp;
    }
      }
      if (doscale ==1) { /* also scale the regression */
    /* this cannot be done per stratum */
    temp =0;
        temp2 =0;
        for (person=0; person<nused; person++) {
          temp += weights[person] * fabs(covar[i][person]);
          temp2 += weights[person];
        }
        if (temp >0) temp = temp2/temp;  /* 1/scale */
    else temp = 1.0;  /* rare case of a constant covariate */
    scale[i] = temp;
    for (person=0; person<nused; person++) {
      covar[i][person] *= temp;
    }
      }
    }

    if (doscale ==1) {
      for (i=0; i<nvar; i++) beta[i] /= scale[i]; /* rescale initial betas */
    }
    else {for (i=0; i<nvar; i++) scale[i] = 1.0;}
    /* keep[] will have number of event times for which this subject is at risk */
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
      if (event[p]) {
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
    /* First iteration, which has different ending criteria */
    for (person=0; person<nused; person++) {
      zbeta = 0;      /* form the term beta*z   (vector mult) */
    for (i=0; i<nvar; i++)
      zbeta += beta[i]*covar[i][person];
    eta[person] = zbeta + offset[person];
    }

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'dtime' is a scratch variable holding the time of current interest
    **  'indx1' walks through the start times.
    */
    newlk =0;
    for (i=0; i<nvar; i++) {
      u[i] =0;
      for (j=0; j<nvar; j++) imat[i][j] =0;
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
      for (j=0; j<nvar; j++) cmat[i][j] =0;
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
            for (j=0; j<nvar; j++) cmat[i][j] =0;
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
                for (j=0; j<=i; j++) cmat[i][j] =0;
              }
            }
            else {
              etasum -= eta[p1];
              risk = exp(eta[p1]) * weights[p1];
              denom -= risk;
              for (i=0; i<nvar; i++) {
                a[i] -= risk*covar[i][p1];
                for (j=0; j<=i; j++)
                  cmat[i][j] -= risk*covar[i][p1]*covar[j][p1];
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
              flag[1]++;  /* a count, for debugging/profiling purposes */
            temp = etasum/nrisk;
            for (i=0; i<nused; i++) eta[i] -= temp;
            temp = exp(-temp);
            denom *= temp;
            for (i=0; i<nvar; i++) {
              a[i] *= temp;
              for (j=0; j<nvar; j++) {
                cmat[i][j]*= temp;
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
                cmat2[i][j]=0;
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
              newlk += weights[p]* eta[p];
              for (i=0; i<nvar; i++) {
                u[i] += weights[p] * covar[i][p];
                a2[i]+= risk*covar[i][p];
                for (j=0; j<=i; j++)
                  cmat2[i][j] += risk*covar[i][p]*covar[j][p];
              }
            }
            else if (keep[p] >0) {
              nrisk++;
              etasum += eta[p];
              denom += risk;
              for (i=0; i<nvar; i++) {
                a[i] += risk*covar[i][p];
                for (j=0; j<=i; j++)
                  cmat[i][j] += risk*covar[i][p]*covar[j][p];
              }
            }
            }
            /*
            ** Add results into u and imat for all events at this time point
            */
            if (method==0 || deaths ==1) { /*Breslow */
            denom += denom2;
              newlk -= meanwt*log(denom);  /* sum of death weights*/
            for (i=0; i<nvar; i++) {
              a[i] += a2[i];
              temp = a[i]/denom;   /*mean covariate at this time */
            u[i] -= meanwt*temp;
            for (j=0; j<=i; j++) {
              cmat[i][j] += cmat2[i][j];
              imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
            }
            }
            }
            else {
              meanwt /= deaths;
              for (k=0; k<deaths; k++) {
                denom += denom2/deaths;
                newlk -= meanwt*log(denom);
                for (i=0; i<nvar; i++) {
                  a[i] += a2[i]/deaths;
                  temp = a[i]/denom;
                  u[i] -= meanwt*temp;
                  for (j=0; j<=i; j++) {
                    cmat[i][j] += cmat2[i][j]/deaths;
                    imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
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
              flag[1]++;  /* a count, for debugging/profiling purposes */
            temp = etasum/nrisk;
            for (i=0; i<nused; i++) eta[i] -= temp;
            temp = exp(-temp);
            denom *= temp;
            for (i=0; i<nvar; i++) {
              a[i] *= temp;
              for (j=0; j<nvar; j++) {
                cmat[i][j]*= temp;
              }
            }
            etasum =0;
            }
          }
    }   /* end  of accumulation loop */
            loglik[0] = newlk;   /* save the loglik for iteration zero  */



