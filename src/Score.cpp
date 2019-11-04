#include <Rcpp.h>
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

/*
** Compute the martingale residual for a Cox model
**
** Input
**      n       number of subjects
**      method  will be ==1 for the Efron method
**      time    vector of times
**      status  vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      strata  is =1 for the last obs of a strata
**      mark    carried forward from the coxfit routine
**
** Output
**      expected the expected number of events for the subject
**
** The martingale residual is more of a nuisance for the Efron method
**
*/
#include <stdio.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

//[[Rcpp::export()]]
Rcpp::List Score(int  n,     int  method,    NumericVector time,
                 IntegerVector status, IntegerVector strata,   NumericVector score,
                 NumericVector wt)
{
  int i,j;
  int lastone;

  double deaths, denom=0, e_denom=0, denom2=0;
  double hazard;
  double temp, wtsum;
  double downwt;

  Rcpp::NumericVector expect(n);
  Rcpp::NumericVector infor(n);

  strata[n-1] =1;  /* Failsafe */

/* Pass 1-- store the risk denominator in 'expect' */
for (i= n -1; i>=0; i--) {
  if (strata[i]==1) denom =0;
  denom += score[i]*wt[i];
  if (i==0 || strata[i-1]==1 ||  time[i-1]!=time[i])
    expect[i] = denom;
  else    {
    expect[i] =0;
    infor[i] = 0;
  }
}

/* Pass 2-- now do the work */
deaths=0;
wtsum =0;
e_denom=0;
hazard =0;
denom2=0;
lastone = 0;
for (i= 0; i<n; i++) {
  if (expect[i]!=0) denom = expect[i];
  expect[i] = status[i]*wt[i];
  deaths += status[i];
  wtsum += status[i]*wt[i];
  e_denom += score[i]*status[i] *wt[i];
  if (strata[i]==1 ||  time[i+1]!=time[i]) {
    /*last subject of a set of tied times */
    if (deaths<2 || method==0) {
      hazard += wtsum/denom;
      denom2 += wtsum/pow(denom, 2);
      for (j=lastone; j<=i; j++) {
        expect[j] -= score[j]*wt[j]*hazard;
        infor[j] += score[j]*wt[j]*hazard - pow(score[j]*wt[j],2)*denom2;
      }
    }
    else {
      temp = hazard;
      wtsum /=deaths;
      for (j=0; j<deaths; j++) {
        downwt = j /deaths;
        hazard +=  wtsum/(denom - e_denom* downwt);
        temp   +=  wtsum*(1-downwt)/(denom - e_denom* downwt);
      }
      for (j=lastone; j<=i; j++) {
        if (status[j]==0) expect[j] = -score[j]*hazard;
        else  expect[j] -=  score[j]* temp;
      }
    }
    lastone =i +1;
    deaths =0;
    wtsum=0;
    e_denom =0;
  }
  if (strata[i]==1) hazard =0;
}

for (j=lastone; j<n; j++)  {
  expect[j] -= score[j] * wt[j]*hazard;
  infor[j] += score[j]*wt[j]*hazard - pow(score[j]*wt[j],2)*denom2;
}

Rcpp::List res = List::create(Named("score") = expect, _["infor"] = infor);
return(res);
}
