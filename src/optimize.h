#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <stdio.h>
#include <float.h>
#include <iostream>
using namespace Rcpp;

double ss2(int j, NumericVector tmpb, arma::mat Q, arma::vec B, int n);
NumericVector wshoot (arma::mat x, arma::vec y, NumericVector init, NumericVector weight, double lambda, int maxiter, double tol, int n);
NumericVector muf(NumericVector b, int M, int d);
NumericVector wshoot1 (arma::mat x, arma::vec y, NumericVector init, int penalty, NumericVector weight, NumericVector wbeta, double lambda, double alpha, int maxiter, double tol, int n);

#endif

