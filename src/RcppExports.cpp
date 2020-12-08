// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fagfit_init
List fagfit_init(NumericMatrix surv2, NumericMatrix covar2, IntegerVector strata2, NumericVector weights, NumericVector offset, NumericVector ibeta, IntegerVector sort1, IntegerVector sort2, int method, double eps);
RcppExport SEXP _fcoxph_fagfit_init(SEXP surv2SEXP, SEXP covar2SEXP, SEXP strata2SEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP ibetaSEXP, SEXP sort1SEXP, SEXP sort2SEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type surv2(surv2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort1(sort1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort2(sort2SEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(fagfit_init(surv2, covar2, strata2, weights, offset, ibeta, sort1, sort2, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// fagfit_loglik
double fagfit_loglik(NumericMatrix surv2, NumericMatrix covar2, IntegerVector strata2, NumericVector weights, NumericVector offset, NumericVector ibeta, IntegerVector sort1, IntegerVector sort2, int method, double eps);
RcppExport SEXP _fcoxph_fagfit_loglik(SEXP surv2SEXP, SEXP covar2SEXP, SEXP strata2SEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP ibetaSEXP, SEXP sort1SEXP, SEXP sort2SEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type surv2(surv2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort1(sort1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort2(sort2SEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(fagfit_loglik(surv2, covar2, strata2, weights, offset, ibeta, sort1, sort2, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// fag_score
NumericMatrix fag_score(int n, int nvar, NumericMatrix y, NumericMatrix covar2, IntegerVector strata, NumericVector score, NumericVector weights, int method);
RcppExport SEXP _fcoxph_fag_score(SEXP nSEXP, SEXP nvarSEXP, SEXP ySEXP, SEXP covar2SEXP, SEXP strataSEXP, SEXP scoreSEXP, SEXP weightsSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(fag_score(n, nvar, y, covar2, strata, score, weights, method));
    return rcpp_result_gen;
END_RCPP
}
// fcoxfit_init
List fcoxfit_init(NumericVector time, IntegerVector status, NumericMatrix covar2, NumericVector offset, NumericVector weights, IntegerVector strata2, double eps, int method, NumericVector ibeta);
RcppExport SEXP _fcoxph_fcoxfit_init(SEXP timeSEXP, SEXP statusSEXP, SEXP covar2SEXP, SEXP offsetSEXP, SEXP weightsSEXP, SEXP strata2SEXP, SEXP epsSEXP, SEXP methodSEXP, SEXP ibetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fcoxfit_init(time, status, covar2, offset, weights, strata2, eps, method, ibeta));
    return rcpp_result_gen;
END_RCPP
}
// fcoxfit_loglik
double fcoxfit_loglik(NumericVector time, IntegerVector status, NumericMatrix covar2, NumericVector offset, NumericVector weights, IntegerVector strata2, double eps, int method, NumericVector ibeta);
RcppExport SEXP _fcoxph_fcoxfit_loglik(SEXP timeSEXP, SEXP statusSEXP, SEXP covar2SEXP, SEXP offsetSEXP, SEXP weightsSEXP, SEXP strata2SEXP, SEXP epsSEXP, SEXP methodSEXP, SEXP ibetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fcoxfit_loglik(time, status, covar2, offset, weights, strata2, eps, method, ibeta));
    return rcpp_result_gen;
END_RCPP
}
// fcox_score
NumericMatrix fcox_score(int n, int nvar, NumericMatrix y, NumericMatrix covar2, IntegerVector strata, NumericVector score, NumericVector weights, int method);
RcppExport SEXP _fcoxph_fcox_score(SEXP nSEXP, SEXP nvarSEXP, SEXP ySEXP, SEXP covar2SEXP, SEXP strataSEXP, SEXP scoreSEXP, SEXP weightsSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(fcox_score(n, nvar, y, covar2, strata, score, weights, method));
    return rcpp_result_gen;
END_RCPP
}
// fagfit_cpp
List fagfit_cpp(NumericMatrix surv2, NumericMatrix covar2, IntegerVector strata2, NumericVector weights, NumericVector offset, NumericVector ibeta, IntegerVector sort1, IntegerVector sort2, int method, int maxiter, double eps, NumericMatrix H, NumericMatrix Dstar, NumericMatrix G, NumericVector wbeta, NumericVector lambda, double alpha, double gamma, int M, int d, int n_npvar, int Dnrow, int pen, IntegerVector penalty_where, Function f, Function df_f);
RcppExport SEXP _fcoxph_fagfit_cpp(SEXP surv2SEXP, SEXP covar2SEXP, SEXP strata2SEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP ibetaSEXP, SEXP sort1SEXP, SEXP sort2SEXP, SEXP methodSEXP, SEXP maxiterSEXP, SEXP epsSEXP, SEXP HSEXP, SEXP DstarSEXP, SEXP GSEXP, SEXP wbetaSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP MSEXP, SEXP dSEXP, SEXP n_npvarSEXP, SEXP DnrowSEXP, SEXP penSEXP, SEXP penalty_whereSEXP, SEXP fSEXP, SEXP df_fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type surv2(surv2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort1(sort1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sort2(sort2SEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H(HSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dstar(DstarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wbeta(wbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n_npvar(n_npvarSEXP);
    Rcpp::traits::input_parameter< int >::type Dnrow(DnrowSEXP);
    Rcpp::traits::input_parameter< int >::type pen(penSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type penalty_where(penalty_whereSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< Function >::type df_f(df_fSEXP);
    rcpp_result_gen = Rcpp::wrap(fagfit_cpp(surv2, covar2, strata2, weights, offset, ibeta, sort1, sort2, method, maxiter, eps, H, Dstar, G, wbeta, lambda, alpha, gamma, M, d, n_npvar, Dnrow, pen, penalty_where, f, df_f));
    return rcpp_result_gen;
END_RCPP
}
// fcoxfit_cpp
List fcoxfit_cpp(NumericVector time, IntegerVector status, NumericMatrix covar2, NumericVector offset, NumericVector weights, IntegerVector strata2, int maxiter, double eps, NumericMatrix H, NumericMatrix Dstar, NumericMatrix G, NumericVector wbeta, int method, NumericVector ibeta, NumericVector lambda, double alpha, double gamma, int M, int d, int n_npvar, int Dnrow, int pen, IntegerVector penalty_where, Function f, Function df_f);
RcppExport SEXP _fcoxph_fcoxfit_cpp(SEXP timeSEXP, SEXP statusSEXP, SEXP covar2SEXP, SEXP offsetSEXP, SEXP weightsSEXP, SEXP strata2SEXP, SEXP maxiterSEXP, SEXP epsSEXP, SEXP HSEXP, SEXP DstarSEXP, SEXP GSEXP, SEXP wbetaSEXP, SEXP methodSEXP, SEXP ibetaSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP MSEXP, SEXP dSEXP, SEXP n_npvarSEXP, SEXP DnrowSEXP, SEXP penSEXP, SEXP penalty_whereSEXP, SEXP fSEXP, SEXP df_fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar2(covar2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata2(strata2SEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H(HSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dstar(DstarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wbeta(wbetaSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ibeta(ibetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n_npvar(n_npvarSEXP);
    Rcpp::traits::input_parameter< int >::type Dnrow(DnrowSEXP);
    Rcpp::traits::input_parameter< int >::type pen(penSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type penalty_where(penalty_whereSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< Function >::type df_f(df_fSEXP);
    rcpp_result_gen = Rcpp::wrap(fcoxfit_cpp(time, status, covar2, offset, weights, strata2, maxiter, eps, H, Dstar, G, wbeta, method, ibeta, lambda, alpha, gamma, M, d, n_npvar, Dnrow, pen, penalty_where, f, df_f));
    return rcpp_result_gen;
END_RCPP
}
// wshoot1
NumericVector wshoot1(arma::mat x, arma::vec y, NumericVector init, int pen, NumericVector weight, NumericVector wbeta, double lambda, double alpha, int maxiter, double tol, int n);
RcppExport SEXP _fcoxph_wshoot1(SEXP xSEXP, SEXP ySEXP, SEXP initSEXP, SEXP penSEXP, SEXP weightSEXP, SEXP wbetaSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type pen(penSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wbeta(wbetaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(wshoot1(x, y, init, pen, weight, wbeta, lambda, alpha, maxiter, tol, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fcoxph_fagfit_init", (DL_FUNC) &_fcoxph_fagfit_init, 10},
    {"_fcoxph_fagfit_loglik", (DL_FUNC) &_fcoxph_fagfit_loglik, 10},
    {"_fcoxph_fag_score", (DL_FUNC) &_fcoxph_fag_score, 8},
    {"_fcoxph_fcoxfit_init", (DL_FUNC) &_fcoxph_fcoxfit_init, 9},
    {"_fcoxph_fcoxfit_loglik", (DL_FUNC) &_fcoxph_fcoxfit_loglik, 9},
    {"_fcoxph_fcox_score", (DL_FUNC) &_fcoxph_fcox_score, 8},
    {"_fcoxph_fagfit_cpp", (DL_FUNC) &_fcoxph_fagfit_cpp, 26},
    {"_fcoxph_fcoxfit_cpp", (DL_FUNC) &_fcoxph_fcoxfit_cpp, 25},
    {"_fcoxph_wshoot1", (DL_FUNC) &_fcoxph_wshoot1, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_fcoxph(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
