#ifndef AUXILIARY_H
#define AUXILIARY_H

int sample1(int n, double *prob);
double log_add_exp(double x, double y);
double log_sum_exp(Rcpp::NumericVector x);
Rcpp::NumericVector elogdiri(Rcpp::NumericVector a);
Rcpp::NumericVector plogdiri(Rcpp::NumericVector a, Rcpp::NumericVector b);

#endif
