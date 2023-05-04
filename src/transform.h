#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <Rcpp.h>
Rcpp::NumericVector logit_pi(
   Rcpp::NumericVector pi, int nclass, int ref
);
Rcpp::NumericMatrix logit_tau(
   Rcpp::NumericMatrix tau, int nk, int nl,
   int *restr, int *ref
);
Rcpp::NumericVector logit_rho(
   Rcpp::NumericVector rho, Rcpp::IntegerVector ncat, int nclass,
   int *restr, int *ref
);
Rcpp::NumericVector log2logit(
   Rcpp::List param, int npar, Rcpp::List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   Rcpp::IntegerVector nclass_root, Rcpp::IntegerVector nclass_leaf,
   Rcpp::IntegerVector nclass_u, Rcpp::IntegerVector nclass_v,
   Rcpp::LogicalVector restr, Rcpp::IntegerVector ref
);
Rcpp::List logit2log(
   Rcpp::NumericVector logit, Rcpp::List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
   Rcpp::IntegerVector nclass_root, Rcpp::IntegerVector nclass_leaf,
   Rcpp::IntegerVector nclass_u, Rcpp::IntegerVector nclass_v,
   Rcpp::LogicalVector restr, Rcpp::IntegerVector ref
);
Rcpp::List splitlogit(
   Rcpp::NumericVector logit, Rcpp::List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
   Rcpp::IntegerVector nclass_root, Rcpp::IntegerVector nclass_u,
   Rcpp::IntegerVector nclass_v, Rcpp::IntegerVector nclass_leaf,
   Rcpp::LogicalVector restr, Rcpp::IntegerVector ref
);
Rcpp::List splitSE(
   Rcpp::NumericVector se, Rcpp::List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   Rcpp::IntegerVector root, Rcpp::IntegerVector ulv, Rcpp::IntegerVector vlv,
   Rcpp::IntegerVector nclass, Rcpp::IntegerVector nclass_u,
   Rcpp::IntegerVector nclass_v, Rcpp::IntegerVector nclass_leaf,
   Rcpp::LogicalVector restr, Rcpp::IntegerVector ref
);

#endif
