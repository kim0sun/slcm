#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pi_gnr(int nk) {
   NumericVector pi(nk);
   pi.fill(-log(nk));
   return pi;
}

// [[Rcpp::export]]
NumericMatrix tau_gnr(int nk, int nl) {
   NumericMatrix tau(nk, nl);
   double *ptau = tau.begin();

   for (int l = 0; l < nl; l ++) {
      NumericVector p = runif(nk, 0, 1);
      for (int k = 0; k < nk; k ++) {
         ptau[k] = log(p[k]) - log(sum(p));
      }
      ptau += nk;
   }

   return tau;
}

// [[Rcpp::export]]
NumericMatrix rho_gnr(int nk, IntegerVector ncat) {
   NumericMatrix pv(sum(ncat), nk);
   double *ppv = pv.begin();

   for (int k = 0; k < nk; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         NumericVector p = runif(ncat[m], 0, 1);
         for (int r = 0; r < ncat[m]; r ++) {
            ppv[r] = log(p[r] / sum(p));
         }
         ppv += ncat[m];
      }
   }
   return pv;
}


// [[Rcpp::export]]
List par_gnr(
      int nobs, IntegerVector nvar, List ncat,
      int nroot, int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v,
      LogicalVector init, List init_param
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   if (init[0]) {
      List piList = init_param["pi"];
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi_init = piList[r];
         NumericVector pi = clone(pi_init);
         lst_pi[r] = pi;
      }
   } else {
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi = pi_gnr(nclass[root[r]]);
         lst_pi[r] = pi;
      }
   }

   if (init[1]) {
      List tauList = init_param["tau"];
      for (int d = 0; d < nlink_unique; d ++) {
         NumericMatrix tau_init = tauList[d];
         NumericMatrix tau = clone(tau_init);
         lst_tau[d] = tau;
      }
   } else {
      for (int d = 0; d < nlink_unique; d ++) {
         NumericMatrix tau = tau_gnr(nclass_u[d], nclass_v[d]);
         lst_tau[d] = tau;
      }
   }

   if (init[2]) {
      List rhoList = init_param["rho"];
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericMatrix rho_init = rhoList[v];
         NumericMatrix lrho = clone(rho_init);
         lst_rho[v] = lrho;
      }
   } else {
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericMatrix lrho = rho_gnr(nclass_leaf[v], ncatv);
         lst_rho[v] = lrho;
      }
   }

   List par;
   par["pi"] = lst_pi;
   par["tau"] = lst_tau;
   par["rho"] = lst_rho;

   return par;
}

