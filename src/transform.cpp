#include <Rcpp.h>
using namespace Rcpp;

void logit_pi(
   double *logit, NumericVector pi, int nclass, int ref
) {
   for (int k = 0; k < nclass; k ++) {
      if (k == ref) logit[k] = R_NaN;
      else logit[k] = pi[k] - pi[ref];
   }
}

void logit_tau(
   double *logit, NumericMatrix tau, int nk, int nl,
   int *restr, int *ref
) {
   double *tau_ = tau.begin();

   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk; k ++) {
         if (k == ref[l]) logit[k] = R_NaN;
         else logit[k] = tau_[k] - tau_[ref[l]];
      }
      logit += nk;
      restr += nk;
      tau_  += nk;
   }
}

void logit_rho(
   double *logit, NumericVector rho, IntegerVector ncat,
   int nclass, int *restr, int *ref
) {
   double *rho_ = rho.begin();

   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         for (int r = 0; r < ncat[m]; r ++) {
            if (r == ref[m]) logit[r] = R_NaN;
            else logit[r] = rho_[r] - rho_[ref[m]];
         }
         logit += ncat[m];
         restr += ncat[m];
         rho_  += ncat[m];
      }
      ref += ncat.length();
   }
}

// [[Rcpp::export]]
List logit2log(
   NumericVector logit, List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   IntegerVector root, IntegerVector ulv, IntegerVector vlv,
   IntegerVector nclass_root, IntegerVector nclass_leaf,
   IntegerVector nclass_u, IntegerVector nclass_v,
   LogicalVector restr, IntegerVector ref
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *logit_ = logit.begin();
   int *ref_ = ref.begin();

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass_root[r];
      NumericVector pi(nk);
      double denom = 1;
      for (int k = 0; k < nk; k ++) {
         if (k == ref_[r]) continue;
         pi[k] = logit_[k];
         denom += exp(logit_[k]);
      }
      for (int k = 0; k < nk; k ++) {
         pi[k] -= log(denom);
      }

      logit_ += nk;
      lst_pi[r] = pi;
   }
   ref_ += nroot;

   double denom = 0;
   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d];
      int nl = nclass_v[d];
      NumericMatrix tau(nk, nl);
      double *tau_ = tau.begin();
      for (int l = 0; l < nl; l ++) {
         denom = 1;
         for (int k = 0; k < nk; k ++) {
            if (k == ref_[l]) continue;
            tau_[k] = logit_[k];
            denom += exp(logit_[k]);
         }
         for (int k = 0; k < nk; k ++) {
            tau_[k] -= log(denom);
         }
         tau_ += nk;
         logit_ += nk;
      }

      lst_tau[d] = tau;
      ref_ += nl;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      int nk = nclass_leaf[v];
      NumericVector rho(nk * sum(ncatv));
      double *rho_ = rho.begin();
      for (int k = 0; k < nk; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            denom = 1;
            for (int r = 0; r < ncatv[m] - 1; r ++) {
               if (r == ref_[m]) continue;
               rho_[r] = logit_[r];
               denom += exp(logit_[r]);
            }
            for (int r = 0; r < ncatv[m]; r ++) {
               rho_[r] -= log(denom);
            }
            rho_  += ncatv[m];
            logit_ += ncatv[m];
         }
         ref_ += ncat.length();
      }
      lst_rho[v] = rho;
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}



// [[Rcpp::export]]
NumericVector log2logit(
   List param, int npar, List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   IntegerVector nclass_root, IntegerVector nclass_leaf,
   IntegerVector nclass_u, IntegerVector nclass_v,
   LogicalVector restr, IntegerVector ref
) {
   List lst_pi = param["pi"];
   List lst_tau = param["tau"];
   List lst_rho = param["rho"];

   NumericVector logit(npar);
   double *logit_ = logit.begin();
   int *restr_ = restr.begin();
   int *ref_ = ref.begin();

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass_root[r];
      logit_pi(logit_, lst_pi[r], nk, ref[r]);

      logit_ += nk;
   }
   ref_ += nroot;

   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d];
      int nl = nclass_v[d];
      logit_tau(logit_, lst_tau[d], nk, nl, restr_, ref_);

      logit_ += nl * nk;
      restr_ += nl * nk;
      ref_   += nl;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      int nk = nclass_leaf[v];
      logit_rho(logit_, lst_rho[v], ncat[v], nk, restr_, ref_);

      logit_ += nk * sum(ncatv);
      restr_ += nk * sum(ncatv);
      ref_   += nk * ncatv.length();
   }

   return logit;
}


// [[Rcpp::export]]
List param2list(
   NumericVector param, List ncat,
   int nroot, int nlink_unique, int nleaf_unique,
   IntegerVector root, IntegerVector ulv, IntegerVector vlv,
   IntegerVector nclass_root, IntegerVector nclass_u,
   IntegerVector nclass_v, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);

   double *param_ = param.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi(nclass_root[r]);
      for (int k = 0; k < nclass_root[r]; k ++) {
         pi[k] = param_[k];
      }
      param_ += nclass_root[r];
      lst_pi[r] = pi;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau(nclass_u[d], nclass_v[d]);
      for (int l = 0; l < nclass_v[d]; l ++) {
         for (int k = 0; k < nclass_u[d]; k ++) {
            tau(k, l) = param_[k];
         }
         param_ += nclass_u[d];
      }
      lst_tau[d] = tau;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericMatrix rho(sum(ncatv), nclass_leaf[v]);
      double *rho_ = rho.begin();
      for (int k = 0; k < nclass_leaf[v]; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            for (int r = 0; r < ncatv[m]; r ++) {
               rho_[r] = param_[r];
            }
            param_ += ncatv[m];
            rho_ += ncatv[m];
         }
      }
      lst_rho[v] = rho;
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}

