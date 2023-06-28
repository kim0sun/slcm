#include <Rcpp.h>
#include "auxiliary.h"
using namespace Rcpp;

// [[Rcpp::export]]
List simModel(
      int nobs, IntegerVector nvar, List nlev, NumericVector par,
      int nlv, int nrl, int nlf, int npi, int ntau, int nrho,
      IntegerVector ul, IntegerVector vl, IntegerVector lf,
      IntegerVector rt, IntegerVector eqrl, IntegerVector eqlf,
      IntegerVector nc, IntegerVector nk, IntegerVector nl, IntegerVector ncl,
      IntegerVector nc_pi, IntegerVector nk_tau, IntegerVector nl_tau,
      IntegerVector nc_rho, IntegerVector nr_rho
) {
   List y(nlf);
   IntegerMatrix cl(nobs, nlv);
   double *_par_ = par.begin();

   std::vector<double*> _pi_(npi);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);

   for (int r = 0; r < npi; r ++) {
      for (int i = 0; i < nobs; i ++) {
         cl(i, rt[r]) = sample1(nc_pi[r], _par_);
      }
      _par_ += nc_pi[r];
   }
   for (int d = 0; d < ntau; d ++) {
      _tau_[d] = _par_;
      _par_ += nk_tau[d] * nl_tau[d];
   }
   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      _par_ += nr_rho[v] * nc_rho[v];
   }

   for (int d = 0; d < nrl; d ++) {
      double* tau = _tau_[eqrl[d]];
      for (int i = 0; i < nobs; i ++) {
         cl(i, ul[d]) = sample1(nk[d], tau + nk[d] * cl(i, vl[d]));
      }
   }

   for (int u = 0 ; u < nlf; u ++) {
      IntegerVector lev = nlev[eqlf[u]];
      NumericMatrix yy(nobs, nvar[eqlf[u]]);
      for (int i = 0; i < nobs; i ++) {
         double* rho = _rho_[eqlf[u]] + sum(lev) * cl(i, lf[u]);
         for (int m = 0; m < nvar[eqlf[u]]; m ++) {
            yy(i, m) = sample1(lev[m], rho);
            rho += lev[m];
         }
      }
      y[u] = yy;
   }

   List res;
   res["class"] = cl;
   res["y"] = y;

   return res;
}
