#include <Rcpp.h>
#include "upward.h"
#include "downward.h"
using namespace Rcpp;

// [[Rcpp::export]]
List calcModel(
      IntegerVector y,
      int nobs, IntegerVector nvar, IntegerVector nlev,
      NumericVector par, LogicalVector fix0, IntegerVector ref,
      int nlv, int nrl, int nlf, int npi, int ntau, int nrho,
      IntegerVector ul, IntegerVector vl, IntegerVector lf,
      IntegerVector tr, IntegerVector rt, IntegerVector eqrl, IntegerVector eqlf,
      IntegerVector nc, IntegerVector nk, IntegerVector nl, IntegerVector ncl,
      IntegerVector nc_pi, IntegerVector nk_tau, IntegerVector nl_tau,
      IntegerVector nc_rho, IntegerVector nr_rho
) {
   int *_y_;
   int *_ref_ = ref.begin();
   NumericVector score(nobs * par.length());
   double *_par_ = par.begin();
   double *_scr_ = score.begin();

   std::vector<double*> _pi_(npi);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);

   std::vector<double*> _scr1_(npi);
   std::vector<double*> _scr2_(ntau);
   std::vector<double*> _scr3_(nrho);

   std::vector<int*> _ref1_(ntau);
   std::vector<int*> _ref2_(nrho);

   std::vector<double*> _ll_(npi);
   std::vector<double*> _a_(nlv), _l_(nlv), _j_(nrl);
   std::vector<double*> _post_(nlv), _joint_(nrl);
   std::vector<int*> _nlev_(nrho);

   for (int r = 0; r < npi; r ++) {
      _pi_[r] = _par_;
      _scr1_[r] = _scr_;

      _par_ += nc_pi[r];
      _scr_ += nc_pi[r] * nobs;
   }
   _ref_ += npi;

   for (int d = 0; d < ntau; d ++) {
      _tau_[d] = _par_;
      _scr2_[d] = _scr_;
      _ref1_[d] = _ref_;
      _par_ += nk_tau[d] * nl_tau[d];
      _scr_ += nk_tau[d] * nl_tau[d] * nobs;
      _ref_ += nl_tau[d];
   }

   int *_lev_ = nlev.begin();
   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      _scr3_[v] = _scr_;
      _ref2_[v] = _ref_;
      _nlev_[v] = _lev_;

      _par_ += nr_rho[v] * nc_rho[v];
      _scr_ += nr_rho[v] * nc_rho[v] * nobs;
      _ref_ += nvar[v] * nc_rho[v];
      _lev_ += nvar[v];
   }

   NumericMatrix ll(nobs, npi);
   double *_tmp1_ = ll.begin();
   for (int r = 0; r < npi; r ++) {
      _ll_[r] = _tmp1_;
      _tmp1_ += nobs;
   }

   NumericVector j(nobs * sum(nl));
   NumericVector joint(nobs * sum(nk * nl));
   _tmp1_ = j.begin();
   double *_tmp2_ = joint.begin();
   for (int d = 0; d < nrl; d ++) {
      _j_[d] = _tmp1_;
      _joint_[d] = _tmp2_;
      _tmp1_ += nl[d] * nobs;
      _tmp2_ += nk[d] * nl[d] * nobs;
   }

   NumericVector alpha(nobs * sum(nc));
   NumericVector lambda(nobs * sum(nc));
   NumericVector post(nobs * sum(nc));
   _tmp1_ = alpha.begin();
   _tmp2_ = lambda.begin();
   double *_tmp3_ = post.begin();
   for (int u = 0; u < nlv; u ++) {
      _a_[u] = _tmp1_;
      _l_[u] = _tmp2_;
      _post_[u] = _tmp3_;
      _tmp1_ += nc[u] * nobs;
      _tmp2_ += nc[u] * nobs;
      _tmp3_ += nc[u] * nobs;
   }

   // (expectation-step)
   // initiate lambda
   _y_ = y.begin();
   for (int v = 0; v < nlf; v ++) {
      upInit(_y_, _rho_[eqlf[v]], _l_[lf[v]], ncl[v],
             nobs, nvar[eqlf[v]], _nlev_[eqlf[v]]);
      _y_ += nobs * nvar[eqlf[v]];
   }

   // upward recursion
   for (int d = nrl - 1; d > -1; d --) {
      int u = ul[d], v = vl[d];
      upRec(_l_[v], _j_[d], _l_[u], _tau_[eqrl[d]],
            nobs, nk[d], nl[d], false);
   }

   // calculate loglik
   for (int r = 0; r < npi; r ++) {
      calclri(_l_[rt[r]], _pi_[r], _ll_[r],
              nobs, nc_pi[r], false);
   }

   // initiate alpha
   for (int r = 0; r < npi; r ++) {
      dnInit(_a_[rt[r]], _l_[rt[r]], _pi_[r], _post_[rt[r]],
             _ll_[r], nobs, nc_pi[r], false);
   }

   // Downward recursion
   for (int d = 0; d < nrl; d ++) {
      int u = ul[d], v = vl[d];
      dnRec(_a_[u], _a_[v], _l_[u], _l_[v], _j_[d], nobs, nk[d], nl[d],
            _tau_[eqrl[d]], _post_[u], _joint_[d], _ll_[tr[d]], false);
   }

   // score function
   _ref_ = ref.begin();
   for (int r = 0; r < npi; r ++) {
      _scr_ = _scr1_[r];
      _tmp1_ = _post_[rt[r]];
      _tmp2_ = _pi_[r];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < nc_pi[r]; k ++) {
            if (k == _ref_[r]) _scr_[k] = R_NaN;
            else _scr_[k] = exp(_tmp1_[k]) - exp(_tmp2_[k]);
         }
         _tmp1_ += nc_pi[r];
         _scr_ += nc_pi[r];
      }
   }

   for (int d = 0; d < nrl; d ++) {
      int u = eqrl[d];
      _scr_ = _scr2_[u];
      _ref_ = _ref1_[u];
      _tmp1_ = _joint_[d];
      _tmp2_ = _post_[vl[d]];
      for (int i = 0; i < nobs; i ++) {
         _tmp3_ = _tau_[u];
         for (int l = 0; l < nl[d]; l ++) {
            for (int k = 0; k < nk[d]; k ++) {
               if (k == _ref_[l]) _scr_[k] = R_NaN;
               else _scr_[k] += exp(_tmp1_[k]) - exp(_tmp2_[l] + _tmp3_[k]);
            }
            _scr_ += nk[d];
            _tmp1_ += nk[d];
            _tmp3_ += nk[d];
         }
         _tmp2_ += nl[d];
      }
   }

   _y_ = y.begin();
   for (int v = 0; v < nlf; v ++) {
      int u = eqlf[v];
      _scr_ = _scr3_[u];
      _lev_ = _nlev_[u];
      _tmp1_ = _post_[lf[v]];
      for (int i = 0; i < nobs; i ++) {
         _tmp2_ = _rho_[u];
         _ref_ = _ref2_[u];
         for (int k = 0; k < nc_rho[u]; k ++) {
            for (int m = 0; m < nvar[u]; m ++) {
               _scr_[_y_[m] - 1] += exp(_tmp1_[k]);
               for (int r = 0; r < _lev_[m]; r ++) {
                  if (r == _ref_[m]) _scr_[r] = R_NaN;
                  else _scr_[r] -= exp(_tmp1_[k] + _tmp2_[r]);
               }
               _scr_ += _lev_[m];
               _tmp2_ += _lev_[m];
            }
            _ref_ += nvar[u];
         }
         _tmp1_ += nc_rho[u];
         _y_ += nvar[u];
      }
   }

   List res;
   res["ll"] = ll;
   res["score"] = score;
   res["post"] = post;
   res["joint"] = joint;

   return res;
}
