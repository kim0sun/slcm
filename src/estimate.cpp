#include <Rcpp.h>
#include "auxiliary.h"
#include "upward.h"
#include "downward.h"
#include "mstep.h"
using namespace Rcpp;

// [[Rcpp::export]]
List em_est(
      IntegerVector y,
      int nobs, IntegerVector nvar, IntegerVector nlev,
      NumericVector par_origin, LogicalVector fix0,
      int nlv, int nrl, int nlf, int npi, int ntau, int nrho,
      IntegerVector ul, IntegerVector vl, IntegerVector lf,
      IntegerVector tr, IntegerVector rt, IntegerVector eqrl, IntegerVector eqlf,
      IntegerVector nc, IntegerVector nk, IntegerVector nl, IntegerVector ncl,
      IntegerVector nc_pi, IntegerVector nk_tau, IntegerVector nl_tau,
      IntegerVector nc_rho, IntegerVector nr_rho,
      int max_iter, double tol, bool verbose, int newiter
) {
   List res;
   int *_y_;
   NumericVector par = clone(par_origin);
   NumericVector ss(par.length());
   NumericVector denom(npi + sum(nl_tau) + sum(nc_rho));

   std::vector<double*> _pi_(npi);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);

   std::vector<double*> _pi_d_(npi), _pi_ss_(npi);
   std::vector<double*> _tau_d_(ntau), _tau_ss_(ntau);
   std::vector<double*> _rho_d_(nrho), _rho_ss_(nrho);

   std::vector<double*> _ll_(npi);
   std::vector<double*> _a_(nlv), _l_(nlv), _j_(nrl);
   std::vector<double*> _post_(nlv), _joint_(nrl);

   int *_fix0_;
   std::vector<int*> _nlev_(nrho);

   double *_par_ = par.begin();
   double *_ss_ = ss.begin();
   double *_denom_ = denom.begin();

   for (int r = 0; r < npi; r ++) {
      _pi_[r] = _par_;
      _pi_ss_[r] = _ss_;
      _pi_d_[r] = _denom_;

      _par_ += nc_pi[r];
      _ss_ += nc_pi[r];
      _denom_ ++;
   }

   for (int d = 0; d < ntau; d ++) {
      _tau_[d] = _par_;
      _tau_ss_[d] = _ss_;
      _tau_d_[d] = _denom_;

      _par_ += nk_tau[d] * nl_tau[d];
      _ss_ += nk_tau[d] * nl_tau[d];
      _denom_ += nl_tau[d];
   }

   int *_tmp_ = nlev.begin();
   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      _rho_ss_[v] = _ss_;
      _rho_d_[v] = _denom_;
      _nlev_[v] = _tmp_;

      _par_ += nr_rho[v] * nc_rho[v];
      _ss_ += nr_rho[v] * nc_rho[v];
      _denom_ += nc_rho[v];
      _tmp_ += nvar[v];
   }

   NumericVector ll(nobs * npi);
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

   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;

   while ( (iter < max_iter) && (dll > tol) ) {
      // ss to parameter
      if (iter > 0) {
         _fix0_ = fix0.begin();
         for (int r = 0; r < npi; r ++) {
            updatePi(_pi_[r], _pi_ss_[r], _pi_d_[r], nc_pi[r]);
            _fix0_ += nc_pi[r];
         }
         for (int d = 0; d < ntau; d ++) {
            updateTau(_tau_[d], _tau_ss_[d], _tau_d_[d],
                      nk_tau[d], nl_tau[d], _fix0_);
            _fix0_ += nk_tau[d] * nl_tau[d];
         }
         for (int v = 0; v < nrho; v ++) {
            updateRho(_rho_[v], _rho_ss_[v], _rho_d_[v], nobs,
                      nc_rho[v], nvar[v], _nlev_[v], _fix0_);
            _fix0_ += nr_rho[v] * nc_rho[v];
         }
      }
      iter ++;
      lastll = currll;

      // alpha, lambda initiate
      alpha.fill(R_NegInf);
      lambda.fill(0);

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

      // build sufficient statistics
      ss.fill(R_NegInf);
      denom.fill(R_NegInf);
      // pi
      for (int r = 0; r < npi; r ++) {
         cumPi(_pi_ss_[r], _pi_d_[r], _post_[rt[r]], nobs, nc_pi[r]);
      }
      // tau
      for (int d = 0; d < nrl; d ++) {
         cumTau(_tau_ss_[eqrl[d]], _tau_d_[eqrl[d]], _joint_[d], nobs, nk[d], nl[d]);
      }
      // rho
      _y_ = y.begin();
      for (int v = 0; v < nlf; v ++) {
         int u = lf[v], w = eqlf[v];
         cumRho(_rho_ss_[w], _rho_d_[w], _y_, nobs, nvar[w], _nlev_[w],
                ncl[v], _post_[u], _rho_[w]);
         _y_ += nobs * nvar[w];
      }
      currll = 0;
      currll = sum(ll);

      if (lastll == R_NegInf) dll = R_PosInf;
      else dll = currll - lastll;
      if (verbose) {
         Rcout << iter << " iterations  logLik: " <<
            std::fixed << std::setprecision(2) << currll << "  diff: ";
         if (dll > 1e-5) Rcout << std::setprecision(5) << dll;
         else Rcout << std::scientific << std::setprecision(1) << dll;
         Rcout << "      \r" << std::flush;

         if (iter % newiter == 0) {
            Rcout << "\n" << std::flush;
         }
      }
   }
   if (verbose) {
      if (iter % newiter != 0) Rcout << "\n";
   }

   // Final estimates
   _fix0_ = fix0.begin();
   for (int r = 0; r < npi; r ++) {
      updatePi(_pi_[r], _pi_ss_[r], _pi_d_[r], nc_pi[r]);
      _fix0_ += nc_pi[r];
   }
   for (int d = 0; d < ntau; d ++) {
      updateTau(_tau_[d], _tau_ss_[d], _tau_d_[d],
                nk_tau[d], nl_tau[d], _fix0_);
      _fix0_ += nk_tau[d] * nl_tau[d];
   }
   for (int v = 0; v < nrho; v ++) {
      updateRho(_rho_[v], _rho_ss_[v], _rho_d_[v], nobs,
                nc_rho[v], nvar[v], _nlev_[v], _fix0_);
      _fix0_ += nc_rho[v] * nr_rho[v];
   }

   res["param"] = par;
   res["converged"] = dll < tol;
   res["niter"] = iter;
   res["ll"] = currll;
   return res;
}


// [[Rcpp::export]]
double fll(
      IntegerVector y, NumericVector par,
      int nobs, IntegerVector nvar, IntegerVector nlev,
      int nlv, int nrl, int nlf, int npi, int ntau, int nrho,
      IntegerVector ul, IntegerVector vl, IntegerVector lf,
      IntegerVector tr, IntegerVector rt, IntegerVector eqrl, IntegerVector eqlf,
      IntegerVector nc, IntegerVector nk, IntegerVector nl, IntegerVector ncl,
      IntegerVector nc_pi, IntegerVector nk_tau, IntegerVector nl_tau,
      IntegerVector nc_rho, IntegerVector nr_rho
) {
   std::vector<double*> _pi_(npi);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);
   std::vector<int*> _nlev_(nrho);

   double *_par_ = par.begin();
   for (int r = 0; r < npi; r ++) {
      _pi_[r] = _par_;
      _par_ += nc_pi[r];
   }
   for (int d = 0; d < ntau; d ++) {
      _tau_[d] = _par_;
      _par_ += nk_tau[d] * nl_tau[d];
   }
   int *_tmp_ = nlev.begin();
   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      _nlev_[v] = _tmp_;
      _par_ += nr_rho[v] * nc_rho[v];
      _tmp_ += nvar[v];
   }

   NumericVector lambda(nobs * sum(nc));
   std::vector<double*> _l_(nlv);
   double *_tmp2_ = lambda.begin();
   for (int v = 0; v < nlv; v ++) {
      _l_[v] = _tmp2_;
      _tmp2_ += nc[v] * nobs;
   }

   double ll = 0.0;
   // initiate lambda
   int *_y_ = y.begin();
   for (int v = 0; v < nlf; v ++) {
      upInit(_y_, _rho_[eqlf[v]], _l_[lf[v]], ncl[v],
             nobs, nvar[eqlf[v]], _nlev_[eqlf[v]]);
      _y_ += nobs * nvar[eqlf[v]];
   }

   // upward recursion
   for (int d = nrl - 1; d > -1; d --) {
      upRec2(_l_[vl[d]], _l_[ul[d]], _tau_[eqrl[d]],
             nobs, nk[d], nl[d], false);
   }

   for (int r = 0; r < npi; r ++) {
      ll += calclr(_l_[rt[r]], _pi_[r], nobs, nc_pi[r], false);
   }

   return -ll;
}

