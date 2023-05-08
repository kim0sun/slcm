#include <Rcpp.h>
#include "src/auxiliary.h"
#include "src/transform.h"
#include "src/param.h"
#include "src/upward.h"
#include "src/downward.h"
#include "src/mstep.h"
using namespace Rcpp;

// [[Rcpp::export]]
List emFit12(
      IntegerVector y, int nobs,
      IntegerVector nvar, List ncat,
      int nlv, int nroot, int nlink, int nleaf,
      int nlink_unique, int nleaf_unique,
      IntegerVector tree_index, IntegerVector root,
      IntegerVector ulv, IntegerVector vlv, IntegerVector leaf,
      IntegerVector cstr_link, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v,
      NumericVector par, LogicalVector restr0
) {
   int max_iter = 1000; double tol = 1e-6;
   bool verbose = true; int newiter = 100;
   bool reg = false;
   int *py;
   NumericVector ss(par.length());
   NumericVector denom(par.length());

   std::vector<double*> _pi_(nroot);
   std::vector<double*> _tau_(nlink_unique);
   std::vector<double*> _rho_(nleaf_unique);

   std::vector<double*> _pi_d_(nroot), _pi_ss_(nroot);
   std::vector<double*> _tau_d_(nlink_unique), _tau_ss_(nlink_unique);
   std::vector<double*> _rho_d_(nleaf_unique), _rho_ss_(nleaf_unique);

   std::vector<double*> _ll_(nroot);
   std::vector<double*> _a_(nlv), _l_(nlv), _j_(nlink);
   std::vector<double*> _post_(nlv), _joint_(nlink);

   int *restr_;


   double *_par_ = par.begin();
   double *_ss_ = ss.begin();
   double *_denom_ = denom.begin();

   for (int r = 0; r < nroot; r ++) {
      _pi_[r] = _par_;
      _pi_ss_[r] = _ss_;
      _pi_d_[r] = _denom_;

      _par_ += nclass[r];
      _ss_ += nclass[r];
      _denom_ ++;
   }

   for (int d = 0; d < nlink_unique; d ++) {
      _tau_[d] = _par_;
      _tau_ss_[d] = _ss_;
      _tau_d_[d] = _denom_;

      _par_ += nclass_u[d] * nclass_v[d];
      _ss_ += nclass_u[d] * nclass_v[d];
      _denom_ += nclass_v[d];
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      _rho_[v] = _par_;
      _rho_ss_[v] = _ss_;
      _rho_d_[v] = _denom_;

      _par_ += sum(ncatv) * nclass_leaf[v];
      _ss_ += sum(ncatv) * nclass_leaf[v];
      _denom_ += nclass_leaf[v];
   }

   NumericMatrix ll(nobs, nroot);
   double* tmp1 = ll.begin();
   for (int r = 0; r < nroot; r ++) {
      _ll_[r] = tmp1;
      tmp1 += nobs;
   }

   NumericVector j(nobs * sum(nclass));
   NumericVector joint(nobs * sum(nclass));
   tmp1 = j.begin();
   double* tmp2 = joint.begin();
   for (int d = 0; d < nlink; d ++) {
      int u = ulv[d]; int v = vlv[d];
      _j_[d] = tmp1;
      _joint_[d] = tmp2;
      tmp1 += nclass[v] * nobs;
      tmp2 += nclass[u] * nclass[v] * nobs;
   }

   NumericVector alpha(nobs * sum(nclass));
   NumericVector lambda(nobs * sum(nclass));
   NumericVector post(nobs * sum(nclass));
   tmp1 = alpha.begin();
   tmp2 = lambda.begin();
   double* tmp3 = post.begin();
   for (int u = 0; u < nlv; u ++) {
      _a_[u] = tmp1;
      _l_[u] = tmp2;
      _post_[u] = tmp3;
      tmp1 += nclass[u] * nobs;
      tmp2 += nclass[u] * nobs;
      tmp3 += nclass[u] * nobs;
   }
   delete tmp1;
   delete tmp2;
   delete tmp3;

   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;

   while ( (iter < max_iter) && (dll > tol) ) {
      // ss to parameter
      if (iter > 0) {

      }
      iter ++;
      lastll = currll;

      // alpha, lambda initiate
      alpha.fill(R_NegInf);
      lambda.fill(0);

      // (expectation-step)
      // initiate lambda
      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         upInit(py, _rho_[cstr_leaf[v]], _l_[leaf[v]], nclass[leaf[v]],
                nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
         py += nobs * nvar[cstr_leaf[v]];
      }

      // upward recursion
      for (int d = nlink - 1; d > -1; d --) {
         int u = ulv[d], v = vlv[d];
         upRec(_l_[v], _j_[d], _l_[u], _tau_[cstr_link[d]],
               nobs, nclass[u], nclass[v], reg);
      }

      // initiate alpha
      for (int r = 0; r < nroot; r ++) {
         dnInit(_a_[root[r]], _l_[root[r]], _pi_[r],
                _post_[root[r]], _ll_[r], nobs,
                nclass[root[r]], reg);
      }

      // Downward recursion
      for (int d = 0; d < nlink; d ++) {
         int u = ulv[d];
         int v = vlv[d];
         dnRec(_a_[u], _a_[v], _l_[u], _l_[v], _j_[d],
               nobs, nclass[u], nclass[v], _tau_[cstr_link[d]],
                                                _post_[u], _joint_[d], _ll_[tree_index[d]], reg);
      }

      // (maximization-step)
      // pi updates
      for (int r = 0; r < nroot; r ++) {
         updatePi(_pi_[r], _post_[root[r]], nobs,
                  nclass[root[r]]);
      }

      restr_ = restr0.begin();
      // tau updates
      ss.fill(R_NegInf);
      denom.fill(R_NegInf);
      for (int d = 0; d < nlink; d ++) {
         int u = ulv[d]; int v = vlv[d];
         cumTau(_tau_ss_[cstr_link[d]], _tau_d_[cstr_link[d]],
                _joint_[d], nobs, nclass[u], nclass[v]);
      }
      for (int d = 0; d < nlink_unique; d ++) {
         updateTau(_tau_[d], _tau_ss_[d], _tau_d_[d],
                   nclass_u[d], nclass_v[d], restr_);
         restr_ += nclass_u[d] * nclass_v[d];
      }

      // rho updates
      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         int u = leaf[v], cu = cstr_leaf[v];
         IntegerVector ncatv = ncat[cu];
         cumRho(_rho_ss_[cu], _rho_d_[cu], py, nobs, nvar[cu],
                ncatv, nclass[u], _post_[u], _rho_[cu]);
         py += nobs * nvar[cu];
      }

      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         updateRho(_rho_[v], _rho_ss_[v], _rho_d_[v],
                   nobs, nclass_leaf[v], nvar[v], ncatv, restr_);
         restr_ += nclass_leaf[v] * sum(ncatv);
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


   List res;
   res["params"] = par;
   res["converged"] = dll < tol;
   res["niter"] = iter;

   return res;
}
