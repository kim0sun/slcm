#include <Rcpp.h>
#include "auxiliary.h"
#include "transform.h"
#include "param.h"
#include "upward.h"
#include "downward.h"
#include "mstep.h"
using namespace Rcpp;

// [[Rcpp::export]]
List emFit1(
   IntegerVector y, List num, List idx, List nclass, List dim,
   NumericVector par, LogicalVector reg, LogicalVector restr0, List control
) {
   int nobs = num["nobs"];
   IntegerVector nvar = num["nvar"];
   IntegerVector nresp = num["nresp"];

   int nlv = num["nlv"];
   int nr = num["nr"];
   int nrl = num["nrl"];
   int nlf = num["nlf"];
   int ntau = num["ntau"];
   int nrho = num["nrho"];

   IntegerVector ul = idx["ul"];
   IntegerVector vl = idx["vl"];
   IntegerVector lf = idx["lf"];
   IntegerVector tr = idx["tr"];
   IntegerVector rt = idx["rt"];
   IntegerVector eqrl = idx["eqrl"];
   IntegerVector eqlf = idx["eqlf"];

   IntegerVector nc = nclass["nc"];
   IntegerVector nk = nclass["nk"];
   IntegerVector nl = nclass["nl"];
   IntegerVector ncl = nclass["ncl"];

   IntegerVector nc_pi = dim["nc_pi"];
   IntegerVector nk_tau = dim["nk_tau"];
   IntegerVector nl_tau = dim["nl_tau"];
   IntegerVector nr_rho = dim["nr_rho"];
   IntegerVector nc_rho = dim["nc_rho"];

   int max_iter = control["max_iter"];
   double tol = control["tol"];
   bool verbose = control["verbose"];
   int newiter = control["newiter"];

   int *_y_;
   NumericVector ss(par.length());
   NumericVector denom(par.length());

   std::vector<double*> _pi_(nr);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);

   std::vector<double*> _pi_d_(nr), _pi_ss_(nr);
   std::vector<double*> _tau_d_(ntau), _tau_ss_(ntau);
   std::vector<double*> _rho_d_(nrho), _rho_ss_(nrho);

   std::vector<double*> _ll_(nr);
   std::vector<double*> _a_(nlv), _l_(nlv), _j_(nrl);
   std::vector<double*> _post_(nlv), _joint_(nrl);

   int *_restr_;
   std::vector<double*> _nresp_(nrho);

   double *_par_ = par.begin();
   double *_ss_ = ss.begin();
   double *_denom_ = denom.begin();

   for (int r = 0; r < nr; r ++) {
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

   double* _tmp1_ = nresp.begin();
   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      _rho_ss_[v] = _ss_;
      _rho_d_[v] = _denom_;
      _nresp_[v] = _tmp1_;

      _par_ += nr_rho[v] * nc_rho[v];
      _ss_ += nr_rho[v] * nc_rho[v];
      _denom_ += nc_rho[v];
      _tmp1_ += nvar[v];
   }

   NumericMatrix ll(nobs, nr);
   _tmp1_ = ll.begin();
   for (int r = 0; r < nr; r ++) {
      _ll_[r] = _tmp1_;
      _tmp1_ += nobs;
   }

   NumericVector j(nobs * sum(nl));
   NumericVector joint(nobs * sum(nk * nl));
   _tmp1_ = j.begin();
   double* _tmp2_ = joint.begin();
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
   double* _tmp3_ = post.begin();
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
         _restr_ = restr0.begin();
         for (int d = 0; d < ntau; d ++) {
            updateTau(_tau_[d], _tau_ss_[d], _tau_d_[d],
                      nk_tau[d], nl_tau[d], _restr_);
            _restr_ += nk_tau[d] * nl_tau[d];
         }
         _nr_ = nr.begin();
         for (int v = 0; v < nrho; v ++) {
            updateRho(_rho_[v], _rho_ss_[v], _rho_d_[v], nobs,
                      nc_rho[v], nvar[v], _nresp_[v], _restr_);
            _restr_ += nr_rho[v] * nc_rho[v];
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
                nobs, nvar[eqlf[v]], nresp[eqlf[v]]);
         _y_ += nobs * nvar[eqlf[v]];
      }

      // upward recursion
      for (int d = nrl - 1; d > -1; d --) {
         int u = ul[d], v = vl[d];
         upRec(_l_[v], _j_[d], _l_[u], _tau_[eqrl[d]],
               nobs, nk[d], nl[d], false);
      }

      // initiate alpha
      for (int r = 0; r < nr; r ++) {
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
      for (int r = 0; r < nr; r ++) {
         updatePi(_pi_[r], _post_[rt[r]], nobs, nc_pi[r]);
      }
      // tau
      for (int d = 0; d < nrl; d ++) {
         cumTau(_tau_ss_[eqrl[d]], _tau_d_[eqrl[d]], _joint_[d], nobs, nk[d], nl[d]);
      }
      // rho
      _y_ = y.begin();
      for (int v = 0; v < nlf; v ++) {
         int u = lf[v], w = eqlf[v];
         IntegerVector nrw = nresp[w];
         cumRho(_rho_ss_[w], _rho_d_[w], _y_, nobs, nvar[w], _nresp_[v],
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
   _restr_ = restr0.begin();
   for (int d = 0; d < ntau; d ++) {
      updateTau(_tau_[d], _tau_ss_[d], _tau_d_[d],
                nk_tau[d], nl_tau[d], _restr_);
      _restr_ += nk_tau[d] * nl_tau[d];
   }
   for (int v = 0; v < nrho; v ++) {
      IntegerVector nr = nresp[v];
      updateRho(_rho_[v], _rho_ss_[v], _rho_d_[v], nobs,
                nc_rho[v], nvar[v], _nresp_[v], _restr_);
      _restr_ += nc_rho[v] * sum(nr);
   }


   List res;
   res["params"] = par;
   res["converged"] = dll < tol;
   res["niter"] = iter;

   return res;
}

// [[Rcpp::export]]
List emFit(
      IntegerVector y, int nobs,
      IntegerVector nvar, List ncat,
      int nlv, int nroot, int nlink, int nleaf,
      int nlink_unique, int nleaf_unique,
      IntegerVector tree_index, IntegerVector root,
      IntegerVector ulv, IntegerVector vlv, IntegerVector leaf,
      IntegerVector cstr_link, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v,
      List init_param, LogicalVector restr0,
      int max_iter, double tol,
      bool verbose, int newiter = 100,
      bool reg = false
) {
   int *py;

   List lst_pi(nroot), lst_tau(nlink_unique), lst_rho(nleaf_unique);
   List lst_ntau_d(nlink_unique), lst_ntau_n(nlink_unique);
   List lst_nrho_d(nleaf_unique), lst_nrho_n(nleaf_unique);
   std::vector<double*> ptr_pi(nroot);
   std::vector<double*> ptr_tau(nlink_unique);
   std::vector<double*> ptr_ntau_d(nlink_unique), ptr_ntau_n(nlink_unique);
   std::vector<double*> ptr_rho(nleaf_unique);
   std::vector<double*> ptr_nrho_d(nleaf_unique), ptr_nrho_n(nleaf_unique);

   List lst_ll(nroot);
   List lst_a(nlv), lst_l(nlv), lst_j(nlink);
   std::vector<double*> ptr_ll(nroot);
   std::vector<double*> ptr_a(nlv), ptr_l(nlv), ptr_j(nlink);

   List lst_post(nlv), lst_joint(nlink);
   std::vector<double*> ptr_post(nlv), ptr_joint(nlink);

   int *restr_;

   List piList = init_param["pi"];
   for (int r = 0; r < nroot; r ++) {
      NumericVector pi_init = piList[r];
      NumericVector pi = clone(pi_init);
      ptr_pi[r] = pi.begin();
      lst_pi[r] = pi;
   }

   List tauList = init_param["tau"];
   for (int d = 0; d < nlink_unique; d ++) {
      NumericMatrix tau_init = tauList[d];
      NumericMatrix tau = clone(tau_init);
      ptr_tau[d] = tau.begin();
      lst_tau[d] = tau;
      NumericMatrix ntau_n(nclass_u[d], nclass_v[d]);
      NumericVector ntau_d(nclass_v[d]);
      ptr_ntau_n[d] = ntau_n.begin();
      ptr_ntau_d[d] = ntau_d.begin();
      lst_ntau_n[d] = ntau_n;
      lst_ntau_d[d] = ntau_d;
   }

   List rhoList = init_param["rho"];
   for (int v = 0; v < nleaf_unique; v ++) {
      NumericMatrix rho_init = rhoList[v];
      NumericMatrix lrho = clone(rho_init);
      ptr_rho[v] = lrho.begin();
      lst_rho[v] = lrho;
   }
   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector denom(nclass_leaf[v]);
      NumericMatrix numer(sum(ncatv), nclass_leaf[v]);
      ptr_nrho_d[v] = denom.begin();
      ptr_nrho_n[v] = numer.begin();
      lst_nrho_d[v] = denom;
      lst_nrho_n[v] = numer;
   }

   for (int r = 0; r < nroot; r ++) {
      NumericVector ll(nobs);
      lst_ll[r] = ll;
      ptr_ll[r] = ll.begin();
   }

   for (int d = 0; d < nlink; d ++) {
      int u = ulv[d]; int v = vlv[d];
      NumericMatrix joint(nclass[u] * nclass[v], nobs);
      ptr_joint[d] = joint.begin();
      lst_joint[d] = joint;
      NumericMatrix jlambda(nclass[v], nobs);
      ptr_j[d] = jlambda.begin();
      lst_j[d] = jlambda;
   }

   for (int u = 0; u < nlv; u ++) {
      NumericMatrix post(nclass[u], nobs);
      NumericMatrix alpha(nclass[u], nobs);
      NumericMatrix lambda(nclass[u], nobs);
      ptr_post[u] = post.begin();
      ptr_a[u] = alpha.begin();
      ptr_l[u] = lambda.begin();
      lst_post[u] = post;
      lst_a[u] = alpha;
      lst_l[u] = lambda;
   }

   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;

   while ( (iter < max_iter) && (dll > tol) ) {
      iter ++;
      lastll = currll;

      // alpha, lambda initiate
      for (int v = 0; v < nlv; v ++) {
         NumericMatrix alpha = lst_a[v];
         NumericMatrix lambda = lst_l[v];
         alpha.fill(R_NegInf);
         lambda.fill(0);
      }

      // (expectation-step)
      // initiate lambda
      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         upInit(py, ptr_rho[cstr_leaf[v]], ptr_l[leaf[v]], nclass[leaf[v]],
                nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
         py += nobs * nvar[cstr_leaf[v]];
      }

      // upward recursion
      for (int d = nlink - 1; d > -1; d --) {
         int u = ulv[d], v = vlv[d];
         upRec(ptr_l[v], ptr_j[d], ptr_l[u], ptr_tau[cstr_link[d]],
               nobs, nclass[u], nclass[v], reg);
      }

      // initiate alpha
      for (int r = 0; r < nroot; r ++) {
         dnInit(ptr_a[root[r]], ptr_l[root[r]], ptr_pi[r],
                ptr_post[root[r]], ptr_ll[r], nobs,
                nclass[root[r]], reg);
      }

      // Downward recursion
      for (int d = 0; d < nlink; d ++) {
         int u = ulv[d];
         int v = vlv[d];
         dnRec(ptr_a[u], ptr_a[v], ptr_l[u], ptr_l[v], ptr_j[d],
               nobs, nclass[u], nclass[v], ptr_tau[cstr_link[d]],
               ptr_post[u], ptr_joint[d], ptr_ll[tree_index[d]], reg);
      }

      // (maximization-step)
      // pi updates
      for (int r = 0; r < nroot; r ++) {
         updatePi(ptr_pi[r], ptr_post[root[r]], nobs,
                  nclass[root[r]]);
      }

      restr_ = restr0.begin();
      // tau updates
      for (int d = 0; d < nlink_unique; d ++) {
         NumericMatrix numer = lst_ntau_n[d];
         NumericVector denom = lst_ntau_d[d];
         numer.fill(R_NegInf);
         denom.fill(R_NegInf);
      }
      for (int d = 0; d < nlink; d ++) {
         int u = ulv[d]; int v = vlv[d];
         cumTau(ptr_ntau_n[cstr_link[d]], ptr_ntau_d[cstr_link[d]],
                ptr_joint[d], nobs, nclass[u], nclass[v]);
      }
      for (int d = 0; d < nlink_unique; d ++) {
         updateTau(ptr_tau[d], ptr_ntau_n[d], ptr_ntau_d[d],
                   nclass_u[d], nclass_v[d], restr_);
         restr_ += nclass_u[d] * nclass_v[d];
      }

      // rho updates
      py = y.begin();
      for (int v = 0; v < nleaf_unique; v ++) {
         NumericMatrix numer = lst_nrho_n[v];
         NumericVector denom = lst_nrho_d[v];
         numer.fill(R_NegInf);
         denom.fill(R_NegInf);
      }
      for (int v = 0; v < nleaf; v ++) {
         int u = leaf[v], cu = cstr_leaf[v];
         IntegerVector ncatv = ncat[cu];
         cumRho(ptr_nrho_n[cu], ptr_nrho_d[cu], py, nobs, nvar[cu],
                ncatv, nclass[u], ptr_post[u], ptr_rho[cu]);
         py += nobs * nvar[cu];
      }

      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         updateRho(ptr_rho[v], ptr_nrho_n[v], ptr_nrho_d[v],
                   nobs, nclass_leaf[v], nvar[v], ncatv, restr_);
         restr_ += nclass_leaf[v] * sum(ncatv);
      }

      currll = 0;
      for (int r = 0; r < nroot; r ++) {
         double *ll = ptr_ll[r];
         for (int i = 0; i < nobs; i ++) {
            currll += ll[i];
         }
      }

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


   List par;
   par["pi"] = lst_pi;
   par["tau"] = lst_tau;
   par["rho"] = lst_rho;

   List res;
   res["params"] = par;
   res["converged"] = dll < tol;
   res["niter"] = iter;

   return res;
}



// [[Rcpp::export]]
double fll(
      NumericVector logit, IntegerVector y, int nobs, IntegerVector nvar, List ncat,
      int nlv, int nr, int nrl, int nlf, int ntau, int nrho, int nprob,
      IntegerVector rl, IntegerVector ul, IntegerVector vl, IntegerVector lf,
      IntegerVector eqrl, IntegerVector eqlf,
      IntegerVector nc, IntegerVector nk, IntegerVector nl, IntegerVector ncl,
      IntegerVector nc_pi, IntegerVector nk_tau, IntegerVector nl_tau, IntegerVector nc_rho,
      LogicalVector restr0, IntegerVector ref, bool reg = false
) {
   int *py;
   NumericVector par(nprob);
   std::vector<double*> _pi_(nr);
   std::vector<double*> _tau_(ntau);
   std::vector<double*> _rho_(nrho);

   double *_par_ = par.begin();
   double *_logit_ = logit.begin();
   int *_restr_ = restr0.begin();
   int *_ref_ = ref.begin();

   double denom = 0.0;
   for (int r = 0; r < nr; r ++) {
      _pi_[r] = _par_;
      denom = 1.0;
      for (int k = 0; k < nc_pi[r]; k ++) {
         if (k == _ref_[r]) continue;
         _par_[k] = *_logit_;
         denom += exp(*_logit_);
         _logit_ ++;
      }
      for (int k = 0; k < nc_pi[r]; k ++) {
         _par_[k] -= log(denom);
      }
   }
   _ref_ += nr;

   for (int d = 0; d < ntau; d ++) {
      _tau_[d] = _par_;
      for (int l = 0; l < nl_tau[d]; l ++) {
         denom = 1.0;
         for (int k = 0; k < nk_tau[d]; k ++) {
            if (_restr_[k] | (k == _ref_[l])) continue;
            _par_[k] = *_logit_;
            denom += exp(*_logit_);
            _logit_ ++;
         }
         for (int k = 0; k < nk_tau[d]; k ++) {
            if (_restr_[k]) _par_[k] = R_NegInf;
            else _par_[k] -= log(denom);
         }
         _par_ += nk_tau[d];
         _restr_ += nk_tau[d];
      }
      _ref_ += nl_tau[d];
   }

   for (int v = 0; v < nrho; v ++) {
      _rho_[v] = _par_;
      IntegerVector ncatv = ncat[v];
      for (int k = 0; k < nc_rho[v]; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            denom = 1;
            for (int r = 0; r < ncatv[m] - 1; r ++) {
               if (_restr_[r] | (r == _ref_[m])) continue;
               _par_[r] = *_logit_;
               denom += exp(*_logit_);
               _logit_ ++;
            }
            for (int r = 0; r < ncatv[m]; r ++) {
               if (_restr_[r]) _par_[r] = R_NegInf;
               else _par_[r] -= log(denom);
            }
            _par_  += ncatv[m];
            _restr_ += ncatv[m];
         }
         _ref_ += nvar[v];
      }
   }

   NumericVector lambda(nobs * sum(nc));
   std::vector<double*> _l_(nlv);
   double * _tmp_ = lambda.begin();
   for (int v = 0; v < nlv; v ++) {
      _l_[v] = _tmp_;
      _tmp_ += nc[v] * nobs;
   }

   double ll = 0.0;
   // initiate lambda
   py = y.begin();
   for (int v = 0; v < nlf; v ++) {
      upInit(py, _rho_[eqlf[v]], _l_[lf[v]], ncl[v],
             nobs, nvar[eqlf[v]], ncat[eqlf[v]]);
      py += nobs * nvar[eqlf[v]];
   }

   // upward recursion
   for (int d = nrl - 1; d > -1; d --) {
      upRec2(_l_[vl[d]], _l_[ul[d]], _tau_[eqrl[d]],
             nobs, nk[d], nl[d], false);
   }

   for (int r = 0; r < nr; r ++) {
      ll += calcll(_l_[rl[r]], _pi_[r], nobs, nc_pi[r], false);
   }

   return -ll;
}




// [[Rcpp::export]]
double floglik(
      NumericVector logit,
      IntegerVector y, int nobs, IntegerVector nvar, List ncat,
      int nlv, int nroot, int nlink, int nlink_unique,
      int nleaf, int nleaf_unique, IntegerVector root,
      IntegerVector ulv, IntegerVector vlv, IntegerVector cstr_link,
      IntegerVector leaf, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v,
      LogicalVector restr0, IntegerVector ref, bool reg = false
) {
   int *py;
   List lst_pi(nroot);
   List lst_tau(nlink_unique);
   List lst_rho(nleaf_unique);
   std::vector<double*> ptr_pi(nroot);
   std::vector<double*> ptr_tau(nlink_unique);
   std::vector<double*> ptr_rho(nleaf_unique);

   List lst_l(nlv);
   std::vector<double*> ptr_l(nlv);

   double *logit_ = logit.begin();
   int *restr_ = restr0.begin();
   int *ref_ = ref.begin();

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass[root[r]];
      NumericVector pi(nk);
      double denom = 1;
      for (int k = 0; k < nk; k ++) {
         if (k == ref_[r]) continue;
         pi[k] = *logit_;
         denom += exp(*logit_);
         logit_ ++;
      }
      for (int k = 0; k < nk; k ++) {
         pi[k] -= log(denom);
      }

      lst_pi[r] = pi;
      ptr_pi[r] = pi.begin();
   }
   ref_ += nroot;

   double denom = 0.0;
   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d];
      int nl = nclass_v[d];
      NumericMatrix tau(nk, nl);
      double *tau_ = tau.begin();
      for (int l = 0; l < nl; l ++) {
         denom = 1;
         for (int k = 0; k < nk; k ++) {
            if (restr_[k] | (k == ref[l])) continue;
            tau_[k] = *logit_;
            denom += exp(*logit_);
            logit_ ++;
         }
         for (int k = 0; k < nk; k ++) {
            if (restr_[k]) tau_[k] = R_NegInf;
            else tau_[k] -= log(denom);
         }
         tau_ += nk;
         restr_ += nk;
      }
      ref_   += nl;

      lst_tau[d] = tau;
      ptr_tau[d] = tau.begin();
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
               if (restr_[r] | (r == ref[m])) continue;
               rho_[r] = *logit_;
               denom += exp(*logit_);
               logit_ ++;
            }
            for (int r = 0; r < ncatv[m]; r ++) {
               if (restr_[r]) rho_[r] = R_NegInf;
               else rho_[r] -= log(denom);
            }
            rho_  += ncatv[m];
            restr_ += ncatv[m];
         }
         ref_ += ncat.length();
      }
      lst_rho[v] = rho;
      ptr_rho[v] = rho.begin();
   }

   for (int v = 0; v < nlv; v ++) {
      NumericMatrix lambda(nclass[v], nobs);
      ptr_l[v] = lambda.begin();
      lst_l[v] = lambda;
   }

   double ll = 0.0;
   // initiate lambda
   py = y.begin();
   for (int v = 0; v < nleaf; v ++) {
      upInit(py, ptr_rho[cstr_leaf[v]], ptr_l[leaf[v]],
             nclass[leaf[v]], nobs, nvar[cstr_leaf[v]],
             ncat[cstr_leaf[v]]);
      py += nobs * nvar[cstr_leaf[v]];
   }

   // upward recursion
   for (int d = nlink - 1; d > -1; d --) {
      int u = ulv[d];
      int v = vlv[d];
      upRec2(ptr_l[v], ptr_l[u], ptr_tau[cstr_link[d]],
             nobs, nclass[u], nclass[v], reg);
   }

   for (int r = 0; r < nroot; r ++) {
      ll += calcll(ptr_l[root[r]], ptr_pi[r],
                   nobs, nclass[root[r]], reg);
   }

   return -ll;
}

