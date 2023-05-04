// #include <Rcpp.h>
// #include "auxiliary.h"
// #include "transform.h"
// #include "param.h"
// #include "upward.h"
// #include "downward.h"
// #include "mstep.h"
// using namespace Rcpp;
//
// double pldiri(NumericVector lx, NumericVector p) {
//    return lgammaf(sum(p)) - sum(lgamma(p)) + sum((p - 1) * lx);
// }
//
// double pldiri2(NumericVector lx, double p) {
//    return lgammaf(lx.length() * p) - lx.length() * lgammaf(p) +
//       sum((p - 1) * lx);
// }
//
// NumericVector eldiri(NumericVector x, bool a) {
//    if (a) return digamma(x) - R::digamma(sum(x));
//    else return log(x) - log(sum(x));
// }
//
// // [[Rcpp::export]]
// List vemfit(
//       IntegerVector y, int nobs, IntegerVector nvar, List ncat,
//       int nlv, int nroot, int nlink, int nleaf,
//       int nlink_unique, int nleaf_unique,
//       IntegerVector tree_index, IntegerVector root,
//       IntegerVector ulv, IntegerVector vlv, IntegerVector leaf,
//       IntegerVector cstr_link, IntegerVector cstr_leaf,
//       IntegerVector nclass, IntegerVector nclass_leaf,
//       IntegerVector nclass_u, IntegerVector nclass_v,
//       List init_param, LogicalVector restr0,
//       int max_iter, double tol,
//       bool verbose, int periter = 100
// ) {
//    int *py;
//
//    List lst_pi(nroot), lst_tau(nlink_unique), lst_rho(nleaf_unique);
//    List lst_ntau(nlink_unique);
//    List lst_nrho_d(nleaf_unique), lst_nrho_n(nleaf_unique);
//    std::vector<double*> ptr_pi(nroot);
//    std::vector<double*> ptr_tau(nlink_unique), ptr_ntau(nlink_unique);
//    std::vector<double*> ptr_rho(nleaf_unique);
//    std::vector<double*> ptr_nrho_d(nleaf_unique), ptr_nrho_n(nleaf_unique);
//
//    List lst_ll(nroot);
//    List lst_a(nlv), lst_l(nlv), lst_j(nlink);
//    std::vector<double*> ptr_ll(nroot);
//    std::vector<double*> ptr_a(nlv), ptr_l(nlv), ptr_j(nlink);
//
//    List lst_post(nlv), lst_joint(nlink);
//    std::vector<double*> ptr_post(nlv), ptr_joint(nlink);
//
//    int *restr_;
//
//    List piList = init_param["pi"];
//    for (int r = 0; r < nroot; r ++) {
//       NumericVector pi_init = piList[r];
//       NumericVector pi = clone(pi_init);
//       ptr_pi[r] = pi.begin();
//       lst_pi[r] = pi;
//    }
//
//    List tauList = init_param["tau"];
//    for (int d = 0; d < nlink_unique; d ++) {
//       NumericMatrix tau_init = tauList[d];
//       NumericMatrix tau = clone(tau_init);
//       ptr_tau[d] = tau.begin();
//       lst_tau[d] = tau;
//    }
//    for (int d = 0; d < nlink_unique; d ++) {
//       NumericMatrix ntau(nclass_u[d], nclass_v[d]);
//       ptr_ntau[d] = ntau.begin();
//       lst_ntau[d] = ntau;
//    }
//
//    List rhoList = init_param["rho"];
//    for (int v = 0; v < nleaf_unique; v ++) {
//       NumericMatrix rho_init = rhoList[v];
//       NumericMatrix lrho = clone(rho_init);
//       ptr_rho[v] = lrho.begin();
//       lst_rho[v] = lrho;
//    }
//    for (int v = 0; v < nleaf_unique; v ++) {
//       IntegerVector ncatv = ncat[v];
//       NumericVector denom(nclass_leaf[v] * nvar[v]);
//       NumericVector numer(nclass_leaf[v] * sum(ncatv));
//       ptr_nrho_d[v] = denom.begin();
//       ptr_nrho_n[v] = numer.begin();
//       lst_nrho_d[v] = denom;
//       lst_nrho_n[v] = numer;
//    }
//
//    for (int r = 0; r < nroot; r ++) {
//       NumericVector ll(nobs);
//       lst_ll[r] = ll;
//       ptr_ll[r] = ll.begin();
//    }
//
//    for (int d = 0; d < nlink; d ++) {
//       int u = ulv[d]; int v = vlv[d];
//       NumericMatrix joint(nclass[u] * nclass[v], nobs);
//       ptr_joint[d] = joint.begin();
//       lst_joint[d] = joint;
//       NumericMatrix jlambda(nclass[v], nobs);
//       ptr_j[d] = jlambda.begin();
//       lst_j[d] = jlambda;
//    }
//
//    for (int u = 0; u < nlv; u ++) {
//       NumericMatrix post(nclass[u], nobs);
//       NumericMatrix alpha(nclass[u], nobs);
//       NumericMatrix lambda(nclass[u], nobs);
//       ptr_post[u] = post.begin();
//       ptr_a[u] = alpha.begin();
//       ptr_l[u] = lambda.begin();
//       lst_post[u] = post;
//       lst_a[u] = alpha;
//       lst_l[u] = lambda;
//    }
//
//    int iter = 0;
//    double currll = R_NegInf;
//    double lastll = R_NegInf;
//    double dll = R_PosInf;
//
//    while ( (iter < max_iter) && (dll > tol) ) {
//       iter ++;
//       lastll = currll;
//
//       // lambda, cleaning
//       for (int v = 0; v < nlv; v ++) {
//          NumericMatrix lambda = lst_l[v];
//          lambda.fill(0);
//       }
//
//       // (expectation-step)
//       // initiate lambda
//       py = y.begin();
//       for (int v = 0; v < nleaf; v ++) {
//          upInit(py, ptr_rho[cstr_leaf[v]],
//                 ptr_l[leaf[v]], nclass[leaf[v]],
//                 nobs, nvar[cstr_leaf[v]],
//                 ncat[cstr_leaf[v]]);
//          py += nobs * nvar[cstr_leaf[v]];
//       }
//
//       // upward recursion
//       for (int d = nlink - 1; d > -1; d --) {
//          int u = ulv[d];
//          int v = vlv[d];
//          upRec(ptr_l[v], ptr_j[d], ptr_l[u], ptr_tau[cstr_link[d]],
//                nobs, nclass[u], nclass[v]);
//       }
//
//       // initiate alpha
//       for (int r = 0; r < nroot; r ++) {
//          dnInit(ptr_a[root[r]], ptr_l[root[r]],
//                 ptr_pi[r], ptr_post[root[r]],
//                 ptr_ll[r], nobs, nclass[root[r]]);
//       }
//
//       // Downward recursion
//       for (int d = 0; d < nlink; d ++) {
//          int u = ulv[d];
//          int v = vlv[d];
//          dnRec(ptr_a[u], ptr_a[v], ptr_l[u], ptr_l[v], ptr_j[d],
//                nobs, nclass[u], nclass[v], ptr_tau[cstr_link[d]],
//                ptr_post[u], ptr_joint[d], ptr_ll[tree_index[d]]);
//       }
//
//       // (maximization-step)
//       // pi updates
//       for (int r = 0; r < nroot; r ++) {
//          updatePi(ptr_pi[r], ptr_post[root[r]], nobs,
//                   nclass[root[r]]);
//       }
//
//       restr_ = restr0.begin();
//       // tau updates
//       for (int d = 0; d < nlink_unique; d ++) {
//          NumericVector ntau = lst_ntau[d];
//          ntau.fill(0);
//       }
//       for (int d = 0; d < nlink; d ++) {
//          int u = ulv[d]; int v = vlv[d];
//          cumTau(ptr_joint[d], ptr_ntau[cstr_link[d]], nobs,
//                 nclass[u], nclass[v]);
//       }
//       for (int d = 0; d < nlink_unique; d ++) {
//          updateTau(ptr_tau[d], ptr_ntau[d],
//                    nclass_u[d], nclass_v[d], restr_);
//          restr_ += nclass_u[d] * nclass_v[d];
//       }
//
//       // rho updates
//       py = y.begin();
//       for (int v = 0; v < nleaf_unique; v ++) {
//          NumericVector denom = lst_nrho_d[v];
//          NumericVector numer = lst_nrho_n[v];
//          denom.fill(0);
//          numer.fill(0);
//       }
//       for (int v = 0; v < nleaf; v ++) {
//          int u = leaf[v], cu = cstr_leaf[v];
//          IntegerVector ncatv = ncat[cu];
//          cumRho(ptr_nrho_d[cu], ptr_nrho_n[cu], py, nobs, nvar[cu],
//                 ncatv, nclass[u], ptr_post[u], ptr_rho[cu]);
//          py += nobs * nvar[cu];
//       }
//       for (int v = 0; v < nleaf_unique; v ++) {
//          IntegerVector ncatv = ncat[v];
//          updateRho(ptr_rho[v], ptr_nrho_n[v], ptr_nrho_d[v],
//                    nobs, nclass_leaf[v], nvar[v], ncatv, restr_);
//          restr_ += nclass_leaf[v] * sum(ncatv);
//       }
//
//       currll = 0;
//       for (int r = 0; r < nroot; r ++) {
//          double *ll = ptr_ll[r];
//          for (int i = 0; i < nobs; i ++) {
//             currll += ll[i];
//          }
//       }
//
//       if (lastll == R_NegInf) dll = R_PosInf;
//       else dll = currll - lastll;
//       if (verbose) {
//          Rcout << iter << " iterations  logLik: " <<
//             std::fixed << std::setprecision(2) << currll << "  diff: ";
//          if (dll > 1e-5) Rcout << std::setprecision(5) << dll;
//          else Rcout << std::scientific << std::setprecision(1) << dll;
//          Rcout << "      \r" << std::flush;
//
//          if (iter % periter == 0) {
//             Rcout << "\n" << std::flush;
//          }
//       }
//    }
//    if (verbose) {
//       if (iter % periter != 0) Rcout << "\n";
//    }
//
//
//    List par;
//    par["pi"] = lst_pi;
//    par["tau"] = lst_tau;
//    par["rho"] = lst_rho;
//
//    List res;
//    res["params"] = par;
//    res["converged"] = dll < tol;
//    res["niter"] = iter;
//
//    return res;
// }
//
//
//
// // [[Rcpp::export]]
// double felbo(
//       NumericVector logit,
//       IntegerVector y, int nobs, IntegerVector nvar, List ncat,
//       int nlv, int nroot, int nlink, int nlink_unique,
//       int nleaf, int nleaf_unique, IntegerVector root,
//       IntegerVector ulv, IntegerVector vlv, IntegerVector cstr_link,
//       IntegerVector leaf, IntegerVector cstr_leaf,
//       IntegerVector nclass, IntegerVector nclass_leaf,
//       IntegerVector nclass_u, IntegerVector nclass_v,
//       LogicalVector restr0, IntegerVector ref
// ) {
//    int *py;
//    List lst_pi(nroot);
//    List lst_tau(nlink_unique);
//    List lst_rho(nleaf_unique);
//    std::vector<double*> ptr_pi(nroot);
//    std::vector<double*> ptr_tau(nlink_unique);
//    std::vector<double*> ptr_rho(nleaf_unique);
//
//    List lst_l(nlv);
//    std::vector<double*> ptr_l(nlv);
//
//    double *logit_ = logit.begin();
//    int *restr_ = restr0.begin();
//    int *ref_ = ref.begin();
//
//    for (int r = 0; r < nroot; r ++) {
//       int nk = nclass[root[r]];
//       NumericVector pi(nk);
//       double denom = 1;
//       for (int k = 0; k < nk; k ++) {
//          if (k == ref_[r]) continue;
//          pi[k] = *logit_;
//          denom += exp(*logit_);
//          logit_ ++;
//       }
//       for (int k = 0; k < nk; k ++) {
//          pi[k] -= log(denom);
//       }
//
//       lst_pi[r] = pi;
//       ptr_pi[r] = pi.begin();
//    }
//    ref_ += nroot;
//
//    double denom = 0;
//    for (int d = 0; d < nlink_unique; d ++) {
//       int nk = nclass_u[d];
//       int nl = nclass_v[d];
//       NumericMatrix tau(nk, nl);
//       double *tau_ = tau.begin();
//       for (int l = 0; l < nl; l ++) {
//          denom = 1;
//          for (int k = 0; k < nk; k ++) {
//             if (restr_[k] | (k == ref[l])) continue;
//             tau_[k] = *logit_;
//             denom += exp(*logit_);
//             logit_ ++;
//          }
//          for (int k = 0; k < nk; k ++) {
//             if (restr_[k]) tau_[k] = R_NegInf;
//             else tau_[k] -= log(denom);
//          }
//          tau_ += nk;
//          restr_ += nk;
//       }
//       ref_   += nl;
//
//       lst_tau[d] = tau;
//       ptr_tau[d] = tau.begin();
//    }
//
//    for (int v = 0; v < nleaf_unique; v ++) {
//       IntegerVector ncatv = ncat[v];
//       int nk = nclass_leaf[v];
//       NumericVector rho(nk * sum(ncatv));
//       double *rho_ = rho.begin();
//       for (int k = 0; k < nk; k ++) {
//          for (int m = 0; m < ncatv.length(); m ++) {
//             denom = 1;
//             for (int r = 0; r < ncatv[m] - 1; r ++) {
//                if (restr_[r] | (r == ref[m])) continue;
//                rho_[r] = *logit_;
//                denom += exp(*logit_);
//                logit_ ++;
//             }
//             for (int r = 0; r < ncatv[m]; r ++) {
//                if (restr_[r]) rho_[r] = R_NegInf;
//                else rho_[r] -= log(denom);
//             }
//             rho_  += ncatv[m];
//             restr_ += ncatv[m];
//          }
//          ref_ += ncat.length();
//       }
//       lst_rho[v] = rho;
//       ptr_rho[v] = rho.begin();
//    }
//
//    for (int v = 0; v < nlv; v ++) {
//       NumericMatrix lambda(nclass[v], nobs);
//       ptr_l[v] = lambda.begin();
//       lst_l[v] = lambda;
//    }
//
//    double ll = 0;
//    // initiate lambda
//    py = y.begin();
//    for (int v = 0; v < nleaf; v ++) {
//       upInit(py, ptr_rho[cstr_leaf[v]], ptr_l[leaf[v]],
//              nclass[leaf[v]], nobs, nvar[cstr_leaf[v]],
//                                         ncat[cstr_leaf[v]]);
//       py += nobs * nvar[cstr_leaf[v]];
//    }
//
//    // upward recursion
//    for (int d = nlink - 1; d > -1; d --) {
//       int u = ulv[d];
//       int v = vlv[d];
//       upRec2(ptr_l[v], ptr_l[u], ptr_tau[cstr_link[d]],
//              nobs, nclass[u], nclass[v]);
//    }
//
//    for (int r = 0; r < nroot; r ++) {
//       ll += calcll(ptr_l[root[r]], ptr_pi[r],
//                    nobs, nclass[root[r]]);
//    }
//
//    return -ll;
// }
//
//
