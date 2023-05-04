// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_add_exp
double log_add_exp(double x, double y);
RcppExport SEXP _slcm_log_add_exp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(log_add_exp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calcModel
List calcModel(List param, IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nlink, int nleaf, int nlink_unique, int nleaf_unique, IntegerVector root, IntegerVector tree_index, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_link, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, List ref, bool reg);
RcppExport SEXP _slcm_calcModel(SEXP paramSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nlinkSEXP, SEXP nleafSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP tree_indexSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_linkSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP refSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink(nlinkSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tree_index(tree_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_link(cstr_linkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< List >::type ref(refSEXP);
    Rcpp::traits::input_parameter< bool >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(calcModel(param, y, nobs, nvar, ncat, nlv, nroot, nlink, nleaf, nlink_unique, nleaf_unique, root, tree_index, ulv, vlv, leaf, cstr_link, cstr_leaf, nclass, nclass_leaf, nclass_u, nclass_v, ref, reg));
    return rcpp_result_gen;
END_RCPP
}
// emFit
List emFit(IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nlink, int nleaf, int nlink_unique, int nleaf_unique, IntegerVector tree_index, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector leaf, IntegerVector cstr_link, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, List init_param, LogicalVector restr0, int max_iter, double tol, bool verbose, int newiter, bool reg);
RcppExport SEXP _slcm_emFit(SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nlinkSEXP, SEXP nleafSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP tree_indexSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP leafSEXP, SEXP cstr_linkSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP init_paramSEXP, SEXP restr0SEXP, SEXP max_iterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP newiterSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink(nlinkSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tree_index(tree_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_link(cstr_linkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< List >::type init_param(init_paramSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type restr0(restr0SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type newiter(newiterSEXP);
    Rcpp::traits::input_parameter< bool >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(emFit(y, nobs, nvar, ncat, nlv, nroot, nlink, nleaf, nlink_unique, nleaf_unique, tree_index, root, ulv, vlv, leaf, cstr_link, cstr_leaf, nclass, nclass_leaf, nclass_u, nclass_v, init_param, restr0, max_iter, tol, verbose, newiter, reg));
    return rcpp_result_gen;
END_RCPP
}
// floglik
double floglik(NumericVector logit, IntegerVector y, int nobs, IntegerVector nvar, List ncat, int nlv, int nroot, int nlink, int nlink_unique, int nleaf, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector cstr_link, IntegerVector leaf, IntegerVector cstr_leaf, IntegerVector nclass, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, LogicalVector restr0, IntegerVector ref, bool reg);
RcppExport SEXP _slcm_floglik(SEXP logitSEXP, SEXP ySEXP, SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP nrootSEXP, SEXP nlinkSEXP, SEXP nlink_uniqueSEXP, SEXP nleafSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP cstr_linkSEXP, SEXP leafSEXP, SEXP cstr_leafSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP restr0SEXP, SEXP refSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logit(logitSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink(nlinkSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_link(cstr_linkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type restr0(restr0SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< bool >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(floglik(logit, y, nobs, nvar, ncat, nlv, nroot, nlink, nlink_unique, nleaf, nleaf_unique, root, ulv, vlv, cstr_link, leaf, cstr_leaf, nclass, nclass_leaf, nclass_u, nclass_v, restr0, ref, reg));
    return rcpp_result_gen;
END_RCPP
}
// pi_gnr
NumericVector pi_gnr(int nk);
RcppExport SEXP _slcm_pi_gnr(SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(pi_gnr(nk));
    return rcpp_result_gen;
END_RCPP
}
// tau_gnr
NumericMatrix tau_gnr(int nk, int nl);
RcppExport SEXP _slcm_tau_gnr(SEXP nkSEXP, SEXP nlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< int >::type nl(nlSEXP);
    rcpp_result_gen = Rcpp::wrap(tau_gnr(nk, nl));
    return rcpp_result_gen;
END_RCPP
}
// rho_gnr
NumericMatrix rho_gnr(int nk, IntegerVector ncat);
RcppExport SEXP _slcm_rho_gnr(SEXP nkSEXP, SEXP ncatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ncat(ncatSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_gnr(nk, ncat));
    return rcpp_result_gen;
END_RCPP
}
// par_gnr
List par_gnr(int nobs, IntegerVector nvar, List ncat, int nroot, int nlink_unique, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, LogicalVector init, List init_param);
RcppExport SEXP _slcm_par_gnr(SEXP nobsSEXP, SEXP nvarSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclassSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP initSEXP, SEXP init_paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< List >::type init_param(init_paramSEXP);
    rcpp_result_gen = Rcpp::wrap(par_gnr(nobs, nvar, ncat, nroot, nlink_unique, nleaf_unique, root, ulv, vlv, nclass, nclass_leaf, nclass_u, nclass_v, init, init_param));
    return rcpp_result_gen;
END_RCPP
}
// calcfreq
List calcfreq(IntegerVector mis, IntegerVector nrep, int nmis, IntegerVector freq, IntegerVector xobs, int nc, int N, double tol, int max_iter);
RcppExport SEXP _slcm_calcfreq(SEXP misSEXP, SEXP nrepSEXP, SEXP nmisSEXP, SEXP freqSEXP, SEXP xobsSEXP, SEXP ncSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type mis(misSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type nmis(nmisSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type xobs(xobsSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(calcfreq(mis, nrep, nmis, freq, xobs, nc, N, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// root_gnr
IntegerVector root_gnr(int nobs, int nk, Nullable<NumericVector> prob);
RcppExport SEXP _slcm_root_gnr(SEXP nobsSEXP, SEXP nkSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(root_gnr(nobs, nk, prob));
    return rcpp_result_gen;
END_RCPP
}
// cls_gnr
IntegerVector cls_gnr(int nobs, int nk, int nl, IntegerVector v, Nullable<NumericMatrix> prob);
RcppExport SEXP _slcm_cls_gnr(SEXP nobsSEXP, SEXP nkSEXP, SEXP nlSEXP, SEXP vSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< int >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(cls_gnr(nobs, nk, nl, v, prob));
    return rcpp_result_gen;
END_RCPP
}
// y_gnr
IntegerMatrix y_gnr(int nobs, int nk, IntegerVector ncat, IntegerVector cls, Nullable<NumericMatrix> prob);
RcppExport SEXP _slcm_y_gnr(SEXP nobsSEXP, SEXP nkSEXP, SEXP ncatSEXP, SEXP clsSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cls(clsSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(y_gnr(nobs, nk, ncat, cls, prob));
    return rcpp_result_gen;
END_RCPP
}
// ysim
List ysim(int nsim, List ncat, int nlv, IntegerVector root, IntegerVector leaf, Nullable<IntegerVector> ulv, Nullable<IntegerVector> vlv, Nullable<IntegerVector> cstr_link, IntegerVector cstr_leaf, int nroot, int nleaf, int nlink, IntegerVector nclass, List pi, List tau, List rho, bool print_class);
RcppExport SEXP _slcm_ysim(SEXP nsimSEXP, SEXP ncatSEXP, SEXP nlvSEXP, SEXP rootSEXP, SEXP leafSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP cstr_linkSEXP, SEXP cstr_leafSEXP, SEXP nrootSEXP, SEXP nleafSEXP, SEXP nlinkSEXP, SEXP nclassSEXP, SEXP piSEXP, SEXP tauSEXP, SEXP rhoSEXP, SEXP print_classSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nlv(nlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type leaf(leafSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type cstr_link(cstr_linkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cstr_leaf(cstr_leafSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf(nleafSEXP);
    Rcpp::traits::input_parameter< int >::type nlink(nlinkSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass(nclassSEXP);
    Rcpp::traits::input_parameter< List >::type pi(piSEXP);
    Rcpp::traits::input_parameter< List >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< List >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type print_class(print_classSEXP);
    rcpp_result_gen = Rcpp::wrap(ysim(nsim, ncat, nlv, root, leaf, ulv, vlv, cstr_link, cstr_leaf, nroot, nleaf, nlink, nclass, pi, tau, rho, print_class));
    return rcpp_result_gen;
END_RCPP
}
// logit2log
List logit2log(NumericVector logit, List ncat, int nroot, int nlink_unique, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass_root, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, LogicalVector restr, IntegerVector ref);
RcppExport SEXP _slcm_logit2log(SEXP logitSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclass_rootSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP restrSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logit(logitSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_root(nclass_rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type restr(restrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(logit2log(logit, ncat, nroot, nlink_unique, nleaf_unique, root, ulv, vlv, nclass_root, nclass_leaf, nclass_u, nclass_v, restr, ref));
    return rcpp_result_gen;
END_RCPP
}
// log2logit
NumericVector log2logit(List param, int npar, List ncat, int nroot, int nlink_unique, int nleaf_unique, IntegerVector nclass_root, IntegerVector nclass_leaf, IntegerVector nclass_u, IntegerVector nclass_v, LogicalVector restr, IntegerVector ref);
RcppExport SEXP _slcm_log2logit(SEXP paramSEXP, SEXP nparSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP nclass_rootSEXP, SEXP nclass_leafSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP restrSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< int >::type npar(nparSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_root(nclass_rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type restr(restrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(log2logit(param, npar, ncat, nroot, nlink_unique, nleaf_unique, nclass_root, nclass_leaf, nclass_u, nclass_v, restr, ref));
    return rcpp_result_gen;
END_RCPP
}
// param2list
List param2list(NumericVector param, List ncat, int nroot, int nlink_unique, int nleaf_unique, IntegerVector root, IntegerVector ulv, IntegerVector vlv, IntegerVector nclass_root, IntegerVector nclass_u, IntegerVector nclass_v, IntegerVector nclass_leaf);
RcppExport SEXP _slcm_param2list(SEXP paramSEXP, SEXP ncatSEXP, SEXP nrootSEXP, SEXP nlink_uniqueSEXP, SEXP nleaf_uniqueSEXP, SEXP rootSEXP, SEXP ulvSEXP, SEXP vlvSEXP, SEXP nclass_rootSEXP, SEXP nclass_uSEXP, SEXP nclass_vSEXP, SEXP nclass_leafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< List >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< int >::type nroot(nrootSEXP);
    Rcpp::traits::input_parameter< int >::type nlink_unique(nlink_uniqueSEXP);
    Rcpp::traits::input_parameter< int >::type nleaf_unique(nleaf_uniqueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type root(rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ulv(ulvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vlv(vlvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_root(nclass_rootSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_u(nclass_uSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_v(nclass_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nclass_leaf(nclass_leafSEXP);
    rcpp_result_gen = Rcpp::wrap(param2list(param, ncat, nroot, nlink_unique, nleaf_unique, root, ulv, vlv, nclass_root, nclass_u, nclass_v, nclass_leaf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_slcm_log_add_exp", (DL_FUNC) &_slcm_log_add_exp, 2},
    {"_slcm_calcModel", (DL_FUNC) &_slcm_calcModel, 24},
    {"_slcm_emFit", (DL_FUNC) &_slcm_emFit, 28},
    {"_slcm_floglik", (DL_FUNC) &_slcm_floglik, 24},
    {"_slcm_pi_gnr", (DL_FUNC) &_slcm_pi_gnr, 1},
    {"_slcm_tau_gnr", (DL_FUNC) &_slcm_tau_gnr, 2},
    {"_slcm_rho_gnr", (DL_FUNC) &_slcm_rho_gnr, 2},
    {"_slcm_par_gnr", (DL_FUNC) &_slcm_par_gnr, 15},
    {"_slcm_calcfreq", (DL_FUNC) &_slcm_calcfreq, 9},
    {"_slcm_root_gnr", (DL_FUNC) &_slcm_root_gnr, 3},
    {"_slcm_cls_gnr", (DL_FUNC) &_slcm_cls_gnr, 5},
    {"_slcm_y_gnr", (DL_FUNC) &_slcm_y_gnr, 5},
    {"_slcm_ysim", (DL_FUNC) &_slcm_ysim, 17},
    {"_slcm_logit2log", (DL_FUNC) &_slcm_logit2log, 14},
    {"_slcm_log2logit", (DL_FUNC) &_slcm_log2logit, 12},
    {"_slcm_param2list", (DL_FUNC) &_slcm_param2list, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_slcm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}