#include <Rcpp.h>
#include "auxiliary.h"
using namespace Rcpp;

void dnInit(
      double *alpha, double *lambda, double *pi, double *post,
      double *ll, int nobs, int nclass, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      ll[i] = R_NegInf;
      for (int k = 0; k < nclass; k ++) {
         alpha[k] = pi[k];
         post[k] = alpha[k] + lambda[k];
         ll[i] = log_add_exp(ll[i], post[k]);
      }
      for (int k = 0; k < nclass; k ++) {
         post[k] -= ll[i];
      }

      alpha += nclass;
      lambda += nclass;
      post += nclass;
      if (reg) pi += nclass;
   }
}

double calcll(
      double *lambda, double *pi,
      int nobs, int nclass, bool reg
) {
   double ll = 0;
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         lik = log_add_exp(lik, pi[k] + lambda[k]);
      }
      ll += lik;

      lambda += nclass;
      if (reg) pi += nclass;
   }
   return ll;
}

void calclli(
      double *lambda, double *pi, double *ll,
      int nobs, int nclass, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++) {
         ll[i] = log_add_exp(ll[i], pi[k] + lambda[k]);
      }
      if (reg) pi += nclass;

      lambda += nclass;
   }
}

void dnRec(
      double *alpha, double *ualpha,
      double *lambda, double *ulambda, double *jlambda,
      int nobs, int nk, int nl, double *tau,
      double *post, double *joint, double *ll, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         for (int l = 0; l < nl; l ++) {
            val = tau[k + l * nk] + ualpha[l] + ulambda[l] - jlambda[l];
            joint[k + l * nk] = val + lambda[k] - ll[i];
            alpha[k] = log_add_exp(alpha[k], val);
         }
         post[k] = alpha[k] + lambda[k] - ll[i];
      }

      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      lambda += nk; ulambda += nl; jlambda += nl;
      if (reg) tau += nk * nl;
   }
}
