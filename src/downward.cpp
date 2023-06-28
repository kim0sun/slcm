#include <Rcpp.h>
#include "auxiliary.h"
using namespace Rcpp;

void dnInit(
      double *alpha, double *lambda, double *pi, double *post,
      double *lr, int nobs, int nclass, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      lr[i] = R_NegInf;
      for (int k = 0; k < nclass; k ++) {
         alpha[k] = pi[k];
         post[k] = alpha[k] + lambda[k];
         lr[i] = log_add_exp(lr[i], post[k]);
      }
      for (int k = 0; k < nclass; k ++) {
         post[k] -= lr[i];
      }

      alpha += nclass;
      lambda += nclass;
      post += nclass;
      if (reg) pi += nclass;
   }
}

double calclr(
      double *lambda, double *pi,
      int nobs, int nclass, bool reg
) {
   double lr = 0;
   for (int i = 0; i < nobs; i ++) {
      double lik = R_NegInf;
      for (int k = 0; k < nclass; k ++) {
         lik = log_add_exp(lik, pi[k] + lambda[k]);
      }
      lr += lik;

      lambda += nclass;
      if (reg) pi += nclass;
   }
   return lr;
}

void calclri(
      double *lambda, double *pi, double *lr,
      int nobs, int nclass, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      lr[i] = R_NegInf;
      for (int k = 0; k < nclass; k ++) {
         lr[i] = log_add_exp(lr[i], pi[k] + lambda[k]);
      }
      if (reg) pi += nclass;

      lambda += nclass;
   }
}

void dnRec(
      double *alpha, double *ualpha,
      double *lambda, double *ulambda, double *jlambda,
      int nobs, int nk, int nl, double *tau,
      double *post, double *joint, double *lr, bool reg
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         alpha[k] = R_NegInf;
         for (int l = 0; l < nl; l ++) {
            val = tau[k + l * nk] + ualpha[l] + ulambda[l] - jlambda[l];
            joint[k + l * nk] = val + lambda[k] - lr[i];
            alpha[k] = log_add_exp(alpha[k], val);
         }
         post[k] = alpha[k] + lambda[k] - lr[i];
      }

      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      lambda += nk; ulambda += nl; jlambda += nl;
      if (reg) tau += nk * nl;
   }
}
