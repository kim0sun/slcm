#include <Rcpp.h>
#include "auxiliary.h"
using namespace Rcpp;

void upInit(
      int *y, const double *ptr_rho, double *lambda,
      int nk, int nobs, int nvar, int* ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *rho = (double *) ptr_rho;
      for (int k = 0; k < nk; k ++) {
         for (int m = 0; m < nvar; m ++) {
            if (y[m] > 0)
               lambda[k] += rho[y[m] - 1];
            rho += ncat[m];
         }
      }
      y += nvar;
      lambda += nk;
   }
}

void upRec(
      double *lambda, double *jlambda, double *llambda,
      const double *tau, int nobs, int nk, int nl, bool reg
) {
   double ml = R_NegInf;
   double *tau_ = (double *) tau;
   for (int i = 0; i < nobs; i ++) {
      if (!reg) tau_ = (double *) tau;
      for (int l = 0; l < nl; l ++) {
         ml = R_NegInf;
         for (int k = 0; k < nk; k ++) {
            ml = log_add_exp(ml, tau_[k] + llambda[k]);
         }

         tau_ += nk;
         jlambda[l] = ml;
         lambda[l] += ml;
      }
      llambda += nk;
      jlambda += nl;
      lambda  += nl;
   }
}

void upRec2(
      double *lambda, double *llambda, const double *tau,
      int nobs, int nk, int nl, bool reg
) {
   double ml = 0;
   double *tau_ = (double *) tau;
   for (int i = 0; i < nobs; i ++) {
      if (!reg) tau_ = (double *) tau;
      for (int l = 0; l < nl; l ++) {
         ml = R_NegInf;
         for (int k = 0; k < nk; k ++)
            ml = log_add_exp(ml, tau_[k] + llambda[k]);

         tau_ += nk;
         lambda[l] += ml;
      }
      llambda += nk;
      lambda  += nl;
   }
}


