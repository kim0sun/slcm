#include <Rcpp.h>
#include "auxiliary.h"
using namespace Rcpp;

void cumPi(
   double *numer, double *denom, double *post,
   int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++) {
         numer[k] = log_add_exp(numer[k], post[k]);
      }
      post += nclass;
   }
   for (int k = 0; k < nclass; k ++) {
      *denom = log_add_exp(*denom, numer[k]);
   }
}

void cumTau(
      double *numer, double *denom, double *joint,
      int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      for (int l = 0; l < nl; l ++) {
         for (int k = 0; k < nk; k ++) {
            numer[k + l * nk] = log_add_exp(numer[k + l * nk], joint[k]);
            denom[l] = log_add_exp(denom[l], joint[k]);
         }
         joint += nk;
      }
   }
}

void cumRho(
      const double *numer, double *denom,
      int *y, int nobs, int nvar, int *ncat,
      int nk, double *post, const double *old_rho
) {
   for (int i = 0; i < nobs; i ++) {
      double *nmr = (double *) numer;
      double *rho = (double *) old_rho;
      for (int k = 0; k < nk; k ++) {
         denom[k] = log_add_exp(denom[k], post[k]);
         for (int m = 0; m < nvar; m ++) {
            if (y[m] > 0) {
               nmr[y[m] - 1] = log_add_exp(nmr[y[m] - 1], post[k]);
            } else {
               for (int r = 0; r < ncat[m]; r ++) {
                  nmr[r] = log_add_exp(nmr[r], post[k] + rho[r]);
               }
            }
            nmr += ncat[m];
            rho += ncat[m];
         }
      }
      post += nk;
      y += nvar;
   }
}


void updatePi(
   double *pi, double *numer, double *denom, int nclass
) {
   for (int k = 0; k < nclass; k ++)
      pi[k] = numer[k] - *denom;
}


void updateTau(
      double *tau, double *numer, double *denom,
      int nk, int nl, int *restr
) {
   for (int l = 0; l < nl; l ++) {
      if (denom[l] == R_NegInf) {
         for (int k = 0; k < nk; k ++) {
            tau[k] = R_NegInf;
         }
      } else {
         for (int k = 0; k < nk; k ++) {
            if (restr[k]) tau[k] = R_NegInf;
            else tau[k] = numer[k] - denom[l];
         }
      }
      tau  += nk;
      numer += nk;
      restr += nk;
   }
}


void updateRho(
      double *rho, double *numer, double *denom,
      int nobs, int nclass, int nvar, int *ncat,
      int *restr
) {
   for (int k = 0; k < nclass; k ++) {
      if (denom[k] == R_NegInf) {
         for (int m = 0; m < nvar; m ++) {
            for (int r = 0; r < ncat[m]; r ++) {
               rho[r] = R_NegInf;
            }
            rho   += ncat[m];
            numer += ncat[m];
            restr += ncat[m];
         }
      } else {
         for (int m = 0; m < nvar; m ++) {
            for (int r = 0; r < ncat[m]; r ++) {
               if (restr[r]) rho[r] = R_NegInf;
               else rho[r] = numer[r] - denom[k];
            }
            rho   += ncat[m];
            numer += ncat[m];
            restr += ncat[m];
         }
      }
   }
}


void updateA(
   double *pi, double *post, int nobs, int nclass
) {
   NumericVector npi(nclass);
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++)
         npi[k] += exp(post[k]);
      post += nclass;
   }
   for (int k = 0; k < nclass; k ++)
      pi[k] = R::digamma(npi[k]) - R::digamma(sum(npi));
}


void updateB(
   double *tau, double *ntau, int nk, int nl,
   int *restr
) {
   for (int l = 0; l < nl; l ++) {
      double stau = 0;
      for (int k = 0; k < nk; k ++) {
         if (restr[k]) continue;
         stau += ntau[k];
      }
      for (int k = 0; k < nk; k ++) {
         if (restr[k]) tau[k] = R_NegInf;
         else if (stau == 0) tau[k] = -log(nk);
         else tau[k] = R::digamma(ntau[k]) - R::digamma(stau);
      }
      ntau += nk;
      tau  += nk;
      restr += nk;
   }
}


void updateC(
      double *rho, double *numer, double *denom,
      int nobs, int nclass, int nvar, int *ncat,
      int *restr
) {
   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < nvar; m ++) {
         for (int r = 0; r < ncat[m]; r ++) {
            if (restr[r]) rho[r] = R_NegInf;
            else
               rho[r] = R::digamma(numer[r]) - R::digamma(denom[m]);
         }
         rho   += ncat[m];
         numer += ncat[m];
         restr += ncat[m];
      }
      denom += nvar;
   }
}
