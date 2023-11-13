#ifndef MSTEP_H
#define MSTEP_H

// M-step
void cumPi(double *numer, double *denom, double *post, int nobs, int nclass);
void cumTau(
   double *numer, double *denom, double *joint, int nobs, int nk, int nl);
void cumRho(
   const double *numer, double *denom,
   int *y, int nobs, int nvar, int *ncat,
   int nk, double *post, const double *old_rho
);
void updatePi(double *pi, double *numer, double *denom, int nclass);
void updateTau(double *tau, double *numer, double *denom, int nk, int nl, int* restr);
void updateRho(
   double *rho, double *numer, double *denom,
   int nobs, int nclass, int nvar, int *ncat,
   int* restr
);
void updateA(double *pi, double *post, int nobs, int nclass);
void updateB(double *tau, double *ntau, int nk, int nl, int* restr);
void updateC(
      double *rho, double *numer, double *denom,
      int nobs, int nclass, int nvar, int *ncat,
      int* restr
);

#endif
