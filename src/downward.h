#ifndef DOWNWARD_H
#define DOWNWARD_H

void dnInit(
   double *alpha, double *lambda, double *pi, double *post,
   double *ll, int nobs, int nclass, bool reg
);

double calclr(
   double *lambda, double *pi,
   int nobs, int nclass, bool reg
);

void calclri(
   double *lambda, double *pi, double *ll,
   int nobs, int nclass, bool reg
);

void dnRec(
   double *alpha, double *ualpha,
   double *lambda, double *ulambda, double *jlambda,
   int nobs, int nk, int nl, double *tau,
   double *post, double *joint, double *ll, bool reg
);

#endif
