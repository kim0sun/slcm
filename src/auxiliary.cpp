#include <Rcpp.h>
using namespace Rcpp;

int sample1(int n, double *prob) {
   double ran = R::runif(0, 1);
   double cp = 0.0;
   for (int i = 0; i < n; i ++) {
      cp += exp(prob[i]);
      if (ran < cp) return i;
   }
   return n - 1;
}

// [[Rcpp::export]]
double log_add_exp(double x, double y) {
   double ans;
   if (x == R_NegInf) return y;
   if (x > y) {
      ans = x + log(1 + exp(y - x));
   } else {
      ans = y + log(1 + exp(x - y));
   }
   return ans;
}

// [[Rcpp::export]]
double log_sum_exp(NumericVector x) {
   double max_x = max(x);
   double sum = 0.0;
   for(int i = 0; i < x.size(); i ++) {
      sum += exp(x[i] - max_x);
   }
   return log(sum) + max_x;
}

NumericVector elogdiri(NumericVector a) {
   return digamma(a) - R::digamma(sum(a));
}

NumericVector plogdiri(
      NumericVector a, NumericVector b
) {
   return lgammaf(sum(a)) - lgamma(a) +
      sum( (a - 1) * elogdiri(b) );;
}
