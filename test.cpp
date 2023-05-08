#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

void inside(double* x, double* y) {
   *y = *x;
}
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
   NumericVector y(x.length());
   double *_x_ = x.begin();
   double *_y_ = y.begin();

   for (int i = 0; i < x.length(); i ++) {
      inside(_x_ ++, _y_ ++);
   }

   return y;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(1:5)
*/
