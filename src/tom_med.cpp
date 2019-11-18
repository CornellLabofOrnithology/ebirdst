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
//' @export
// [[Rcpp::export]]
double tom_med(const Rcpp::NumericVector& x) {
  Rcpp::NumericVector xx = Rcpp::clone(x);
  Rcpp::NumericVector::iterator last = std::remove_if(xx.begin(),xx.end(),R_IsNA);
  Rcpp::NumericVector::iterator start = xx.begin();
  double F;
  const int sz=last-start,middle=sz/2-1;
  if(sz%2==0){
    std::nth_element(start,start+middle,last);
    F=(xx[middle]+*(std::min_element(start+middle+1,last)))/2.0;
  }else{
    std::nth_element(start,start+middle+1,last);
    F=xx[middle+1];
  }
  return F;
}
