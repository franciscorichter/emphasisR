#include <Rcpp.h>
#include <cassert>
#include <numeric>
#include <vector>
#include <tuple>
#include "random_thijs.h"
#include <string>
#include <cmath>

#include "augment_tree.h"

using namespace Rcpp;

// please note that I use -1e6 as an indicator for Inf (Infinity).
// This should be parsed in R later on.






// [[Rcpp::export]]
NumericMatrix augment_cpp(NumericVector brts_in,
                          NumericVector pars,
                          std::string model,
                          int soc) {
  
  double mu = std::max(0.0, pars(0));
  
  std::vector<double> brts(brts_in.begin(), brts_in.end());
  std::vector<double> parameters(pars.begin(), pars.end());
  
  std::vector< std::vector< double > > tree = augment(brts,
                                                      parameters,
                                                      mu,
                                                      soc,
                                                      model);
  
  NumericMatrix output(tree[0].size(), 5);
  for(int i = 0; i < 5; ++i) {
    NumericVector temp(tree[i].begin(), tree[i].end());
    output(_, i) = temp;
  }
  return output;
}
  
