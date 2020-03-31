#include <Rcpp.h>
#include <string>
#include <numeric>
#include <vector>
#include <math.h>  

#include "likelihood.h"

using namespace Rcpp;

void cout_vec(const std::vector<double>& v, std::string vec_name) {
  Rcout << vec_name << ": ";
  for(auto i : v) {
    Rcout << i << " ";
  }
  Rcout << "\n";
  return;
}

// [[Rcpp::export]]
double loglik_tree_cpp(std::string model,
                       NumericVector input_pars,
                       NumericMatrix input_tree) {
  
  std::vector< double > pars(input_pars.begin(), input_pars.end());
  
  std::vector< double > to(input_tree.nrow());
  std::vector< double > brts(input_tree.nrow());
  std::vector< double > n(input_tree.nrow());
  std::vector< double > pd(input_tree.nrow());
  
  // input_tree =   [brts, t_ext, to_, pd, n];
  for(int i = 0; i < input_tree.nrow(); ++i) {
    brts[i] = input_tree(i, 0);
    to[i] = input_tree(i, 2);
    n[i]  = input_tree(i, 4);
    pd[i] = input_tree(i, 3);
  }

  if(model == "rpd1") {
    return(calc_ll_rpd1(pars,to, brts, n));
  }
  if(model == "rpd5c") {
   return(calc_ll_rpd5c(pars, to, brts, n, pd));
  }
  if(model == "rpd5") {
    return(calc_ll_rpd5(pars,to, brts, n, pd));
  }
  return -1;
}
