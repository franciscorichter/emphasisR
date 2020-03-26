#include "augment_tree.hpp"
#include <cmath>
#include "random_thijs.hpp"
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;


std::vector<double> calc_sum_diff_std(std::vector<double> input) {
  //   brts = cumsum(-diff(c(brts,0)))
  
  std::vector<double> output(input.size());
  input.push_back(0);
  double sum = 0;
  for(int i = 1; i < input.size(); ++i) {
    sum += (-1 * (input[i] - input[i-1]));
    output[i-1] = sum;
  }
  return(output);
  // VERIFIED TO BE CORRECT IN R //
}

std::vector< node > create_tree(const std::vector<double>& brts) {
  std::vector< node > tree;
  for(int i = 0; i < brts.size(); ++i) {
    tree.push_back(std::make_tuple(brts[i], 2, -1e6)); 
  }
  return(tree);
}

double get_next_bt(const std::vector<node>& tree, 
                   double cbt) {
  // TODO: use find_if with lambda expression
  // somehow this didn't work for now (maybe because only c++11 ?)
  for(const auto& i : tree) {
    if(std::get<0>(i) > cbt) {
      return std::get<0>(i);
    }
  }  
  return std::get<0>(tree.back());
}

int n_from_time_single_std(double bt,
                           const std::vector<node> tree,
                           double soc) {
  // we want the species diversity at t = bt
  int n = soc;
  for(auto i : tree) {
    if (std::get<0>(i) >= bt) break;
    
    if(std::get<1>(i) == 0) {
      n--;
    } else {
      n++;
    }
  }
  return n;
}

bool check_i3(double t_ext,
              double tm,
              const std::vector< node >& tree) {
  // tree$t_ext %in% tree$brts[i2]
  for(auto i : tree) {
    double brts = std::get<0>(i);
    
    if(t_ext == brts) {
      if(brts <= tm && std::get<1>(i) == 0) return true;
      return false;
    }
  }
  
  return false;
}




double calc_pd_std_single(double tm,
                          const std::vector< node >& tree,
                          double soc) {
  // return 0;
  /*  R code:
   i1<-tree$brts  <= tm 
   i2<-tree$to==0 & i1
   i3<-tree$t_ext %in% tree$brts[i2]
   dt<-diff(c(0, tree$brts [i1& !i2& !i3], tm))
   
   return(sum(dt*(soc:(length(dt)+soc-1))))
   */
  
  std::vector<double> dt;
  dt.push_back(0);
  for(auto i : tree) { 
    double brts = std::get<0>(i);
    if(brts > tm) break;  // assuming a tree sorted on brts, all next entries are false
    
    bool i1 = true;
    
    int to = std::get<1>(i);  // tree$to
    bool i2 = to == 0; // i2 = tree$to == 0 & i1
    
    bool i3 = check_i3(std::get<2>(i), tm, tree);
    
    if(i1 && !i2 && !i3) {
      dt.push_back(brts);
    }
  }
  dt.push_back(tm);
  
  double sum = 0.0;
  for(int i = 1; i < dt.size(); ++i) {
    sum += ((dt[i] - dt[i-1])) * (soc + i - 1);
  }
  
  return sum;
}

double speciation_r_single_std(double bt,
                               const std::vector<node>& tree,
                               const std::vector<double>& pars,
                               double soc,
                               double N,
                               std::string model) {
  
  // pd = sapply(tm, phylodiversity,tree=tree,soc=soc)
  //  N = sapply(tm, n_from_time,tree=tree,soc=soc)
  //  lambda = max(0, pars[2] + pars[3]*N + pars[4] * (pd/N) )
  //  return(lambda)
  
  if(model == "rpd1") {
    double lambda = pars[1] + pars[3] * N;
    return lambda;
  } 
  
  double pd = calc_pd_std_single(bt, tree, soc);
  if(model == "rpd5c") pd -= bt;
  
  double lambda = pars[1] + pars[2] * N + pars[3] * pd / N;
  if(lambda < 0) lambda = 0;
  return lambda;
}


double sum_speciation_rate_std(double bt,
                               const std::vector<node>& tree,
                               const std::vector<double>& pars,
                               double soc,
                               std::string model) {
  
  /*
   N = sapply(x, n_from_time,tree=tree,soc=soc)
   speciation_r = get(paste0("lambda.", model))
   lambda = speciation_r(x,tree,pars,soc=soc)
   return(N*lambda)
   */
  
  double N = n_from_time_single_std(bt, tree, soc);
  double lambda = speciation_r_single_std(bt, tree, pars, soc, 
                                          N, model);
  
  return(lambda * N);
}

double calc_lambda_std(double bt,
                       const std::vector<node>& tree,
                       const std::vector<double>& pars,
                       double soc,
                       double mu,
                       double b,
                       std::string model) {
  
  double sum_specrate = sum_speciation_rate_std(bt, tree, pars, soc, model); 
  double mult = 1 - exp(-mu * (b - bt));
  return sum_specrate * mult;
}

std::vector< std::vector< double > > augment(std::vector< double > brts,
                                             const std::vector< double > parameters,
                                             double  mu,
                                             double soc,
                                             std::string model) {
  
  rnd_t rndgen;
  
  brts = calc_sum_diff_std(brts);
  
  double b = *std::max_element(brts.begin(), brts.end());
  
  double cbt = 0.0;
  
  auto tree = create_tree(brts);
  std::sort(tree.begin(), tree.end()); // sort based on first, which is brts.
  
  int num_missing_branches = 0;
  
  while(cbt < b) {
    double next_bt = get_next_bt(tree, cbt);
    
    double lambda1 = calc_lambda_std(cbt, tree, parameters, soc, mu, b, model);
    double lambda2 = calc_lambda_std(next_bt, tree, parameters, soc, mu, b, model);
    double lambda_max = lambda1;
    if(lambda2 > lambda_max) lambda_max = lambda2;
    
    double u1 = rndgen.uniform();
    double next_speciation_time = cbt - std::log(u1) / lambda_max;
    if(next_speciation_time < next_bt) {
      double u2 = rndgen.uniform();
      double pt = calc_lambda_std(next_speciation_time, tree, 
                                  parameters, soc, mu, b, model) / lambda_max;
      if(u2 < pt) {
        double extinction_time = next_speciation_time + rndgen.trunc_exp(0, 
                                                                         b - next_speciation_time,
                                                                         mu);
        node to_add = std::make_tuple(next_speciation_time, 1, extinction_time);
        tree.push_back(to_add);
        node to_add2 = std::make_tuple(extinction_time, 0, -1e6);
        tree.push_back(to_add2);
        
        std::sort(tree.begin(), tree.end()); // sort based on first, which is brts.
        num_missing_branches++;
        if(num_missing_branches > 1000) {
          stop("Current parameters leds to a large number of species");
        }
      }
    }
    cbt = std::min(next_speciation_time,next_bt);
  }
  
  std::vector< std::vector < double > > output(5);
  
  for(auto i : tree) {
    double brts  = std::get<0>(i);   // brts, to, t_ext
    int    to    = std::get<1>(i);
    double t_ext = std::get<2>(i);
    
    double n     = n_from_time_single_std(brts, tree, soc);
    double pd    = calc_pd_std_single(brts, tree, soc);
    
    // brts t_ext            to  n          pd
    output[0].push_back(brts);
    output[1].push_back(t_ext);
    output[2].push_back(to);
    output[3].push_back(n);
    output[4].push_back(pd);
  }
  return output;
}




