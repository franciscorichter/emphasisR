#include <Rcpp.h>
#include <cassert>
#include <numeric>
#include "random_thijs.hpp"

using namespace Rcpp;

enum tree_type {brts_ = 0, 
                to_ext_ = 1, 
                to_ = 2};

void  update_console() {
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}

void cout_vec(const NumericVector& v, std::string vec_name) {
  Rcout << vec_name << ": ";
  for(auto i : v) {
    Rcout << i << " ";
  }
  Rcout << "\n";
  return;
}

void cout_tree(const NumericMatrix& v, std::string vec_name) {
  Rcout << vec_name << ": \n";
  for(int i = 0; i < v.nrow(); ++i) {
    Rcout << v(i, 0) << " " << v(i, 1) << " " <<  v(i, 2) << "\n"; 
  }

  return;
}

void cout_mb(const NumericMatrix& v, std::string vec_name) {
  Rcout << vec_name << ": \n";
  for(int i = 0; i < v.nrow(); ++i) {
    Rcout << v(i, 0) << " " << v(i, 1) << "\n"; 
  }
  
  return;
}

NumericVector n_from_time_vec(NumericVector t_input,
                               NumericMatrix tree,
                               int soc) {
  
  NumericVector to = tree(_, 2);
  to = Rcpp::head(to, -1);
  
 to[ to == 2] = 1;

 NumericVector cs1 = Rcpp::cumsum(to);
 NumericVector cs2 = Rcpp::cumsum(to - 1);

 
 NumericVector n = soc + cs1 + cs2;
 n.push_front(soc);

 NumericVector output(t_input.size(), R_NaReal);
 
 for(int i = 0; i < t_input.size(); ++i) {
   int index = tree.nrow();
   bool found = false;
   for(; index >= 0; --index) {
     if(tree(index, 0) < t_input[i]) {
       found = true;
       break;
     }
   }
   if(found) {
     index++;
     if(index < n.size()) {
      output[i] = n[index];
     }
   } else {
     if(-1 < t_input[i]) {
       output[i] = n[0];
     }
   }
 }
 return(output);
}




double phylodiv_cpp(double tm,
                    NumericMatrix tree,
                    double soc) {
  
  NumericVector i1(tree.nrow(), FALSE);
  NumericVector i2(tree.nrow(), FALSE);
  NumericVector i3(tree.nrow(), FALSE);
  
  i1[ tree(_, brts_) <= tm] = TRUE;
  for(int i = 0; i < i1.size(); ++i) {
    if(tree(i, to_) == 0 && i1[i] == TRUE) {
      i2[i] = TRUE;
    }
  }
  
  NumericVector preselect;
  for(int i = 0; i < i2.size(); ++i) {
    if(i2[i]) {
      preselect.push_back(tree(i, brts_));
    }
  }

  for(int i = 0; i < i3.size(); ++i) {
    auto local_val = tree(i, to_ext_);
    for(auto j : preselect) {
      if(j == local_val) {
        i3[i] = TRUE;
        break;
      }
    }
  }
  
  NumericVector middle;
  middle.push_back(0);
  for(int i = 0; i < tree.nrow(); ++i) {
    if(i1[i] && !i2[i] && !i3[i]) {
      middle.push_back(tree(i, brts_));
    }
  }
  middle.push_back(tm);
 
  auto dt = Rcpp::diff(middle);
  
  NumericVector mult1(dt.size()); // = Rcpp::seq(soc, dt.size() + soc);
  std::iota(mult1.begin(), mult1.end(), soc);
 
  double sum = Rcpp::sum(dt * mult1);
  
  
  return(sum);
}

NumericVector phylodiv_vec_cpp(NumericVector t,
                               NumericMatrix tree,
                               double soc) {
  
  NumericVector output(t.size());
  for(int i = 0;i < t.size(); ++i) {
    output[i] = phylodiv_cpp(t[i], tree, soc);
  }
  return(output);
}



// [[Rcpp::export]]
double speciation_r_cpp(const NumericVector& tm,
                        const NumericMatrix& tree,
                        const NumericVector& pars,
                        double soc) {
  
  NumericVector pd(tm.size(), 0);
  for(int i = 0; i < tm.size(); ++i) {
    pd[i] = phylodiv_cpp(tm[i], tree, soc) - tm[i];
  }
  NumericVector N = n_from_time_vec(tm, tree, soc);
  
  double lambda = 0;
  for(int i = 0; i < pd.size(); ++i) {
    if(!isnan(N[i])) {
      double temp1 = pars[1];
      double temp2 = pars[2] * N[i];
      double temp3 = pars[3] * pd[i] / N[i];
      
      double temp = temp1 + temp2 + temp3;
      if(temp > lambda) lambda = temp;
      //  lambda = max(0, pars[2] + pars[3]*N + pars[4] * pd/N ) 
    }
  }
  
  return(lambda);
}

 
// [[Rcpp::export]]
double sum_speciation_rate(double cbt, 
                           const NumericMatrix& tree, 
                           const NumericVector& pars,
                           double soc) {
  NumericVector input_cbt(1, cbt);
  NumericVector N = n_from_time_vec(input_cbt, tree, soc);  // this seems correct
  
  double lambda = speciation_r_cpp(input_cbt, tree, pars, soc); // this is wrong...
  
  return(N[0] * lambda);
}

double calc_lambda(double cbt, 
                   const NumericMatrix& tree, 
                   const NumericVector& pars, 
                   double soc, 
                   double mu, 
                   double b) {
  
  double sum_spec_rate = sum_speciation_rate(cbt, tree, 
                                             pars, soc);
  
  double mult = 1.0 - exp(-1 * mu * (b - cbt));
  return(sum_spec_rate * mult);
}

NumericVector cum_sum_diff(NumericVector input) {
  //brts = cumsum(-diff(c(brts,0)))
  //auto diff_input = Rcpp::diff(input);
  NumericVector diff_input(input.size());
  input.push_back(0);
  
  for(int i = 1; i < input.size(); ++i) {
    diff_input[i-1] = -1 * (input[i] - input[i-1]);
  }
  return(Rcpp::cumsum(diff_input));
}

NumericMatrix create_tree(const NumericMatrix& missing_branches,
                          const NumericVector& brts) {
  
  if(missing_branches(0, 0) != -1e6) { // -1e6 indicates an empty matrix
    NumericVector bb(missing_branches(_, 0));
    for(int i = 0; i < brts.size(); ++i) {
      bb.push_back(brts[i]);
    }
    for(int i = 0; i < missing_branches.nrow(); ++i) {
      bb.push_back(missing_branches(i,1));
    }
    
    NumericVector t_ext = missing_branches(_, 1);
    for(int i = 0; i < brts.size() + missing_branches.nrow(); ++i) {
      t_ext.push_back(std::numeric_limits<double>::max());
    }
    
    NumericVector to(1, missing_branches.size()); 
    for(int i = 0; i < brts.size(); ++i) {
      to.push_back(2);
    }                 
    for(int i = 0; i < missing_branches.nrow(); ++i) {
      to.push_back(0);
    }
    
    NumericMatrix output(brts.size(), 3);
    output(_, 0) = bb;
    output(_, 1) = t_ext;
    output(_, 2) = to;
  
    return(output);
  } else {
    NumericVector bb = brts;
    NumericVector t_ext(brts.size(), std::numeric_limits<double>::max());
    NumericVector to(brts.size(), 2);
    NumericMatrix output(brts.size(), 3);
    output(_, 0) = bb;
    output(_, 1) = t_ext;
    output(_, 2) = to;
    
    return(output);
  }
}

NumericMatrix sort_tree(NumericMatrix input) {
  
  NumericVector y = input(_, 0);
 
  IntegerVector idx = seq_along(y) - 1;
 
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  
  NumericMatrix output(idx.size(), 3);
  
  for(int i = 0; i < 3; ++i) {
    NumericVector temp = input(_, i);
    temp = temp[idx];
    output(_, i) = temp;
  }
  
  return output;
}

double get_next_bt(const NumericVector& v, double cbt) {
  double m = 1e10;
  for(const auto& i : v) {
    if(i > cbt) {
      if(i < m) m = i;
    }
  }
  return(m);
}

NumericMatrix extend_missing_branches(const NumericMatrix& old,
                                      double s,
                                      double e) {
  
  
  if(old(0, 0) != -1e6) {
    NumericMatrix new_matrix(old.nrow() + 1, 2);
    for(int i = 0; i < old.nrow(); ++i) {
        new_matrix(i,_) = old(i, _);
    }
    
    int j = new_matrix.nrow() - 1;
    new_matrix(j, 0) = s;
    new_matrix(j, 1) = e;
    return(new_matrix);
  } else {
    NumericMatrix new_matrix(1, 2);
    new_matrix(0, 0) = s;
    new_matrix(0, 1) = e;
    return(new_matrix);
  }
}

// [[Rcpp::export]]
List augment_cpp(NumericVector brts_in,
                 NumericVector pars,
                 int soc) {
  
  rnd_t rndgen;
  
  double mu = std::max(0.0, pars(0));
  
  auto brts = cum_sum_diff(brts_in); //  brts = cumsum(-diff(c(brts,0)))
 
   //cout_vec(brts, "brts"); 
  
  double b = *std::max_element(brts.begin(), brts.end()); //   b = max(brts)
 
 
 
   double cbt = 0;
  
  //   missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  // [0] = speciation
  // [1] = extinction
  NumericMatrix missing_branches(1, 2);
  missing_branches(0, 0) = -1e6;
  missing_branches(0, 1) = -1e6; // initialize with fake number to indicate empty.
  
  while(cbt < b) {
    auto tree = create_tree(missing_branches, brts);
   
    // tree = matrix with [ brts, to_ext, to];
    tree = sort_tree(tree);
    double next_bt = get_next_bt(tree(_, brts_), cbt);
    

    
    double lambda_1 = calc_lambda(cbt,     tree, pars, soc, mu, b);
    double lambda_2 = calc_lambda(next_bt, tree, pars, soc, mu, b);
    
    double lambda_max = std::max(lambda_1, lambda_2);
    
    
    
    double u1 = rndgen.uniform();
    double next_speciation_time = cbt - log(u1) / lambda_max;
      
      if(next_speciation_time < next_bt){
        // do something 
        double mult = (1 - exp( -1 * mu * (b - next_speciation_time))) / lambda_max;
        
        double mult2 = sum_speciation_rate(next_speciation_time, tree, pars, soc);
        double pt =  mult2 * mult; 
        
        double u2 = rndgen.uniform();
        
        if (u2 < pt) {
            double extinction_time = next_speciation_time + 
                              rndgen.trunc_exp(0, b - next_speciation_time, mu);
            
            //// CORRECT UP TO HERE /// 
            
            missing_branches = extend_missing_branches(missing_branches, 
                                                       next_speciation_time,
                                                       extinction_time);
            
            if(missing_branches.nrow() > 1000) {
              stop("Current parameters leds to a large number of species");
            }
          }
      }
      cbt = std::min(next_speciation_time,next_bt);
      
      cbt = next_bt;
      if(next_speciation_time < next_bt) cbt = next_speciation_time;
  }
  
  auto tree = create_tree(missing_branches, brts);
  tree = sort_tree(tree);

  NumericVector tree_n = n_from_time_vec(tree(_, brts_),
                                         tree,
                                         soc);
  
  NumericVector tree_pd = phylodiv_vec_cpp(tree(_, brts_),
                                       tree, soc);
    
    
  return List::create( Named("tree") = tree,
                       Named("tree_n") = tree_n,
                       Named("tree_pd") = tree_pd);
}

/*

 
 std::vector<double> n_from_time_vec1(const std::vector<double>& t,
 const std::vector< std::vector< double > >& tree,
 double soc) {
 
 std::vector<double> local_to = tree[to_];
 local_to.pop_back();
 for(auto i = local_to.begin(); i != local_to.end(); ++i) {
 if((*i) == 2) (*i) = 1;
 }
 
 int initspec = soc;
 Rcout << "starting cumsum part\n";
 // this part can potentially be optimized, by only doing one 
 // movement along local_to:
 std::vector<double> cs1(local_to.size());
 std::vector<double> cs2(local_to.size());
 std::partial_sum(local_to.begin(), local_to.end(), cs1.begin());
 for(auto i = local_to.begin(); i != local_to.end(); ++i) {
 (*i)--;
 }
 std::partial_sum(local_to.begin(), local_to.end(), cs2.begin());
 // end of part to be optimized
 
 std::vector< double > n(cs1.size() + 1, 0);
 n[0] = initspec;
 for(int i = 0; i < cs1.size(); ++i) {
 n[i+1] = initspec + cs1[i] + cs2[i];
 }
 
 std::vector<double> output;
 for(auto local_t = t.begin(); local_t != t.end(); ++local_t) {
 int index = static_cast<int>(tree[brts_].size() - 1);
 for(; index >= 0; index--) {
 double val = tree[brts_][index];
 if(val < (*local_t)) {
 break;
 }
 }
 index++; // N = n[max(which(c(-1,tree$brts) < tm))]
 if(index >= n.size()) {
 output.push_back(n.back());
 } else {
 output.push_back(n[index]);
 }
 }
 
 return(output);
 } 
 
 
double n_from_time_single(double t,
                          const std::vector< std::vector< double > >& tree,
                          double soc) {
  
  std::vector<double> local_to = tree[to_];
  local_to.pop_back();
  for(auto& i : local_to) {
    if(i == 2) i = 1;
  }
  
  int initspec = soc;
  std::vector<double> cs1, cs2;
  std::partial_sum(local_to.begin(), local_to.end(), cs1);
  for(auto& i : local_to) {
    i--;
  }
  std::partial_sum(local_to.begin(), local_to.end(), cs2);
  std::vector< double > n;
  n.push_back(initspec);
  for(int i = 0; i < cs1.size(); ++i) {
    double add = initspec + cs1[i] + cs2[i];
    n.push_back(add);
  }
  int index;
  for(index = tree[brts_].size(); index >= 0; index--) {
    if(tree[brts_][index] < t) {
      break;
    }
  }
  index++; // N = n[max(which(c(-1,tree$brts) < tm))]
  
  double N = n[index];
  return N;
}

 // [[Rcpp::export]]
 NumericVector n_from_time_cpp1(NumericVector t,
 NumericMatrix tree,
 int soc) {
 
 std::vector< double > temp;
 std::vector< std::vector< double > > tree_2(3, temp);
 for(int i = 0; i < 3; ++i) {
 std::vector< double > temp;
 for(int j = 0; j < tree.nrow(); ++j) {
 temp.push_back(tree[j, i]);
 }
 tree_2[i] = temp;
 }
 
 std::vector< double > t_input(t.begin(), t.end());
 
 std::vector<double> output = n_from_time_vec1(t_input, tree_2, soc);
 NumericVector output_v(output.begin(), output.end());
 return(output_v);                                     
 }
 
 
 */