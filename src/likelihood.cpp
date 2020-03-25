#include <Rcpp.h>
#include <string>
#include <numeric>
#include <vector>
#include <math.h>  


using namespace Rcpp;

double calc_ll_rpd1(const std::vector<double>& pars,
                    const std::vector<double>& to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n);

void cout_vec(const std::vector<double>& v, std::string vec_name) {
  Rcout << vec_name << ": ";
  for(auto i : v) {
    Rcout << i << " ";
  }
  Rcout << "\n";
  return;
}

double calc_sum_inte(const std::vector<double>& pars,
                     const std::vector<double>& Pt,
                     const std::vector<double>& brts,
                     const std::vector<double>& n,
                     const std::vector<double>& wt,
                     double mu) {
  
  
  return sum_inte;
}




double calc_ll_rpd5(const std::vector<double>& pars,
                    std::vector<double> to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n,
                    const std::vector<double>& pd) {
  
  if(pars[3] == 0) {
    return(calc_ll_rpd1(pars, to, brts, n));
  }
  
  double mu = std::max(pars[0], 0.0);
  
  std::vector<double> wt(brts.size());
  wt[0] = brts[0] - 0;
  
  std::vector< double > Pt(pd.size() + 1, 0);
  for(int i = 0; i < pd.size(); ++i) {
    Pt[i+1] = pd[i];
    wt[i] = brts[i] - brts[i-1];
  }

  double sum_rho = 0.0;
  for(int i = 0; i < (to.size() - 1); ++i) {
    double to_ = to[i];
    if(to_ == 2) to_ = 1;
    
    double pd2 = Pt[i] + n[i] * wt[i];
    
    double lambda = pars[1] + pars[2] * n[i] + pars[3] * pd2 / n[i];
    if(lambda < 0) lambda = 0;
    
    double rho = lambda * to_ + mu * (1 - to_);
    if(rho > 0) sum_rho += std::log(rho);
  }
  
  double sum_inte = 0.0;
  for(int i = 0; i < brts.size(); ++i)  {
    double brts_i = brts[i];
    double brts_im1 = 0;
    if(i > 0) brts_im1 = brts[i-1];
    
    double c1 = pars[1] + pars[2] * n[i] + (pars[3] / n[i]) * (Pt[i] - brts_im1 * n[i]);
    
    
    if( (brts_im1 > (-c1/pars[3])) & (brts_i < (-c1/pars[3])) ){
      if(pars[3]>0){
        brts_im1 = -c1/pars[3];
      }else{
        brts_i = -c1/pars[3];
      }
    }
    
    sum_inte += n[i]*(mu * wt[i] + c1*(brts_i-brts_im1) + pars[3]*(brts_i * brts_i - 
      brts_im1*brts_im1)/2);
  }
  
  double log_lik = -1 * sum_inte + sum_rho;
  return log_lik;
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
   // return(calc_ll_rpd5c(pars, tree));
  }
  if(model == "rpd5") {
    return(calc_ll_rpd5(pars,to, brts, n, pd));
  }
  return -1;
}




double calc_ll_rpd1(const std::vector<double>& pars,
                    const std::vector<double>& to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n) {
  
  double mu = std::max(pars[0], 0.0);
  
  double sum_rho = 0.0;
  double sum_sigma_over_tree = 0.0;
  
  for (int i = 0; i < n.size(); ++i) {
    double lambda = pars[1] + pars[2] * n[i];
    if(lambda < 0) lambda = 0;
    
    double wt;
    if(i == 0) {
      wt = brts[i] - 0;
    } else {
      wt = brts[i] - brts[i-1];
    }
    
    sum_sigma_over_tree += n[i] * (mu + lambda) * wt;
    
    
    if(i < (n.size() - 1))  {
      int to_ = to[i];
      if(to_ == 2) to_ = 1;
      
      double rho = lambda * to_ + mu * (1 - to_);
      if(rho < 0) rho = 0;
      sum_rho += std::log(rho);
    }
  }
  
  double loglik = -1.0 * sum_sigma_over_tree + sum_rho;
  
  return(loglik);
}


double calc_ll_rpd5_verbose(const std::vector<double>& pars,
                    std::vector<double> to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n,
                    const std::vector<double>& pd) {
  
  if(pars[3] == 0) {
    return(calc_ll_rpd1(pars, to, brts, n));
  }
  
  double mu = std::max(pars[0], 0.0);
  
  std::vector<double> wt(brts.size());
  wt[0] = brts[0] - 0;
  for(int i = 1; i < brts.size(); ++i) {
    wt[i] = brts[i] - brts[i-1];
  }
  
  std::vector< double > Pt(pd.size() + 1, 0);
  for(int i = 0; i < pd.size(); ++i) {
    Pt[i+1] = pd[i];
    if(to[i] == 2) to[i] = 1;
  }
  Pt.pop_back();
  
  std::vector<double> pd2(Pt.size());
  for(int i = 0; i < pd2.size(); ++i) {
    pd2[i] = Pt[i] + n[i] * wt[i];
  }
  
  std::vector<double> brts_i = brts;
  std::vector<double> brts_im1(1, 0);
  for(auto i : brts) {
    brts_im1.push_back(i);
  }
  brts_im1.pop_back();
  
  std::vector<double> lambda(pd2.size());
  for(int i = 0; i < pd2.size(); ++i) {
    lambda[i] = pars[1] + pars[2] * n[i] + pars[3] * pd2[i] / n[i];
    if(lambda[i] < 0)  lambda[i] = 0;
  }
  
  std::vector<double> rho(to.size());
  for(int i = 0; i < rho.size(); ++i) {
    rho[i] = lambda[i] * to[i] + mu * (1 - to[i]);
    if(rho[i] < 0) rho[i] = 0;
    rho[i] = std::log(rho[i]);
  }
  
  std::vector<double> c1(n.size());
  for(int i = 0; i < c1.size(); ++i) {
    c1[i] = pars[1] + pars[2] * n[i] + (pars[3] / n[i]) * (Pt[i] - brts_im1[i] * n[i]);
  }
  
  std::vector<double> inte(brts.size()); 
  for(int i = 0; i < inte.size(); ++i)  {
    if( (brts_im1[i] > (-c1[i]/pars[3])) & (brts_i[i] < (-c1[i]/pars[3])) ){
      if(pars[3]>0){
        brts_im1[i] = -c1[i]/pars[3];
      }else{
        brts_i[i] = -c1[i]/pars[3];
      }
    }
    
    inte[i] = n[i]*(mu*wt[i] + c1[i]*(brts_i[i]-brts_im1[i]) + pars[3]*(brts_i[i] * brts_i[i] - 
      brts_im1[i]*brts_im1[i])/2);
  }
  
  double sum_rho = std::accumulate(rho.begin(), rho.end(), 0.0);
  double sum_inte = std::accumulate(inte.begin(), inte.end(), 0.0);
  double log_lik = -1 * sum_inte + sum_rho;
  return log_lik;
}