#include <string>
#include <numeric>
#include <cmath>  

#include "likelihood.h"




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
    if(i > 0) wt[i] = brts[i] - brts[i-1];
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
    
    
    if( (brts_im1 > (-c1/pars[3])) && (brts_i < (-c1/pars[3])) ){
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

double calc_ll_rpd5c(const std::vector<double>& pars,
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
    if(i > 0) wt[i] = brts[i] - brts[i-1];
  }
  
  double sum_rho = 0.0;
  for(int i = 0; i < (to.size() - 1); ++i) {
    double to_ = to[i];
    if(to_ == 2) to_ = 1;
    
    double pd2 = Pt[i] + n[i] * wt[i];
    
    double lambda = pars[1] + pars[2] * n[i] + pars[3] * (pd2 - brts[i]) / n[i];
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
    double c2 = pars[3] * ((n[i] - 1) / n[i]);
    
    
    if( (brts_im1 > (-c1/c2)) && (brts_i < (-c1/c2)) ){
      if(c2 > 0){
        brts_im1 = -c1 / c2;
      }else{
        brts_i   = -c1 / c2;
      }
    }
    
    sum_inte += n[i]*(mu * wt[i] + c1*(brts_i-brts_im1) + c2*(brts_i * brts_i - 
      brts_im1 * brts_im1) / 2);
  }
  
  double log_lik = -1 * sum_inte + sum_rho;
  return log_lik;
}
