#include <vector>
#include <string>
#include <numeric>
#include <cmath>  

double calc_ll_rpd1(const std::vector<double>& pars,
                    const std::vector<double>& to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n);

double calc_ll_rpd5(const std::vector<double>& pars,
                    std::vector<double> to,
                    const std::vector<double>& brts,
                    const std::vector<double>& n,
                    const std::vector<double>& pd);

double calc_ll_rpd5c(const std::vector<double>& pars,
                     std::vector<double> to,
                     const std::vector<double>& brts,
                     const std::vector<double>& n,
                     const std::vector<double>& pd);