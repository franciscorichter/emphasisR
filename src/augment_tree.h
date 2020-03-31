#include <vector>


using node = std::tuple< double, int, double >; // brts, to, t_ext

std::vector<double> calc_sum_diff_std(std::vector<double> input);
std::vector< node > create_tree(const std::vector<double>& brts);
double get_next_bt(const std::vector<node>& tree, 
                   double cbt);

int n_from_time_single_std(double bt,
                           const std::vector<node> tree,
                           double soc);

bool check_i3(double t_ext,
              double tm,
              const std::vector< node >& tree);

double calc_pd_std_single(double tm,
                          const std::vector< node >& tree,
                          double soc);

double speciation_r_single_std(double bt,
                               const std::vector<node>& tree,
                               const std::vector<double>& pars,
                               double soc,
                               double N,
                               std::string model);

double sum_speciation_rate_std(double bt,
                               const std::vector<node>& tree,
                               const std::vector<double>& pars,
                               double soc,
                               std::string model);

double calc_lambda_std(double bt,
                       const std::vector<node>& tree,
                       const std::vector<double>& pars,
                       double soc,
                       double mu,
                       double b,
                       std::string model);

std::vector< std::vector< double > > augment(std::vector< double > brts,
                                             const std::vector< double > parameters,
                                             double  mu,
                                             double soc,
                                             std::string model);