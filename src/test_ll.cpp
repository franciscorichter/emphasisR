//
//  main.cpp
//  test_emphasis_ll
//
//  Created by Thijs Janzen on 25/03/2020.
//  Copyright © 2020 Thijs Janzen. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "likelihood.hpp"

void get_from_file(std::string file_name,
                   std::vector<double>& to,
                   std::vector<double>& brts,
                   std::vector<double>& n,
                   std::vector<double>& pd);


int main(int argc, const char * argv[]) {

  std::vector<double> pars = {0.0625999, 29.9000000, -0.8603454 , 0.1000000};

  std::string model = "rpd5";

  std::vector<double> to;
  std::vector<double> brts;
  std::vector<double> n;
  std::vector<double> pd;

  get_from_file("tree.txt", to, brts, n, pd);



  const double rpd5_ll_verified = -7631.669;
  const double rpd5c_ll_verified = -7593.219;
  const double rpd1_ll_verified = -7143.491;

  double rpd5_ll = calc_ll_rpd5(pars, to, brts, n, pd);
  if(rpd5_ll - rpd5_ll_verified > 1e-3) {
    std::cout << "ERROR ERROR ERROR with rpd5\n";
  }

  double rpd5c_ll = calc_ll_rpd5c(pars, to, brts, n, pd);
  if(rpd5c_ll - rpd5c_ll_verified > 1e-3) {
    std::cout << "ERROR ERROR ERROR with rpd5c\n";
  }

  double rpd1_ll = calc_ll_rpd1(pars, to, brts, n);
  if(rpd1_ll - rpd1_ll_verified > 1e-3) {
    std::cout << "ERROR ERROR ERROR with rpd1\n";
  }

  std::cout << "done\n";

  int num_repl = 1e7;
  for(int i = 0; i < num_repl; ++i) {
    int number = std::rand() % 3;
    switch(number) {
      case 0: {
        rpd5_ll = calc_ll_rpd5(pars, to, brts, n, pd);
        break;
      }
      case 1: {
        rpd5c_ll = calc_ll_rpd5c(pars, to, brts, n, pd);
        break;
      }
      case 2: {
        rpd1_ll = calc_ll_rpd1(pars, to, brts, n);
        break;
      }
    }
  }

  std::cout << "done again\n";

  return 0;
}

void get_from_file(std::string file_name,
std::vector<double>& v_to,
std::vector<double>& v_brts,
std::vector<double>& v_n,
std::vector<double>& v_pd) {
  std::ifstream input(file_name.c_str());
  while(!input.eof()) {
    double brts, t_ext, to, pd, n;
    input >> brts;
    input >> t_ext;
    input >> to;
    input >> pd;
    input >> n;

    v_brts.push_back(brts);
    v_to.push_back(to);
    v_n.push_back(n);
    v_pd.push_back(pd);
  }
  v_brts.pop_back();
  v_to.pop_back();
  v_n.pop_back();
  v_pd.pop_back();
}
