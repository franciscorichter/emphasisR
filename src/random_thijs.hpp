#ifndef RANDOM_THIJS_H
#define RANDOM_THIJS_H

  // Copyright 2019 Thijs Janzen

#include <random>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

struct rnd_t {
  std::mt19937 rndgen_;

  rnd_t() {
    std::random_device rd;
    std::mt19937 rndgen_t(rd());
    rndgen_ = rndgen_t;
  }

  std::uniform_real_distribution<float> unif_dist =
    std::uniform_real_distribution<float>(0, 1.0);
  
  
  int random_number(size_t n)    {
    if(n <= 1) return 0;
    return std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(rndgen_);
  }

  float uniform()    {
    return unif_dist(rndgen_);
  }
  
  double trunc_exp(double lower, double upper, double rate) {
    std::exponential_distribution<double> exp_dist(rate);
    double result = exp_dist(rndgen_);
    while(result < lower || result > upper) {
      result = exp_dist(rndgen_);
    }
    return result;
  }
  
  
};

#endif  // RANDOM_THIJS_H
