#ifndef RQMC_HPP
#define RQMC_HPP
#include <iostream>
#include <functional>
#include "Normal_distribution_LD.hpp"
#include <random>
#include <numeric>
#include <cmath>        // std::abs

class RQMC_Gaussian {
  int I; // How many different shifted
  int N; // How many samples
  int dimension;
  std::function<double(std::vector<double>)> f;
  std::uniform_real_distribution<> U;
  std::mt19937 gen;
  //std::vector<double> _x; //shift for shifted sequence
  //Gaussian_SH _gaussian;
public:
  std::vector<double> estimator;
  std::vector<double> half_CI;
  RQMC_Gaussian() {};
  ~RQMC_Gaussian(){};
  friend std::ostream& operator<<(std::ostream& os, const RQMC_Gaussian rqmc);
  RQMC_Gaussian(int I, int N, int dimension,std::function<double(std::vector<double>)> f,std::mt19937 gen);
  void operator ()();
};














#endif
